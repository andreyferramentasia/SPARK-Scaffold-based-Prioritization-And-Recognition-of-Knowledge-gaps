# Part IV: bioactivity mining and global occurrence context

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(httr)
  library(jsonlite)
  library(writexl)
  library(readxl)
  library(progress)
  library(mongolite)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(stringr)
  library(grid)
})

run_date <- cfg$run_tag_date %||% Sys.Date()
tag_base <- paste(cfg$taxon_mode, paste(cfg$taxon_values, collapse="-"),
                  format(run_date, "%Y%m%d"), sep="_")
base_tag <- paste0("lotus_", gsub("[^A-Za-z0-9._-]+", "_", tag_base))

stopifnot(dir.exists(OUT_DIR))

cat("Target:", paste(cfg$taxon_values, collapse=", "), "\n")
cat("Output:", OUT_DIR, "\n\n")

if (exists("uni_enriched")) {

  vec_primary_term <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "NA" | x == ""] <- NA_character_
    x <- trimws(sub(";.*", "", gsub("\\|", ";", x)))
    x[is.na(x) | !nzchar(x) | x == "NA"] <- NA_character_
    x
  }

  vec_standardize_class <- function(x) {
    dplyr::case_when(
      is.na(x)                                             ~ NA_character_,
      grepl("Fatty", x, ignore.case = TRUE)                ~ "Fatty acyls",
      tolower(x) == "prenol lipids"                        ~ "Prenol lipids",
      grepl("Glycerophospholipids", x, ignore.case = TRUE) ~ "Glycerophospholipids",
      grepl("Tannin",    x, ignore.case = TRUE)            ~ "Tannins",
      grepl("Coumarin",  x, ignore.case = TRUE)            ~ "Coumarins",
      grepl("Flavonoid", x, ignore.case = TRUE)            ~ "Flavonoids",
      TRUE                                                 ~ x
    )
  }

  uni_enriched <- uni_enriched %>%
    dplyr::mutate(
      .temp_np   = vec_primary_term(chemicalTaxonomyNPclassifierSuperclass),
      .temp_cf   = vec_primary_term(chemicalTaxonomyClassyfireClass),
      .raw_class = dplyr::coalesce(.temp_np, .temp_cf),
      hybrid_class = vec_standardize_class(.raw_class)
    ) %>%
    dplyr::select(-any_of(c(".temp_np", ".temp_cf", ".raw_class")))

  cat("Hybrid class applied. Top classes:\n")
  print(head(sort(table(uni_enriched$hybrid_class), decreasing=TRUE), 5))

} else {
  stop("'uni_enriched' not found. Load data first.")
}


## ChEMBL bioactivity

ACTIVITY_CUTOFF_NM <- 10000
VALID_TYPES <- c("IC50", "Ki", "EC50", "Kd", "MIC", "AC50", "GI50")

inchis_to_map <- uni_enriched %>%
  dplyr::filter(!is.na(inchikey), nzchar(inchikey), inchikey != "NA") %>%
  dplyr::distinct(inchikey) %>%
  dplyr::pull(inchikey)

cat("Unique compounds to query ChEMBL:", length(inchis_to_map), "\n")

get_chembl_ids_batch <- function(inchikeys_chunk) {
  resp <- tryCatch(
    httr::GET(
      "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
      query = list(
        molecule_structures__standard_inchi_key__in = paste(inchikeys_chunk, collapse = ","),
        limit = length(inchikeys_chunk)
      ),
      httr::config(ssl_verifypeer = 0),
      httr::timeout(30)
    ),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)

  json <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE),
    error = function(e) NULL
  )
  if (is.null(json) || !is.data.frame(json$molecules) || nrow(json$molecules) == 0) return(NULL)

  mols   <- json$molecules
  ik_col <- grep("standard_inchi_key", names(mols), value = TRUE, ignore.case = TRUE)[1]
  id_col <- grep("molecule_chembl_id", names(mols), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(ik_col) || is.na(id_col)) return(NULL)

  result <- data.frame(inchikey = mols[[ik_col]], chembl_id = mols[[id_col]],
                       stringsAsFactors = FALSE)
  result[!is.na(result$inchikey) & !is.na(result$chembl_id), ]
}

CHEMBL_BATCH <- 50L
ik_chunks <- split(inchis_to_map, ceiling(seq_along(inchis_to_map) / CHEMBL_BATCH))
cat(sprintf("Mapping %d InChIKeys via ChEMBL (%d requests)... ",
            length(inchis_to_map), length(ik_chunks)))

map_list <- vector("list", length(ik_chunks))
for (i in seq_along(ik_chunks)) {
  map_list[[i]] <- get_chembl_ids_batch(ik_chunks[[i]])
  if (i %% 10 == 0) Sys.sleep(0.3)
}

if (any(!vapply(map_list, is.null, logical(1)))) {
  chembl_map <- dplyr::bind_rows(map_list[!vapply(map_list, is.null, logical(1))]) %>%
    dplyr::distinct(inchikey, chembl_id)
  cat("Done. Mapped:", nrow(chembl_map), "IDs.\n")
} else {
  cat("No ChEMBL IDs found. Proceeding without bioactivity data.\n")
  chembl_map <- data.frame(inchikey = character(0), chembl_id = character(0))
}

activities_list <- list()

if (nrow(chembl_map) > 0) {
  cat("Fetching bioactivities... ")

  get_chembl_bioactivity <- function(chembl_ids) {
    url     <- "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    ids_str <- paste(chembl_ids, collapse = ",")
    resp <- tryCatch({
      httr::GET(url, query = list(molecule_chembl_id__in = ids_str,
                                  standard_type__in = paste(VALID_TYPES, collapse=","),
                                  limit = 1000),
                httr::config(ssl_verifypeer = 0))
    }, error = function(e) NULL)

    if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
    cont <- httr::content(resp, as="text", encoding="UTF-8")
    json <- tryCatch(jsonlite::fromJSON(cont, flatten = TRUE), error = function(e) NULL)
    if (is.null(json) || length(json$activities) == 0) return(NULL)

    json$activities %>%
      dplyr::select(molecule_chembl_id, standard_type, standard_relation,
                    standard_value, standard_units, target_pref_name, target_organism) %>%
      as.data.frame()
  }

  ids_vec <- unique(chembl_map$chembl_id)
  batches <- split(ids_vec, ceiling(seq_along(ids_vec) / 50))

  for (k in seq_along(batches)) {
    act_df <- get_chembl_bioactivity(batches[[k]])
    if (!is.null(act_df)) activities_list[[length(activities_list) + 1]] <- act_df
    if (k %% 5 == 0) Sys.sleep(0.3)
  }
}

if (length(activities_list) > 0) {
  raw_activities <- dplyr::bind_rows(activities_list) %>%
    dplyr::rename(chembl_id = molecule_chembl_id)

  clean_activities <- raw_activities %>%
    dplyr::mutate(standard_value = as.numeric(standard_value)) %>%
    dplyr::filter(!is.na(standard_value), standard_units == "nM",
                  standard_value <= ACTIVITY_CUTOFF_NM) %>%
    dplyr::left_join(chembl_map, by = "chembl_id")

  bio_compound_summary <- clean_activities %>%
    dplyr::group_by(inchikey) %>%
    dplyr::arrange(standard_value) %>%
    dplyr::summarise(
      Best_Potency_nM = first(standard_value),
      Best_Target     = first(target_pref_name),
      Best_Organism   = first(target_organism),
      N_Assays        = n(),
      Evidence_Flag_A = "E1 (Exp. Active)"
    )

  outfile_sum <- file.path(OUT_DIR, paste0(base_tag, "_BIO_A_Summary.xlsx"))
  writexl::write_xlsx(list(Summary = bio_compound_summary, Raw = clean_activities),
                      path = outfile_sum)

  bio_evidence_A <- bio_compound_summary
  cat("Done. Saved:", basename(outfile_sum), "\n")
} else {
  cat("No valid bioactivities (<10uM) found.\n")
  bio_evidence_A <- data.frame(inchikey = character(0))
}

rm(list=intersect(ls(), c("res_list","activities_list","raw_activities","chembl_map","batches")))
gc()


## global occurrence context (MongoDB)

MONGO_URL <- cfg$mongo_url %||% "mongodb://127.0.0.1:27017"
lotus_db  <- mongo(collection = "lotusUniqueNaturalProduct", db = "lotus", url = MONGO_URL)

target_inchikeys <- unique(na.omit(uni_enriched$inchikey))
cat("Querying global occurrence for", length(target_inchikeys), "compounds.\n")

extract_family_nuclear <- function(tax_obj) {
  if (is.null(tax_obj)) return(character(0))
  flat <- unlist(tax_obj)
  idx  <- grep("family$", names(flat), ignore.case = TRUE)
  if (length(idx) > 0) {
    fams <- as.character(flat[idx])
    return(unique(fams[nzchar(fams) & fams != "NA"]))
  }
  return(character(0))
}

BATCH_SIZE   <- 2000
chunks       <- split(target_inchikeys, ceiling(seq_along(target_inchikeys)/BATCH_SIZE))
results_list <- list()

pb <- progress::progress_bar$new(
  format = "   Processing [:bar] :percent | Batch :current/:total | ETA: :eta",
  total  = length(chunks), width = 70
)

for (i in seq_along(chunks)) {
  pb$tick()

  ids_vec  <- chunks[[i]]
  ids_json <- if (length(ids_vec) == 1) paste0('["', ids_vec, '"]') else
    jsonlite::toJSON(ids_vec, auto_unbox = TRUE)
  q_json   <- sprintf('{"inchikey": {"$in": %s}}', ids_json)

  tryCatch({
    b_data <- lotus_db$find(
      query  = q_json,
      fields = '{"inchikey": 1, "taxonomyReferenceObjects": 1, "_id": 0}'
    )

    if (nrow(b_data) > 0) {
      row_indices <- seq_len(nrow(b_data))

      families_list <- lapply(row_indices, function(idx) {
        if (!is.null(b_data$taxonomyReferenceObjects)) {
          if (is.data.frame(b_data$taxonomyReferenceObjects)) {
            return(extract_family_nuclear(b_data$taxonomyReferenceObjects[idx, , drop=FALSE]))
          } else if (is.list(b_data$taxonomyReferenceObjects)) {
            return(extract_family_nuclear(b_data$taxonomyReferenceObjects[[idx]]))
          }
        }
        return(character(0))
      })

      counts  <- sapply(families_list, length)
      strings <- sapply(families_list, function(x) {
        if (length(x) == 0) return("Unknown")
        paste(head(as.character(x), 5), collapse = "; ")
      })

      results_list[[i]] <- data.frame(
        inchikey            = b_data$inchikey,
        Global_Family_Count = counts,
        Family_List_String  = strings,
        stringsAsFactors    = FALSE
      )
    }
  }, error = function(e) {
    warning(paste("Batch", i, "error:", e$message))
  })

  if (i %% 10 == 0) gc()
}

df_pass1 <- dplyr::bind_rows(results_list)

missing_keys <- setdiff(target_inchikeys, df_pass1$inchikey)

if (length(missing_keys) > 0) {
  cat("\n  Filling", length(missing_keys), "missing compounds with 0 global families...\n")
  df_missing <- data.frame(
    inchikey            = missing_keys,
    Global_Family_Count = 0,
    Family_List_String  = "NOT_IN_DB_OR_NO_FAMILY",
    stringsAsFactors    = FALSE
  )
  df_pass1 <- dplyr::bind_rows(df_pass1, df_missing)
}

global_context_C <- df_pass1 %>%
  dplyr::distinct(inchikey, .keep_all = TRUE) %>%
  dplyr::mutate(Biogeography_Status = dplyr::case_when(
    Global_Family_Count == 0  ~ "Not Found / Error",
    Global_Family_Count == 1  ~ "Exclusive (1 Family)",
    Global_Family_Count <= 10 ~ "Restricted",
    TRUE                      ~ "Ubiquitous"
  ))

outfile <- file.path(OUT_DIR, paste0(base_tag, "_BIO_C_Global_Context.xlsx"))
writexl::write_xlsx(global_context_C, outfile)
cat("Global context saved:", basename(outfile), "\n")


## integration and prioritization

TOP_N_LABELS <- 10

blacklist_terms <- c(
  "oleic", "linoleic", "palmitic", "stearic", "myristic",
  "ergosterol", "sitosterol", "stigmasterol", "campesterol", "cholesterol",
  "quercetin", "rutin", "kaempferol", "catechin", "epicatechin",
  "gallic acid", "ellagic acid", "chlorogenic", "caffeic",
  "lupeol", "amyrin", "friedelin", "betulin",
  "chlorophyll", "carotene", "squalene", "tocopherol", "sucrose", "glucose"
)
BLACKLIST_REGEX <- paste(blacklist_terms, collapse = "|")

scaffold_counts <- uni_enriched %>%
  filter(!is.na(murko_framework), murko_framework != "") %>%
  count(murko_framework, name = "Local_Scaffold_Freq")

master <- uni_enriched %>%
  select(inchikey, smiles, iupac_name, molecular_formula, molecular_weight,
         xlogp, heavy_atom_number,
         Primary_Class  = hybrid_class,
         Scaffold       = murko_framework,
         Taxon_Families = family, Taxon_Genera = genus, Taxon_Species = species) %>%
  left_join(scaffold_counts, by=c("Scaffold"="murko_framework")) %>%
  left_join(bio_evidence_A,  by="inchikey") %>%
  left_join(global_context_C, by="inchikey") %>%
  mutate(
    Global_Family_Count = ifelse(is.na(Global_Family_Count), 0, Global_Family_Count),
    Best_Potency_nM     = suppressWarnings(as.numeric(Best_Potency_nM)),
    Has_Activity        = !is.na(Best_Potency_nM),
    pIC50               = ifelse(Has_Activity & Best_Potency_nM > 0,
                                 -log10(Best_Potency_nM * 1e-9), 0),
    xlogp_val           = suppressWarnings(as.numeric(xlogp)),
    LLE_Score           = ifelse(!is.na(xlogp_val) & Has_Activity, pIC50 - xlogp_val, -99),
    Is_Primary          = grepl(BLACKLIST_REGEX, iupac_name, ignore.case = TRUE),
    Score_Local         = case_when(replace_na(Local_Scaffold_Freq, 1) == 1 ~ 3,
                                    Local_Scaffold_Freq <= 5 ~ 1, TRUE ~ 0),
    Score_Global        = case_when(Global_Family_Count == 0  ~ -10,
                                    Global_Family_Count == 1  ~  3,
                                    Global_Family_Count <= 3  ~  1, TRUE ~ -5),
    Score_Potency       = case_when(Has_Activity & (LLE_Score > 2 | Best_Potency_nM < 200) ~ 3,
                                    Has_Activity ~ 1, TRUE ~ 0),
    PRIORITY_SCORE      = Score_Local + Score_Global + (Score_Potency * 2),
    Class = case_when(
      Is_Primary                                                              ~ "DISCARD (Primary Metabolite)",
      Global_Family_Count == 0                                                ~ "DISCARD (Unverified Data)",
      Global_Family_Count > 10                                                ~ "DISCARD (Ubiquitous)",
      Has_Activity & Score_Global > 0 & (LLE_Score > 0 | Best_Potency_nM < 500) ~ "STAR (Exclusive & Active)",
      Has_Activity & Score_Global > 0                                         ~ "HIT (Low Efficiency)",
      Score_Global > 0                                                        ~ "HIDDEN GEM (High Novelty)",
      TRUE                                                                    ~ "BASELINE"
    )
  ) %>%
  arrange(desc(PRIORITY_SCORE))

outfile <- file.path(OUT_DIR, paste0(base_tag, "_MASTER_LIST_GLOBAL.xlsx"))
writexl::write_xlsx(
  list(
    All_Ranked = master,
    Top_STARS  = filter(master, grepl("STAR", Class)) %>% head(50),
    Top_GEMS   = filter(master, grepl("GEM", Class))  %>% head(50)
  ),
  path = outfile
)

cat("Integration complete. Master list saved.\n")
cat("  STARs identified:", sum(grepl("STAR", master$Class)), "\n")


## prioritization matrix plot

if (exists("master") && nrow(master) > 0) {
  suppressPackageStartupMessages({
    library(ggplot2); library(ggrepel); library(dplyr); library(scales); library(grid)
  })

  set.seed(cfg$seed %||% 42)

  plot_data <- master %>%
    dplyr::mutate(
      Novelty_Score_Raw  = Score_Local + Score_Global,
      Activity_Score_Raw = Score_Potency,
      Primary_Class      = ifelse(is.na(Primary_Class) | Primary_Class == "",
                                  "Unclassified", Primary_Class),
      X_Jitter           = Novelty_Score_Raw  + runif(dplyr::n(), -0.3, 0.3),
      Y_Jitter           = Activity_Score_Raw + runif(dplyr::n(), -0.3, 0.3)
    )

  top_classes <- plot_data %>%
    dplyr::count(Primary_Class, sort=TRUE) %>% head(5) %>% dplyr::pull(Primary_Class)
  plot_data <- plot_data %>%
    dplyr::mutate(Visual_Class = factor(
      ifelse(Primary_Class %in% top_classes, Primary_Class, "Other Classes"),
      levels = c(top_classes, "Other Classes")
    ))

  labels_star <- plot_data %>%
    filter(grepl("STAR", Class)) %>%
    arrange(desc(PRIORITY_SCORE), inchikey) %>%
    head(TOP_N_LABELS)

  labels_gem <- plot_data %>%
    filter(grepl("GEM", Class)) %>%
    arrange(desc(PRIORITY_SCORE), inchikey) %>%
    head(TOP_N_LABELS)

  cat("  Labels selected:", nrow(labels_star), "STARs and", nrow(labels_gem), "GEMs.\n")

  cols_npg <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#7E6148","#B09C85")
  if (length(levels(plot_data$Visual_Class)) <= length(cols_npg))
    names(cols_npg) <- levels(plot_data$Visual_Class)

  p_final <- ggplot(plot_data, aes(x = X_Jitter, y = Y_Jitter)) +
    annotate("rect", xmin=2.5, xmax=Inf, ymin=0.5,  ymax=Inf,  fill="#00A087", alpha=0.05) +
    annotate("rect", xmin=2.5, xmax=Inf, ymin=-Inf, ymax=0.5,  fill="#E64B35", alpha=0.05) +
    geom_vline(xintercept=2.5, linetype="longdash", color="grey60") +
    geom_hline(yintercept=0.5, linetype="longdash", color="grey60") +
    geom_point(aes(fill=Visual_Class, size=molecular_weight),
               shape=21, color="white", stroke=0.3, alpha=0.85) +
    scale_fill_manual(values=cols_npg, name="Chemical Class") +
    scale_size_continuous(range=c(2.5, 7), name="Molecular Weight") +
    guides(size=guide_legend(override.aes=list(fill="grey60")),
           fill=guide_legend(override.aes=list(size=4))) +
    annotate("text", x=5.5, y=3.3,  label="Q1: HIGH PRIORITY",
             fontface="bold", color="#2F4F4F", size=4) +
    annotate("text", x=5.5, y=-0.2, label="Q3: DARK MATTER",
             fontface="bold", color="#A52A2A", size=4) +
    geom_label_repel(data=labels_star, aes(label=inchikey),
                     size=2.0, fill="white", box.padding=0.3, point.padding=0.2,
                     segment.color="grey50", force=20, max.overlaps=Inf) +
    geom_label_repel(data=labels_gem, aes(label=inchikey),
                     size=2.0, fill="white", box.padding=0.3, point.padding=0.2,
                     segment.color="grey50", force=20, max.overlaps=Inf) +
    theme_classic(base_size=14) + coord_cartesian(clip="off") +
    labs(x="Chemotaxonomic Novelty Score", y="Pharmacological Evidence Score") +
    theme(panel.border=element_rect(colour="black", fill=NA, size=0.8),
          axis.line=element_blank(), legend.position="right")

  ggsave(file.path(OUT_DIR, paste0(base_tag, "_Matrix_NatureStyle_V4.pdf")),
         plot=p_final, width=11, height=8)
  cat("Matrix plot saved.\n")
}


## extraction profiler

if (nrow(master) > 0) {
  props_ref <- uni_enriched %>%
    dplyr::select(inchikey, xlogp, number_of_carbons, number_of_oxygens,
                  contains_sugar, iupac_name) %>%
    dplyr::distinct(inchikey, .keep_all = TRUE)

  prof_data <- master %>%
    dplyr::select(inchikey, Class, PRIORITY_SCORE, Primary_Class) %>%
    dplyr::inner_join(props_ref, by = "inchikey") %>%
    dplyr::filter(grepl("STAR|HIDDEN GEM", Class)) %>%
    dplyr::mutate(
      xlogp_raw  = suppressWarnings(as.numeric(xlogp)),
      nC         = as.numeric(number_of_carbons),
      nO         = as.numeric(number_of_oxygens),
      OC_Ratio   = ifelse(nC > 0, nO / nC, 0),
      name_lower = tolower(coalesce(iupac_name, "")),
      Is_Glycoside = (contains_sugar == TRUE) |
        grepl("glucoside|rhamnoside|arabinoside|rutinoside|glycoside",
              name_lower, fixed = FALSE),
      Label_Name = inchikey
    ) %>%
    dplyr::filter(!is.na(xlogp_raw))

  prof_data_smart <- prof_data %>%
    mutate(
      xlogp_corrected = case_when(
        Is_Glycoside & xlogp_raw > 1.5   ~ 0.5,
        OC_Ratio > 0.6 & xlogp_raw > 2.0 ~ 1.0,
        TRUE                              ~ xlogp_raw
      )
    )

  top_cp    <- prof_data_smart %>% count(Primary_Class, sort=TRUE) %>% head(6) %>% pull(Primary_Class)
  plot_prof <- prof_data_smart %>% filter(Primary_Class %in% top_cp)

  connectors <- bind_rows(
    prof_data_smart %>% filter(grepl("STAR", Class)) %>%
      arrange(desc(PRIORITY_SCORE)) %>% head(TOP_N_LABELS),
    prof_data_smart %>% filter(grepl("GEM",  Class)) %>%
      arrange(desc(PRIORITY_SCORE)) %>% head(TOP_N_LABELS)
  ) %>% mutate(Type = ifelse(grepl("STAR", Class), "STAR", "GEM"))

  zones <- data.frame(
    xmin   = c(-Inf, 0.5, 2.5, 4.5),
    xmax   = c(0.5,  2.5, 4.5, Inf),
    Solvent = c("Water/Hydroalc.", "Ethanol", "EtOAc", "Hexane"),
    Fill    = c("#E3F2FD", "#BBDEFB", "#C8E6C9", "#FFF9C4")
  )
  zones$Solvent <- factor(zones$Solvent,
                          levels=c("Water/Hydroalc.", "Ethanol", "EtOAc", "Hexane"))

  p_prof <- ggplot() +
    geom_rect(data=zones,
              aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=Solvent), alpha=0.5) +
    geom_density(data=plot_prof, aes(x=xlogp_corrected), fill="grey40", alpha=0.4) +
    geom_point(data=connectors, aes(x=xlogp_corrected, y=0.02, color=Type),
               size=4, shape=18) +
    geom_label_repel(
      data=connectors,
      aes(x=xlogp_corrected, y=0.02, label=Label_Name, color=Type),
      size=2.2, nudge_y=0.15, direction="y", max.overlaps=50, max.time=2, max.iter=2000
    ) +
    facet_wrap(~Primary_Class, ncol=2, scales="free_y") +
    scale_fill_manual(values=setNames(zones$Fill, zones$Solvent),
                      name="Extraction Solvent") +
    scale_color_manual(values=c("STAR"="#2E7D32", "GEM"="#C62828"), name="Priority") +
    coord_cartesian(xlim=c(-3, 8)) + theme_bw(base_size=12) +
    labs(title="Target-Oriented Extraction Profiler",
         x="Polarity Index", y="Relative Density") +
    theme(legend.position="bottom", panel.grid=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())

  ggsave(file.path(OUT_DIR, paste0(base_tag, "_Profiler_V16_Hydroalcoholic.pdf")),
         plot=p_prof, width=10, height=11)
  cat("Profiler plot saved.\n")
}


## top 25 discussion export

if (exists("master") && exists("uni_enriched")) {
  enrich_top_list <- function(ranked_df) {
    ranked_df %>%
      dplyr::select(inchikey, PRIORITY_SCORE, Class, Score_Global, Score_Local,
                    Score_Potency, Global_Family_Count, Family_List_String,
                    Best_Potency_nM, Best_Target, LLE_Score, Primary_Class) %>%
      dplyr::left_join(uni_enriched %>% select(-hybrid_class), by = "inchikey") %>%
      dplyr::select(inchikey, iupac_name, PRIORITY_SCORE, Class, Primary_Class,
                    Global_Family_Count, everything())
  }

  top25_stars <- master %>%
    filter(grepl("STAR", Class)) %>%
    arrange(desc(PRIORITY_SCORE)) %>% head(25) %>% enrich_top_list()
  top25_gems  <- master %>%
    filter(grepl("GEM",  Class)) %>%
    arrange(desc(PRIORITY_SCORE)) %>% head(25) %>% enrich_top_list()

  out_file_disc <- file.path(OUT_DIR, paste0(base_tag, "_Top25_Discussion_FullData.xlsx"))
  writexl::write_xlsx(list(Top25_STARS = top25_stars, Top25_GEMS = top25_gems),
                      path = out_file_disc)

  cat("Discussion data saved:", basename(out_file_disc), "\n")
}

cat("\nPipeline completed.\n")
