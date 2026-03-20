# Part I: extract LOTUS compounds for target taxon, enrich with properties, optionally normalize via WFO

suppressPackageStartupMessages({
  library(mongolite)
  library(jsonlite)
  library(dplyr)
  library(progress)
  library(data.table)
  library(stringr)
  library(writexl)
  library(readr)
  library(stringi)
})

options(OutDec = ".", scipen = 999)

## helpers

fix_ref_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[!nzchar(x)] <- NA_character_
  x <- gsub("\\$x\\$x\\$", ".", x, perl = TRUE)
  x
}

norm_ascii   <- function(x) stringi::stri_trans_general(x, "Latin-ASCII")
tidy_space   <- function(x) trimws(gsub("\\s+", " ", x))

title_case_1 <- function(x){
  x <- as.character(x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

canon_genus <- function(x){
  x <- tidy_space(norm_ascii(x))
  x <- sub("\\s+.*$", "", x)
  x <- gsub("[^A-Za-z-]", "", x)
  ifelse(nzchar(x), title_case_1(x), NA_character_)
}

is_binomial <- function(x){
  x <- tidy_space(as.character(x))
  grepl("^\\S+\\s+\\S+", x)
}

# "Ocoteaodorifera" -> "Ocotea odorifera"
fix_glued_species <- function(genus, species){
  genus   <- as.character(genus)
  species <- as.character(species)
  needs   <- !is.na(genus) & !is.na(species) &
    mapply(function(g, s) grepl(paste0("^", g, "[A-Za-z]"), s), genus, species)
  if (any(needs)) {
    species[needs] <- mapply(function(g, s) sub(paste0("^(", g, ")([A-Za-z])"), "\\1 \\2", s),
                             genus[needs], species[needs])
  }
  species
}

## config

normalize_taxon_mode <- function(x){
  x0 <- tolower(trimws(as.character(x %||% "")))
  map <- c(
    "family"="family","families"="family","fam"="family",
    "genus"="genus","genero"="genus","gênero"="genus","gen"="genus",
    "species"="species","specie"="species","sp"="species","especie"="species","espécie"="species",
    "kingdom"="kingdom","reino"="kingdom","kingdoms"="kingdom"
  )
  x1 <- map[[x0]]
  if (is.null(x1)) stop(
    sprintf("Invalid TAXON_MODE: '%s'. Use: family | genus | species | kingdom.", x0),
    call. = FALSE
  )
  x1
}

TAXON_MODE   <- normalize_taxon_mode(cfg$taxon_mode %||% "family")
TAXON_VALUES <- cfg$taxon_values %||% character(0)
TAXON_VALUES <- unique(na.omit(trimws(as.character(TAXON_VALUES))))
if (length(TAXON_VALUES) == 0L) stop("TAXON_VALUES is empty after normalization.", call.=FALSE)

MONGO_URL <- cfg$mongo_url %||%
  "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000"
DB_NAME    <- cfg$db_name    %||% "lotus"
COLL_NAME  <- cfg$coll_name  %||% "lotusUniqueNaturalProduct"
opts       <- cfg$mongo_opts %||% '{"allowDiskUse": true, "batchSize": 5000}'

PAGE       <- cfg$page_size_lines     %||% 50000L
chunk_size <- cfg$chunk_size_inchikey %||% 1000L

safe_tag <- function(mode, values, run_date, suffix = NULL) {
  tag <- paste(mode,
               paste(values, collapse = "-"),
               format(run_date, "%Y%m%d"),
               sep = "_")
  if (!is.null(suffix) && nzchar(suffix)) {
    tag <- paste0(tag, "_", gsub("[^A-Za-z0-9._-]+","_",suffix))
  }
  gsub("[^A-Za-z0-9._-]+", "_", tag)
}

tag_base <- if (!is.null(cfg$prefix_base_tag) && nzchar(cfg$prefix_base_tag)) {
  cfg$prefix_base_tag
} else {
  safe_tag(
    TAXON_MODE,
    TAXON_VALUES,
    cfg$run_tag_date %||% Sys.Date(),
    cfg$custom_tag_suffix %||% NULL
  )
}

out_dir_base <- cfg$out_dir_base %||% getwd()
OUT_DIR <- file.path(out_dir_base, paste0("lotus_", tag_base))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

safe_file <- function(base, ext) file.path(OUT_DIR, paste0(base, ext))
base_tag  <- paste0("lotus_", tag_base)

PROPS_CORE_FIELDS <- cfg$props_core_fields %||% c(
  "lotus_id","wikidata_id","inchikey","smiles","iupac_name",
  "molecular_formula","molecular_weight",
  "xlogp","alogp","amralogp","manholdlogp",
  "topoPSA","tpsaEfficiency","fsp3",
  "hBondAcceptorCount","hBondDonorCount","LipinskiRuleOf5Failures",
  "contains_sugar","contains_ring_sugars","contains_linear_sugars",
  "number_of_carbons","number_of_oxygens","number_of_nitrogens",
  "total_atom_number","heavy_atom_number","max_number_of_rings","min_number_of_rings",
  "murko_framework","ertlFunctionalFragmentsPseudoSmiles",
  "chemicalTaxonomyNPclassifierSuperclass","chemicalTaxonomyNPclassifierClass",
  "chemicalTaxonomyClassyfireSuperclass","chemicalTaxonomyClassyfireClass",
  "traditional_name","allWikidataIds"
)

collapse <- function(x) {
  x <- unique(x); x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_character_)
  paste(x, collapse=";")
}
safe_first <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_character_)
  x[1]
}

## connect

lotus <- mongo(collection = COLL_NAME, db = DB_NAME, url = MONGO_URL)
cat("Part I: connected, total docs:", lotus$count(), "\n")
cat(sprintf("TAXON_MODE = %s | TAXON_VALUES = %s\n",
            TAXON_MODE, paste(TAXON_VALUES, collapse=", ")))

## taxon filter

# regex special chars in taxon names need escaping before using as MongoDB $regex
regex_escape <- function(x) {
  x <- as.character(x); x <- stringr::str_squish(x)
  gsub("([\\.^$|()?*+\\[\\]{}-])", "\\\\\\1", x, perl = TRUE)
}

build_taxon_match <- function(prefix, mode, values) {
  vals <- unique(na.omit(trimws(values)))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("Empty/invalid TAXON_VALUES.")

  or_list <- list()
  for (v0 in vals) {
    v <- regex_escape(v0)
    if (mode == "genus") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")), paste0(prefix, "genus"))),
        list(setNames(list(list("$regex" = paste0("^", v, "\\s"), "$options" = "i")), paste0(prefix, "species")))
      )
    } else if (mode == "family") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")),
                      paste0(prefix, "family")))
      )
    } else if (mode == "species") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "species"))),
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "organism_value"))),
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "cleaned_organism_id")))
      )
    } else if (mode == "kingdom") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")),
                      paste0(prefix, "kingdom")))
      )
    } else {
      stop("Invalid TAXON_MODE.")
    }
  }
  list("$or" = or_list)
}

## count matching rows

pipe_count <- list(
  list("$project" = list(tx1 = list("$objectToArray" = "$taxonomyReferenceObjects"))),
  list("$unwind"  = "$tx1"),
  list("$project" = list(tx2 = list("$objectToArray" = "$tx1.v"))),
  list("$unwind"  = "$tx2"),
  list("$unwind"  = "$tx2.v"),
  list("$match"   = build_taxon_match("tx2.v.", TAXON_MODE, TAXON_VALUES)),
  list("$count"   = "n")
)
cnt <- lotus$aggregate(jsonlite::toJSON(pipe_count, auto_unbox = TRUE), options = opts)
total_lines <- if (nrow(cnt)) cnt$n[1] else 0L
cat("Total matched rows (compound x ref x source x species):", total_lines, "\n")
if (total_lines == 0) stop("Nothing to extract for this taxonomic criterion.")

## extract lin table

n_batches <- ceiling(total_lines / PAGE)
pb <- progress::progress_bar$new(
  format = "Batch :current/:total [:bar] :percent | rows=:rows | :elapsed (ETA :eta)",
  total = n_batches, clear = FALSE, width = 80
)
processed <- 0L; pages <- vector("list", n_batches)

for (i in seq_len(n_batches)) {
  skip_rows <- (i - 1L) * PAGE
  pipe_page <- list(
    list("$project" = list(`_id`=0, lotus_id=1, smiles=1, inchikey=1, iupac_name=1, molecular_formula=1,
                           tx1 = list("$objectToArray" = "$taxonomyReferenceObjects"))),
    list("$unwind"  = "$tx1"),
    list("$project" = list(lotus_id=1, smiles=1, inchikey=1, iupac_name=1, molecular_formula=1,
                           ref_id = "$tx1.k", srcMap = "$tx1.v")),
    list("$project" = list(lotus_id=1, smiles=1, inchikey=1, iupac_name=1, molecular_formula=1, ref_id=1,
                           src = list("$objectToArray" = "$srcMap"))),
    list("$unwind"  = "$src"),
    list("$unwind"  = "$src.v"),
    list("$match"   = build_taxon_match("src.v.", TAXON_MODE, TAXON_VALUES)),
    list("$project" = list(lotus_id=1, smiles=1, inchikey=1, iupac_name=1, molecular_formula=1,
                           ref_id=1, source="$src.k",
                           family="$src.v.family", genus="$src.v.genus", species="$src.v.species")),
    list("$sort"    = list("lotus_id"=1L, "ref_id"=1L, "source"=1L, "species"=1L)),
    list("$skip"    = skip_rows),
    list("$limit"   = PAGE)
  )
  page <- lotus$aggregate(jsonlite::toJSON(pipe_page, auto_unbox = TRUE), options = opts)
  if (nrow(page)) {
    pages[[i]] <- data.table::as.data.table(page)
    processed <- processed + nrow(page)
  } else {
    pages[[i]] <- NULL
  }
  pb$tick(tokens = list(rows = format(processed, big.mark = ",", decimal.mark = ".", scientific = FALSE)))
}

lin <- data.table::rbindlist(pages, use.names = TRUE, fill = TRUE)
data.table::setDF(lin)
lin <- dplyr::distinct(lin, lotus_id, ref_id, source, family, genus, species,
                       inchikey, smiles, iupac_name, molecular_formula, .keep_all = TRUE)

cat("\nRows loaded into memory:", nrow(lin), "\n")
lin <- lin %>%
  dplyr::mutate(ref_id = fix_ref_id(ref_id))

## genus-mode consistency filter
# ensures retrieved species actually belong to the target genus
if (identical(TAXON_MODE, "genus")) {
  target_gen <- unique(tolower(trimws(TAXON_VALUES)))
  target_gen <- target_gen[nzchar(target_gen)]
  esc <- function(x) gsub("([\\.^$|()?*+\\[\\]{}-])","\\\\\\1", x, perl=TRUE)
  pat_gen  <- paste0("^(", paste(esc(target_gen), collapse="|"), ")$")
  pat_pref <- paste0("^(", paste(esc(target_gen), collapse="|"), ")\\s")

  lin <- lin %>%
    dplyr::mutate(
      .g0   = tolower(trimws(genus)),
      .s0   = tolower(trimws(species)),
      .s_ok = !is.na(.s0) & grepl(pat_pref, .s0, perl = TRUE),
      .g_ok = !is.na(.g0) & grepl(pat_gen,  .g0, perl = TRUE)
    )

  keep_mask <- with(lin, (.s_ok) | (.g_ok & (is.na(.s0) | .s_ok)))

  drop_df <- lin[!keep_mask, c("lotus_id","ref_id","source","family","genus","species")]
  if (nrow(drop_df)) {
    diag_dir <- file.path(OUT_DIR, "diagnostics")
    if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      drop_df,
      file = file.path(diag_dir, paste0(base_tag, "_genus_mode_inconsistent_pairs.tsv")),
      sep = "\t", quote = TRUE
    )
    message(sprintf("[FILTER] genus-mode: %d rows removed due to genus-species inconsistency (see diagnostics/).",
                    nrow(drop_df)))
  }

  lin <- lin[keep_mask, , drop = FALSE]
  lin$.g0 <- lin$.s0 <- lin$.s_ok <- lin$.g_ok <- NULL
}

## retrieve compound properties

fields_props <- sprintf(
  '{"_id":0,%s}',
  paste(sprintf('"%s":1', PROPS_CORE_FIELDS), collapse = ",")
)

scalarize_field <- function(doc, nm) {
  val <- doc[[nm]]
  if (is.null(val)) return(NA_character_)
  if (is.atomic(val) && length(val) == 1) return(as.character(val))
  v <- tryCatch(unlist(val, use.names = FALSE), error = function(e) val)
  v <- v[!is.na(v)]
  v <- v[is.atomic(v)]
  if (!length(v)) return(NA_character_)
  paste(unique(as.character(v)), collapse = ";")
}

inchis <- unique(stats::na.omit(lin$inchikey))

if ((!length(inchis) || all(is.na(inchis))) &&
    isTRUE(cfg$stop_if_no_inchikey %||% TRUE)) {
  stop("No valid InChIKeys found in 'lin'. Cannot enrich properties.")
}

props_list <- vector("list", ceiling(length(inchis)/chunk_size))
for (j in seq_along(props_list)) {
  idx <- ((j-1)*chunk_size + 1) : min(j*chunk_size, length(inchis))
  q <- jsonlite::toJSON(list(inchikey = list(`$in` = unname(inchis[idx]))), auto_unbox = TRUE)
  it <- lotus$iterate(query = q, fields = fields_props)
  rows <- list(); k <- 0L
  repeat {
    doc <- it$one(); if (is.null(doc)) break
    k <- k + 1L
    if (!is.list(doc)) doc <- as.list(doc)
    rows[[k]] <- setNames(
      lapply(PROPS_CORE_FIELDS, function(nm) scalarize_field(doc, nm)),
      PROPS_CORE_FIELDS
    )
  }
  props_list[[j]] <- if (k > 0L) data.table::rbindlist(lapply(rows, as.list), use.names = TRUE, fill = TRUE) else NULL
}
props_core <- data.table::rbindlist(props_list, use.names = TRUE, fill = TRUE)
props_core <- unique(as.data.frame(props_core))
cat("Chemical properties retrieved for", nrow(props_core), "unique InChIKeys.\n")

## enrich and clean

unify_dupes <- function(df, bases = c("lotus_id","smiles","iupac_name",
                                      "molecular_formula","molecular_weight")) {
  for (b in bases) {
    x <- paste0(b, ".x"); y <- paste0(b, ".y")
    has_x <- x %in% names(df); has_y <- y %in% names(df)
    if (has_x && has_y)      { df[[b]] <- dplyr::coalesce(df[[x]], df[[y]]); df[[x]] <- NULL; df[[y]] <- NULL }
    else if (has_x && !has_y){ df[[b]] <- df[[x]]; df[[x]] <- NULL }
    else if (!has_x && has_y){ df[[b]] <- df[[y]]; df[[y]] <- NULL }
  }
  df
}

lin_enriched <- lin %>%
  dplyr::left_join(props_core, by = "inchikey") %>%
  unify_dupes()

NUMERIC_PROPS <- cfg$numeric_props %||% c(
  "molecular_weight","xlogp","alogp","amralogp","manholdlogp",
  "topoPSA","tpsaEfficiency","fsp3",
  "hBondAcceptorCount","hBondDonorCount",
  "number_of_carbons","number_of_oxygens","number_of_nitrogens",
  "total_atom_number","heavy_atom_number",
  "max_number_of_rings","min_number_of_rings",
  "LipinskiRuleOf5Failures"
)
LOGICAL_PROPS <- cfg$logical_props %||% c(
  "contains_ring_sugars","contains_linear_sugars","contains_sugar"
)

lin_enriched <- lin_enriched %>%
  dplyr::mutate(dplyr::across(all_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::mutate(dplyr::across(all_of(LOGICAL_PROPS), ~ {
    if (is.logical(.)) . else tolower(as.character(.)) %in% c("true","t","1")
  }))

## WFO normalization (optional)

USE_WFO <- isTRUE(cfg$use_WFO_normalization %||% TRUE)

if (USE_WFO) {
  cat("WFO: normalizing taxonomy...\n")

  canon_name <- function(x) {
    x <- tolower(stringr::str_squish(as.character(x)))
    x <- gsub("\\s+\\(.*?\\)", "", x)
    x <- gsub("\\b(ex|sensu|auct\\.|non)\\b.*$", "", x)
    x <- gsub("\\b(subsp\\.|ssp\\.|var\\.|subvar\\.|f\\.|forma|cv\\.|group)\\b", "", x)
    x <- stringr::str_squish(x)
    toks <- strsplit(x, "\\s+")
    vapply(toks, function(tt) {
      if (length(tt) >= 3) paste(tt[1:3], collapse = " ")
      else if (length(tt) >= 2) paste(tt[1:2], collapse = " ")
      else tt[1] %||% NA_character_
    }, FUN.VALUE = character(1))
  }

  mk_join_key_vec <- function(genus, species) {
    cand <- ifelse(!is.na(species) & nzchar(species), species, paste(genus, species))
    tolower(stringr::str_squish(cand))
  }

  first_non_na_chr <- function(x) {
    x <- as.character(x); x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }

  split_scientific2 <- function(nm){
    nm <- stringr::str_squish(as.character(nm))
    parts <- strsplit(nm, "\\s+")
    data.frame(
      corrected_genus             = vapply(parts, function(p) if (length(p)>=1) p[1] else NA_character_, ""),
      corrected_specific_epithet  = vapply(parts, function(p) if (length(p)>=2) p[2] else NA_character_, ""),
      corrected_infraspecific     = vapply(parts, function(p) if (length(p)>=3) paste(p[3:length(p)], collapse=" ") else NA_character_, ""),
      stringsAsFactors = FALSE
    )
  }

  WFO_CSV_PATH <- cfg$wfo_csv_path %||%
    file.path(dirname(normalizePath(sys.frame(1)$ofile,
                                    mustWork = FALSE)),
              "..", "DBs", "classification.tsv")
  TODAY_STR    <- format(Sys.Date())
  col_taxonID <- cfg$wfo_cols$taxonID             %||% "taxonID"
  col_sciName <- cfg$wfo_cols$scientificName      %||% "scientificName"
  col_status  <- cfg$wfo_cols$taxonomicStatus     %||% "taxonomicStatus"
  col_accID   <- cfg$wfo_cols$acceptedNameUsageID %||% "acceptedNameUsageID"
  col_family  <- cfg$wfo_cols$family              %||% "family"
  col_genus   <- cfg$wfo_cols$genus               %||% "genus"

  wfo_raw <- readr::read_tsv(
    file      = WFO_CSV_PATH,
    col_types = readr::cols(.default = readr::col_character()),
    progress  = TRUE,
    locale    = readr::locale(encoding = "UTF-8"),
    na = c("", "NA", "NULL")
  )
  if (ncol(wfo_raw) == 1L) stop("WFO file invalid (check delimiters).")

  wfo <- wfo_raw %>%
    transmute(
      taxonID  = .data[[col_taxonID]],
      name     = .data[[col_sciName]],
      status   = tolower(stringr::str_squish(.data[[col_status]])),
      accID    = .data[[col_accID]],
      family   = .data[[col_family]],
      genus    = .data[[col_genus]]
    )

  accepted <- wfo %>%
    filter(status == "accepted") %>%
    transmute(accepted_id = taxonID, accepted_name = name)

  synonyms <- wfo %>%
    filter(status != "accepted", !is.na(accID), nzchar(accID))

  syn_map <- synonyms %>%
    left_join(accepted, by = c("accID" = "accepted_id")) %>%
    transmute(
      synonym_name  = name,
      synonym_id    = taxonID,
      accepted_name = accepted_name,
      accepted_id   = accID
    )

  acc_key <- accepted %>%
    transmute(key = canon_name(accepted_name), accepted_name, accepted_id) %>%
    filter(!is.na(key) & nzchar(key)) %>%
    distinct(key, .keep_all = TRUE)

  syn_key <- syn_map %>%
    transmute(key = canon_name(synonym_name), accepted_name, accepted_id) %>%
    filter(!is.na(key) & nzchar(key)) %>%
    distinct(key, .keep_all = TRUE)

  dict <- bind_rows(acc_key, syn_key) %>%
    distinct(key, .keep_all = TRUE)

  cat("WFO: dictionary loaded:", nrow(dict), "unique keys.\n")

  lin_join <- mk_join_key_vec(lin$genus, lin$species)
  resolved <- data.frame(
    original_name = lin_join,
    key           = canon_name(lin_join),
    stringsAsFactors = FALSE
  ) %>%
    left_join(dict, by = "key") %>%
    mutate(
      tax_provider   = ifelse(!is.na(accepted_name), "WFO_offline", NA_character_),
      tax_status     = dplyr::case_when(
        !is.na(accepted_name) & canon_name(original_name) != canon_name(accepted_name) ~ "synonym",
        !is.na(accepted_name)                                                          ~ "accepted",
        TRUE                                                                           ~ NA_character_
      ),
      tax_checked_at = TODAY_STR
    )

  crosswalk_wfo <- resolved %>%
    transmute(original_name, accepted_name, accepted_id, tax_status,
              tax_provider, tax_checked_at) %>%
    group_by(original_name) %>%
    summarise(
      accepted_name  = first_non_na_chr(accepted_name),
      accepted_id    = first_non_na_chr(accepted_id),
      tax_status     = first_non_na_chr(tax_status),
      tax_provider   = first_non_na_chr(tax_provider),
      tax_checked_at = first_non_na_chr(tax_checked_at),
      .groups = "drop"
    )

  lin_wfo <- lin_enriched %>%
    mutate(.original_name = mk_join_key_vec(genus, species)) %>%
    left_join(crosswalk_wfo, by = c(".original_name" = "original_name"))

  accepted_meta <- wfo %>%
    filter(status == "accepted") %>%
    transmute(
      accepted_id     = taxonID,
      accepted_name_w = name,
      accepted_family = family,
      accepted_genus  = genus
    ) %>%
    distinct(accepted_id, .keep_all = TRUE)

  for (col in c("genus","species","family")) {
    if (!col %in% names(lin_wfo)) lin_wfo[[col]] <- NA_character_
  }
  lin_wfo <- lin_wfo %>% mutate(
    genus   = as.character(genus),
    species = as.character(species),
    family  = as.character(family)
  )

  lin_applied <- lin_wfo %>%
    mutate(
      .corr_scientific = dplyr::coalesce(accepted_name, .original_name),
      .corr_scientific = stringr::str_squish(.corr_scientific)
    ) %>%
    left_join(accepted_meta, by = "accepted_id") %>%
    { cbind(., split_scientific2(.$.corr_scientific)) } %>%
    mutate(
      corrected_family = dplyr::coalesce(accepted_family, family),
      taxonomy_action  = dplyr::case_when(
        is.na(accepted_name) ~ "unresolved",
        canon_name(.original_name) == canon_name(accepted_name) ~ "accepted_confirmed",
        TRUE ~ "synonym_replaced"
      ),
      tax_checked_at   = TODAY_STR
    ) %>%
    mutate(
      species_original = species,
      genus_original   = genus,
      family_original  = family,
      species = dplyr::coalesce(.corr_scientific, species),
      genus   = dplyr::coalesce(corrected_genus, genus),
      family  = dplyr::coalesce(corrected_family, family)
    )

  lin_enriched <- lin_applied
  cat("WFO: normalization applied.\n")

  # if WFO reclassified a species out of the target genus, drop it
  genus_or_lower <- function(x) tolower(tidy_space(as.character(x)))
  if (nrow(lin_enriched) > 0) {
    genus_original <- genus_or_lower(coalesce(lin_enriched$genus_original, lin_enriched$genus))
    genus_final    <- genus_or_lower(coalesce(lin_enriched$accepted_genus,
                                              lin_enriched$corrected_genus,
                                              lin_enriched$genus))
    target_genera_norm <- genus_or_lower(cfg$taxon_values)

    unresolved_flag   <- grepl("unresolved",
                               tolower(coalesce(lin_enriched$tax_status, lin_enriched$taxonomy_action, "")),
                               fixed = TRUE)
    reclassified_flag <- (nzchar(genus_original) & nzchar(genus_final) & (genus_original != genus_final))
    outside_target    <- !(genus_final %in% target_genera_norm)

    remove_vec <- rep(FALSE, nrow(lin_enriched))
    if (identical(tolower(cfg$taxon_mode), "genus")) {
      remove_vec <- (outside_target | reclassified_flag)
    }
    keep_row <- (!unresolved_flag) & (!remove_vec)

    lin_enriched <- lin_enriched[keep_row, , drop = FALSE]
    cat(sprintf("WFO: after genus filter: %d rows kept.\n", nrow(lin_enriched)))
  }
}

## canonicalize

lin_enriched <- lin_enriched %>%
  mutate(
    genus   = canon_genus(genus),
    species = gsub("\\bsp\\.?\\b", "", species, ignore.case = TRUE),
    species = stringr::str_squish(species),
    species = fix_glued_species(genus, species)
  )

## final deduplication

dedup_before <- nrow(lin_enriched)

lin_enriched <- lin_enriched %>%
  dplyr::arrange(inchikey, family, genus, species, ref_id, source) %>%
  dplyr::distinct(
    inchikey, family, genus, species, ref_id,
    .keep_all = TRUE
  )

dedup_after <- nrow(lin_enriched)
cat(sprintf(
  "[DEDUP] Compound x Species x Reference: %d -> %d rows (removed %d duplicates).\n",
  dedup_before, dedup_after, dedup_before - dedup_after
))

## build tables

build_map_tax_inchi <- function(lin_tbl, tax_col = c("family","genus","species")){
  tax_col <- match.arg(tax_col)

  out <- lin_tbl %>%
    dplyr::transmute(
      inchikey = as.character(inchikey),
      family   = as.character(family),
      genus    = canon_genus(genus),
      species  = species
    ) %>%
    dplyr::mutate(
      species = gsub("\\bsp\\.?\\b", "", species, ignore.case = TRUE),
      species = stringr::str_squish(species),
      species = fix_glued_species(genus, species)
    ) %>%
    dplyr::filter(!is.na(inchikey), nzchar(inchikey))

  if (identical(tax_col, "genus")) {
    out <- out %>%
      dplyr::filter(!is.na(genus), nzchar(genus))
  } else if (identical(tax_col, "species")) {
    out <- out %>%
      dplyr::filter(is_binomial(species))
  } else if (identical(tax_col, "family")) {
    out <- out %>%
      dplyr::filter(!is.na(family), nzchar(family))
  }

  out <- out %>%
    dplyr::mutate(
      taxon = dplyr::case_when(
        identical(tax_col, "family")  ~ family,
        identical(tax_col, "genus")   ~ genus,
        identical(tax_col, "species") ~ species,
        TRUE                          ~ species
      )
    ) %>%
    dplyr::filter(!is.na(taxon), nzchar(taxon)) %>%
    dplyr::distinct(inchikey, taxon, family, genus, species)

  if (!tax_col %in% names(out)) out[[tax_col]] <- out$taxon
  out
}

tax_col_for_map <- cfg$analysis_tax_level %||% "genus"
map_tax_inchi <- build_map_tax_inchi(lin_enriched, tax_col = tax_col_for_map)
cat(sprintf("[Map] %s-level map: %d rows.\n", tax_col_for_map, nrow(map_tax_inchi)))

lin_dt2 <- data.table::as.data.table(lin_enriched)
uni <- lin_dt2[, .(
  lotus_id          = safe_first(lotus_id),
  smiles            = safe_first(smiles),
  iupac_name        = safe_first(iupac_name),
  molecular_formula = safe_first(molecular_formula),
  genus             = collapse(genus),
  family            = collapse(family),
  species           = collapse(species),
  ref_ids           = collapse(ref_id)
), by = .(inchikey)] |> as.data.frame()

uni_enriched <- uni %>%
  dplyr::left_join(props_core, by = "inchikey") %>%
  unify_dupes() %>%
  dplyr::mutate(
    dplyr::across(all_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.))),
    dplyr::across(all_of(LOGICAL_PROPS), ~ {
      if (is.logical(.)) . else tolower(as.character(.)) %in% c("true","t","1")
    })
  )

cat("UNI:", nrow(uni), "rows | UNI_ENRICHED:", nrow(uni_enriched), "rows.\n")

dedup_compound_species <- function(df) {
  num_cols  <- names(df)[vapply(df, is.numeric,  logical(1))]
  logi_cols <- names(df)[vapply(df, is.logical, logical(1))]
  chr_cols  <- names(df)[vapply(df, is.character, logical(1))]
  chr_cols  <- setdiff(chr_cols, c("inchikey","family","genus","species","ref_id","source"))

  df %>%
    dplyr::group_by(inchikey, species) %>%
    dplyr::summarise(
      family = collapse(family),
      genus  = collapse(genus),
      ref_id = collapse(ref_id),
      source = collapse(source),
      dplyr::across(dplyr::all_of(num_cols), ~ { v <- .[!is.na(.)]; if (length(v)) v[1] else NA_real_ }, .names = "{.col}"),
      dplyr::across(dplyr::all_of(logi_cols), ~ any(., na.rm = TRUE), .names = "{.col}"),
      dplyr::across(dplyr::all_of(chr_cols), ~ safe_first(.), .names = "{.col}"),
      .groups = "drop"
    )
}

lin_cs <- dedup_compound_species(lin_enriched)
cat(sprintf("[Dedup CS] Compound x Species table: %d rows.\n", nrow(lin_cs)))

## export

EXPORT_EXCEL   <- isTRUE(cfg$export_excel   %||% TRUE)
EXPORT_PARQUET <- isTRUE(cfg$export_parquet %||% TRUE)

if (EXPORT_EXCEL) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    if (isTRUE(cfg$auto_install_missing_packages %||% TRUE)) install.packages("writexl")
  }
  MAX_XLSX <- cfg$max_xlsx_cell_chars %||% 32767L
  BIG_COLS <- cfg$big_cols %||% c(
    "allWikidataIds","xrefs","pubchemBitsString","pubchemBits",
    "circularFingerprint","extendedFingerprint","pfCounts",
    "ertlFunctionalFragmentsPseudoSmiles"
  )
  collapse_list <- function(x){
    sapply(x, function(v){
      if (is.null(v)) return(NA_character_)
      v <- tryCatch(unlist(v, use.names = FALSE), error = function(e) v)
      v <- v[!is.na(v)]
      if (!length(v)) return(NA_character_)
      paste(unique(as.character(v)), collapse = ";")
    })
  }
  trim_cell <- function(x, max_len = MAX_XLSX){
    if (is.null(x)) return(x)
    if (!is.character(x)) x <- as.character(x)
    n <- nchar(x, allowNA = TRUE); too_long <- !is.na(n) & n > max_len
    x[too_long] <- paste0(substr(x[too_long], 1, max_len - 3), "...")
    x
  }
  sanitize_for_excel <- function(df, drop_cols = NULL){
    df <- dplyr::mutate(df, dplyr::across(where(is.factor), as.character))
    df <- dplyr::mutate(df, dplyr::across(where(is.list), collapse_list))
    if (!is.null(drop_cols)) {
      keep <- setdiff(names(df), drop_cols)
      df <- df[keep]
    }
    df <- dplyr::mutate(df, dplyr::across(where(is.character), trim_cell))
    df
  }

  lin_x          <- sanitize_for_excel(lin,           drop_cols = intersect(BIG_COLS, names(lin)))
  lin_enriched_x <- sanitize_for_excel(lin_enriched,  drop_cols = intersect(BIG_COLS, names(lin_enriched)))
  uni_x          <- sanitize_for_excel(uni,           drop_cols = intersect(BIG_COLS, names(uni)))
  uni_enriched_x <- sanitize_for_excel(uni_enriched,  drop_cols = intersect(BIG_COLS, names(uni_enriched)))
  lin_cs_x       <- sanitize_for_excel(lin_cs,        drop_cols = intersect(BIG_COLS, names(lin_cs)))

  xlsx_path <- safe_file(base_tag, ".xlsx")
  writexl::write_xlsx(
    list(
      lin                  = lin_x,
      lin_enriched         = lin_enriched_x,
      lin_compound_species = lin_cs_x,
      uni                  = uni_x,
      uni_enriched         = uni_enriched_x
    ),
    path = xlsx_path
  )
  cat("Excel saved:", normalizePath(xlsx_path), "\n")
}

if (EXPORT_PARQUET) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    if (isTRUE(cfg$auto_install_missing_packages %||% TRUE)) install.packages("arrow")
  }
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(lin_enriched, safe_file(paste0(base_tag, "_lin_enriched"),         ".parquet"))
    arrow::write_parquet(uni_enriched, safe_file(paste0(base_tag, "_uni_enriched"),         ".parquet"))
    arrow::write_parquet(lin_cs,       safe_file(paste0(base_tag, "_lin_compound_species"), ".parquet"))
    cat("Parquet saved:", normalizePath(OUT_DIR), "\n")
  } else {
    warning("Package 'arrow' not available; skipping Parquet export.")
  }
}
