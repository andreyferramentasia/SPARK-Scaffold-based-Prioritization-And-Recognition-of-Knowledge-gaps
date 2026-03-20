## Part III — figures and statistics

## packages
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(rlang)
  library(ggplot2); library(ggrepel); library(scales); library(forcats)
  library(vegan); library(writexl); library(tibble)
})

## helpers
`%||%` <- function(a,b) if (is.null(a)) b else a
to_num  <- function(x) suppressWarnings(as.numeric(x))
nzchar0 <- function(x){ x <- as.character(x); ifelse(is.na(x), FALSE, nzchar(trimws(x))) }

normalize_taxcol <- function(df, tax_col) {
  if (!is.data.frame(df)) return(df)
  if (!tax_col %in% names(df)) stop("normalize_taxcol: missing column '", tax_col,
                                    "'. Names: ", paste(names(df), collapse=", "))
  df %>%
    dplyr::mutate(!!tax_col := as.character(.data[[tax_col]])) %>%
    dplyr::filter(!is.na(.data[[tax_col]]), nzchar(.data[[tax_col]]))
}

safe_expose_species <- function(df, tax_col, prefer = NULL) {
  if (!is.data.frame(df)) return(df)
  if (!"species" %in% names(df)) {
    src <- prefer
    if (is.null(src) || !src %in% names(df)) src <- if (tax_col %in% names(df)) tax_col else NULL
    if (!is.null(src) && src %in% names(df)) df$species <- as.character(df[[src]])
  }
  df
}

ensure_tax_col_present <- function(df, tax_col,
                                   fallbacks = c("species","species.x","species.y",
                                                 "species_std",".species_std","accepted_name",
                                                 "genus","family","taxon")) {
  if (!is.data.frame(df)) return(df)
  if (tax_col %in% names(df)) return(df)
  hit <- intersect(fallbacks, names(df))
  if (length(hit)) {
    df[[tax_col]] <- as.character(df[[hit[1]]])
    return(df)
  }
  stop("ensure_tax_col_present: '", tax_col,
       "' not found and no fallback available. Names: ",
       paste(names(df), collapse=", "))
}


## export flags
SHOW_PLOTS       <- isTRUE(cfg$show_plots %||% TRUE)
SAVE_STATIC_PNGS <- FALSE
SAVE_STATIC_PDFS <- isTRUE(cfg$save_pdfs %||% TRUE)
OPEN_AFTER_SAVE  <- isTRUE(cfg$open_after_save %||% FALSE)

save_plot_multi <- function(plot_obj, file_base, width = 8, height = 12, bg = "white") {
  if (isTRUE(SAVE_STATIC_PDFS)) {
    grDevices::cairo_pdf(paste0(file_base, ".pdf"),
                         width = width, height = height,
                         onefile = FALSE, family = "sans", bg = bg)
    print(plot_obj); grDevices::dev.off()
  }
  if (isTRUE(SHOW_PLOTS)) print(plot_obj)
}

## validate env
need_objs <- c("lin_enriched","uni_enriched","OUT_DIR","base_tag",
               "TAXON_VALUES","TAXON_MODE","cfg")
miss <- need_objs[!vapply(need_objs, exists, logical(1))]
if (length(miss)) stop("Missing required objects in Global Env: ",
                       paste(miss, collapse=", "))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cfg_level_raw  <- cfg$analysis_level %||% cfg$analysis_tax_level
analysis_level <- tolower(cfg_level_raw %||% TAXON_MODE)

if (!analysis_level %in% c("species","genus","family")) {
  analysis_level <- if (identical(TAXON_MODE, "species")) "species" else "genus"
}

# fall back to genus if no usable species column exists
species_candidates_global <- c("species","accepted_name","species_std","scientificName")
if (identical(analysis_level, "species") &&
    !any(species_candidates_global %in% names(lin_enriched))) {
  analysis_level <- "genus"
}

level_internal <- analysis_level
tax_col        <- analysis_level   # "species" | "genus" | "family"

sp_priority <- cfg$species_col_priority %||%
  c("accepted_name","species","species_std","scientificName")
sp_priority <- unique(c(sp_priority,
                        "accepted_name","species","species_std","scientificName"))

lin_scope <- lin_enriched %>%
  dplyr::filter(!is.na(inchikey), nzchar(inchikey))

.bad_pat <- "(^|\\s)(sp\\.?|spp\\.|cf\\.|aff\\.|nr\\.|gr\\.|indet\\.|\\?)($|\\s|-)"
.is_binomial <- function(x) {
  x <- trimws(as.character(x))
  pat_binom <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+$"
  pat_infra <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+\\s+(subsp\\.|ssp\\.|var\\.|f\\.)\\s+[a-z][a-z\\-]+$"
  grepl(pat_binom, x) | grepl(pat_infra, x)
}

species_col <- NULL
for (cand in sp_priority) {
  if (cand %in% names(lin_enriched)) { species_col <- cand; break }
}

if (!is.null(species_col)) {
  if (identical(species_col, "scientificName")) {
    species_std <- vapply(
      strsplit(as.character(lin_enriched$scientificName), "\\s+"),
      function(tok) if (length(tok) >= 2) paste(tok[1], tok[2]) else NA_character_,
      character(1)
    )
  } else {
    species_std <- as.character(lin_enriched[[species_col]])
  }
} else {
  species_std <- rep(NA_character_, nrow(lin_enriched))
}

# Aplica filtro de escopo com base em TAXON_MODE/TAXON_VALUES
lin_scope <- lin_scope %>%
  dplyr::mutate(.species_std = species_std) %>%
  dplyr::filter(
    dplyr::case_when(
      identical(TAXON_MODE, "family")  ~ family %in% TAXON_VALUES,
      identical(TAXON_MODE, "genus")   ~ genus  %in% TAXON_VALUES,
      identical(TAXON_MODE, "species") ~ .species_std %in% TAXON_VALUES,
      TRUE ~ TRUE
    )
  )

face_for_taxon <- function(tax_col) {
  tc <- tolower(as.character(tax_col %||% "species"))
  if (tc %in% c("genus","species")) "italic" else "plain"
}

## --- Optional flags in cfg (ADVANCED CONFIG) ---
# cfg$allow_qualifiers_species    <- FALSE  # TRUE to allow "sp.", "cf.", etc.
# cfg$accept_non_strict_binomial  <- FALSE  # TRUE to allow unmarked trinomials

## --- Species cleaning (only if analyzing by species) ---
if (identical(tax_col, "species")) {
  species_src <- species_col %||% "species"
  lin_scope <- lin_scope %>%
    dplyr::mutate(species_std = trimws(as.character(.data[[species_src]]))) %>%
    dplyr::filter(!is.na(species_std), nzchar(species_std))
  
  if (!isTRUE(cfg$allow_qualifiers_species %||% FALSE)) {
    lin_scope <- lin_scope %>%
      dplyr::filter(!grepl(.bad_pat, species_std, ignore.case = TRUE))
  }
  
  if (!isTRUE(cfg$accept_non_strict_binomial %||% FALSE)) {
    lin_scope <- lin_scope %>% dplyr::filter(.is_binomial(species_std))
  } else {
    .is_binomial_relaxed <- function(x) {
      x <- trimws(as.character(x))
      pat_binom   <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+$"
      pat_infra   <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+\\s+(subsp\\.|ssp\\.|var\\.|f\\.)\\s+[a-z][a-z\\-]+$"
      pat_trinom  <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+\\s+[a-z][a-z\\-]+$"
      grepl(pat_binom, x) | grepl(pat_infra, x) | grepl(pat_trinom, x)
    }
    lin_scope <- lin_scope %>% dplyr::filter(.is_binomial_relaxed(species_std))
  }
  
  lin_scope <- lin_scope %>% dplyr::mutate(species = species_std)
}

# Taxon usado internamente para rótulos/títulos posteriores
lin_scope <- lin_scope %>%
  dplyr::mutate(
    taxon = dplyr::case_when(
      tax_col == "species" ~ as.character(species),
      tax_col == "genus"   ~ as.character(genus),
      tax_col == "family"  ~ as.character(family),
      TRUE                 ~ as.character(species)
    )
  )

if (!nrow(lin_scope)) {
  cat("Diagnostics — top 30 by accepted_name/species:\n")
  if ("accepted_name" %in% names(lin_enriched)) {
    print(
      lin_enriched %>%
        dplyr::count(accepted_name, sort = TRUE) %>%
        head(30)
    )
  }
  if (!is.null(species_col)) {
    print(
      lin_enriched %>%
        dplyr::count(!!rlang::sym(species_col), sort = TRUE) %>%
        head(30)
    )
  }
  stop("Scope empty after applying TAXON_MODE/TAXON_VALUES and species cleaning.")
}

fam_alvo <- lin_scope$family |> unique() |> sort() |> as.character()
fam_alvo <- fam_alvo[nzchar(fam_alvo) & !is.na(fam_alvo)]
if (!length(fam_alvo)) {
  fam_alvo <- unique(na.omit(as.character(lin_enriched$family)))
}
fam_label <- paste(fam_alvo, collapse = ", ")

scope_label <- (function() {
  mode <- TAXON_MODE
  if (identical(mode, "family")) {
    lbl <- fam_label
    if (!nzchar(lbl)) lbl <- "All families in set"
    return(lbl)
  }
  if (identical(mode, "genus")) {
    vals <- TAXON_VALUES
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) vals <- sort(unique(na.omit(as.character(lin_scope$genus))))
    if (!length(vals)) return("Genus set")
    return(paste0("Genus: ", paste(vals, collapse = ", ")))
  }
  if (identical(mode, "species")) {
    vals <- TAXON_VALUES
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) vals <- sort(unique(na.omit(as.character(lin_scope$species))))
    if (!length(vals)) return("Species set")
    if (length(vals) <= 5) {
      return(paste0("Species: ", paste(vals, collapse = ", ")))
    } else {
      return(paste0("Species: ", paste(vals[1:5], collapse = ", "),
                    " + ", length(vals)-5, " more"))
    }
  }
  "Selected set"
})()

## Back-compat: anything still using fam_label_safe will get the scope-aware label
fam_label_safe <- scope_label

min_n_taxon <- as.integer(cfg$analysis_min_compounds_per_taxon %||% 10L)
TOP_N_TAXA  <- as.integer(cfg$analysis_top_taxa %||% 40L)

use_superclass   <- isTRUE(cfg$analysis_heatmap_use_superclass %||%
                             cfg$chem_use_superclass %||% FALSE)
occ_pct_min      <- as.numeric(cfg$analysis_heatmap_min_occ_pct %||%
                                 cfg$chem_class_min_presence_pct %||% 0.00)
n_mols_min       <- as.integer(cfg$analysis_heatmap_min_molecules %||%
                                 cfg$chem_class_min_molecules %||% 5L)
heatmap_max_cols <- as.integer(cfg$analysis_heatmap_max_cols %||%
                                 cfg$chem_heatmap_max_classes %||% 40L)

CAP_TPSA   <- as.numeric(cfg$props_cap_tpsa %||% 500)
CAP_RINGS  <- as.integer(cfg$props_cap_max_rings %||% 10)
USE_MEDIAN <- isTRUE(cfg$props_summary_use_median %||%
                       cfg$taxon_property_use_median %||% FALSE)

DO_RICHNESS_DIVERSITY <- isTRUE(cfg$do_richness_diversity %||% TRUE)
DO_HEATMAP_QUIMICO    <- isTRUE(cfg$do_heatmap_classes %||%
                                  cfg$do_heatmap_chem_classes %||% TRUE)
DO_BOX_PROPS          <- isTRUE(cfg$do_boxplots_props %||%
                                  cfg$do_boxplots_physchem %||% TRUE)
DO_LIPINSKI_SUGAR     <- isTRUE(cfg$do_lipinski_sugars %||%
                                  cfg$do_lipinski_and_sugars %||% TRUE)
DO_ELEM_RATIOS        <- isTRUE(cfg$do_elemental_ratios %||% TRUE)
DO_MURCKO             <- isTRUE(cfg$do_murcko_frameworks %||% TRUE)
DO_SHARED_COMPOUNDS   <- isTRUE(cfg$do_shared_compounds %||% TRUE)
DO_BIBLIOMETRICS      <- isTRUE(cfg$do_bibliometrics %||% TRUE)
DO_PCA_PROPS          <- isTRUE(cfg$do_pca_props %||% TRUE)
DO_NMDS_PERMANOVA     <- isTRUE(cfg$do_nmds_permanova %||% TRUE)
DO_PHYSCHEM_HEATMAP   <- isTRUE(cfg$do_physchem_heatmap %||% TRUE)
DO_PHYSCHEM_VIOLIN    <- isTRUE(cfg$do_physchem_violin %||% TRUE)

## Italic face for taxa on axes
y_face <- if (tax_col %in% c("genus","species")) "italic" else "plain"
x_face <- y_face

make_label_context <- function(cfg, block_name, distance_used = NULL) {
  strat <- cfg$count_strategy_per_block[[block_name]] %||%
    cfg$count_strategy %||% "unique_per_taxon"
  dist_info <- if (!is.null(distance_used))
    paste0(" | distance: ", distance_used) else ""
  label <- dplyr::case_when(
    strat == "unique_per_taxon" ~ "(counts = unique compounds per taxon)",
    strat == "by_occurrence"    ~ "(counts = all occurrences per taxon)",
    TRUE                        ~ "(counts = unspecified)"
  )
  paste0(label, dist_info)
}

## core tables
df_lin <- lin_scope %>%
  dplyr::filter(!is.na(.data[[tax_col]]),
                nzchar(as.character(.data[[tax_col]]))) %>%
  dplyr::distinct(.data[[tax_col]], inchikey)

df_lin <- ensure_tax_col_present(df_lin, tax_col)
df_lin <- safe_expose_species(df_lin, tax_col)

# ---- contar compostos por táxon e limitar pelo pipeline ----
taxa_counts <- df_lin %>%
  dplyr::count(.data[[tax_col]], name = "n") %>%
  dplyr::arrange(dplyr::desc(n))

taxa_keep <- taxa_counts %>%
  dplyr::filter(n >= min_n_taxon) %>%
  dplyr::slice_head(n = TOP_N_TAXA) %>% 
  dplyr::pull(!!rlang::sym(tax_col))

if (!length(taxa_keep)) {
  print(utils::head(taxa_counts, 40))
  stop("No taxon has ≥ ", min_n_taxon,
       " compounds after applying TAXON_MODE/TAXON_VALUES.")
}


# ---- map_tax_inchi: reaproveita, se existir; senão, constrói a partir de lin_scope ----
base_map <- NULL
if (exists("map_tax_inchi") && is.data.frame(map_tax_inchi)) {
  base_map <- map_tax_inchi
} else {
  base_map <- lin_scope %>%
    dplyr::transmute(
      inchikey = as.character(inchikey),
      family   = as.character(family),
      genus    = as.character(genus),
      species  = as.character(species)
    ) %>%
    dplyr::distinct()
}

base_map <- ensure_tax_col_present(base_map, tax_col)
base_map <- safe_expose_species(base_map, tax_col)

map_tax_inchi <- base_map %>%
  dplyr::filter(
    !is.na(inchikey), nzchar0(inchikey),
    !is.na(.data[[tax_col]]), nzchar0(.data[[tax_col]]),
    .data[[tax_col]] %in% taxa_keep
  ) %>%
  dplyr::distinct(inchikey, .data[[tax_col]], family, genus, species,
                  .keep_all = TRUE)

if (!"taxon" %in% names(map_tax_inchi)) {
  map_tax_inchi <- map_tax_inchi %>%
    dplyr::mutate(taxon = .data[[tax_col]])
}

## hybrid class creation

if ("chemicalTaxonomyNPclassifierSuperclass" %in% names(uni_enriched) &&
    "chemicalTaxonomyClassyfireClass" %in% names(uni_enriched)) {

  # vectorized helpers — takes first term before any pipe/semicolon separator
  get_primary_term_v <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "NA" | x == ""] <- NA_character_
    x <- trimws(sub(";.*", "", gsub("\\|", ";", x)))
    x[!nzchar(x) | x == "NA"] <- NA_character_
    x
  }

  standardize_class_v <- function(x) {
    dplyr::case_when(
      is.na(x)                                    ~ NA_character_,
      grepl("Fatty", x, ignore.case = TRUE)       ~ "Fatty acyls",
      x == "Prenol Lipids"                        ~ "Prenol lipids",
      grepl("Glycerophospholipids", x)            ~ "Glycerophospholipids",
      TRUE                                        ~ x
    )
  }

  uni_enriched <- uni_enriched %>%
    dplyr::mutate(
      .np          = get_primary_term_v(chemicalTaxonomyNPclassifierSuperclass),
      .cf          = get_primary_term_v(chemicalTaxonomyClassyfireClass),
      hybrid_class = standardize_class_v(dplyr::coalesce(.np, .cf))
    ) %>%
    dplyr::select(-.np, -.cf)

  CHEM_COL_ALVO <- "hybrid_class"

} else {
  warning("NPClassifier or ClassyFire column missing — falling back to ClassyFire class.")
  CHEM_COL_ALVO <- "chemicalTaxonomyClassyfireClass"
}

## column selection
want_cols <- c(
  "inchikey","molecular_weight","xlogp","topoPSA","fsp3",
  "hBondAcceptorCount","hBondDonorCount",
  "number_of_carbons","number_of_oxygens","number_of_nitrogens",
  "total_atom_number","heavy_atom_number",
  "max_number_of_rings","min_number_of_rings",
  "chemicalTaxonomyClassyfireSuperclass","chemicalTaxonomyClassyfireClass",
  "chemicalTaxonomyNPclassifierSuperclass","chemicalTaxonomyNPclassifierClass",
  "iupac_name","traditional_name",
  "murko_framework","ertlFunctionalFragmentsPseudoSmiles",
  "alogp","amralogp","manholdlogp","tpsaEfficiency"
)

# Dynamically add the target column to ensure it is included in the list
want_cols <- unique(c(want_cols, CHEM_COL_ALVO))

# Now perform intersection (hybrid_class now exists in uni_enriched!)
keep_cols <- intersect(want_cols, names(uni_enriched))

df_props <- uni_enriched %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(
    map_tax_inchi %>%
      dplyr::select(inchikey, !!rlang::sym(tax_col)),
    by = "inchikey",
    relationship = "many-to-many"
  )

df_props <- ensure_tax_col_present(df_props, tax_col)
df_props <- df_props %>%
  dplyr::filter(!is.na(.data[[tax_col]]),
                nzchar(as.character(.data[[tax_col]])))
df_props <- safe_expose_species(df_props, tax_col)

# ---- pa_mat: matriz presença/ausência taxon × composto (flexível) ----
pa_mat <- map_tax_inchi %>%
  dplyr::filter(
    !is.na(.data[[tax_col]]),
    nzchar0(.data[[tax_col]]),
    !is.na(inchikey),
    nzchar0(inchikey)
  ) %>%
  # 1 linha por táxon × inchikey
  dplyr::distinct(!!rlang::sym(tax_col), inchikey) %>%
  dplyr::mutate(val = 1L) %>%
  tidyr::pivot_wider(
    names_from  = inchikey,
    values_from = val,
    values_fill = 0L
  ) %>%
  as.data.frame()

# guarda vetor de táxons
tax_vec <- as.character(pa_mat[[tax_col]])

# remove coluna do nível taxonômico da matriz
pa_mat[[tax_col]] <- NULL

# converte para inteiro/matriz
pa_mat[] <- lapply(pa_mat, function(x) as.integer(x))
pa_mat  <- as.matrix(pa_mat)

# garante nomes de linha únicos (caso extremo)
rownames(pa_mat) <- make.unique(tax_vec)

# ordena táxons por riqueza (nº de compostos)
ord_rich <- order(rowSums(pa_mat > 0, na.rm = TRUE), decreasing = TRUE)
pa_mat   <- pa_mat[ord_rich, , drop = FALSE]

storage.mode(pa_mat) <- "double"
pa_mat[!is.finite(pa_mat)] <- 0


num_cols_all <- c(
  "molecular_weight","xlogp","alogp","amralogp","manholdlogp",
  "topoPSA","tpsaEfficiency","fsp3",
  "hBondAcceptorCount","hBondDonorCount",
  "number_of_carbons","number_of_oxygens","number_of_nitrogens",
  "total_atom_number","heavy_atom_number",
  "max_number_of_rings","min_number_of_rings"
)


## richness & diversity
if (isTRUE(DO_RICHNESS_DIVERSITY) && nrow(pa_mat) > 0 && ncol(pa_mat) > 0) {
  # --- Metrics ---
  richness <- rowSums(pa_mat > 0, na.rm = TRUE)
  div_shan <- vegan::diversity(pa_mat, index = "shannon")
  div_simp <- vegan::diversity(pa_mat, index = "simpson")
  
  stats_df <- tibble::tibble(
    taxon    = names(richness),
    richness = as.integer(richness),
    shannon  = as.numeric(div_shan),
    simpson  = as.numeric(div_simp)
  ) %>%
    dplyr::arrange(dplyr::desc(richness))
  
  # --- Top N selection ---
  n_show <- min(TOP_N_TAXA, nrow(stats_df))
  df_top <- dplyr::slice_head(stats_df, n = n_show)
  
  # --- Subtitles / labels (generalized for family/genus/species) ---
  core_lab <- make_label_context(cfg, "richness_diversity")
  sub_lab  <- paste0(
    core_lab,
    " | taxon level: ",
    tax_col
  )
  
  x_lab <- dplyr::case_when(
    tax_col == "genus"   ~ "Genus",
    tax_col == "species" ~ "Species",
    tax_col == "family"  ~ "Family",
    TRUE                 ~ "Taxon"
  )
  
  # --- Font style for taxon names on y-axis ---
  y_face <- face_for_taxon(tax_col)
  
  # --- Plot ---
  p_rich <- ggplot(df_top, aes(x = reorder(taxon, richness), y = richness)) +
    geom_col(width = 0.72, fill = "grey40") +
    geom_text(aes(label = richness),
              hjust = -0.15, size = 3.4) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = x_lab,
      y = "Richness (unique InChIKeys)",
      title    = sprintf("Top %d by richness — %s", n_show, scope_label),
      subtitle = sub_lab
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = .5, face = "bold"),
      plot.subtitle = element_text(hjust = .5, margin = margin(t = 2, b = 6)),
      axis.title.y  = element_text(margin = margin(r = 6)),
      axis.title.x  = element_text(margin = margin(t = 6)),
      axis.text.y   = element_text(face = y_face),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin = margin(t = 6, r = 18, b = 6, l = 6)
    )
  
  # --- Save figure (PDF only; pipeline-integrated) ---
  save_plot_multi(
    p_rich,
    file.path(OUT_DIR, paste0(base_tag, "_richness_top"))
  )
  
  # --- Export summary table ---
  out_stats <- file.path(OUT_DIR, paste0(base_tag, "_richness_diversity.xlsx"))
  writexl::write_xlsx(list(richness_diversity = stats_df), path = out_stats)
  cat("✔ Saved:", normalizePath(out_stats, winslash = "/"), "\n")
}

## chemical class heatmap

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(ComplexHeatmap); library(circlize); library(grid)
  library(viridisLite)
})

make_class_taxon_heatmap <- function(
    df_props,
    tax_col,
    use_superclass = FALSE,    # TRUE → Superclass; FALSE → Class
    log_transform  = TRUE,     # log1p nos counts
    drop_zeros     = TRUE,     # remove linhas/colunas com soma 0
    top_n_rows     = NULL,     # mantém top-N classes por total
    top_n_cols     = NULL,     # mantém top-N taxa por total
    occ_pct_min    = 0,        # fração mínima de taxa em que a classe aparece
    n_mols_min     = 0,        # nº mínimo de compostos por classe (soma total)
    scope_label    = NULL,     # rótulo de escopo para título (ex.: família/gênero alvo)
    outdir         = getwd(),
    file_prefix    = "class_taxon_heatmap",
    width_in       = 13,
    height_in      = 7.5
){
  if (!is.data.frame(df_props)) stop("df_props must be a data.frame")
  if (!nzchar(tax_col) || !tax_col %in% names(df_props)) {
    stop("tax_col '", tax_col, "' not found in df_props. Names: ",
         paste(names(df_props), collapse = ", "))
  }
  
  axis_col <- if (isTRUE(use_superclass))
    "chemicalTaxonomyClassyfireSuperclass" else "chemicalTaxonomyClassyfireClass"
  if (!axis_col %in% names(df_props)) {
    stop("Column '", axis_col, "' not found in df_props.")
  }
  if (!"inchikey" %in% names(df_props)) {
    stop("Column 'inchikey' not found in df_props.")
  }
  
  # Se scope_label não for passado, tenta pegar do ambiente global (Part III)
  if (is.null(scope_label) && exists("scope_label", envir = .GlobalEnv)) {
    scope_label <- get("scope_label", envir = .GlobalEnv)
  }
  
  # ---- 1) Tabela CLASS × TAXON (n_distinct(inchikey)) ----
  df_hm <- df_props %>%
    filter(
      !is.na(.data[[axis_col]]) & nzchar(trimws(as.character(.data[[axis_col]]))),
      !is.na(.data[[tax_col]])  & nzchar(trimws(as.character(.data[[tax_col]])))
    ) %>%
    mutate(
      !!axis_col := trimws(as.character(.data[[axis_col]])),
      !!tax_col  := trimws(as.character(.data[[tax_col]]))
    ) %>%
    group_by(!!rlang::sym(axis_col), !!rlang::sym(tax_col)) %>%
    summarise(n = dplyr::n_distinct(inchikey), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = !!rlang::sym(tax_col),
      values_from = n,
      values_fill = 0
    )
  
  if (nrow(df_hm) == 0) stop("No data after class × taxon grouping.")
  
  mat <- as.data.frame(df_hm)
  rownames(mat) <- mat[[axis_col]]
  mat[[axis_col]] <- NULL
  mat <- as.matrix(mat)
  
  # ---- 2) Filtros básicos: zeros, thresholds e top-N ----
  # Remove classes/colunas com tudo zero (se solicitado)
  if (drop_zeros) {
    if (nrow(mat) > 1) mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    if (ncol(mat) > 1) mat <- mat[, colSums(mat) > 0, drop = FALSE]
  }
  
  # Filtro por nº mínimo de moléculas por classe (soma total na linha)
  if (!is.null(n_mols_min) && n_mols_min > 0 && nrow(mat) > 0) {
    keep_rows <- rowSums(mat, na.rm = TRUE) >= n_mols_min
    mat <- mat[keep_rows, , drop = FALSE]
  }
  
  # Filtro por % mínima de taxa em que a classe aparece (linha)
  if (!is.null(occ_pct_min) && occ_pct_min > 0 && nrow(mat) > 0 && ncol(mat) > 0) {
    frac_taxa <- rowSums(mat > 0, na.rm = TRUE) / ncol(mat)
    keep_rows <- frac_taxa >= occ_pct_min
    mat <- mat[keep_rows, , drop = FALSE]
  }
  
  # Top-N classes (linhas)
  if (!is.null(top_n_rows) && nrow(mat) > top_n_rows) {
    idx <- order(rowSums(mat), decreasing = TRUE)[seq_len(top_n_rows)]
    mat <- mat[idx, , drop = FALSE]
  }
  # Top-N taxa (colunas)
  if (!is.null(top_n_cols) && ncol(mat) > top_n_cols) {
    idx <- order(colSums(mat), decreasing = TRUE)[seq_len(top_n_cols)]
    mat <- mat[, idx, drop = FALSE]
  }
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Matrix is empty after filters (zeros / thresholds / top-N).")
  }
  
  # ---- 3) Transformação (log1p por default) ----
  mat_plot <- if (isTRUE(log_transform)) log1p(mat) else mat
  
  # ---- 4) Paleta de cores (30 tons, sequencial) ----
  brks <- seq(min(mat_plot, na.rm = TRUE),
              max(mat_plot, na.rm = TRUE), length.out = 30)
  cols <- colorRampPalette(
    c("#FEFFF7","#FFF7BC","#FEE391","#FEC44F",
      "#FE9929","#EC7014","#CC4C02","#993404","#662506")
  )(30)
  col_fun <- circlize::colorRamp2(breaks = brks, colors = cols)
  
  # ---- 5) Títulos / rótulos ----
  main_title_core <- paste0(
    if (use_superclass) "Superclass" else "Class",
    " × ", tools::toTitleCase(tax_col),
    " (distinct compounds)"
  )
  if (!is.null(scope_label) && nzchar(scope_label)) {
    main_title <- paste0(main_title_core, " — ", scope_label)
  } else {
    main_title <- main_title_core
  }
  
  row_title  <- if (use_superclass) "ClassyFire Superclass" else "ClassyFire Class"
  legend_nm  <- if (log_transform) "log1p(n)" else "n"
  
  # Fonte para nomes de coluna (taxa): itálico só para genus/species
  col_face <- if (tolower(tax_col) %in% c("genus","species")) "italic" else "plain"
  
  # ---- 6) Heatmap (clustering average) ----
  ht <- Heatmap(
    mat_plot,
    name = legend_nm,
    col  = col_fun,
    cluster_rows    = (nrow(mat_plot) > 1),
    cluster_columns = (ncol(mat_plot) > 1),
    clustering_method_rows    = "average",
    clustering_method_columns = "average",
    show_row_dend    = (nrow(mat_plot) > 1),
    show_column_dend = (ncol(mat_plot) > 1),
    row_title        = row_title,
    column_title     = main_title,
    row_names_gp     = grid::gpar(fontsize = 9),
    column_names_gp  = grid::gpar(fontsize = 9, fontface = col_face),
    column_names_rot = 55,
    column_names_centered = FALSE
  )
  
  # ---- 7) Desenhar no painel de Plots ----
  grid.newpage()
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  
  # ---- 8) Export PDF ----
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
  pdf_path <- file.path(outdir, paste0(file_prefix, "_", stamp, ".pdf"))
  grDevices::pdf(pdf_path, width = width_in, height = height_in, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()
  
  message("Saved: ", pdf_path)
  invisible(list(pdf = pdf_path, matrix = mat, breaks = brks))
}

## chemical class heatmap — call

OUT_DIR_SAFE  <- if (exists("OUT_DIR")) OUT_DIR else getwd()
BASE_TAG_SAFE <- if (exists("base_tag")) base_tag else "lotus_chem"

if (isTRUE(DO_HEATMAP_QUIMICO)) {
  if (!exists("df_props") || !is.data.frame(df_props) || !nrow(df_props)) {
    message("DO_HEATMAP_QUIMICO = TRUE, but df_props is empty or missing — skipping heatmap.")
  } else {
    message("Running chemical class heatmap (Hybrid NP/ClassyFire × ", tax_col, ")...")
    
    df_plot_input <- df_props
    if (exists("CHEM_COL_ALVO") && CHEM_COL_ALVO %in% names(df_plot_input)) {
      df_plot_input$chemicalTaxonomyClassyfireClass <- df_plot_input[[CHEM_COL_ALVO]]
    }
    
    try({
      hm_res <- make_class_taxon_heatmap(
        df_props       = df_plot_input,
        tax_col        = tax_col,
        use_superclass = FALSE,
        log_transform  = TRUE,
        drop_zeros     = TRUE,
        top_n_rows     = NULL,
        top_n_cols     = heatmap_max_cols,
        occ_pct_min    = occ_pct_min,
        n_mols_min     = n_mols_min,
        scope_label    = scope_label,
        outdir         = OUT_DIR_SAFE,
        file_prefix    = paste0(BASE_TAG_SAFE, "_chem_heatmap_Hybrid"), # Nome atualizado
        width_in       = 13,
        height_in      = 7.5
      )
    }, silent = FALSE)
  }
} 

## distribution plots

if (exists("df_props") && exists("taxa_keep")) {
  df_visuals_ready <- df_props %>%
    dplyr::filter(.data[[tax_col]] %in% taxa_keep) %>%
    dplyr::filter(nzchar(as.character(.data[[tax_col]])), !is.na(.data[[tax_col]])) %>%
    dplyr::mutate(
      C_num = to_num(number_of_carbons),
      O_num = to_num(number_of_oxygens),
      OC_Ratio = dplyr::if_else(C_num > 0, O_num / C_num, 0),
      topoPSA_cap = pmin(to_num(topoPSA), CAP_TPSA),
      n_rings_cap = pmin(to_num(max_number_of_rings), CAP_RINGS)
    )
} else {
  stop("Critical Error: 'df_props' or 'taxa_keep' missing. Run Section 4 first.")
}

plot_box <- function(df_source, var, label = NULL, use_cap = FALSE) {
  v <- if (use_cap) paste0(var, "_cap") else var
  if (!v %in% names(df_source)) { return(invisible(NULL)) }
  
  d <- df_source %>%
    dplyr::filter(!is.na(.data[[v]])) %>%
    dplyr::group_by(.data[[tax_col]]) %>%
    dplyr::filter(dplyr::n() >= min_n_taxon) %>%
    dplyr::ungroup()
  
  if (!nrow(d)) return(invisible(NULL))
  
  # Labels
  lbl_x <- tools::toTitleCase(tax_col)
  lbl_y <- label %||% var
  ttl   <- sprintf("%s by %s – %s", lbl_y, lbl_x, scope_label)
  
  # Sort by MEAN
  ord_fun  <- mean
  stat_fun <- mean
  y_face   <- face_for_taxon(tax_col)
  
  p <- ggplot(d, aes(x = reorder(.data[[tax_col]], .data[[v]], FUN = ord_fun), y = .data[[v]])) +
    geom_boxplot(fill = "white", color = "black", outlier.alpha = 0.2, outlier.size = 1) +
    stat_summary(fun = stat_fun, geom = "point", size = 1.8, color = "red") + 
    coord_flip() +
    labs(x = lbl_x, y = lbl_y, title = ttl, subtitle = "(Red dot = Mean)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(face = y_face, size = 10, hjust = 0),
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank()
    )
  
  fn <- file.path(OUT_DIR, paste0(base_tag, "_box_", var))
  save_plot_multi(p, fn, width = 8, height = max(6, length(unique(d[[tax_col]]))*0.2))
}

if (isTRUE(DO_BOX_PROPS)) {
  plot_box(df_visuals_ready, "molecular_weight", "Molecular Weight")
  plot_box(df_visuals_ready, "xlogp", "xLogP")
  plot_box(df_visuals_ready, "topoPSA", "Topological PSA", use_cap = TRUE)
  plot_box(df_visuals_ready, "fsp3", "Fsp3 (Complexity)")
  plot_box(df_visuals_ready, "OC_Ratio", "O/C Ratio (Oxidation)")
}


## violin plots

if (isTRUE(DO_PHYSCHEM_VIOLIN)) {
  v_list <- c("molecular_weight", "xlogp", "topoPSA", "fsp3", "OC_Ratio", "hBondDonorCount")
  v_labs <- c("MW", "xLogP", "TPSA", "Fsp3", "O/C Ratio", "HBD")
  
  # Filter Data
  df_vio_long <- df_visuals_ready %>%
    dplyr::select(dplyr::all_of(c(tax_col, "inchikey", v_list))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(v_list), names_to = "Var", values_to = "Val") %>%
    dplyr::filter(!is.na(Val))
  
  keep_vio <- df_vio_long %>%
    dplyr::count(.data[[tax_col]], Var) %>%
    dplyr::filter(n >= min_n_taxon) %>%
    dplyr::group_by(.data[[tax_col]]) %>%
    dplyr::filter(dplyr::n_distinct(Var) == length(v_list)) %>% 
    dplyr::pull(!!rlang::sym(tax_col)) %>%
    unique()
  
  df_vio_final <- df_vio_long %>% dplyr::filter(.data[[tax_col]] %in% keep_vio)
  
  if (nrow(df_vio_final) > 0) {
    # Fix Labels
    df_vio_final$Var <- factor(df_vio_final$Var, levels = v_list, labels = v_labs)
    
    # Sorting: Sort taxa by MEAN MW
    ord_levels <- df_visuals_ready %>%
      dplyr::group_by(.data[[tax_col]]) %>%
      dplyr::summarise(mean_mw = mean(to_num(molecular_weight), na.rm=TRUE)) %>%
      dplyr::arrange(mean_mw) %>%
      dplyr::pull(!!rlang::sym(tax_col))
    
    df_vio_final[[tax_col]] <- factor(df_vio_final[[tax_col]], levels = ord_levels)
    
    # Plot
    p_vio <- ggplot(df_vio_final, aes(x = .data[[tax_col]], y = Val)) +
      geom_violin(aes(fill = Var), scale = "width", trim = TRUE, alpha = 0.6, size = 0.3) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.8) +
      stat_summary(fun = mean, geom = "point", size = 0.8, color = "black") +
      facet_wrap(~Var, scales = "free_x", nrow = 1) + 
      coord_flip() +
      scale_fill_viridis_d(option = "D", guide = "none") +
      labs(
        title = paste0("Physicochemical Distributions (Inc. O/C) — ", scope_label),
        x = tools::toTitleCase(tax_col), y = "Value"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(face = face_for_taxon(tax_col), size = 9),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    save_plot_multi(p_vio, file.path(OUT_DIR, paste0(base_tag, "_physchem_violins_new")), width = 16, height = 10)
  } else {
    message("Insufficient data for Violin plots.")
  }
}

## physicochemical heatmap
if (isTRUE(DO_PHYSCHEM_HEATMAP)) {
  needed <- c("df_props", "tax_col", "OUT_DIR", "base_tag")
  miss   <- needed[!vapply(needed, exists, logical(1))]
  if (length(miss)) stop("[Physchem heatmap] Missing: ", paste(miss, collapse = ", "))
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
    library(ComplexHeatmap); library(circlize); library(grid)
  })
  
  # ---- Candidate columns ----
  phys_candidates <- c(
    "molecular_weight","xlogp","topoPSA","fsp3",
    "hBondAcceptorCount","hBondDonorCount",
    "min_number_of_rings","max_number_of_rings","number_of_rings"
  )
  phys_cols0 <- intersect(phys_candidates, names(df_props))
  if (!length(phys_cols0)) {
    message("[Physchem heatmap] No numeric columns found; skipped.")
  } else {
    df_phys0 <- df_props
    if (!"number_of_rings" %in% names(df_phys0)) {
      if ("min_number_of_rings" %in% names(df_phys0))
        df_phys0$number_of_rings <- df_phys0$min_number_of_rings
      else if ("max_number_of_rings" %in% names(df_phys0))
        df_phys0$number_of_rings <- df_phys0$max_number_of_rings
    }
    
    phys_cols <- intersect(
      c("molecular_weight","xlogp","topoPSA","fsp3",
        "hBondAcceptorCount","hBondDonorCount","number_of_rings"),
      names(df_phys0)
    )
    if (!length(phys_cols)) {
      message("[Physchem heatmap] No valid physicochemical columns; skipped.")
    } else {
      # ---- Long table ----
      df_phys_long <- df_phys0 %>%
        select(inchikey, !!rlang::sym(tax_col), all_of(phys_cols)) %>%
        mutate(across(all_of(phys_cols), ~ suppressWarnings(as.numeric(.)))) %>%
        pivot_longer(cols = all_of(phys_cols),
                     names_to = "property", values_to = "value") %>%
        filter(is.finite(value), !is.na(.data[[tax_col]]), nzchar0(.data[[tax_col]]))
      
      if (!nrow(df_phys_long)) {
        message("[Physchem heatmap] No numeric values; skipped.")
      } else {
        # ---- Select top taxa ----
        tmp <- df_phys_long %>%
          distinct(!!rlang::sym(tax_col), inchikey) %>%
          count(!!rlang::sym(tax_col), name = "n") %>%
          arrange(desc(n))
        n_show <- min(TOP_N_TAXA, nrow(tmp))
        top_taxa <- tmp %>% slice_head(n = n_show) %>% pull(!!rlang::sym(tax_col))
        
        # ---- Median + z-score ----
        df_heat <- df_phys_long %>%
          filter(.data[[tax_col]] %in% top_taxa) %>%
          group_by(!!rlang::sym(tax_col), property) %>%
          summarise(median_value = median(value, na.rm = TRUE), .groups = "drop") %>%
          group_by(property) %>%
          mutate(median_scaled = {
            m  <- mean(median_value, na.rm = TRUE)
            sdv <- sd(median_value, na.rm = TRUE)
            if (is.finite(sdv) && sdv > 0) (median_value - m)/sdv else 0
          }) %>%
          ungroup()
        
        # ---- Build matrix taxon × property ----
        mat_phys <- df_heat %>%
          select(!!rlang::sym(tax_col), property, median_scaled) %>%
          pivot_wider(names_from = property, values_from = median_scaled) %>%
          tibble::column_to_rownames(var = tax_col) %>%
          as.matrix()
        mat_phys[is.na(mat_phys)] <- 0
        
        # ---- Colors ----
        brks <- seq(min(mat_phys, na.rm = TRUE),
                    max(mat_phys, na.rm = TRUE), length.out = 30)
        col_fun <- circlize::colorRamp2(
          breaks = brks,
          colors = colorRampPalette(c("#2166AC","white","#B2182B"))(30)
        )
        
        # ---- Titles ----
        main_title <- paste0(
          "Physicochemical profile by ", tax_col,
          " — ", scope_label
        )
        row_title <- tools::toTitleCase(tax_col)
        legend_nm <- "Scaled median"
        
        # ---- Heatmap ----
        ht_phys <- Heatmap(
          mat_phys,
          name = legend_nm,
          col  = col_fun,
          cluster_rows    = TRUE,
          cluster_columns = TRUE,
          clustering_method_rows    = "average",
          clustering_method_columns = "average",
          show_row_dend    = TRUE,
          show_column_dend = TRUE,
          row_names_gp     = grid::gpar(fontsize = 9,
                                        fontface = if (tax_col %in% c("genus","species")) "italic" else "plain"),
          column_names_gp  = grid::gpar(fontsize = 9),
          column_names_rot = 35,
          column_names_centered = FALSE,
          row_title        = row_title,
          column_title     = main_title,
          heatmap_legend_param = list(title_position = "topcenter")
        )
        
        # ---- Draw ----
        grid.newpage()
        draw(ht_phys, heatmap_legend_side = "right",
             annotation_legend_side = "right")
        
        # ---- Export PDF ----
        pdf_path <- file.path(
          OUT_DIR,
          paste0(base_tag, "_physchem_heatmap_Complex.pdf")
        )
        grDevices::pdf(pdf_path, width = 8.5, height = 6.0, useDingbats = FALSE)
        draw(ht_phys, heatmap_legend_side = "right",
             annotation_legend_side = "right")
        grDevices::dev.off()
        message("Saved: ", pdf_path)
      }
    }
  }
}




## lipinski & sugars

if (isTRUE(DO_LIPINSKI_SUGAR)) {
  as_bool <- function(x) {
    if (is.null(x)) return(logical(0))
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) {
      xl <- tolower(trimws(x))
      return(xl %in% c("true","t","1","yes","y","sim","s"))
    }
    rep(NA, length(x))
  }
  
  derive_any_sugar <- function(base_df) {
    # 1) se já existir any_sugar razoavelmente preenchido, usa
    if ("any_sugar" %in% names(base_df) && any(!is.na(base_df$any_sugar))) {
      base_df$any_sugar <- as_bool(base_df$any_sugar)
      return(base_df)
    }
    # 2) se existir contains_sugar, usa como base
    if ("contains_sugar" %in% names(base_df) && any(!is.na(base_df$contains_sugar))) {
      base_df$any_sugar <- as_bool(base_df$contains_sugar)
      if (any(!is.na(base_df$any_sugar))) return(base_df)
    }
    # 3) fallback: heurística no nome
    n <- nrow(base_df)
    any_sugar_name <- rep(FALSE, n)
    nm <- rep("", n)
    for (cand in c("iupac_name","name","traditional_name")) {
      if (cand %in% names(base_df)) {
        nm <- paste(nm, base_df[[cand]])
      }
    }
    nm <- tolower(nm)
    kw <- c(
      "glycoside","glucoside","rhamnoside","arabinoside","galactoside",
      "mannoside","xyloside","rutinoside","neohesperidoside","sophoroside",
      "diglycoside","triglycoside","dirhamnoside","rutinosyl","glucosyl"
    )
    if (n > 0 && any(nzchar(nm))) {
      any_sugar_name <- grepl(paste(kw, collapse = "|"), nm)
    }
    base_df$any_sugar <- any_sugar_name
    base_df
  }
  
  base0 <- df_props %>%
    dplyr::inner_join(
      dplyr::select(map_tax_inchi, inchikey, !!rlang::sym(tax_col)),
      by = "inchikey",
      relationship = "many-to-many"
    ) %>%
    { ensure_tax_col_present(., tax_col) } %>%
    normalize_taxcol(tax_col = tax_col) %>%
    dplyr::mutate(
      molecular_weight = to_num(molecular_weight),
      xlogp            = to_num(xlogp),
      HBD = to_num(if ("hBondDonorCount"    %in% names(.)) hBondDonorCount    else NA),
      HBA = to_num(if ("hBondAcceptorCount" %in% names(.)) hBondAcceptorCount else NA)
    ) %>%
    dplyr::mutate(
      LipinskiRuleOf5Failures = dplyr::case_when(
        is.na(molecular_weight) | is.na(xlogp) | is.na(HBD) | is.na(HBA) ~ NA_real_,
        TRUE ~ (as.numeric(molecular_weight > 500) +
                  as.numeric(xlogp            > 5) +
                  as.numeric(HBD              > 5) +
                  as.numeric(HBA              > 10))
      ),
      Lipinski = dplyr::case_when(
        is.na(LipinskiRuleOf5Failures) ~ NA_character_,
        LipinskiRuleOf5Failures > 0    ~ "Violates",
        TRUE                           ~ "OK"
      )
    ) %>%
    derive_any_sugar() %>%
    { safe_expose_species(., tax_col) }
  
  base0_filt <- base0 %>%
    dplyr::filter(!is.na(.data[[tax_col]]), nzchar(.data[[tax_col]]))
  
  tab <- base0_filt %>%
    dplyr::group_by(.data[[tax_col]]) %>%
    dplyr::summarise(
      n                   = dplyr::n(),
      n_Lipinski_measured = sum(!is.na(Lipinski)),
      pct_Lipinski_OK     = mean(Lipinski == "OK",    na.rm = TRUE),
      pct_with_sugar      = mean(any_sugar %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n))
  
  print(tab)
  
  lip_long <- base0_filt %>%
    dplyr::filter(!is.na(Lipinski)) %>%
    dplyr::count(!!rlang::sym(tax_col), Lipinski, name = "n") %>%
    tidyr::complete(
      !!rlang::sym(tax_col),
      Lipinski = c("OK","Violates"),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(!!rlang::sym(tax_col)) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  if (nrow(lip_long)) {
    taxa_ord <- tab %>% dplyr::pull(!!rlang::sym(tax_col))
    sub_lab  <- make_label_context(cfg, "lipinski")
    
    tax_label <- dplyr::case_when(
      tax_col == "genus"   ~ "genus",
      tax_col == "species" ~ "species",
      tax_col == "family"  ~ "family",
      TRUE                 ~ "taxon"
    )
    x_lab <- dplyr::case_when(
      tax_col == "genus"   ~ "Genus",
      tax_col == "species" ~ "Species",
      tax_col == "family"  ~ "Family",
      TRUE                 ~ "Taxon"
    )
    
    p_lip <- ggplot(
      lip_long %>%
        dplyr::mutate(taxon = factor(.data[[tax_col]], levels = taxa_ord)),
      aes(x = taxon, y = prop, fill = Lipinski)
    ) +
      geom_col(position = position_fill(reverse = TRUE),
               width = 0.8, color = NA) +
      coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_fill_manual(
        values = c("OK" = "#1E88E5", "Violates" = "#9E9E9E"),
        name   = NULL
      ) +
      labs(
        x = x_lab,
        y = "Proportion",
        title    = paste0("Lipinski Rule by ", tax_label, " — ", scope_label),
        subtitle = sub_lab
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title        = element_text(hjust = .5, face = "bold"),
        axis.text.y       = element_text(face = face_for_taxon(tax_col)),
        panel.grid.major.y= element_blank(),
        legend.position   = "top"
      )
    
    save_plot_multi(
      p_lip,
      file.path(OUT_DIR, paste0(base_tag, "_lipinski_rule"))
    )
  }
  
  sugar_long <- base0_filt %>%
    dplyr::mutate(
      sugar_lbl = dplyr::case_when(
        any_sugar %in% TRUE  ~ "contains sugar",
        any_sugar %in% FALSE ~ "no sugar",
        TRUE                 ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(sugar_lbl)) %>%
    dplyr::count(.data[[tax_col]], sugar_lbl, name = "n") %>%
    tidyr::complete(
      !!rlang::sym(tax_col),
      sugar_lbl = c("no sugar","contains sugar"),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(!!rlang::sym(tax_col)) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  if (nrow(sugar_long) > 0) {
    taxa_ord2 <- sugar_long %>%
      dplyr::count(.data[[tax_col]], wt = n, name = "ntot") %>%
      dplyr::arrange(dplyr::desc(ntot)) %>%
      dplyr::pull(!!rlang::sym(tax_col))
    
    sub_lab  <- make_label_context(cfg, "sugars")
    
    tax_label <- dplyr::case_when(
      tax_col == "genus"   ~ "genus",
      tax_col == "species" ~ "species",
      tax_col == "family"  ~ "family",
      TRUE                 ~ "taxon"
    )
    x_lab <- dplyr::case_when(
      tax_col == "genus"   ~ "Genus",
      tax_col == "species" ~ "Species",
      tax_col == "family"  ~ "Family",
      TRUE                 ~ "Taxon"
    )
    
    p_sugar <- ggplot(
      sugar_long %>%
        dplyr::mutate(taxon = factor(.data[[tax_col]], levels = taxa_ord2)),
      aes(x = taxon, y = prop, fill = sugar_lbl)
    ) +
      geom_col(position = "fill", width = 0.8, color = NA) +
      coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_fill_manual(
        values = c("no sugar" = "grey70", "contains sugar" = "steelblue"),
        name   = NULL
      ) +
      labs(
        x = x_lab,
        y = "Proportion",
        title    = paste0("Glycosides by ", tax_label, " — ", scope_label),
        subtitle = sub_lab
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title        = element_text(hjust = .5, face = "bold"),
        axis.text.y       = element_text(face = face_for_taxon(tax_col)),
        panel.grid.major.y= element_blank(),
        legend.position   = "top"
      )
    
    save_plot_multi(
      p_sugar,
      file.path(OUT_DIR, paste0(base_tag, "_glycosides"))
    )
  }
  
  out_lip <- file.path(OUT_DIR, paste0(base_tag, "_lipinski_sugars_summary.xlsx"))
  writexl::write_xlsx(
    list(summary_by_taxon = tab),
    path = out_lip
  )
  cat("✔ Saved:", normalizePath(out_lip, winslash = "/"), "\n")
}

## elemental ratios
if (isTRUE(DO_ELEM_RATIOS)) {
  if (exists("CHEM_COL_ALVO") && CHEM_COL_ALVO %in% names(uni_enriched)) {
    eixo_tax <- CHEM_COL_ALVO
    eixo_label <- "Chemical Class (Hybrid: NPClassifier > ClassyFire)"
  } else {
    eixo_tax <- "chemicalTaxonomyClassyfireClass"
    eixo_label <- "Class (ClassyFire)"
  }
  
  elem <- uni_enriched %>%
    dplyr::select(inchikey, number_of_carbons, number_of_oxygens) %>%
    dplyr::mutate(
      C = to_num(number_of_carbons),
      O = to_num(number_of_oxygens),
      C = dplyr::if_else(is.na(C) | C <= 0, NA_real_, C),
      OC = O / C,
      OC_cap = pmin(OC, 1.20)
    ) %>%
    dplyr::filter(is.finite(OC_cap)) %>%
    dplyr::inner_join(
      dplyr::select(map_tax_inchi, inchikey, !!rlang::sym(tax_col)),
      by = "inchikey",
      relationship = "many-to-many"
    )
  
  if (nrow(elem)) {
    d_tax <- elem %>%
      dplyr::group_by(.data[[tax_col]]) %>%
      dplyr::filter(dplyr::n() >= min_n_taxon) %>%
      dplyr::ungroup()
    
    if (nrow(d_tax)) {
      sub_lab <- make_label_context(cfg, "oc_taxon")
      y_face  <- face_for_taxon(tax_col)
      
      p_oc_taxon <- ggplot(d_tax, aes(x = reorder(.data[[tax_col]], OC_cap, FUN = median, na.rm = TRUE), y = OC_cap)) +
        geom_boxplot(outlier.alpha = .25, width = .6, color = "black", fill = "white", linewidth = .5) +
        stat_summary(fun = median, geom = "point", size = 1.8, color = "black") +
        geom_jitter(width = .12, alpha = .25, size = .9, color = "grey30") +
        coord_flip() +
        labs(x = tools::toTitleCase(tax_col), y = "O/C (cap=1.20)", 
             title = paste0("O/C by ", tax_col, " — ", scope_label), subtitle = sub_lab) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = .5, face = "bold"), axis.text.y = element_text(face = y_face))
      
      save_plot_multi(p_oc_taxon, file.path(OUT_DIR, paste0(base_tag, "_OC_taxon")))
    }
  }
  
  if (eixo_tax %in% names(uni_enriched)) {
    
    chem_elem <- uni_enriched %>%
      dplyr::filter(!is.na(.data[[eixo_tax]]), nzchar(as.character(.data[[eixo_tax]]))) %>%
      dplyr::select(
        inchikey, !!rlang::sym(eixo_tax),
        number_of_carbons, number_of_oxygens
      ) %>%
      dplyr::mutate(
        C = to_num(number_of_carbons),
        O = to_num(number_of_oxygens),
        C = dplyr::if_else(is.na(C) | C <= 0, NA_real_, C),
        OC = O / C
      ) %>%
      dplyr::filter(is.finite(OC)) %>%
      dplyr::inner_join(
        dplyr::select(map_tax_inchi, inchikey, !!rlang::sym(tax_col)),
        by = "inchikey",
        relationship = "many-to-many"
      )
    
    if (nrow(chem_elem)) {
      present_tab <- chem_elem %>%
        dplyr::distinct(.data[[tax_col]], .data[[eixo_tax]]) %>%
        dplyr::count(.data[[eixo_tax]], name = "n_tax_present")
      
      n_taxa_all <- chem_elem %>% dplyr::distinct(.data[[tax_col]]) %>% nrow()
      
      keep_occ  <- present_tab %>%
        dplyr::filter(n_tax_present / n_taxa_all >= occ_pct_min) %>%
        dplyr::pull(!!rlang::sym(eixo_tax))
      
      keep_mols <- chem_elem %>%
        dplyr::distinct(inchikey, .data[[eixo_tax]]) %>%
        dplyr::count(.data[[eixo_tax]], name = "n_mols") %>%
        dplyr::filter(n_mols >= n_mols_min) %>%
        dplyr::pull(!!rlang::sym(eixo_tax))
      
      chem_elem_f <- chem_elem %>%
        dplyr::filter(.data[[eixo_tax]] %in% intersect(keep_occ, keep_mols))
      
      if (nrow(chem_elem_f)) {
        order_classes <- chem_elem_f %>%
          dplyr::group_by(.data[[eixo_tax]]) %>%
          dplyr::summarise(OC_med = median(OC, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(OC_med)) %>%
          dplyr::pull(!!rlang::sym(eixo_tax))
        
        wrap_lab <- function(x, w = 36) stringr::str_wrap(x, width = w)
        
        d_plot <- chem_elem_f %>%
          dplyr::mutate(
            class  = factor(.data[[eixo_tax]], levels = order_classes),
            OC_cap = pmin(OC, 1.20)
          )
        
        sub_lab <- make_label_context(cfg, "oc_class")
        
        p_oc_class <- ggplot(d_plot, aes(x = class, y = OC_cap)) +
          geom_boxplot(outlier.alpha = .25, width = .6, color = "black", fill = "white") +
          stat_summary(fun = median, geom = "point", size = 1.6) +
          geom_jitter(width = .12, alpha = .20, size = .8, color = "grey35") +
          coord_flip() +
          labs(
            x = eixo_label,
            y = "O/C (cap=1.20)",
            title    = paste0("O/C by class — ", scope_label),
            subtitle = sub_lab
          ) +
          scale_x_discrete(labels = function(x) wrap_lab(x, 36)) +
          theme_classic(base_size = 11) +
          theme(
            plot.title   = element_text(hjust = .5, face = "bold"),
            axis.ticks.y = element_blank(),
            panel.grid   = element_blank()
          )
        
        save_plot_multi(p_oc_class, file.path(OUT_DIR, paste0(base_tag, "_OC_class")))
      }
    }
  }
}

## murcko frameworks
if (isTRUE(DO_MURCKO)) {
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    message("Murcko skipped (ChemmineR not installed).")
  } else if (!"murko_framework" %in% names(df_props)) {
    message("Murcko skipped: column 'murko_framework' not found in df_props.")
  } else {
    library(ChemmineR)
    
    MW_MIN <- as.numeric(cfg$murcko_min_fragment_mw %||% 200)
    TOPK   <- as.integer(cfg$murcko_top_k %||% 15)
    
    safe_sdf_from_smiles <- function(sm) {
      if (is.null(sm) || !nzchar(sm)) return(NULL)
      
      s1 <- try(
        suppressWarnings(suppressMessages(ChemmineR::smiles2sdf(sm))),
        silent = TRUE
      )
      if (!inherits(s1, "try-error") && length(s1) >= 1) {
        return(s1)
      }
      
      sm2 <- gsub("@@?", "", sm)
      s2 <- try(
        suppressWarnings(suppressMessages(ChemmineR::smiles2sdf(sm2))),
        silent = TRUE
      )
      if (!inherits(s2, "try-error") && length(s2) >= 1) {
        return(s2)
      }
      
      NULL
    }
    
    mw_from_smiles <- function(sm) {
      sdfset <- safe_sdf_from_smiles(sm)
      if (is.null(sdfset)) return(NA_real_)
      as.numeric(ChemmineR::MW(sdfset))[1]
    }
    
    base_fw <- df_props %>%
      dplyr::inner_join(
        dplyr::select(
          map_tax_inchi,
          inchikey,
          taxon_plot = !!rlang::sym(tax_col)
        ),
        by = "inchikey",
        relationship = "many-to-many"
      ) %>%
      dplyr::filter(
        !is.na(murko_framework), nzchar(murko_framework),
        !is.na(taxon_plot),      nzchar(taxon_plot)
      ) %>%
      dplyr::transmute(taxon_plot, murko_framework)
    
    if (!nrow(base_fw)) {
      message("Murcko: no frameworks after filters.")
    } else {
      uniq_fw <- dplyr::distinct(base_fw, murko_framework)
      
      uniq_fw$mw_fragment <- vapply(
        uniq_fw$murko_framework,
        mw_from_smiles,
        numeric(1)
      )
      
      uniq_fw_f <- uniq_fw %>%
        dplyr::filter(
          !is.na(mw_fragment),
          is.finite(mw_fragment),
          mw_fragment >= MW_MIN
        )
      
      if (!nrow(uniq_fw_f)) {
        message(sprintf(
          "Murcko: no framework with MW(fragment) ≥ %g.",
          MW_MIN
        ))
      } else {
        mf_long <- base_fw %>%
          dplyr::semi_join(uniq_fw_f, by = "murko_framework") %>%
          dplyr::count(murko_framework, taxon_plot, name = "n")
        
        mf_total <- mf_long %>%
          dplyr::group_by(murko_framework) %>%
          dplyr::summarise(
            total = sum(n),
            .groups = "drop"
          ) %>%
          dplyr::left_join(uniq_fw_f, by = "murko_framework") %>%
          dplyr::arrange(dplyr::desc(total))
        
        TOPK_EFF <- min(TOPK, nrow(mf_total))
        
        if (exists("mf_total") && nrow(mf_total) > 0) {
          n_avail_rank <- nrow(mf_total)
          n_show_rank  <- min(TOPK_EFF, n_avail_rank)
          
          rank_tbl <- mf_total %>%
            dplyr::slice_head(n = n_show_rank) %>%
            dplyr::mutate(rank = dplyr::row_number())
        } else {
          rank_tbl <- tibble::tibble()
        }
        
        xls_fw <- rank_tbl %>%
          dplyr::select(
            rank, murko_framework, mw_fragment, total
          ) %>%
          dplyr::arrange(rank)
        
        xls_by_taxon <- mf_long %>%
          dplyr::semi_join(rank_tbl, by = "murko_framework") %>%
          dplyr::arrange(murko_framework, dplyr::desc(n))
        
        xls_path <- file.path(
          OUT_DIR,
          sprintf(
            "%s_murcko_top_%d_MWfrag_ge_%g.xlsx",
            base_tag, TOPK_EFF, MW_MIN
          )
        )
        
        writexl::write_xlsx(
          list(
            frameworks      = xls_fw,
            counts_by_taxon = xls_by_taxon
          ),
          path = xls_path
        )
        message("✔ Murcko Excel saved at: ",
                normalizePath(xls_path, winslash = "/"))
        
        plot_df <- rank_tbl %>%
          dplyr::left_join(mf_long, by = "murko_framework") %>%
          dplyr::mutate(
            rank_fac = factor(
              sprintf("%02d", rank),
              levels = sprintf("%02d", rev(rank_tbl$rank))
            )
          )
        
        sub_lab <- make_label_context(cfg, "murcko")
        
        fill_lab <- dplyr::case_when(
          tax_col == "genus"   ~ "Genus",
          tax_col == "species" ~ "Species",
          tax_col == "family"  ~ "Family",
          TRUE                 ~ "Taxon"
        )
        
        p_murcko <- ggplot(
          plot_df,
          aes(x = rank_fac, y = n, fill = taxon_plot)
        ) +
          geom_col(position = "stack", width = 0.75) +
          coord_flip() +
          labs(
            x = "Framework rank (1 = most frequent)",
            y = "Count",
            title = sprintf(
              "Murcko frameworks (Top %d) — %s | MW(fragment) ≥ %g",
              TOPK_EFF, scope_label, MW_MIN
            ),
            subtitle = sub_lab,
            fill = fill_lab
          ) +
          theme_classic(base_size = 12)
        
        save_plot_multi(
          p_murcko,
          file.path(OUT_DIR, paste0(base_tag, "_murcko_top"))
        )
      }
    }
  }
}

## shared compounds
if (isTRUE(DO_SHARED_COMPOUNDS)) {
  req_objs <- c("map_tax_inchi", "uni_enriched", "OUT_DIR", "base_tag", "tax_col")
  miss     <- req_objs[!vapply(req_objs, exists, logical(1))]
  if (length(miss)) {
    stop("Block 11 (shared compounds): missing required objects: ",
         paste(miss, collapse = ", "))
  }
  if (!is.data.frame(map_tax_inchi) || !"inchikey" %in% names(map_tax_inchi)) {
    stop("Block 11: 'map_tax_inchi' must be a data.frame with 'inchikey' and tax_col.")
  }
  if (!tax_col %in% names(map_tax_inchi)) {
    stop("Block 11: tax_col = '", tax_col,
         "' not found in map_tax_inchi. Names: ",
         paste(names(map_tax_inchi), collapse = ", "))
  }
  
  MAX_XLSX <- 32767L
  trim_cell <- function(x, max_len = MAX_XLSX) {
    if (is.null(x)) return(x)
    if (!is.character(x)) x <- as.character(x)
    n <- nchar(x, allowNA = TRUE)
    too_long <- !is.na(n) & n > max_len
    x[too_long] <- paste0(substr(x[too_long], 1, max_len - 3), "...")
    enc2utf8(x)
  }
  
  cross <- map_tax_inchi %>%
    dplyr::filter(!is.na(.data[[tax_col]]), nzchar(.data[[tax_col]])) %>%
    dplyr::group_by(inchikey) %>%
    dplyr::summarise(
      taxa   = list(unique(.data[[tax_col]])),
      n_taxa = dplyr::n_distinct(.data[[tax_col]]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_taxa)) %>%
    dplyr::left_join(
      uni_enriched %>%
        dplyr::distinct(
          inchikey,
          iupac_name,
          traditional_name,
          .keep_all = FALSE
        ),
      by = "inchikey"
    )
  
  cross_shared <- cross %>%
    dplyr::filter(n_taxa >= 2)
  
  if (!nrow(cross_shared)) {
    message("Shared compounds: no compound found in ≥2 taxa for tax_col = '",
            tax_col, "'.")
  } else {
    
    resumo_wide <- cross_shared %>%
      dplyr::transmute(
        inchikey,
        iupac_name       = dplyr::coalesce(iupac_name, NA_character_),
        traditional_name = dplyr::coalesce(traditional_name, NA_character_),
        n_taxa,
        taxa = vapply(
          taxa,
          function(v)
            paste(
              sort(unique(as.character(v))),
              collapse = "; "
            ),
          FUN.VALUE = character(1),
          USE.NAMES = FALSE
        ),
        taxon_level = tax_col
      ) %>%
      dplyr::mutate(taxa = trim_cell(taxa))
    
    resumo_long <- cross_shared %>%
      tidyr::unnest_longer(taxa) %>%
      dplyr::rename(taxon = taxa) %>%
      dplyr::transmute(
        inchikey,
        iupac_name       = dplyr::coalesce(iupac_name, NA_character_),
        traditional_name = dplyr::coalesce(traditional_name, NA_character_),
        taxon            = as.character(taxon),
        taxon_level      = tax_col
      )
    
    out_xlsx <- file.path(OUT_DIR, paste0(base_tag, "_shared_compounds.xlsx"))
    writexl::write_xlsx(
      list(
        shared_compounds      = resumo_wide,
        shared_compounds_long = resumo_long
      ),
      path = out_xlsx
    )
    cat("✔ Saved:", normalizePath(out_xlsx, winslash = "/"), "\n")
    
    print(utils::head(resumo_wide, 20))
  }
}

## bibliometrics
if (isTRUE(DO_BIBLIOMETRICS)) {
  req_objs <- c("lin_enriched", "uni_enriched", "map_tax_inchi",
                "OUT_DIR", "base_tag", "tax_col")
  miss <- req_objs[!vapply(req_objs, exists, logical(1))]
  if (length(miss))
    stop("Block 12 (bibliometrics): missing required objects: ",
         paste(miss, collapse = ", "))
  
  if (!"ref_id" %in% names(lin_enriched))
    stop("Block 12: column 'ref_id' (DOI) not found in lin_enriched.")
  
  clean_doi <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("^(https?://(dx\\.)?doi\\.org/|doi:)", "", x, ignore.case = TRUE)
    x <- gsub("\\s+", "", x)
    x[nzchar(x)]
  }
  
  doi_to_url <- function(x) {
    x <- clean_doi(x)
    ifelse(nzchar(x), paste0("https://doi.org/", x), NA_character_)
  }
  
  make_hyperlink <- function(url, label = NULL) {
    label <- label %||% url
    ifelse(!is.na(url) & nzchar(url),
           paste0('=HYPERLINK("', url, '","', label, '")'),
           NA_character_)
  }
  
  comp_tax_total <- map_tax_inchi %>%
    normalize_taxcol(tax_col = tax_col) %>%
    dplyr::filter(!is.na(.data[[tax_col]]), nzchar(.data[[tax_col]])) %>%
    dplyr::distinct(!!rlang::sym(tax_col), inchikey)
  
  biblio_tax_comp_ref <- lin_enriched %>%
    dplyr::filter(!is.na(inchikey), nzchar(inchikey),
                  !is.na(ref_id),    nzchar0(ref_id)) %>%
    dplyr::mutate(ref_id = clean_doi(ref_id)) %>%
    dplyr::distinct(inchikey, ref_id) %>%
    dplyr::inner_join(
      dplyr::select(map_tax_inchi,
                    inchikey,
                    taxon = !!rlang::sym(tax_col)),
      by = "inchikey",
      relationship = "many-to-many"
    ) %>%
    dplyr::rename(!!rlang::sym(tax_col) := taxon) %>%
    normalize_taxcol(tax_col = tax_col) %>%
    dplyr::filter(!is.na(.data[[tax_col]]), nzchar(.data[[tax_col]])) %>%
    dplyr::distinct(!!rlang::sym(tax_col), inchikey, ref_id) %>%
    dplyr::mutate(doi_url = doi_to_url(ref_id))
  
  resumo_taxon <- comp_tax_total %>%
    dplyr::group_by(!!rlang::sym(tax_col)) %>%
    dplyr::summarise(n_compounds = dplyr::n_distinct(inchikey), .groups = "drop") %>%
    dplyr::left_join(
      biblio_tax_comp_ref %>%
        dplyr::group_by(!!rlang::sym(tax_col)) %>%
        dplyr::summarise(
          n_compounds_with_ref = dplyr::n_distinct(inchikey),
          n_refs_unique_total  = dplyr::n_distinct(ref_id),
          .groups = "drop"
        ),
      by = tax_col
    ) %>%
    dplyr::mutate(
      n_compounds_with_ref = coalesce(n_compounds_with_ref, 0L),
      n_refs_unique_total  = coalesce(n_refs_unique_total, 0L),
      pct_with_reference   = ifelse(n_compounds > 0,
                                    n_compounds_with_ref / n_compounds,
                                    NA_real_),
      taxon_level = tax_col
    ) %>%
    dplyr::arrange(dplyr::desc(n_compounds))
  
  comp_with_ref <- biblio_tax_comp_ref %>%
    dplyr::group_by(!!rlang::sym(tax_col), inchikey) %>%
    dplyr::summarise(
      n_refs  = dplyr::n_distinct(ref_id),
      ref_ids = paste(sort(unique(ref_id)), collapse = "; "),
      doi_urls = paste(sort(unique(doi_url)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      uni_enriched %>%
        dplyr::select(inchikey, iupac_name, traditional_name) %>%
        dplyr::distinct(),
      by = "inchikey"
    ) %>%
    dplyr::mutate(
      hyperlink = make_hyperlink(strsplit(doi_urls, ";")[[1]][1],
                                 strsplit(ref_ids, ";")[[1]][1]),
      taxon_level = tax_col
    ) %>%
    dplyr::arrange(!!rlang::sym(tax_col), inchikey)
  
  comp_without_ref <- comp_tax_total %>%
    dplyr::anti_join(
      biblio_tax_comp_ref %>% dplyr::distinct(!!rlang::sym(tax_col), inchikey),
      by = c(tax_col, "inchikey")
    ) %>%
    dplyr::left_join(
      uni_enriched %>%
        dplyr::select(inchikey, iupac_name, traditional_name) %>%
        dplyr::distinct(),
      by = "inchikey"
    ) %>%
    dplyr::mutate(taxon_level = tax_col)
  
  out_xlsx_bib <- file.path(OUT_DIR, paste0(base_tag, "_bibliometrics.xlsx"))
  writexl::write_xlsx(
    list(
      summary_by_taxon      = resumo_taxon,
      compounds_with_refs   = comp_with_ref,
      compounds_without_ref = comp_without_ref
    ),
    path = out_xlsx_bib
  )
  
  cat("✔ Saved:", normalizePath(out_xlsx_bib, winslash = "/"), "\n")
  print(utils::head(resumo_taxon, 20))
}


## pca biplot
if (!exists("DO_PCA_PROPS")) DO_PCA_PROPS <- TRUE

if (isTRUE(DO_PCA_PROPS)) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
    library(stringr); library(tibble); library(grid)
  })
  
  req_objs <- c("lin_enriched", "uni_enriched", "cfg", "OUT_DIR", "base_tag", "scope_label")
  miss <- req_objs[!vapply(req_objs, exists, logical(1))]
  if (length(miss)) stop("PCA classes: missing objects: ", paste(miss, collapse = ", "))
  
  pdf_device <- grDevices::pdf
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  if (!exists("tax_col")) {
    tax_col <- tolower(cfg$analysis_tax_level %||% cfg$taxon_mode %||% "species")
  }
  tax_col <- tolower(tax_col)
  
  group_by_var <- tax_col
  
  label_face <- face_for_taxon(tax_col)
  spec_label <- rep(NA_character_, nrow(lin_enriched))
  if ("accepted_name" %in% names(lin_enriched)) spec_label <- lin_enriched$accepted_name
  if (all(is.na(spec_label)) && "species" %in% names(lin_enriched)) spec_label <- lin_enriched$species
  
  clean_pca_class <- function(x) {
    x <- as.character(x)
    x <- sub("\\|.*", "", x) 
    trimws(x)
  }
  
  if (exists("CHEM_COL_ALVO") && CHEM_COL_ALVO %in% names(uni_enriched)) {
    if (CHEM_COL_ALVO %in% names(lin_enriched)) lin_enriched[[CHEM_COL_ALVO]] <- NULL
    lin_enriched <- lin_enriched %>%
      dplyr::left_join(uni_enriched %>% dplyr::select(inchikey, !!rlang::sym(CHEM_COL_ALVO)), by = "inchikey")
    axis_col  <- CHEM_COL_ALVO
    title_ctx <- "Chemical Classes (Hybrid)"
  } else {
    axis_col  <- "chemicalTaxonomyClassyfireClass"
    title_ctx <- "ClassyFire Classes"
  }
  lin_enriched[[axis_col]] <- clean_pca_class(lin_enriched[[axis_col]])
  
  MIN_N_COMPOUNDS <- as.integer(cfg$pca_classes_min_compounds_per_taxon %||% 30L)
  TOP_LOADINGS    <- as.integer(cfg$pca_classes_top_loadings %||% 12L)
  POINT_ALPHA     <- 0.8
  POINT_SIZE      <- 2.5 
  LABEL_SIZE_TAX  <- 3.0
  LABEL_SIZE_VAR  <- 3.5
  
  lin_use <- lin_enriched %>%
    dplyr::mutate(
      taxon = dplyr::case_when(
        tax_col == "species" ~ spec_label,
        tax_col == "genus"   ~ .data$genus,
        tax_col == "family"  ~ .data$family,
        TRUE                 ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(taxon), nzchar(taxon), !is.na(.data[[axis_col]]), nzchar(.data[[axis_col]]))
  
  keep_taxa <- lin_use %>%
    dplyr::count(taxon, name="n") %>%
    dplyr::filter(n >= MIN_N_COMPOUNDS) %>%
    dplyr::pull(taxon)
  
  lin_use <- lin_use %>% dplyr::filter(taxon %in% keep_taxa)
  
  if (nrow(lin_use) > 0) {
    tbl <- lin_use %>%
      dplyr::group_by(taxon, !!rlang::sym(axis_col)) %>%
      dplyr::summarise(n = dplyr::n_distinct(inchikey), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(axis_col), values_from = n, values_fill = 0)
    
    X <- as.data.frame(tbl); rownames(X) <- X$taxon; X$taxon <- NULL
    X_prop <- sweep(X, 1, pmax(1, rowSums(X)), "/")
    keep_cols <- vapply(X_prop, function(v) sd(v, na.rm=TRUE) > 0, logical(1))
    X_prop <- X_prop[, keep_cols, drop=FALSE]
    
    if (ncol(X_prop) >= 2) {
      pc <- prcomp(X_prop, center = TRUE, scale. = TRUE)
      scores_df <- as.data.frame(pc$x[, 1:2])
      names(scores_df) <- c("PC1", "PC2")
      scores_df$taxon <- rownames(scores_df)
      vexp <- round(100 * (pc$sdev^2) / sum(pc$sdev^2), 1)
      
      loadings <- as.data.frame(pc$rotation[, 1:2])
      loadings$class <- clean_pca_class(rownames(loadings))
      loadings$mag   <- sqrt(loadings$PC1^2 + loadings$PC2^2)
      
      arr_df <- loadings %>% arrange(desc(mag)) %>% slice_head(n = TOP_LOADINGS)
      
      mult <- min(
        (max(scores_df$PC1) - min(scores_df$PC1))/(max(arr_df$PC1)-min(arr_df$PC1)),
        (max(scores_df$PC2) - min(scores_df$PC2))/(max(arr_df$PC2)-min(arr_df$PC2))
      ) * 0.85
      arr_df <- arr_df %>% mutate(xend = PC1 * mult, yend = PC2 * mult)
      
      scores_grp <- scores_df %>%
        dplyr::left_join(lin_use %>% dplyr::distinct(taxon, !!rlang::sym(group_by_var)), by="taxon")
      
      if(!group_by_var %in% names(scores_grp)) {
        if(group_by_var == tax_col) scores_grp[[group_by_var]] <- scores_grp$taxon
        else scores_grp[[group_by_var]] <- "Unknown"
      }
      
      n_groups <- length(unique(scores_grp[[group_by_var]]))
      if (n_groups <= 12) {
        pal_ok <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                    "#E6AB02", "#A6761D", "#666666", "#8DA0CB", "#FC8D62", "#FFD92F", "#A6D854")
        scale_col_layer <- scale_color_manual(values = pal_ok, na.value = "gray60")
      } else {
        scale_col_layer <- scale_color_discrete() 
      }
      
      limx <- max(abs(range(scores_df$PC1))) * 1.15
      limy <- max(abs(range(scores_df$PC2))) * 1.15
      
      p_clean <- ggplot(scores_grp, aes(PC1, PC2, color = .data[[group_by_var]])) +
        geom_hline(yintercept=0, linetype="dotted", color="gray") +
        geom_vline(xintercept=0, linetype="dotted", color="gray") +
        geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
        
        geom_segment(data = arr_df, aes(x=0, y=0, xend=xend, yend=yend), 
                     inherit.aes=FALSE, arrow=arrow(length=unit(0.2,"cm")), color="black", alpha=0.7) +
        
        ggrepel::geom_text_repel(
          data = arr_df, 
          aes(x=xend, y=yend, label=class), 
          inherit.aes=FALSE, 
          size=LABEL_SIZE_VAR, 
          fontface="bold", 
          color="black",
          max.overlaps = Inf,
          min.segment.length = 0,
          box.padding = 0.5
        ) +
        
        scale_col_layer +
        labs(
          title = paste0("PCA (Clean): ", title_ctx, " — ", scope_label),
          subtitle = paste0("Points/Color = ", tools::toTitleCase(group_by_var)),
          x = paste0("PC1 (", vexp[1], "%)"), y = paste0("PC2 (", vexp[2], "%)"),
          color = tools::toTitleCase(group_by_var)
        ) +
        theme_classic(base_size = 12) +
        coord_cartesian(xlim = c(-limx, limx), ylim = c(-limy, limy))
      
      p_labeled <- p_clean +
        labs(title = paste0("PCA (Labeled): ", title_ctx, " — ", scope_label)) +
        ggrepel::geom_text_repel(aes(label=taxon), size=LABEL_SIZE_TAX, fontface=label_face, 
                                 max.overlaps=25, show.legend=FALSE)
      
      fn_clean <- file.path(OUT_DIR, paste0(base_tag, "_PCA_CLEAN.pdf"))
      fn_label <- file.path(OUT_DIR, paste0(base_tag, "_PCA_LABELED.pdf"))
      
      pdf_device(fn_clean, width=11, height=9); print(p_clean); dev.off()
      pdf_device(fn_label, width=11, height=9); print(p_labeled); dev.off()
      
      message("Saved PCA (Clean): ", fn_clean)
      message("Saved PCA (Labeled): ", fn_label)
    }
  }
}

## pcoa

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel); library(vegan); library(writexl)
})

DO_PCOA <- isTRUE(cfg$do_pcoa %||% TRUE)

if (isTRUE(DO_PCOA)) {
  
  trySuppress <- function(expr) suppressWarnings(suppressMessages(try(expr, silent = TRUE)))
  
  target_level <- tolower(cfg$analysis_tax_level %||% cfg$taxon_mode %||% "genus")
  
  if (!target_level %in% names(lin_enriched)) {
    if ("species" %in% names(lin_enriched)) target_level <- "species"
    else if ("genus" %in% names(lin_enriched)) target_level <- "genus"
    else target_level <- "family"
    warning("PCoA: Configured level not found. Falling back to: ", target_level)
  }
  
  tax_col <- target_level

  clean_pcoa_class <- function(x) {
    x <- as.character(x); x <- sub("\\|.*", "", x); trimws(x)
  }
  
  if (exists("CHEM_COL_ALVO") && CHEM_COL_ALVO %in% names(uni_enriched)) {
    if (!CHEM_COL_ALVO %in% names(lin_enriched)) {
      lin_enriched <- lin_enriched %>%
        dplyr::left_join(uni_enriched %>% dplyr::select(inchikey, !!rlang::sym(CHEM_COL_ALVO)), by = "inchikey")
    }
    CHEM_AXIS_PCOA <- CHEM_COL_ALVO
  } else {
    CHEM_AXIS_PCOA <- "chemicalTaxonomyClassyfireClass"
  }
  
  if(CHEM_AXIS_PCOA %in% names(lin_enriched)) {
    lin_enriched[[CHEM_AXIS_PCOA]] <- clean_pcoa_class(lin_enriched[[CHEM_AXIS_PCOA]])
  }

  PCOA_DIST       <- tolower(cfg$pcoa_dist %||% "jaccard")  
  PCOA_BINARY     <- isTRUE(cfg$pcoa_binary %||% TRUE)
  PCOA_MIN_TAX_N  <- as.integer(cfg$pcoa_min_tax_n %||% 3L) 
  PCOA_MIN_FEAT_N <- as.integer(cfg$pcoa_min_feat_n %||% 1L)
  PCOA_POINT_SZ   <- 2.5; PCOA_ALPHA <- 0.8; PCOA_LABEL_SZ <- 3.0; PCOA_ENVFIT_SZ <- 3.5
  KEEP_G_ONLY     <- isTRUE(cfg$keep_genus_only_in_species %||% FALSE)
  pdf_device      <- grDevices::pdf 
  
  if (!exists("OUT_DIR")) OUT_DIR <- getwd()
  
  build_map_pcoa <- function(lin_tbl, tax_c, keep_g){
    if (!tax_c %in% names(lin_tbl)) stop("[PCoA Error] Column '", tax_c, "' missing.")
    
    out <- lin_tbl %>%
      dplyr::select(inchikey, family, genus, species, !!rlang::sym(tax_c), 
                    chem_col = dplyr::any_of(CHEM_AXIS_PCOA)) %>%
      dplyr::mutate(inchikey=as.character(inchikey), chem_class=as.character(chem_col)) %>%
      dplyr::filter(nzchar0(inchikey)) %>%
      dplyr::mutate(chem_class = clean_pcoa_class(chem_class))
    
    if (identical(tax_c, "species")) {
      out <- out %>% 
        dplyr::filter(!is.na(species), species != "NA", !grepl("sp\\.$", species))
    }
    
    out %>% 
      dplyr::mutate(tax_key = as.character(.data[[tax_c]])) %>% 
      dplyr::filter(nzchar0(tax_key), !is.na(tax_key)) %>% 
      dplyr::distinct(inchikey, tax_key, .keep_all=TRUE)
  }
  
  map_tax_inchi_pcoa <- build_map_pcoa(lin_enriched, tax_col, KEEP_G_ONLY)
  
  if (nrow(map_tax_inchi_pcoa) == 0) {
    warning("[PCoA] No data available for level: ", tax_col)
  } else {
    
    pa_mat_raw <- tryCatch({
      xtabs(~ tax_key + inchikey, data = map_tax_inchi_pcoa)
    }, error = function(e) NULL)
    
    if (is.null(pa_mat_raw)) {
      warning("[PCoA] Failed to build crosstab matrix.")
    } else {
      
      pa_mat <- as.matrix(pa_mat_raw)
      pa_mat[pa_mat > 0] <- 1
      storage.mode(pa_mat) <- "double"
      
      if (nrow(pa_mat) > 0) {
        keep_rows <- rowSums(pa_mat) >= PCOA_MIN_TAX_N
        pa_mat <- pa_mat[keep_rows, , drop=FALSE]
        if(nrow(pa_mat) > 0) pa_mat <- pa_mat[, colSums(pa_mat) > 0, drop=FALSE]
      }
      
      if (nrow(pa_mat) < 3 || ncol(pa_mat) < 2) {
        warning(paste0("[PCoA] Matrix too small (Rows=", nrow(pa_mat), "). Need >= 3 taxa. Reduce 'analysis_min_compounds'."))
      } else {
        
        d_obj <- trySuppress(vegan::vegdist(pa_mat, method=PCOA_DIST, binary=PCOA_BINARY))
        
        if (inherits(d_obj, "try-error") || is.null(d_obj) || any(is.na(d_obj))) {
          message("[PCoA] Jaccard failed/NA. Trying Bray-Curtis...")
          d_obj <- trySuppress(vegan::vegdist(pa_mat, method="bray", binary=PCOA_BINARY))
        }
        
        if (!inherits(d_obj, "try-error") && !is.null(d_obj)) {
          cmd <- stats::cmdscale(d_obj, k=2, eig=TRUE)
          vexp <- round(100 * pmax(cmd$eig, 0)[1:2] / sum(pmax(cmd$eig, 0)), 1)
          
          scores <- as.data.frame(cmd$points); colnames(scores) <- c("PCoA1","PCoA2")
          scores$tax_key <- rownames(scores)
          
          meta_cols <- intersect(c("family", "genus", "species"), names(map_tax_inchi_pcoa))
          ref_table <- map_tax_inchi_pcoa %>% 
            dplyr::select(tax_key, dplyr::all_of(meta_cols)) %>% 
            dplyr::distinct(tax_key, .keep_all = TRUE)
          
          df_plot <- scores %>% dplyr::inner_join(ref_table, by="tax_key")
          df_plot$tax_lab <- df_plot$tax_key
          
          # Color Logic
          lbl_col <- if (tax_col == "species" && "genus" %in% names(df_plot)) "genus" else tax_col
          
          # Classes
          sig_vecs_chem <- NULL
          if ("chem_class" %in% names(map_tax_inchi_pcoa)) {
            mat_chem <- map_tax_inchi_pcoa %>%
              dplyr::filter(tax_key %in% rownames(cmd$points), !is.na(chem_class), nzchar(trimws(chem_class))) %>%
              dplyr::count(tax_key, chem_class) %>%
              tidyr::pivot_wider(names_from=chem_class, values_from=n, values_fill=0) %>%
              as.data.frame()
            
            if(nrow(mat_chem) > 0) {
              rownames(mat_chem) <- mat_chem$tax_key; mat_chem$tax_key <- NULL
              mat_chem <- mat_chem[rownames(cmd$points), , drop=FALSE]
              mat_chem <- mat_chem[, colSums(mat_chem) > 0, drop=FALSE]
              
              if(ncol(mat_chem) > 1) {
                ef_c <- trySuppress(vegan::envfit(cmd, mat_chem, perm=999, na.rm=TRUE))
                if(!inherits(ef_c, "try-error") && !is.null(ef_c)) {
                  sc_c <- as.data.frame(vegan::scores(ef_c, display="vectors"))
                  sc_c$pval <- ef_c$vectors$pvals
                  sc_c$varname <- clean_pcoa_class(rownames(sc_c))
                  if(ncol(sc_c) >= 2) { 
                    colnames(sc_c)[1:2] <- c("PCoA1","PCoA2")
                    sig_vecs_chem <- sc_c[sc_c$pval <= 0.05, ] 
                  }
                }
              }
            }
          }
          
          make_pcoa_plot <- function(data_pts, vectors, title_suffix, vec_color, file_suffix) {
            
            lim_xy <- max(max(abs(data_pts$PCoA1)), max(abs(data_pts$PCoA2))) * 1.2
            
            n_groups <- length(unique(data_pts[[lbl_col]]))
            scale_col <- if(n_groups <= 12) scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DA0CB", "#FC8D62", "#FFD92F", "#A6D854")) else scale_color_discrete()
            
            # CLEAN PLOT
            p_clean <- ggplot(data_pts, aes(x=PCoA1, y=PCoA2)) +
              geom_point(aes(color=.data[[lbl_col]]), size=PCOA_POINT_SZ, alpha=PCOA_ALPHA) +
              coord_cartesian(xlim=c(-lim_xy, lim_xy), ylim=c(-lim_xy, lim_xy)) +
              labs(
                title = paste0("PCoA (Clean): ", title_suffix), 
                x = paste0("PCoA1 (", vexp[1], "%)"), y = paste0("PCoA2 (", vexp[2], "%)"), 
                color=tools::toTitleCase(lbl_col)
              ) +
              theme_minimal() + theme(plot.margin=margin(10,10,10,10)) + scale_col
            
            if (!is.null(vectors) && nrow(vectors) > 0) {
              mult <- 0.7 * lim_xy
              p_clean <- p_clean + 
                geom_segment(data=vectors, aes(x=0,y=0,xend=PCoA1*mult, yend=PCoA2*mult), arrow=arrow(length=unit(0.2,"cm")), color=vec_color, lwd=0.6) +
                ggrepel::geom_text_repel(data=vectors, aes(x=PCoA1*mult, y=PCoA2*mult, label=varname), color=vec_color, fontface="bold", size=PCOA_ENVFIT_SZ, max.overlaps=Inf)
            }
            
            # LABELED PLOT
            p_labeled <- p_clean + 
              labs(title = paste0("PCoA (Labeled): ", title_suffix)) +
              ggrepel::geom_text_repel(aes(label=tax_lab), size=PCOA_LABEL_SZ, max.overlaps= 20, min.segment.length=0, fontface=face_for_taxon(tax_col))
            
            FIG_DIR <- file.path(OUT_DIR, "fig_pcoa")
            dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)
            
            fn_clean <- file.path(FIG_DIR, paste0(base_tag, "_pcoa_", file_suffix, "_CLEAN.pdf"))
            fn_label <- file.path(FIG_DIR, paste0(base_tag, "_pcoa_", file_suffix, "_LABELED.pdf"))
            
            pdf_device(fn_clean, 10, 10); print(p_clean); dev.off()
            pdf_device(fn_label, 10, 10); print(p_labeled); dev.off()
            message("Saved PCoA set: ", file_suffix)
          }
          
          if (!is.null(sig_vecs_chem)) {
            make_pcoa_plot(df_plot, sig_vecs_chem, "(Hybrid Class Vectors)", "black", "ChemClassVectors")
          } else {
            make_pcoa_plot(df_plot, NULL, "(Clean)", "black", "Basic")
          }
          
          # Export Excel
          writexl::write_xlsx(
            df_plot %>% dplyr::select(tax_key, PCoA1, PCoA2, any_of(c("family", "genus", "species"))) %>% dplyr::arrange(desc(PCoA1)), 
            file.path(OUT_DIR, paste0(base_tag, "_PCoA_Coordinates.xlsx"))
          )
        }
      }
    }
  }
}


## statistical suite
if (isTRUE(cfg$run_module3)) {
  cat("Running statistical suite...\n")

  if (!exists("df_props") || nrow(df_props) == 0) {
    stop("df_props not found. Please run Part III sections before Part IV.")
  }
  
  df_stat_base <- df_props %>%
    dplyr::mutate(
      number_of_carbons   = to_num(number_of_carbons),
      number_of_oxygens   = to_num(number_of_oxygens),
      OC_Ratio = dplyr::if_else(number_of_carbons > 0,
                                number_of_oxygens / number_of_carbons, NA_real_),
      molecular_weight = to_num(molecular_weight),
      xlogp            = to_num(xlogp),
      topoPSA          = to_num(topoPSA),
      fsp3             = to_num(fsp3),
      n_rings          = to_num(max_number_of_rings)
    )
  
  vars_to_test <- c("molecular_weight", "xlogp", "topoPSA", "fsp3", 
                    "n_rings", "OC_Ratio", "hBondDonorCount", "hBondAcceptorCount")
  vars_to_test <- intersect(vars_to_test, names(df_stat_base))
  
  taxa_validos <- if (exists("taxa_keep")) taxa_keep else unique(df_stat_base[[tax_col]])
  
  df_stat_clean <- df_stat_base %>%
    dplyr::filter(.data[[tax_col]] %in% taxa_validos) %>%
    dplyr::mutate(Group = as.factor(.data[[tax_col]]))
  
  
  stats_summary_list <- list()
  stats_pairwise_list <- list()
  
  for (var in vars_to_test) {
    # Remove NAs for the current test
    curr_data <- df_stat_clean %>% 
      dplyr::select(Group, Value = !!rlang::sym(var)) %>% 
      tidyr::drop_na()
    
    if (nrow(curr_data) > 0 && length(unique(curr_data$Group)) > 1) {
      kw_test <- kruskal.test(Value ~ Group, data = curr_data)
      
      desc <- curr_data %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(
          Mean = mean(Value),
          Median = median(Value),
          SD = sd(Value),
          N = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(Median)) %>%
        dplyr::mutate(
          Variable = var,
          Global_P_Value = kw_test$p.value,
          Significance = dplyr::case_when(
            kw_test$p.value < 0.001 ~ "***",
            kw_test$p.value < 0.01  ~ "**",
            kw_test$p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
        )
      stats_summary_list[[var]] <- desc
      
      if (!is.na(kw_test$p.value) && kw_test$p.value < 0.05) {
        pair_test <- pairwise.wilcox.test(curr_data$Value, curr_data$Group,
                                          p.adjust.method = "BH")
        if (!is.null(pair_test$p.value)) {
          pair_df <- as.data.frame.table(pair_test$p.value) %>%
            stats::na.omit() %>%
            dplyr::rename(Group1 = Var1, Group2 = Var2, P_Adj = Freq) %>%
            dplyr::mutate(
              Variable = var,
              Is_Significant = P_Adj < 0.05
            ) %>%
            dplyr::filter(Is_Significant == TRUE)
          
          stats_pairwise_list[[var]] <- pair_df
        }
      }
    }
  }
  
  if (length(stats_summary_list) > 0) {
    full_summary <- dplyr::bind_rows(stats_summary_list) %>% dplyr::select(Variable, everything())
    full_pairwise <- dplyr::bind_rows(stats_pairwise_list)
    
    out_phys <- file.path(OUT_DIR, paste0(base_tag, "_STATS_PhysChem.xlsx"))
    writexl::write_xlsx(list(
      Descriptive_Global = full_summary,
      Pairwise_Significant = full_pairwise
    ), path = out_phys)
    message("      [Saved] PhysChem Stats: ", basename(out_phys))
  }
  
  
  chem_col_stat <- if (exists("CHEM_COL_ALVO") && CHEM_COL_ALVO %in% names(df_props)) {
    CHEM_COL_ALVO 
  } else { "chemicalTaxonomyClassyfireClass" }
  
  if (chem_col_stat %in% names(df_props)) {
    
    df_chem <- df_props %>%
      dplyr::filter(!is.na(.data[[chem_col_stat]]), 
                    nzchar(as.character(.data[[chem_col_stat]])),
                    .data[[tax_col]] %in% taxa_validos)
    
    top_classes <- df_chem %>%
      dplyr::count(Class = .data[[chem_col_stat]], sort=TRUE) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(Class)
    
    enrichment_list <- list()
    
    for (taxon in taxa_validos) {
      for (chem_cls in top_classes) {
        a <- sum(df_chem[[tax_col]] == taxon & df_chem[[chem_col_stat]] == chem_cls)
        b <- sum(df_chem[[tax_col]] == taxon & df_chem[[chem_col_stat]] != chem_cls)
        c <- sum(df_chem[[tax_col]] != taxon & df_chem[[chem_col_stat]] == chem_cls)
        d <- sum(df_chem[[tax_col]] != taxon & df_chem[[chem_col_stat]] != chem_cls)
        
        if ((a + b) > 0 && (a + c) > 0) {
          fisher_res <- fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")
          
          if (fisher_res$p.value < 0.05) {
            enrichment_list[[paste(taxon, chem_cls)]] <- data.frame(
              Taxon = taxon,
              Chemical_Class = chem_cls,
              Count_In_Taxon = a,
              Odds_Ratio = as.numeric(fisher_res$estimate),
              P_Value = fisher_res$p.value
            )
          }
        }
      }
    }
    
    if (length(enrichment_list) > 0) {
      df_enrich <- dplyr::bind_rows(enrichment_list) %>%
        dplyr::arrange(Taxon, P_Value) %>%
        dplyr::mutate(FDR_Adj_P = p.adjust(P_Value, method = "BH"))
      
      out_enrich <- file.path(OUT_DIR, paste0(base_tag, "_STATS_Chem_Enrichment.xlsx"))
      writexl::write_xlsx(list(Enrichment_Results = df_enrich), path = out_enrich)
      message("      [Saved] Chemical Enrichment: ", basename(out_enrich))
    }
  } else {
    message("      [Skipped] Chemical column not found for enrichment.")
  }
  
  cat("Statistical suite done.\n")
}

## scaffold innovation
if (isTRUE(cfg$run_module3) && exists("df_props")) {
  cat("Running scaffold innovation analysis...\n")

  df_scaff <- df_props %>%
    dplyr::filter(
      !is.na(murko_framework), 
      nzchar(murko_framework),
      .data[[tax_col]] %in% taxa_keep
    ) %>%
    dplyr::select(inchikey, murko_framework, Taxon = !!rlang::sym(tax_col))
  
  if (nrow(df_scaff) > 0) {
    scaffold_stats <- df_scaff %>%
      dplyr::group_by(Taxon) %>%
      dplyr::summarise(
        N_Compounds = dplyr::n_distinct(inchikey),
        N_Scaffolds = dplyr::n_distinct(murko_framework),
        Innovation_Ratio = round(N_Scaffolds / N_Compounds, 3),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Innovation_Ratio))
    
    scaff_distribution <- df_scaff %>%
      dplyr::distinct(murko_framework, Taxon) %>%
      dplyr::count(murko_framework, name = "Species_Count")
    
    unique_scaffolds <- scaff_distribution %>% 
      dplyr::filter(Species_Count == 1) %>%
      dplyr::pull(murko_framework)
    
    exclusivity_stats <- df_scaff %>%
      dplyr::filter(murko_framework %in% unique_scaffolds) %>%
      dplyr::group_by(Taxon) %>%
      dplyr::summarise(
        N_Exclusive_Scaffolds = dplyr::n_distinct(murko_framework),
        Example_Exclusive_SMILES = dplyr::first(murko_framework),
        .groups = "drop"
      )
    
    final_stats <- scaffold_stats %>%
      dplyr::left_join(exclusivity_stats, by = "Taxon") %>%
      dplyr::mutate(
        N_Exclusive_Scaffolds = tidyr::replace_na(N_Exclusive_Scaffolds, 0),
        Pct_Exclusive = round((N_Exclusive_Scaffolds / N_Scaffolds) * 100, 1)
      ) %>%
      dplyr::arrange(dplyr::desc(N_Exclusive_Scaffolds))
    
    out_scaff <- file.path(OUT_DIR, paste0(base_tag, "_STATS_Scaffold_Innovation.xlsx"))
    writexl::write_xlsx(list(Scaffold_Metrics = final_stats), path = out_scaff)
    
    message("      [Saved] Scaffold Innovation Stats: ", basename(out_scaff))
    
    top_innovator <- final_stats$Taxon[1]
    cat("Top exclusive taxon:", top_innovator,
        "(", final_stats$N_Exclusive_Scaffolds[1], "unique scaffolds)\n")
    
  } else {
    message("      [Skipped] No Murcko frameworks found in dataset.")
  }
  cat("Scaffold analysis done.\n")
}