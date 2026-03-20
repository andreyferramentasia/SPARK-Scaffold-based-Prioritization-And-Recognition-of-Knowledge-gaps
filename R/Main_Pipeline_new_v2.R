`%||%` <- function(a, b) if (is.null(a)) b else a

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
})

# resolve script location whether called via Rscript or source()
script_dir <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  tryCatch({
    dirname(normalizePath(commandArgs(trailingOnly = FALSE) |>
      (\(a) a[grepl("--file=", a)])() |>
      (\(a) sub("--file=", "", a))()))
  }, error = function(e2) getwd())
})

args        <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1 && file.exists(args[1])) args[1] else NULL

if (!is.null(config_file)) {
  cat("Lendo configuração de:", config_file, "\n")
  cfg_json <- jsonlite::read_json(config_file, simplifyVector = TRUE)

  cfg_mode   <- cfg_json$taxon_mode         %||% "family"
  cfg_values <- cfg_json$taxon_values       %||% c("Phyllanthaceae")
  cfg_level  <- cfg_json$analysis_tax_level %||% "genus"

  cfg <- list(
    taxon_mode                       = cfg_mode,
    taxon_values                     = cfg_values,
    analysis_tax_level               = cfg_level,
    run_module1                      = isTRUE(cfg_json$run_module1),
    run_module2                      = isTRUE(cfg_json$run_module2),
    run_module3                      = isTRUE(cfg_json$run_module3),
    run_module4                      = isTRUE(cfg_json$run_module4),
    analysis_top_taxa                = as.integer(cfg_json$analysis_top_taxa                %||% 40L),
    analysis_min_compounds_per_taxon = as.integer(cfg_json$analysis_min_compounds_per_taxon %||% 10L),
    out_dir_base                     = cfg_json$out_dir_base %||% file.path(script_dir, "..", "results"),
    run_tag_date                     = Sys.Date(),
    mongo_url                        = cfg_json$mongo_url   %||% "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000",
    db_name                          = cfg_json$db_name     %||% "lotus",
    coll_name                        = cfg_json$coll_name   %||% "lotusUniqueNaturalProduct",
    use_WFO_normalization            = isTRUE(cfg_json$use_WFO_normalization),
    wfo_csv_path                     = cfg_json$wfo_csv_path %||% file.path(script_dir, "..", "DBs", "classification.tsv"),
    export_excel                     = isTRUE(cfg_json$export_excel  %||% TRUE),
    export_parquet                   = TRUE,
    verbose                          = isTRUE(cfg_json$verbose %||% TRUE)
  )

} else {
  cfg_mode   <- "family"
  cfg_values <- c("Phyllanthaceae")
  cfg_level  <- "genus"

  cfg <- list(
    taxon_mode   = cfg_mode,
    taxon_values = cfg_values,

    analysis_tax_level = cfg_level,

    run_module1 = FALSE,  # Part I:   extração MongoDB
    run_module2 = FALSE,  # Part II:  catálogo de estruturas (opcional)
    run_module3 = FALSE,  # Part III: estatísticas & figuras
    run_module4 = TRUE,   # Part IV:  contexto global & raridade

    analysis_top_taxa                = 40L,
    analysis_min_compounds_per_taxon = 10L,

    out_dir_base = file.path(script_dir, "..", "results"),
    run_tag_date = Sys.Date(),

    mongo_url = "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000",
    db_name   = "lotus",
    coll_name = "lotusUniqueNaturalProduct",

    use_WFO_normalization = TRUE,
    wfo_csv_path          = file.path(script_dir, "..", "DBs", "classification.tsv"),
    export_excel          = TRUE,
    export_parquet        = TRUE,
    verbose               = TRUE
  )
}

cfg$taxon_mode         <- tolower(cfg$taxon_mode)
cfg$analysis_tax_level <- tolower(cfg$analysis_tax_level)

cat("SPARK pipeline iniciado\n")
cat("Modo:", cfg$taxon_mode, "-> nível:", cfg$analysis_tax_level, "\n")
cat("Alvo:", paste(cfg$taxon_values, collapse = ", "), "\n")
cat("Saída:", normalizePath(cfg$out_dir_base, mustWork = FALSE), "\n")

dir.create(cfg$out_dir_base, showWarnings = FALSE, recursive = TRUE)

safe_tag_fn <- function(mode, values, run_date, analysis_level) {
  tag <- paste(mode, analysis_level, paste(values, collapse = "-"),
               format(run_date, "%Y%m%d"), sep = "_")
  gsub("[^A-Za-z0-9._-]+", "_", tag)
}

tag_base_load <- safe_tag_fn(cfg$taxon_mode, cfg$taxon_values,
                              cfg$run_tag_date, cfg$analysis_tax_level)
load_dir      <- file.path(cfg$out_dir_base, paste0("lotus_", tag_base_load))

OUT_DIR      <<- load_dir
base_tag     <<- paste0("lotus_", tag_base_load)
TAXON_MODE   <<- cfg$taxon_mode
TAXON_VALUES <<- cfg$taxon_values

if (isTRUE(cfg$run_module1)) {
  cat("\n[1/4] Part I: extração e limpeza...\n")
  assign("cfg", cfg, envir = .GlobalEnv)
  source(file.path(script_dir, "Part I - Extraction_new_v2.R"))
  cat("Part I concluída.\n")
}

if (isTRUE(cfg$run_module2)) {
  cat("\n[2/4] Part II: catálogo de estruturas...\n")
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    cat("Part II ignorada: pacote 'ChemmineR' não instalado.\n")
    cat("Para instalar: BiocManager::install('ChemmineR')\n")
  } else {
    assign("cfg", cfg, envir = .GlobalEnv)
    tryCatch(
      { source(file.path(script_dir, "Part II - Structures_new_v2.R")); cat("Part II concluída.\n") },
      error = function(e) cat("Part II falhou (continuando):", conditionMessage(e), "\n")
    )
  }
}

## auto-load: tenta recuperar dados de parquet se os objetos não estão em memória
if ((isTRUE(cfg$run_module3) || isTRUE(cfg$run_module4)) &&
    (!exists("lin_enriched") || !exists("uni_enriched"))) {

  cat("\nObjetos não encontrados em memória. Buscando parquet...\n")

  pq_lin <- file.path(load_dir, paste0("lotus_", tag_base_load, "_lin_enriched.parquet"))
  pq_uni <- file.path(load_dir, paste0("lotus_", tag_base_load, "_uni_enriched.parquet"))

  if (file.exists(pq_lin) && file.exists(pq_uni)) {
    cat("Carregando:", basename(pq_lin), "\n")
    lin_enriched <<- arrow::read_parquet(pq_lin)
    uni_enriched <<- arrow::read_parquet(pq_uni)
    cat("Dados carregados.\n")
  } else {
    stop(paste(
      "Arquivos parquet não encontrados em:\n", load_dir,
      "\n\nExecute 'run_module1 = TRUE' para gerar os dados primeiro.",
      "\nCertifique-se de que o MongoDB está ativo e o banco LOTUS importado."
    ))
  }
}

if (isTRUE(cfg$run_module3)) {
  cat("\n[3/4] Part III: estatísticas e figuras...\n")
  assign("cfg", cfg, envir = .GlobalEnv)
  source(file.path(script_dir, "Part III - Figure_Statistic_new_v2.R"))
  cat("Part III concluída.\n")
}

if (isTRUE(cfg$run_module4)) {
  cat("\n[4/4] Part IV: bioatividade e contexto global...\n")
  assign("cfg", cfg, envir = .GlobalEnv)

  if (!exists("uni_enriched")) {
    stop("Part IV requer 'uni_enriched'. Verifique o auto-loader.")
  }

  source(file.path(script_dir, "Part IV - actives_occurence_v3.R"))
  cat("Part IV concluída.\n")
}

cat("\nPipeline finalizado.\nResultados em:", normalizePath(cfg$out_dir_base, mustWork = FALSE), "\n")
