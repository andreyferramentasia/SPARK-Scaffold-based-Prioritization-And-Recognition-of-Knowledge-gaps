# ==========================================================
# MAIN PIPELINE - SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps
# Controller for:
#   Part I   - Extraction & Normalization
#   Part II  - Structure Catalog (Optional)
#   Part III - Statistics & Figures
#   Part IV  - Global Context & Rarity
#
# Focus: Family->Genus or Genus->Species analysis levels.
#
# Usage (interativo):
#   Edite a seção CONFIG abaixo e execute.
#
# Usage (via Streamlit / linha de comando):
#   Rscript Main_Pipeline_new_v2.R /caminho/para/config.json
# ==========================================================

# ── Null-coalesce operator ─────────────────────────────────────────────────────
`%||%` <- function(a, b) if (is.null(a)) b else a

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
})

# ── Determine script directory (funciona em Rscript e source()) ───────────────
script_dir <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  tryCatch({
    dirname(normalizePath(commandArgs(trailingOnly = FALSE) |>
      (\(a) a[grepl("--file=", a)])() |>
      (\(a) sub("--file=", "", a))()))
  }, error = function(e2) getwd())
})

# ── CONFIG ────────────────────────────────────────────────────────────────────
# Modo 1: JSON passado pela linha de comando (ex: via Streamlit)
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1 && file.exists(args[1])) args[1] else NULL

if (!is.null(config_file)) {
  cat("[Config] Lendo configuração de:", config_file, "\n")
  cfg_json <- jsonlite::read_json(config_file, simplifyVector = TRUE)

  cfg_mode   <- cfg_json$taxon_mode   %||% "family"
  cfg_values <- cfg_json$taxon_values %||% c("Phyllanthaceae")
  cfg_level  <- cfg_json$analysis_tax_level %||% "genus"

  cfg <- list(
    taxon_mode                    = cfg_mode,
    taxon_values                  = cfg_values,
    analysis_tax_level            = cfg_level,
    run_module1                   = isTRUE(cfg_json$run_module1),
    run_module2                   = isTRUE(cfg_json$run_module2),
    run_module3                   = isTRUE(cfg_json$run_module3),
    run_module4                   = isTRUE(cfg_json$run_module4),
    analysis_top_taxa             = as.integer(cfg_json$analysis_top_taxa             %||% 40L),
    analysis_min_compounds_per_taxon = as.integer(cfg_json$analysis_min_compounds_per_taxon %||% 10L),
    out_dir_base                  = cfg_json$out_dir_base   %||% file.path(script_dir, "..", "results"),
    run_tag_date                  = Sys.Date(),
    mongo_url                     = cfg_json$mongo_url      %||% "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000",
    db_name                       = cfg_json$db_name        %||% "lotus",
    coll_name                     = cfg_json$coll_name      %||% "lotusUniqueNaturalProduct",
    use_WFO_normalization         = isTRUE(cfg_json$use_WFO_normalization),
    wfo_csv_path                  = cfg_json$wfo_csv_path   %||% file.path(script_dir, "..", "DBs", "classification.tsv"),
    export_excel                  = isTRUE(cfg_json$export_excel  %||% TRUE),
    export_parquet                = TRUE,
    verbose                       = isTRUE(cfg_json$verbose %||% TRUE)
  )

} else {
  # Modo 2: Configuração manual (edite aqui para uso interativo)

  # --- CENÁRIO A: Uma FAMÍLIA, comparando GÊNEROS ---
  # cfg_mode   <- "family"
  # cfg_values <- c("Lauraceae")
  # cfg_level  <- "genus"

  # --- CENÁRIO B: Um GÊNERO, comparando ESPÉCIES ---
  cfg_mode   <- "family"
  cfg_values <- c("Phyllanthaceae")
  cfg_level  <- "genus"

  cfg <- list(
    ## --- GRUPO TAXONÔMICO ALVO ─────────────────────
    taxon_mode   = cfg_mode,
    taxon_values = cfg_values,

    ## --- NÍVEL DE ANÁLISE ──────────────────────────
    analysis_tax_level = cfg_level,

    ## --- MÓDULOS ───────────────────────────────────
    run_module1 = FALSE,  # Part I:   Extração MongoDB (demorado)
    run_module2 = FALSE,  # Part II:  Catálogo de estruturas (opcional)
    run_module3 = FALSE,  # Part III: Estatísticas & Figuras
    run_module4 = TRUE,   # Part IV:  Contexto global & raridade

    ## --- FILTROS PARA VISUALIZAÇÕES ────────────────
    analysis_top_taxa                = 40L,
    analysis_min_compounds_per_taxon = 10L,

    ## --- SAÍDA / CAMINHOS ──────────────────────────
    out_dir_base = file.path(script_dir, "..", "results"),
    run_tag_date = Sys.Date(),

    ## --- CONEXÃO LOTUS (MongoDB) ───────────────────
    mongo_url = "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000",
    db_name   = "lotus",
    coll_name = "lotusUniqueNaturalProduct",

    ## --- EXTRAS ────────────────────────────────────
    use_WFO_normalization = TRUE,
    wfo_csv_path          = file.path(script_dir, "..", "DBs", "classification.tsv"),
    export_excel          = TRUE,
    export_parquet        = TRUE,
    verbose               = TRUE
  )
}

# ── EXECUÇÃO AUTOMATIZADA ─────────────────────────────────────────────────────
# Normalizar entradas
cfg$taxon_mode         <- tolower(cfg$taxon_mode)
cfg$analysis_tax_level <- tolower(cfg$analysis_tax_level)

cat("=====================================================\n")
cat(" SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps INICIADO\n")
cat(" Modo: ", cfg$taxon_mode, " -> Analisando por: ", cfg$analysis_tax_level, "\n")
cat(" Alvo: ", paste(cfg$taxon_values, collapse = ", "), "\n")
cat(" Diretório de saída:", normalizePath(cfg$out_dir_base, mustWork = FALSE), "\n")
cat("=====================================================\n")

# Garantir que o diretório de saída existe
dir.create(cfg$out_dir_base, showWarnings = FALSE, recursive = TRUE)

# Helper de tag (duplicado aqui para o Auto-Loader funcionar sem Part I)
safe_tag_fn <- function(mode, values, run_date, analysis_level) {
  tag <- paste(mode, analysis_level, paste(values, collapse = "-"),
               format(run_date, "%Y%m%d"), sep = "_")
  gsub("[^A-Za-z0-9._-]+", "_", tag)
}

tag_base_load <- safe_tag_fn(cfg$taxon_mode, cfg$taxon_values,
                              cfg$run_tag_date, cfg$analysis_tax_level)
load_dir      <- file.path(cfg$out_dir_base, paste0("lotus_", tag_base_load))

# Variáveis globais usadas pelos scripts filhos
OUT_DIR      <<- load_dir
base_tag     <<- paste0("lotus_", tag_base_load)
TAXON_MODE   <<- cfg$taxon_mode
TAXON_VALUES <<- cfg$taxon_values

# ── PART I: EXTRAÇÃO ─────────────────────────────────────────────────────────
if (isTRUE(cfg$run_module1)) {
  cat("\n[1/4] Executando Part I (Extração & Limpeza)...\n")
  assign("cfg", cfg, envir = .GlobalEnv)
  source(file.path(script_dir, "Part I - Extraction_new_v2.R"))
  cat("✔ Part I concluída.\n")
}

# ── PART II: ESTRUTURAS ───────────────────────────────────────────────────────
if (isTRUE(cfg$run_module2)) {
  cat("\n[2/4] Executando Part II (Catálogo de Estruturas)...\n")
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    cat("⚠ Part II ignorada: pacote 'ChemmineR' não instalado.\n")
    cat("  Para instalar: BiocManager::install('ChemmineR')\n")
  } else {
    assign("cfg", cfg, envir = .GlobalEnv)
    tryCatch(
      { source(file.path(script_dir, "Part II - Structures_new_v2.R")); cat("✔ Part II concluída.\n") },
      error = function(e) cat("⚠ Part II falhou (continuando):", conditionMessage(e), "\n")
    )
  }
}

# ── AUTO-LOADER (carrega Parquet se objetos ausentes) ─────────────────────────
if ((isTRUE(cfg$run_module3) || isTRUE(cfg$run_module4)) &&
    (!exists("lin_enriched") || !exists("uni_enriched"))) {

  cat("\n[Auto-Load] Objetos não encontrados na memória. Buscando Parquet...\n")

  pq_lin <- file.path(load_dir, paste0("lotus_", tag_base_load, "_lin_enriched.parquet"))
  pq_uni <- file.path(load_dir, paste0("lotus_", tag_base_load, "_uni_enriched.parquet"))

  if (file.exists(pq_lin) && file.exists(pq_uni)) {
    cat("  -> Carregando:", basename(pq_lin), "\n")
    lin_enriched <<- arrow::read_parquet(pq_lin)
    uni_enriched <<- arrow::read_parquet(pq_uni)
    cat("✔ Dados carregados! Prosseguindo para análise.\n")
  } else {
    stop(paste(
      "ERRO CRÍTICO: Arquivos Parquet não encontrados em:\n", load_dir,
      "\n\nExecute 'run_module1 = TRUE' para gerar os dados primeiro.",
      "\nCertifique-se de que o MongoDB está ativo e o banco LOTUS importado."
    ))
  }
}

# ── PART III: ESTATÍSTICAS & FIGURAS ──────────────────────────────────────────
if (isTRUE(cfg$run_module3)) {
  cat("\n[3/4] Executando Part III (Estatísticas & Figuras)...\n")
  assign("cfg", cfg, envir = .GlobalEnv)
  source(file.path(script_dir, "Part III - Figure_Statistic_new_v2.R"))
  cat("✔ Part III concluída.\n")
}

# ── PART IV: CONTEXTO GLOBAL ──────────────────────────────────────────────────
if (isTRUE(cfg$run_module4)) {
  cat("\n[4/4] Executando Part IV (Bioatividade & Contexto Global)...\n")
  assign("cfg", cfg, envir = .GlobalEnv)

  if (!exists("uni_enriched")) {
    stop("Part IV requer 'uni_enriched'. Verifique o Auto-Loader.")
  }

  # CORRIGIDO: era 'v2', arquivo correto é 'v3'
  source(file.path(script_dir, "Part IV - actives_occurence_v3.R"))
  cat("✔ Part IV concluída.\n")
}

cat("\nPIPELINE FINALIZADO.\nResultados em:", normalizePath(cfg$out_dir_base, mustWork = FALSE), "\n")
