# ==========================================================
# SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps - Instalador de Pacotes R
# Instala todos os pacotes CRAN e Bioconductor necessários.
# Execute uma vez antes de rodar o pipeline.
# ==========================================================

REPOS <- "https://cloud.r-project.org/"

# ── Pacotes CRAN ───────────────────────────────────────────
CRAN_PKGS <- c(
  # Manipulação de dados
  "dplyr", "tidyr", "stringr", "stringi", "data.table",
  "readr", "readxl", "tibble", "forcats", "rlang",
  # I/O
  "writexl", "jsonlite", "arrow",
  # Banco de dados
  "mongolite",
  # HTTP / APIs
  "httr",
  # Progresso
  "progress",
  # Visualização
  "ggplot2", "ggrepel", "scales", "viridisLite",
  # Estatísticas
  "vegan",
  # Gráficos / imagens
  "png", "ragg", "grid",
  # Heatmaps (dependência do ComplexHeatmap)
  "circlize"
)

# ── Pacotes Bioconductor ───────────────────────────────────
BIOC_PKGS <- c(
  "ComplexHeatmap",
  "ChemmineR"
)

cat("=======================================================\n")
cat(" SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps - Instalacao de Pacotes R\n")
cat("=======================================================\n\n")

# ── 1. CRAN ────────────────────────────────────────────────
cat("--- Pacotes CRAN ---\n")
installed <- rownames(installed.packages())
missing_cran <- setdiff(CRAN_PKGS, installed)

if (length(missing_cran) == 0) {
  cat("OK - todos os pacotes CRAN ja instalados.\n\n")
} else {
  cat(sprintf("Instalando %d pacote(s): %s\n\n", length(missing_cran), paste(missing_cran, collapse = ", ")))
  for (pkg in missing_cran) {
    cat(sprintf("  -> %s ... ", pkg))
    tryCatch({
      install.packages(pkg, repos = REPOS, dependencies = TRUE, quiet = TRUE)
      cat("OK\n")
    }, error = function(e) {
      cat("ERRO:", conditionMessage(e), "\n")
    })
  }
}

# ── 2. Bioconductor ────────────────────────────────────────
cat("\n--- Pacotes Bioconductor ---\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Instalando BiocManager...\n")
  install.packages("BiocManager", repos = REPOS, quiet = TRUE)
}

installed <- rownames(installed.packages())
missing_bioc <- setdiff(BIOC_PKGS, installed)

if (length(missing_bioc) == 0) {
  cat("OK - todos os pacotes Bioconductor ja instalados.\n\n")
} else {
  cat(sprintf("Instalando %d pacote(s): %s\n\n", length(missing_bioc), paste(missing_bioc, collapse = ", ")))
  for (pkg in missing_bioc) {
    cat(sprintf("  -> %s ... ", pkg))
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      cat("OK\n")
    }, error = function(e) {
      cat("ERRO:", conditionMessage(e), "\n")
    })
  }
}

# ── 3. Verificação final ───────────────────────────────────
cat("\n=======================================================\n")
cat(" Verificacao Final\n")
cat("=======================================================\n")

all_pkgs <- c(CRAN_PKGS, BIOC_PKGS)
installed <- rownames(installed.packages())
ok  <- intersect(all_pkgs, installed)
nok <- setdiff(all_pkgs, installed)

for (p in ok)  cat(sprintf("  [OK] %s\n", p))
for (p in nok) cat(sprintf("  [FALTANDO] %s\n", p))

cat("\n")
if (length(nok) == 0) {
  cat("SUCESSO - todos os pacotes instalados!\n")
} else {
  cat(sprintf("ATENCAO - %d pacote(s) nao instalados: %s\n", length(nok), paste(nok, collapse = ", ")))
  if ("ChemmineR" %in% nok) {
    cat("NOTA: ChemmineR requer Rtools no Windows.\n")
    cat("      Download: https://cran.r-project.org/bin/windows/Rtools/\n")
  }
}
cat("=======================================================\n")
