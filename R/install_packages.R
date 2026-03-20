repos <- "https://cloud.r-project.org/"

cran_pkgs <- c(
  "dplyr", "tidyr", "stringr", "stringi", "data.table",
  "readr", "readxl", "tibble", "forcats", "rlang",
  "writexl", "jsonlite", "arrow",
  "mongolite",
  "httr",
  "progress",
  "ggplot2", "ggrepel", "scales", "viridisLite",
  "vegan",
  "png", "ragg", "grid",
  "circlize"  # requerido pelo ComplexHeatmap
)

bioc_pkgs <- c(
  "ComplexHeatmap",
  "ChemmineR"
)

cat("Instalação de pacotes SPARK\n\n")

## CRAN
cat("Pacotes CRAN:\n")
installed    <- rownames(installed.packages())
missing_cran <- setdiff(cran_pkgs, installed)

if (length(missing_cran) == 0) {
  cat("Todos os pacotes CRAN já instalados.\n\n")
} else {
  cat(sprintf("Instalando %d pacote(s): %s\n\n", length(missing_cran), paste(missing_cran, collapse = ", ")))
  for (pkg in missing_cran) {
    cat(sprintf("  %s ... ", pkg))
    tryCatch({
      install.packages(pkg, repos = repos, dependencies = TRUE, quiet = TRUE)
      cat("OK\n")
    }, error = function(e) {
      cat("ERRO:", conditionMessage(e), "\n")
    })
  }
}

## Bioconductor
cat("\nPacotes Bioconductor:\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Instalando BiocManager...\n")
  install.packages("BiocManager", repos = repos, quiet = TRUE)
}

installed   <- rownames(installed.packages())
missing_bioc <- setdiff(bioc_pkgs, installed)

if (length(missing_bioc) == 0) {
  cat("Todos os pacotes Bioconductor já instalados.\n\n")
} else {
  cat(sprintf("Instalando %d pacote(s): %s\n\n", length(missing_bioc), paste(missing_bioc, collapse = ", ")))
  for (pkg in missing_bioc) {
    cat(sprintf("  %s ... ", pkg))
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      cat("OK\n")
    }, error = function(e) {
      cat("ERRO:", conditionMessage(e), "\n")
    })
  }
}

## verificação final
cat("\nVerificação final:\n")
all_pkgs  <- c(cran_pkgs, bioc_pkgs)
installed <- rownames(installed.packages())
ok  <- intersect(all_pkgs, installed)
nok <- setdiff(all_pkgs, installed)

for (p in ok)  cat(sprintf("  [OK] %s\n", p))
for (p in nok) cat(sprintf("  [FALTANDO] %s\n", p))

cat("\n")
if (length(nok) == 0) {
  cat("Todos os pacotes instalados com sucesso.\n")
} else {
  cat(sprintf("Atenção: %d pacote(s) não instalados: %s\n", length(nok), paste(nok, collapse = ", ")))
  if ("ChemmineR" %in% nok) {
    cat("ChemmineR requer Rtools no Windows: https://cran.r-project.org/bin/windows/Rtools/\n")
  }
}
