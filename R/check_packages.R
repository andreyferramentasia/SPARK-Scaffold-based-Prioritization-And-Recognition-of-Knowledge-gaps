# Verifica quais pacotes estao instalados e imprime JSON.
# Usado pela interface Streamlit.
pkgs <- c(
  "dplyr","tidyr","stringr","stringi","data.table","readr","readxl",
  "tibble","forcats","rlang","writexl","jsonlite","arrow","mongolite",
  "httr","progress","ggplot2","ggrepel","scales","viridisLite",
  "vegan","png","ragg","grid","circlize","ComplexHeatmap","ChemmineR"
)
inst <- rownames(installed.packages())
result <- setNames(as.list(pkgs %in% inst), pkgs)
cat(jsonlite::toJSON(result, auto_unbox = TRUE))
