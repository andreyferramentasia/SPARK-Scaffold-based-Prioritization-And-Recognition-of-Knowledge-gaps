# ==========================================================
# MODULE 2 - STRUCTURE CATALOG GENERATION (Part II)
#
# Purpose:
#   1. Generate high-resolution structure images (PNG) for each unique compound.
#   2. Compile a paginated PDF structure catalog (Grid Layout).
#   3. Handle stereochemistry rendering fallbacks.
#   4. Report rendering failures for quality control.
#
# Inputs:
#   - 'cfg' list (global configuration).
#   - 'uni_enriched' (unique compounds table from Part I).
#
# Outputs:
#   - PNG images in <OUT_DIR>/png/
#   - PDF Catalog: <OUT_DIR>/<base_tag>_catalog.pdf
#   - Render reports (TSV) for audit.
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(png)
  library(ragg)
  library(ChemmineR)
})

`%||%`  <- function(a,b) if (is.null(a)) b else a
`%|||%` <- function(a,b) {
  if (is.null(a)) return(b)
  if (length(a) == 1 && (is.na(a) || (is.character(a) && !nzchar(a)))) return(b)
  a
}

# ----------------------------------------------------------
# 0. CONFIGURATION & SANITY CHECKS
# ----------------------------------------------------------
# Read settings from global 'cfg' or apply defaults
png_subdir            <- cfg$png_subdir            %||% "png"
pdf_filename_suffix   <- cfg$pdf_filename_suffix   %||% "_catalog_2x3_by_mass"
pdf_file_ext          <- cfg$pdf_file_ext          %||% ".pdf"

img_width_px          <- cfg$img_width_px          %||% 1200L
img_height_px         <- cfg$img_height_px         %||% 1200L
img_res_dpi           <- cfg$img_res_dpi           %||% 150L

gen_pngs              <- cfg$gen_pngs              %||% TRUE

pdf_page_width_in     <- cfg$pdf_page_width_in     %||% 8.27
pdf_page_height_in    <- cfg$pdf_page_height_in    %||% 11.69

pdf_n_cols            <- cfg$pdf_n_cols            %||% 2L
pdf_n_rows            <- cfg$pdf_n_rows            %||% 3L

order_desc            <- cfg$order_desc            %||% FALSE
max_structures        <- cfg$max_structures        %||% NA_integer_

allow_stereo_relax    <- cfg$allow_stereo_relax    %||% TRUE
export_render_reports <- cfg$export_render_reports %||% TRUE

verbose_local         <- cfg$verbose               %||% TRUE

# Ensure required objects exist (loaded from Part I)
stopifnot(
  exists("uni_enriched"),
  exists("OUT_DIR"),
  exists("base_tag")
)

# ----------------------------------------------------------
# 1. PATHS AND HELPER FUNCTIONS
# ----------------------------------------------------------
PNG_DIR  <- file.path(OUT_DIR, png_subdir)
dir.create(PNG_DIR, showWarnings = FALSE, recursive = TRUE)

safe_file <- function(base, ext) file.path(OUT_DIR, paste0(base, ext))
PDF_FILE  <- safe_file(paste0(base_tag, pdf_filename_suffix), pdf_file_ext)

sanitize_fn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Helper: Suppress chatter from ChemmineR/OpenBabel
quiet_ob <- function(expr) {
  tf_out <- tempfile()
  tf_msg <- tempfile()
  con_out <- file(tf_out, open = "wt")
  con_msg <- file(tf_msg, open = "wt")
  on.exit({
    try(sink(NULL), silent = TRUE)
    try(sink(NULL, type = "message"), silent = TRUE)
    try(close(con_out), silent = TRUE)
    try(close(con_msg), silent = TRUE)
  }, add = TRUE)
  sink(con_out)
  sink(con_msg, type = "message")
  suppressWarnings(suppressMessages(force(expr)))
  invisible(NULL)
}

# Helper: Generate SDF from SMILES with stereochemistry fallback
# 1. Try raw SMILES.
# 2. If valid but fails rendering, and fallback is enabled, strip stereochem markers (@).
safe_sdf_from_smiles <- function(sm) {
  if (is.null(sm) || !nzchar(sm)) return(NULL)
  
  # Attempt 1: Direct conversion
  sdf <- try(quiet_ob(ChemmineR::smiles2sdf(sm)), silent = TRUE)
  if (!(inherits(sdf, "try-error")) && length(sdf) >= 1) {
    return(sdf[[1]])
  }
  
  # Attempt 2: Relaxed stereochemistry
  if (isTRUE(allow_stereo_relax)) {
    sm2 <- gsub("@@?", "", sm)
    sdf <- try(quiet_ob(ChemmineR::smiles2sdf(sm2)), silent = TRUE)
    if (!(inherits(sdf, "try-error")) && length(sdf) >= 1) {
      return(sdf[[1]])
    }
  }
  
  NULL
}

# Helper: Ensure column exists using candidates
ensure_col <- function(df, target, candidates) {
  if (!target %in% names(df)) {
    cand <- intersect(candidates, names(df))
    if (length(cand)) {
      df[[target]] <- df[[cand[1]]]
    }
  }
  df
}

# Helper: Auto-fit title text size for catalog panels
fit_title <- function(txt, max_cex = 0.95, min_cex = 0.70, max_width_frac = 0.96, wrap_width = 56) {
  if (is.na(txt) || !nzchar(txt)) {
    return(list(lines = character(), cex = max_cex))
  }
  lines <- strwrap(txt, width = wrap_width)
  cex_try <- max_cex
  
  too_wide <- function(ls, cexv) any(strwidth(ls, cex = cexv) > max_width_frac)
  
  it <- 0
  while (too_wide(lines, cex_try) && cex_try > min_cex && it < 25) {
    cex_try <- cex_try - 0.03
    it <- it + 1
  }
  
  list(lines = lines, cex = max(cex_try, min_cex))
}

# Helper: Draw a single molecule cell in the PDF grid
draw_cell <- function(sm, iupac, ik, formula, mw,
                      PNG_DIR_local = PNG_DIR,
                      img_w = img_width_px,
                      img_h = img_height_px,
                      img_res = img_res_dpi) {
  
  par(mar = c(0.8,0.8,0.6,0.8))
  plot.new()
  usr <- par("usr")
  x0 <- usr[1]; x1 <- usr[2]; y0 <- usr[3]; y1 <- usr[4]
  y_img_top <- y1 - 0.02
  y_img_bot <- y0 + 0.35*(y1 - y0)
  
  # 1) Try loading pre-rendered PNG from disk (fastest)
  png_pre <- file.path(PNG_DIR_local, paste0("STRUCT_", sanitize_fn(ik), ".png"))
  img <- NULL
  if (file.exists(png_pre) &&
      is.finite(file.info(png_pre)$size) &&
      file.info(png_pre)$size > 0) {
    img <- tryCatch(png::readPNG(png_pre), error = function(e) NULL)
  }
  
  # 2) If missing, render on the fly
  if (is.null(img)) {
    sdf <- safe_sdf_from_smiles(sm)
    if (!is.null(sdf)) {
      tf <- tempfile(fileext = ".png")
      ragg::agg_png(
        tf, width  = img_w, height = img_h,
        units  = "px", res    = img_res, background = "white"
      )
      par(mar = c(1,1,1,1))
      plot.new()
      try(quiet_ob(ChemmineR::plotStruc(sdf)), silent = TRUE)
      dev.off()
      img <- tryCatch(png::readPNG(tf), error = function(e) NULL)
      unlink(tf, force = TRUE)
    }
  }
  
  # 3) Draw image or error placeholder
  if (!is.null(img)) {
    rasterImage(
      img,
      x0+0.02*(x1-x0), y_img_bot,
      x1-0.02*(x1-x0), y_img_top
    )
  } else {
    text((x0+x1)/2, (y_img_bot+y_img_top)/2, "Failed to render SMILES", cex = 0.9)
  }
  
  # 4) Draw Metadata Labels
  title_txt <- if (!is.na(iupac) && nzchar(iupac)) iupac else ik
  fit <- fit_title(title_txt)
  
  y_text <- y_img_bot - 0.05*(y1-y0)
  step   <- 0.09*(y1-y0)
  
  if (length(fit$lines)) {
    for (k in seq_along(fit$lines)) {
      text((x0+x1)/2, y_text - (k-1)*step, fit$lines[k], cex = fit$cex, font = 2)
    }
    y_text <- y_text - length(fit$lines)*step
  }
  
  info <- paste0(
    if (!is.na(formula) && nzchar(formula)) paste0("Formula: ", formula) else "",
    if (suppressWarnings(is.finite(as.numeric(mw)))) sprintf("   |   Mass (calc): %.3f", as.numeric(mw)) else ""
  )
  
  if (nzchar(gsub("\\s+|\\|", "", info))) {
    text((x0+x1)/2, y_text - 0.06*(y1-y0), info, cex = 0.80)
  }
  
  text((x0+x1)/2, y_text - 0.14*(y1-y0), ik, cex = 0.78)
}


# ----------------------------------------------------------
# 2. PREPARE DATASET (Unique Structures)
# ----------------------------------------------------------
if (verbose_local) {
  cat("[Module 2] Preparing unique structure list from uni_enriched...\n")
}

ue <- uni_enriched

# Standardize column names
ue <- ensure_col(ue, "smiles",            c("smiles","smiles.x","smiles.y","sugar_free_smiles","smiles2D","smiles2d"))
ue <- ensure_col(ue, "iupac_name",        c("iupac_name","iupac_name.x","iupac_name.y","traditional_name"))
ue <- ensure_col(ue, "molecular_formula", c("molecular_formula","molecular_formula.x","molecular_formula.y"))
ue <- ensure_col(ue, "molecular_weight",  c("molecular_weight","molecular_weight.x","molecular_weight.y"))

ue$molecular_weight <- suppressWarnings(as.numeric(ue$molecular_weight))

df_base <- ue %>%
  dplyr::select(inchikey, smiles, iupac_name, molecular_formula, molecular_weight) %>%
  dplyr::filter(!is.na(inchikey), nzchar(inchikey), !is.na(smiles), nzchar(smiles)) %>%
  dplyr::distinct(inchikey, .keep_all = TRUE) %>%
  {
    if (isTRUE(order_desc)) dplyr::arrange(., dplyr::desc(molecular_weight), inchikey)
    else                    dplyr::arrange(., molecular_weight, inchikey)
  }

if (!is.na(max_structures)) {
  df_base <- head(df_base, max_structures)
}

if (verbose_local) {
  cat("[Module 2] Unique structures to process: ", nrow(df_base), "\n")
}


# ----------------------------------------------------------
# 3. GENERATE PNG FILES (Per-Structure)
# ----------------------------------------------------------
render_fallback <- list()
render_fail     <- list()

if (isTRUE(gen_pngs)) {
  
  if (verbose_local) {
    cat("[Module 2] Generating PNG images in: ", normalizePath(PNG_DIR, winslash="/"), "\n", sep = "")
  }
  
  ok   <- 0L
  fail <- 0L
  
  for (i in seq_len(nrow(df_base))) {
    
    ik <- df_base$inchikey[i]
    sm <- df_base$smiles[i]
    nm <- df_base$iupac_name[i]
    fm <- df_base$molecular_formula[i]
    mw <- df_base$molecular_weight[i]
    
    fn <- file.path(PNG_DIR, paste0("STRUCT_", sanitize_fn(ik), ".png"))
    
    # Skip if valid PNG already exists
    if (file.exists(fn) && is.finite(file.info(fn)$size) && file.info(fn)$size > 0) {
      ok <- ok + 1L
      next
    }
    
    # Initialize PNG device
    ragg::agg_png(
      fn, width = img_width_px, height = img_height_px,
      units = "px", res = img_res_dpi, background = "white"
    )
    
    # 2-Panel Layout: Structure (Top), Metadata (Bottom)
    layout(matrix(c(1,2), nrow = 2), heights = c(0.70, 0.30))
    par(mar = c(1,1,1,1))
    
    # Try standard rendering first
    sdf1 <- try(ChemmineR::smiles2sdf(sm), silent = TRUE)
    used_fallback <- FALSE
    
    if (!(inherits(sdf1, "try-error")) && length(sdf1) >= 1) {
      quiet_ob(ChemmineR::plotStruc(sdf1[[1]]))
    } else {
      # Fallback: Relax stereochemistry
      sdf2 <- safe_sdf_from_smiles(sm)
      if (!is.null(sdf2)) {
        quiet_ob(ChemmineR::plotStruc(sdf2))
        used_fallback <- TRUE
      } else {
        plot.new()
        text(0.5, 0.6, "Failed to render SMILES.", cex = 1.2)
        text(0.5, 0.48, substr(sm, 1, 140), cex = 0.9)
      }
    }
    
    # Bottom Panel: Labels
    par(mar = c(2,2,1,2))
    plot.new()
    
    title_txt <- if (!is.na(nm) && nzchar(nm)) nm else ik
    fit <- fit_title(title_txt)
    
    if (length(fit$lines)) {
      ys <- 0.78 - (0:(length(fit$lines)-1)) * 0.18
      for (k in seq_along(fit$lines)) {
        text(0.5, ys[k], fit$lines[k], cex = fit$cex, font = 2)
      }
      
      info <- paste0(
        if (!is.na(fm) && nzchar(fm)) fm else "",
        if (is.finite(mw)) sprintf(" • MW: %.3f", mw) else ""
      )
      if (nzchar(info)) {
        text(0.5, min(ys) - 0.20, info, cex = 0.92)
      }
      text(0.5, min(ys) - 0.38, ik, cex = 0.92)
    } else {
      text(0.5, 0.70, ik, cex = 0.95, font = 2)
    }
    
    dev.off()
    
    # Validate Output
    if (file.exists(fn) && file.info(fn)$size > 0) {
      ok <- ok + 1L
      if (used_fallback) render_fallback[[length(render_fallback)+1]] <- ik
    } else {
      fail <- fail + 1L
      render_fail[[length(render_fail)+1]] <- ik
    }
  }
  
  if (verbose_local) {
    cat("[Module 2] PNG Generation Complete: ", ok, " Success; ", fail, " Failed.\n", sep = "")
  }
}


# ----------------------------------------------------------
# 4. GENERATE PDF CATALOG (Grid Layout)
# ----------------------------------------------------------
if (verbose_local) {
  cat("[Module 2] Compiling PDF catalog: ", normalizePath(PDF_FILE, winslash = "/"), "\n", sep = "")
}

per_page <- as.integer(pdf_n_cols * pdf_n_rows)

grDevices::pdf(
  PDF_FILE,
  width   = pdf_page_width_in,
  height  = pdf_page_height_in,
  onefile = TRUE,
  paper   = "special"
)

n        <- nrow(df_base)
n_pages  <- ceiling(n / per_page)
idx      <- 1L

for (p in seq_len(n_pages)) {
  
  layout(
    matrix(seq_len(pdf_n_rows * pdf_n_cols),
           nrow = pdf_n_rows,
           byrow = TRUE)
  )
  par(oma = c(0.2,0.2,0.2,0.2))
  
  for (cell in seq_len(per_page)) {
    if (idx > n) {
      plot.new()
      next
    }
    
    r <- df_base[idx, ]
    
    draw_cell(
      sm      = r$smiles,
      iupac   = r$iupac_name  %|||% "",
      ik      = r$inchikey,
      formula = r$molecular_formula %|||% "",
      mw      = r$molecular_weight
    )
    
    idx <- idx + 1L
  }
}

grDevices::dev.off()
if (verbose_local) cat("[Module 2] PDF saved successfully.\n")


# ----------------------------------------------------------
# 5. EXPORT RENDER REPORTS
# ----------------------------------------------------------
if (isTRUE(export_render_reports)) {
  
  if (length(render_fallback)) {
    data.table::fwrite(
      data.frame(inchikey = unlist(render_fallback)),
      safe_file(paste0(base_tag, "_render_fallback"), ".tsv"),
      sep = "\t"
    )
  }
  
  if (length(render_fail)) {
    data.table::fwrite(
      data.frame(inchikey = unlist(render_fail)),
      safe_file(paste0(base_tag, "_render_fail"), ".tsv"),
      sep = "\t"
    )
  }
  
  if (verbose_local) {
    cat("[Module 2] Report: ", length(render_fallback), " stereo-relaxed | ", length(render_fail), " failed to render.\n", sep = "")
  }
}

if (verbose_local) cat("[Module 2] Finished.\n")