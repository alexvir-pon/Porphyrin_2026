#!/usr/bin/env Rscript
# extract_gepia_values.R
# Reads all GEPIA correlation PDFs in a directory,
# extracts Spearman R and p-value, writes a tidy CSV.
#
# Compatible with run_gepia_pdf_batch.py filename format:
#   ph{phase}__{analysis}__{gene1}_vs_{gene2}__{dataset}__{method}.pdf
#
# Usage:
#   Rscript extract_gepia_values.R <pdf_dir> <output_csv>

suppressPackageStartupMessages({
  library(pdftools)
  library(stringr)
  library(dplyr)
})

# Null-coalescing helper (R < 4.4 compatible)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a) > 0) a else b

# ── Arguments ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript extract_gepia_values.R <pdf_dir> <output_csv>\n")
  quit(status = 1)
}

pdf_dir    <- args[1]
output_csv <- args[2]

if (!dir.exists(pdf_dir)) stop(paste("Directory not found:", pdf_dir))

# ── Find PDFs ────────────────────────────────────────────────────────────────
pdfs <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE)
cat(sprintf("Found %d PDFs in %s\n", length(pdfs), pdf_dir))
if (length(pdfs) == 0) stop("No PDF files found.")

# ── Extraction helpers ────────────────────────────────────────────────────────
extract_one <- function(pdf_path) {
  text <- tryCatch(
    paste(pdf_text(pdf_path), collapse = " "),
    error = function(e) NA_character_
  )
  if (is.na(text)) return(list(R = NA_real_, p_value = NA_character_, status = "read_error"))

  # Normalise unicode minus to hyphen
  text <- gsub("\u2212", "-", text)

  r_match <- str_match(text, "R\\s*=\\s*([+-]?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?)")
  r_val   <- if (!is.na(r_match[1, 1])) as.numeric(r_match[1, 2]) else NA_real_

  p_match <- str_match(text, "p[-\u2212]value\\s*=\\s*([^\\s]+)")
  p_str   <- if (!is.na(p_match[1, 1])) p_match[1, 2] else NA_character_

  list(R      = r_val,
       p_value = p_str,
       status  = if (!is.na(r_val)) "ok" else "extraction_failed")
}

parse_filename <- function(fname) {
  stem  <- tools::file_path_sans_ext(basename(fname))
  parts <- str_split(stem, "__")[[1]]

  # Expected format: ph{phase}__{analysis}__{gene1}_vs_{gene2}__{dataset}__{method}
  # parts[1] = phase
  # parts[2] = analysis
  # parts[3] = gene1_vs_gene2   <- _vs_ is here
  # parts[4] = dataset
  # parts[5] = method  (optional, not used)

  if (length(parts) < 4 || !str_detect(parts[3], "_vs_"))
    return(list(phase = NA_character_, analysis = NA_character_,
                dataset = NA_character_, gene1 = NA_character_, gene2 = NA_character_))

  genes <- str_split(parts[3], "_vs_")[[1]]

  list(phase    = parts[1],
       analysis = parts[2],
       gene1    = genes[1],
       gene2    = if (length(genes) >= 2) genes[2] else NA_character_,
       dataset  = parts[4])
}

# ── Main loop ─────────────────────────────────────────────────────────────────
cat("Extracting R and p-values...\n")
results <- vector("list", length(pdfs))

for (i in seq_along(pdfs)) {
  meta <- parse_filename(pdfs[i])
  vals <- extract_one(pdfs[i])

  results[[i]] <- data.frame(
    phase    = meta$phase    %||% NA_character_,
    analysis = meta$analysis %||% NA_character_,
    gene1    = meta$gene1    %||% NA_character_,
    gene2    = meta$gene2    %||% NA_character_,
    dataset  = meta$dataset  %||% NA_character_,
    R        = vals$R,
    p_value  = vals$p_value  %||% NA_character_,
    status   = vals$status,
    path     = pdfs[i],
    stringsAsFactors = FALSE
  )

  if (vals$status != "ok") cat(sprintf("  FAIL: %s\n", basename(pdfs[i])))
}

# ── Write output ──────────────────────────────────────────────────────────────
df     <- bind_rows(results)
ok_n   <- sum(df$status == "ok",  na.rm = TRUE)
fail_n <- sum(df$status != "ok",  na.rm = TRUE)

cat(sprintf("\nExtracted: %d OK  |  %d failed\n", ok_n, fail_n))
cat(sprintf("Phases: %s\n", paste(sort(unique(df$phase[!is.na(df$phase)])), collapse = ", ")))

write.csv(df, output_csv, row.names = FALSE, quote = TRUE)
cat(sprintf("CSV written to: %s\n", output_csv))

cat("\nCounts per phase:\n")
df %>%
  filter(status == "ok") %>%
  count(phase) %>%
  arrange(phase) %>%
  as.data.frame() %>%
  print()
