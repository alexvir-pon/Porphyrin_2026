#!/usr/bin/env Rscript
# geo_riken_heme_extract.R
# -------------------------
# Maps RIKEN cDNA clone IDs (e.g. 0610005C13Rik) from mouse GEO data
# to mouse gene symbols, then to human orthologs, and extracts
# heme biosynthesis and iron metabolism genes.
#
# Usage:
#   Rscript geo_riken_heme_extract.R
#   Rscript geo_riken_heme_extract.R --geo GSE328851 --out results/
#
# Dependencies (install once):
#   install.packages(c("dplyr","tidyr","stringr","readr"))
#   BiocManager::install(c("GEOquery","biomaRt","org.Mm.eg.db",
#                           "org.Hs.eg.db","AnnotationDbi"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(biomaRt)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ── Configuration ─────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
GEO_ID  <- if ("--geo" %in% args) args[which(args=="--geo")+1] else "GSE328851"
OUT_DIR <- if ("--out" %in% args) args[which(args=="--out")+1] else "."
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Human heme/porphyrin gene list (for ortholog matching) ───────────────────
HEME_HUMAN <- toupper(c(
  # Core biosynthesis — all 8 enzymes
  "ALAS1","ALAS2","ALAD","HMBS","UROS","UROD","CPOX","PPOX","FECH",
  # Porphyrin transporters
  "ABCB6","ABCG2","FLVCR1","FLVCR2","SLC25A38",
  # Iron metabolism
  "FTH1","FTL","TFRC","TFR2","SLC40A1","HEPH","HEPHL1",
  # IRP system
  "ACO1","IREB2",
  # Fe-S cluster assembly (required for FECH)
  "FXN","ISCU","NFS1","GLRX5","BOLA1","BOLA3",
  # Heme oxygenase
  "HMOX1","HMOX2",
  # ALAS regulators
  "CLPX","HSPA9",
  # Heme chaperones / binding
  "HEBP1","HEBP2",
  # Heme receptor / signalling
  "PGRMC1","PGRMC2"
))

# Mouse orthologs of the above (same names, lowercase first letter convention
# doesn't matter — we match on symbol after uppercasing)
# biomaRt will give us the definitive list, but add known mouse names here
# so we can do a quick pre-filter
HEME_MOUSE_APPROX <- c(
  "Alas1","Alas2","Alad","Hmbs","Uros","Urod","Cpox","Ppox","Fech",
  "Abcb6","Abcg2","Flvcr1","Flvcr2","Slc25a38",
  "Fth1","Ftl1","Ftl2","Tfrc","Tfr2","Slc40a1","Heph","Hephl1",
  "Aco1","Ireb2",
  "Fxn","Iscu","Nfs1","Glrx5","Bola1","Bola3",
  "Hmox1","Hmox2",
  "Clpx","Hspa9",
  "Hebp1","Hebp2",
  "Pgrmc1","Pgrmc2"
)

# ── Step 1: Download GEO data ─────────────────────────────────────────────────
cat(sprintf("Downloading %s from GEO...\n", GEO_ID))
gse <- getGEO(GEO_ID, GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
if (is.list(gse)) gse <- gse[[1]]   # take first if multiple platforms

cat(sprintf("  %d samples, %d features\n", ncol(gse), nrow(gse)))
cat(sprintf("  Platform: %s\n", annotation(gse)))

# ── Step 2: Inspect feature data to find identifier column ───────────────────
fd <- fData(gse)
cat("\nFeature data columns:\n")
print(names(fd))
cat("\nFirst 5 rows of feature data:\n")
print(head(fd, 5))

# ── Step 3: Identify the RIKEN / gene ID column ───────────────────────────────
# Look for the column that contains RIKEN-style IDs or gene symbols
identify_id_column <- function(fd) {
  for (col in names(fd)) {
    vals <- as.character(fd[[col]])
    vals <- vals[!is.na(vals) & vals != ""]
    # RIKEN pattern: digits + letters + Rik
    n_rik <- sum(grepl("\\d.*Rik$", vals, ignore.case=TRUE))
    if (n_rik > nrow(fd) * 0.1) {
      cat(sprintf("  RIKEN IDs found in column: '%s' (%d/%d)\n",
                  col, n_rik, nrow(fd)))
      return(list(col=col, type="riken"))
    }
    # Ensembl mouse
    n_ens <- sum(grepl("^ENSMUSG", vals))
    if (n_ens > nrow(fd) * 0.1) {
      cat(sprintf("  Ensembl mouse IDs in column: '%s'\n", col))
      return(list(col=col, type="ensembl_mouse"))
    }
    # Gene symbol (mouse: starts with capital, rest lower)
    n_sym <- sum(grepl("^[A-Z][a-z]", vals))
    if (n_sym > nrow(fd) * 0.3) {
      cat(sprintf("  Gene symbols found in column: '%s'\n", col))
      return(list(col=col, type="symbol"))
    }
  }
  # Fallback: use rownames (probe IDs)
  rn <- rownames(fd)
  if (sum(grepl("\\d.*Rik$", rn)) > length(rn)*0.1) {
    cat("  RIKEN IDs found in rownames\n")
    return(list(col="__rownames__", type="riken"))
  }
  cat("  WARNING: could not identify ID column automatically\n")
  return(list(col=names(fd)[1], type="unknown"))
}

cat("\nDetecting identifier type:\n")
id_info <- identify_id_column(fd)

# Extract the ID vector
if (id_info$col == "__rownames__") {
  id_vec <- rownames(fd)
} else {
  id_vec <- as.character(fd[[id_info$col]])
}

# ── Step 4: Map RIKEN IDs → mouse gene symbols via biomaRt ───────────────────
cat(sprintf("\nMapping %d identifiers (%s) to gene symbols...\n",
            length(id_vec), id_info$type))

mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

if (id_info$type == "riken") {
  # RIKEN IDs are stored in biomaRt as 'mgi_id' or can be searched via
  # 'external_synonym' or 'wikigene_name'. Most reliable: mgi_symbol lookup
  # First try: use mgi_symbol filter (RIKEN IDs are often the mgi_symbol)
  bm_riken <- getBM(
    attributes = c("mgi_symbol","ensembl_gene_id","external_synonym"),
    filters    = "mgi_symbol",
    values     = unique(id_vec),
    mart       = mouse_mart
  )
  # Map: RIKEN id -> MGI symbol
  id_to_symbol <- bm_riken %>%
    select(query_id = mgi_symbol, mouse_symbol = mgi_symbol) %>%
    distinct()

  # Also try external_synonym in case some didn't map
  unmapped <- setdiff(unique(id_vec), id_to_symbol$query_id)
  if (length(unmapped) > 0 && length(unmapped) < 5000) {
    cat(sprintf("  Trying synonym search for %d unmapped IDs...\n",
                length(unmapped)))
    bm_syn <- getBM(
      attributes = c("external_synonym","mgi_symbol"),
      filters    = "external_synonym",
      values     = unmapped,
      mart       = mouse_mart
    )
    if (nrow(bm_syn) > 0) {
      syn_map <- bm_syn %>%
        select(query_id = external_synonym, mouse_symbol = mgi_symbol) %>%
        distinct()
      id_to_symbol <- bind_rows(id_to_symbol, syn_map)
    }
  }

} else if (id_info$type == "ensembl_mouse") {
  bm_ens <- getBM(
    attributes = c("ensembl_gene_id","mgi_symbol"),
    filters    = "ensembl_gene_id",
    values     = unique(str_extract(id_vec, "ENSMUSG[0-9]+")),
    mart       = mouse_mart
  )
  id_to_symbol <- bm_ens %>%
    rename(query_id = ensembl_gene_id, mouse_symbol = mgi_symbol)

} else if (id_info$type == "symbol") {
  # Already gene symbols — just use directly
  id_to_symbol <- data.frame(query_id    = unique(id_vec),
                              mouse_symbol = unique(id_vec),
                              stringsAsFactors = FALSE)
} else {
  # Unknown — attempt mgi_symbol lookup as best guess
  bm_unk <- getBM(
    attributes = c("mgi_symbol"),
    filters    = "mgi_symbol",
    values     = unique(id_vec),
    mart       = mouse_mart
  )
  id_to_symbol <- data.frame(query_id     = bm_unk$mgi_symbol,
                              mouse_symbol = bm_unk$mgi_symbol)
}

cat(sprintf("  Mapped %d unique identifiers to gene symbols\n",
            nrow(id_to_symbol)))

# ── Step 5: Map mouse symbols → human orthologs ───────────────────────────────
cat("\nMapping mouse → human orthologs...\n")

# Use biomaRt with human mart
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

orthologs <- getLDS(
  attributes  = c("mgi_symbol"),
  filters     = "mgi_symbol",
  values      = unique(id_to_symbol$mouse_symbol),
  mart        = mouse_mart,
  attributesL = c("hgnc_symbol"),
  martL       = human_mart
)
names(orthologs) <- c("mouse_symbol","human_symbol")
orthologs$human_symbol <- toupper(orthologs$human_symbol)

cat(sprintf("  Got %d mouse→human ortholog mappings\n", nrow(orthologs)))

# ── Step 6: Build complete annotation table ────────────────────────────────────
# Attach to feature data
fd$query_id <- id_vec

fd_annotated <- fd %>%
  left_join(id_to_symbol, by = "query_id", relationship="many-to-many") %>%
  left_join(orthologs,    by = "mouse_symbol", relationship="many-to-many") %>%
  mutate(
    is_heme = human_symbol %in% HEME_HUMAN |
              toupper(mouse_symbol) %in% toupper(HEME_MOUSE_APPROX)
  )

cat(sprintf("\nHeme pathway genes identified: %d probes / features\n",
            sum(fd_annotated$is_heme, na.rm=TRUE)))

# ── Step 7: Extract expression for heme genes ─────────────────────────────────
heme_idx <- which(fd_annotated$is_heme)

if (length(heme_idx) == 0) {
  cat("\nNo heme genes found. Printing all mapped gene symbols for inspection:\n")
  print(sort(unique(fd_annotated$mouse_symbol[!is.na(fd_annotated$mouse_symbol)])))
  stop("No heme pathway genes found — check mapping above")
}

expr_mat  <- exprs(gse)[heme_idx, , drop=FALSE]
ann_heme  <- fd_annotated[heme_idx, ]

# Build tidy long-format result
result_long <- as.data.frame(expr_mat) %>%
  mutate(
    probe_id     = rownames(expr_mat),
    mouse_symbol = ann_heme$mouse_symbol,
    human_symbol = ann_heme$human_symbol
  ) %>%
  pivot_longer(
    cols      = -c(probe_id, mouse_symbol, human_symbol),
    names_to  = "sample_id",
    values_to = "expression"
  )

# Add sample metadata
pd <- pData(gse) %>%
  rownames_to_column("sample_id") %>%
  select(sample_id, title,
         any_of(c("source_name_ch1","characteristics_ch1",
                  "treatment protocol","condition","genotype")))

result_long <- left_join(result_long, pd, by="sample_id")

# Wide format — average across probes mapping to same gene
result_wide <- result_long %>%
  group_by(human_symbol, mouse_symbol, sample_id) %>%
  summarise(expression = mean(expression, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=sample_id, values_from=expression) %>%
  arrange(human_symbol)

# ── Step 8: Print and save ────────────────────────────────────────────────────
cat("\n── Heme pathway genes — wide format ────────────────────────────────────\n")
print(as.data.frame(result_wide), digits=3, row.names=FALSE)

out_long <- file.path(OUT_DIR, paste0(GEO_ID, "_heme_long.csv"))
out_wide <- file.path(OUT_DIR, paste0(GEO_ID, "_heme_wide.csv"))
out_ann  <- file.path(OUT_DIR, paste0(GEO_ID, "_full_annotation.csv"))

write_csv(result_long, out_long)
write_csv(result_wide, out_wide)
write_csv(fd_annotated, out_ann)

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n", out_long, out_wide, out_ann))
cat(sprintf("\nGenes found: %s\n",
            paste(sort(unique(na.omit(result_wide$human_symbol))),
                  collapse=", ")))
