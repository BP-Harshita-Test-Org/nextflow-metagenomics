#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Maaslin2)
  library(readr)
  library(dplyr)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--abundance",       type="character", help="Merged abundance matrix TSV"),
  make_option("--metadata",        type="character", help="Sample metadata CSV"),
  make_option("--outdir",          type="character", help="Output directory"),
  make_option("--fixed_effects",   type="character", default="diagnosis,age,sex",
              help="Comma-separated fixed effects"),
  make_option("--reference_level", type="character", default="neurotypical",
              help="Reference group level")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ───────────────────────────────────────────────────────────────
# Merged matrix: first column = feature, remaining columns = samples
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)
metadata  <- read_csv(opts$metadata, show_col_types=FALSE)

sample_ids <- colnames(abundance)[-1]

cat("Loaded abundance matrix:", nrow(abundance), "features x", length(sample_ids), "samples\n")
cat("Loaded metadata:", nrow(metadata), "samples\n")

fixed_effects <- strsplit(opts$fixed_effects, ",")[[1]]
cat("Fixed effects:", paste(fixed_effects, collapse=", "), "\n")
cat("Reference level:", opts$reference_level, "\n")

# ── Build MaAsLin2 input ──────────────────────────────────────────────────
# MaAsLin2 expects: rows = samples, columns = features
abund_mat <- abundance %>% select(-feature) %>% as.data.frame()
rownames(abund_mat) <- abundance$feature
abund_mat <- t(abund_mat) %>% as.data.frame()  # samples x features

metadata_df <- metadata %>%
  filter(sample_id %in% rownames(abund_mat)) %>%
  tibble::column_to_rownames("sample_id")

# ── Run MaAsLin2 ──────────────────────────────────────────────────────────
if (nrow(abund_mat) >= 4) {
  fit <- Maaslin2(
    input_data      = abund_mat,
    input_metadata  = metadata_df,
    output          = outdir,
    fixed_effects   = fixed_effects,
    reference       = paste0("diagnosis,", opts$reference_level),
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM",
    max_significance = 0.25,
    plot_heatmap    = TRUE,
    plot_scatter    = TRUE
  )

  all_results <- read_tsv(file.path(outdir, "all_results.tsv"), show_col_types=FALSE)
  sig_results <- all_results %>% filter(qval < 0.25)
  write_tsv(sig_results, file.path(outdir, "significant_results.tsv"))
  cat("MaAsLin2 complete:", nrow(sig_results), "significant associations (q < 0.25)\n")
} else {
  cat("NOTE: Need >= 4 samples for MaAsLin2. Writing formatted tables only.\n")
  write_csv(abund_mat, file.path(outdir, "formatted_abundance.csv"))
  write_csv(metadata, file.path(outdir, "formatted_metadata.csv"))
}
