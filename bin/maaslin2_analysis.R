#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Maaslin2)
  library(readr)
  library(dplyr)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--abundance",      type="character", help="Bracken species abundance TSV"),
  make_option("--metadata",       type="character", help="Sample metadata CSV"),
  make_option("--outdir",         type="character", help="Output directory"),
  make_option("--fixed_effects",  type="character", default="diagnosis,age,sex",
              help="Comma-separated fixed effects"),
  make_option("--reference_level", type="character", default="neurotypical",
              help="Reference group level"),
  make_option("--sample_id",      type="character", default="sample", help="Sample ID")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ───────────────────────────────────────────────────────────────
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)
metadata  <- read_csv(opts$metadata, show_col_types=FALSE)

cat("Loaded abundance table:", nrow(abundance), "taxa\n")
cat("Loaded metadata:", nrow(metadata), "samples\n")

# ── Prepare MaAsLin2 input ──────────────────────────────────────────────────
# MaAsLin2 expects:
#   - abundance: rows=samples, columns=features (species)
#   - metadata:  rows=samples, columns=variables
#
# For single-sample mode, write the formatted tables for later cohort analysis.
# For multi-sample (cohort) mode, the merged abundance matrix is used directly.

fixed_effects <- strsplit(opts$fixed_effects, ",")[[1]]
cat("Fixed effects:", paste(fixed_effects, collapse=", "), "\n")
cat("Reference level:", opts$reference_level, "\n")

# Write formatted data for cohort-level analysis
write_csv(abundance, file.path(outdir, "formatted_abundance.csv"))
write_csv(metadata,  file.path(outdir, "formatted_metadata.csv"))

# ── Run MaAsLin2 (requires multi-sample data) ──────────────────────────────
# When abundance matrix has multiple samples (columns), run:
tryCatch({
  # Reshape: species as columns, samples as rows
  if (ncol(abundance) > 3 && nrow(metadata) > 1) {
    # Assumes merged abundance table with sample columns
    abundance_mat <- abundance %>%
      select(-taxonomy_id, -taxonomy_lvl, -kraken_assigned_reads,
             -added_reads, -new_est_reads, -fraction_total_reads) %>%
      tibble::column_to_rownames("name") %>%
      t() %>%
      as.data.frame()

    metadata_df <- metadata %>%
      tibble::column_to_rownames("sample_id")

    fit <- Maaslin2(
      input_data     = abundance_mat,
      input_metadata = metadata_df,
      output         = outdir,
      fixed_effects  = fixed_effects,
      reference      = paste0("diagnosis,", opts$reference_level),
      normalization  = "TSS",
      transform      = "LOG",
      analysis_method = "LM",
      max_significance = 0.25,
      plot_heatmap   = TRUE,
      plot_scatter   = TRUE
    )
    cat("MaAsLin2 analysis complete.\n")
  } else {
    cat("NOTE: Single-sample detected. MaAsLin2 requires a multi-sample abundance matrix.\n")
    cat("      Merge Bracken outputs across all samples, then re-run.\n")
    writeLines(
      "MaAsLin2 requires multi-sample input. Merge Bracken outputs first.",
      file.path(outdir, "README_run_with_cohort.txt")
    )
  }
}, error = function(e) {
  cat("MaAsLin2 skipped:", conditionMessage(e), "\n")
  writeLines(
    paste("MaAsLin2 error:", conditionMessage(e)),
    file.path(outdir, "maaslin2_error.txt")
  )
})
