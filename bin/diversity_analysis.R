#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--abundance", type="character", help="Bracken species abundance TSV"),
  make_option("--metadata",  type="character", help="Sample metadata CSV"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--fdr_alpha", type="double", default=0.05, help="FDR cutoff for diversity tests"),
  make_option("--fdr_ancombc", type="double", default=0.2, help="FDR cutoff for ANCOM-BC"),
  make_option("--sample_id", type="character", default="sample", help="Sample identifier")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ───────────────────────────────────────────────────────────────
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)
metadata  <- read_csv(opts$metadata, show_col_types=FALSE)

cat("Loaded abundance table:", nrow(abundance), "taxa\n")
cat("Loaded metadata:", nrow(metadata), "samples\n")

# ── Build phyloseq object ───────────────────────────────────────────────────
# Bracken output has columns: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
# added_reads, new_est_reads, fraction_total_reads
# For multi-sample analysis, abundance tables need to be merged first.
# This script handles single-sample; for cohort analysis, merge tables upstream.

otu_data <- abundance %>%
  select(name, new_est_reads) %>%
  tibble::column_to_rownames("name")

otu_mat <- as.matrix(otu_data)
OTU <- otu_table(otu_mat, taxa_are_rows=TRUE)

cat("Created phyloseq OTU table\n")

# ── Alpha diversity ─────────────────────────────────────────────────────────
# Compute on the abundance vector
counts <- abundance$new_est_reads
names(counts) <- abundance$name

shannon  <- diversity(counts, index="shannon")
simpson  <- diversity(counts, index="simpson")
observed <- sum(counts > 0)

alpha_df <- data.frame(
  sample_id = opts$sample_id,
  shannon   = shannon,
  simpson   = simpson,
  observed_species = observed
)

write_csv(alpha_df, file.path(outdir, "alpha_diversity.csv"))
cat("Alpha diversity — Shannon:", round(shannon, 4),
    "Simpson:", round(simpson, 4),
    "Observed:", observed, "\n")

# ── Beta diversity (Bray-Curtis matrix for single sample — placeholder) ─────
# Full beta diversity requires multiple samples. Write the abundance profile
# in a format ready for multi-sample merging.
write_csv(
  abundance %>% select(name, new_est_reads, fraction_total_reads),
  file.path(outdir, "abundance_for_beta_diversity.csv")
)

# ── Taxonomic summary ───────────────────────────────────────────────────────
# Top 20 species by relative abundance
top_species <- abundance %>%
  arrange(desc(fraction_total_reads)) %>%
  head(20)

write_csv(top_species, file.path(outdir, "top20_species.csv"))

# ── Bar plot of top species ─────────────────────────────────────────────────
p <- ggplot(top_species, aes(x=reorder(name, fraction_total_reads), y=fraction_total_reads)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title=paste0("Top 20 Species — ", opts$sample_id),
       x="Species", y="Relative Abundance") +
  theme_minimal()

ggsave(file.path(outdir, "top20_species_barplot.pdf"), p, width=10, height=7)

cat("Diversity analysis complete. Results in:", outdir, "\n")

# ── ANCOM-BC placeholder ───────────────────────────────────────────────────
# ANCOM-BC requires multiple samples with group labels. When running on the
# full cohort (merged abundance table), uncomment:
#
# library(ANCOMBC)
# ps <- phyloseq(OTU, sample_data(metadata_df))
# ancom_out <- ancombc2(
#   data = ps,
#   fix_formula = "diagnosis",
#   p_adj_method = "BH",
#   alpha = opts$fdr_ancombc
# )
# write_csv(ancom_out$res, file.path(outdir, "ancombc_results.csv"))

cat("NOTE: ANCOM-BC and PERMANOVA require multi-sample merged tables.\n")
cat("      Merge Bracken outputs across samples, then re-run with full cohort.\n")
