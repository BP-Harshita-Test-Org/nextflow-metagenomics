#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(reshape2)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--abundance",     type="character", help="Bracken species abundance TSV"),
  make_option("--metadata",      type="character", help="Sample metadata CSV"),
  make_option("--diversity_dir", type="character", help="Diversity results directory"),
  make_option("--maaslin2_dir",  type="character", help="MaAsLin2 results directory"),
  make_option("--outdir",        type="character", help="Output directory for figures"),
  make_option("--sample_id",     type="character", default="sample", help="Sample ID")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load abundance data ─────────────────────────────────────────────────────
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)

# ── 1. Taxonomic composition bar plot (top 20 species) ──────────────────────
top20 <- abundance %>%
  arrange(desc(fraction_total_reads)) %>%
  head(20) %>%
  mutate(name = factor(name, levels=rev(name)))

p_taxa <- ggplot(top20, aes(x=name, y=fraction_total_reads, fill=name)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) +
  labs(
    title = paste0("Taxonomic Composition — ", opts$sample_id),
    subtitle = "Top 20 species by relative abundance",
    x = NULL,
    y = "Relative Abundance"
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.y = element_text(face="italic"))

ggsave(file.path(outdir, "taxonomy_barplot.pdf"), p_taxa, width=12, height=8)
cat("Generated: taxonomy_barplot.pdf\n")

# ── 2. Species abundance distribution ──────────────────────────────────────
p_dist <- ggplot(abundance, aes(x=log10(new_est_reads + 1))) +
  geom_histogram(bins=50, fill="steelblue", color="white") +
  labs(
    title = paste0("Species Abundance Distribution — ", opts$sample_id),
    x = "log10(Estimated Reads + 1)",
    y = "Number of Species"
  ) +
  theme_minimal(base_size=12)

ggsave(file.path(outdir, "abundance_distribution.pdf"), p_dist, width=8, height=5)
cat("Generated: abundance_distribution.pdf\n")

# ── 3. Rank abundance curve ────────────────────────────────────────────────
ranked <- abundance %>%
  arrange(desc(fraction_total_reads)) %>%
  mutate(rank = row_number())

p_rank <- ggplot(ranked, aes(x=rank, y=fraction_total_reads)) +
  geom_line(color="darkred", linewidth=0.8) +
  geom_point(color="darkred", size=0.5) +
  scale_y_log10(labels=scales::percent_format()) +
  labs(
    title = paste0("Rank Abundance Curve — ", opts$sample_id),
    x = "Species Rank",
    y = "Relative Abundance (log scale)"
  ) +
  theme_minimal(base_size=12)

ggsave(file.path(outdir, "rank_abundance_curve.pdf"), p_rank, width=8, height=5)
cat("Generated: rank_abundance_curve.pdf\n")

# ── 4. Read diversity from alpha metrics (if available) ────────────────────
alpha_file <- file.path(opts$diversity_dir, "alpha_diversity.csv")
if (file.exists(alpha_file)) {
  alpha <- read_csv(alpha_file, show_col_types=FALSE)
  cat("Alpha diversity metrics loaded for visualization.\n")

  # Write summary table
  write_csv(alpha, file.path(outdir, "alpha_diversity_summary.csv"))
} else {
  cat("No alpha diversity file found; skipping alpha plot.\n")
}

cat("Visualization complete. Figures in:", outdir, "\n")
