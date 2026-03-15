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
  make_option("--abundance",     type="character", help="Merged abundance matrix TSV"),
  make_option("--metadata",      type="character", help="Sample metadata CSV"),
  make_option("--diversity_dir", type="character", help="Diversity results directory"),
  make_option("--maaslin2_dir",  type="character", help="MaAsLin2 results directory"),
  make_option("--outdir",        type="character", help="Output directory for figures")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ───────────────────────────────────────────────────────────────
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)
metadata  <- read_csv(opts$metadata, show_col_types=FALSE)

sample_ids <- colnames(abundance)[-1]
feature_names <- abundance$feature

cat("Loaded", length(feature_names), "features x", length(sample_ids), "samples\n")

# Detect group column
group_col <- if ("group" %in% names(metadata)) "group" else if ("diagnosis" %in% names(metadata)) "diagnosis" else NULL

# ── 1. Taxonomy barplot (top 20 species, grouped) ──────────────────────────
mean_abund <- data.frame(
  species = feature_names,
  mean_abundance = rowMeans(abundance %>% select(-feature))
)
top20_names <- mean_abund %>% arrange(desc(mean_abundance)) %>% head(20) %>% pull(species)

# Reshape for grouped plotting
long_df <- abundance %>%
  filter(feature %in% top20_names) %>%
  melt(id.vars="feature", variable.name="sample_id", value.name="abundance") %>%
  merge(metadata, by="sample_id", all.x=TRUE)

if (!is.null(group_col)) {
  group_means <- long_df %>%
    group_by(feature, .data[[group_col]]) %>%
    summarise(mean_abundance = mean(abundance, na.rm=TRUE), .groups="drop")

  p_taxa <- ggplot(group_means, aes(x=.data[[group_col]], y=mean_abundance, fill=feature)) +
    geom_bar(stat="identity", position="stack") +
    labs(
      title="Taxonomic Composition by Group",
      subtitle="Top 20 species (mean relative abundance)",
      x=NULL, y="Relative Abundance", fill="Species"
    ) +
    theme_minimal(base_size=12) +
    theme(legend.text = element_text(face="italic", size=8))
} else {
  p_taxa <- ggplot(mean_abund %>% filter(species %in% top20_names),
                   aes(x=reorder(species, mean_abundance), y=mean_abundance)) +
    geom_bar(stat="identity", fill="steelblue") +
    coord_flip() +
    labs(title="Top 20 Species", x=NULL, y="Mean Relative Abundance") +
    theme_minimal(base_size=12) +
    theme(axis.text.y = element_text(face="italic"))
}
ggsave(file.path(outdir, "taxonomy_barplot.pdf"), p_taxa, width=12, height=8)
cat("Generated: taxonomy_barplot.pdf\n")

# ── 2. Abundance distribution (all samples) ────────────────────────────────
abund_mat <- abundance %>% select(-feature) %>% as.matrix()
all_values <- as.vector(abund_mat)
all_values <- all_values[all_values > 0]

p_dist <- ggplot(data.frame(x=log10(all_values)), aes(x=x)) +
  geom_histogram(bins=50, fill="steelblue", color="white") +
  labs(
    title="Species Abundance Distribution (All Samples)",
    x="log10(Relative Abundance)",
    y="Count"
  ) +
  theme_minimal(base_size=12)
ggsave(file.path(outdir, "abundance_distribution.pdf"), p_dist, width=8, height=5)
cat("Generated: abundance_distribution.pdf\n")

# ── 3. Rank abundance curve (mean across samples) ──────────────────────────
ranked <- mean_abund %>%
  arrange(desc(mean_abundance)) %>%
  mutate(rank = row_number())

p_rank <- ggplot(ranked, aes(x=rank, y=mean_abundance)) +
  geom_line(color="darkred", linewidth=0.8) +
  geom_point(color="darkred", size=0.5) +
  scale_y_log10() +
  labs(
    title="Rank Abundance Curve (Mean Across Samples)",
    x="Species Rank",
    y="Mean Relative Abundance (log scale)"
  ) +
  theme_minimal(base_size=12)
ggsave(file.path(outdir, "rank_abundance_curve.pdf"), p_rank, width=8, height=5)
cat("Generated: rank_abundance_curve.pdf\n")

# ── 4. Alpha diversity summary from upstream ────────────────────────────────
alpha_file <- file.path(opts$diversity_dir, "alpha_diversity.csv")
if (file.exists(alpha_file)) {
  alpha <- read_csv(alpha_file, show_col_types=FALSE)
  write_csv(alpha, file.path(outdir, "alpha_diversity_summary.csv"))
  cat("Alpha diversity summary copied.\n")
} else {
  cat("No alpha diversity file found; skipping.\n")
}

# ── 5. MaAsLin2 significant results heatmap ─────────────────────────────────
sig_file <- file.path(opts$maaslin2_dir, "significant_results.tsv")
if (file.exists(sig_file)) {
  sig <- read_tsv(sig_file, show_col_types=FALSE)
  if (nrow(sig) > 0) {
    p_sig <- ggplot(sig %>% head(30),
                    aes(x=metadata, y=feature, fill=coef)) +
      geom_tile() +
      scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
      labs(title="Top Significant Associations (MaAsLin2)",
           x="Metadata Variable", y="Feature") +
      theme_minimal(base_size=10) +
      theme(axis.text.y = element_text(face="italic", size=8))
    ggsave(file.path(outdir, "maaslin2_heatmap.pdf"), p_sig, width=10, height=8)
    cat("Generated: maaslin2_heatmap.pdf\n")
  }
} else {
  cat("No MaAsLin2 results found; skipping heatmap.\n")
}

cat("Visualization complete. Figures in:", outdir, "\n")
