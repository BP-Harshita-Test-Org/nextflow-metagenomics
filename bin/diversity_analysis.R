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
  make_option("--abundance", type="character", help="Merged abundance matrix TSV (features x samples)"),
  make_option("--metadata",  type="character", help="Sample metadata CSV"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--fdr_alpha", type="double", default=0.05, help="FDR cutoff for diversity tests"),
  make_option("--fdr_ancombc", type="double", default=0.2, help="FDR cutoff for ANCOM-BC")
)
opts <- parse_args(OptionParser(option_list=option_list))

outdir <- opts$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Load data ───────────────────────────────────────────────────────────────
# Merged matrix: first column = feature, remaining columns = samples
abundance <- read_tsv(opts$abundance, show_col_types=FALSE)
metadata  <- read_csv(opts$metadata, show_col_types=FALSE)

sample_ids <- colnames(abundance)[-1]
feature_names <- abundance$feature

cat("Loaded abundance matrix:", length(feature_names), "features x", length(sample_ids), "samples\n")
cat("Loaded metadata:", nrow(metadata), "samples\n")

# ── Build abundance matrix (samples as rows, species as columns) ────────────
abund_mat <- abundance %>% select(-feature) %>% as.data.frame()
rownames(abund_mat) <- feature_names
abund_mat <- t(abund_mat)  # now: samples x species

# ── Alpha diversity per sample ──────────────────────────────────────────────
alpha_df <- data.frame(
  sample_id = rownames(abund_mat),
  shannon   = diversity(abund_mat, index="shannon"),
  simpson   = diversity(abund_mat, index="simpson"),
  observed_species = apply(abund_mat > 0, 1, sum)
)

# Merge with metadata for group info
group_col <- if ("group" %in% names(metadata)) "group" else if ("diagnosis" %in% names(metadata)) "diagnosis" else NULL
alpha_df <- merge(alpha_df, metadata, by="sample_id", all.x=TRUE)

write_csv(alpha_df, file.path(outdir, "alpha_diversity.csv"))
cat("Alpha diversity computed for", nrow(alpha_df), "samples\n")

# ── Alpha diversity group comparison ────────────────────────────────────────
if (!is.null(group_col)) {
  p_alpha <- ggplot(alpha_df, aes(x=.data[[group_col]], y=shannon, fill=.data[[group_col]])) +
    geom_boxplot(outlier.shape=NA, alpha=0.7) +
    geom_jitter(width=0.2, size=2, alpha=0.6) +
    labs(title="Shannon Diversity by Group", x=NULL, y="Shannon Index") +
    theme_minimal(base_size=12) +
    theme(legend.position="none")
  ggsave(file.path(outdir, "alpha_diversity_boxplot.pdf"), p_alpha, width=8, height=6)

  # Wilcoxon test between groups
  groups <- unique(alpha_df[[group_col]])
  if (length(groups) == 2) {
    test_result <- wilcox.test(
      alpha_df$shannon[alpha_df[[group_col]] == groups[1]],
      alpha_df$shannon[alpha_df[[group_col]] == groups[2]]
    )
    cat("Shannon diversity Wilcoxon p-value:", test_result$p.value, "\n")
    writeLines(
      paste("Wilcoxon p-value (Shannon):", test_result$p.value),
      file.path(outdir, "alpha_diversity_test.txt")
    )
  }
}

# ── Beta diversity (Bray-Curtis PCoA) ───────────────────────────────────────
if (nrow(abund_mat) >= 3) {
  bc_dist <- vegdist(abund_mat, method="bray")
  pcoa_res <- cmdscale(bc_dist, k=2, eig=TRUE)

  pcoa_df <- data.frame(
    sample_id = rownames(abund_mat),
    PC1 = pcoa_res$points[,1],
    PC2 = pcoa_res$points[,2]
  )
  pcoa_df <- merge(pcoa_df, metadata, by="sample_id", all.x=TRUE)

  var_explained <- round(100 * pcoa_res$eig[1:2] / sum(pcoa_res$eig[pcoa_res$eig > 0]), 1)

  if (!is.null(group_col)) {
    p_pcoa <- ggplot(pcoa_df, aes(x=PC1, y=PC2, color=.data[[group_col]])) +
      geom_point(size=3, alpha=0.8) +
      stat_ellipse(level=0.95, linetype=2) +
      labs(
        title="Beta Diversity — Bray-Curtis PCoA",
        x=paste0("PC1 (", var_explained[1], "%)"),
        y=paste0("PC2 (", var_explained[2], "%)")
      ) +
      theme_minimal(base_size=12)
    ggsave(file.path(outdir, "beta_diversity_pcoa.pdf"), p_pcoa, width=10, height=7)
  }

  write_csv(pcoa_df, file.path(outdir, "pcoa_coordinates.csv"))

  # PERMANOVA
  if (!is.null(group_col)) {
    perm_meta <- metadata[match(rownames(abund_mat), metadata$sample_id), ]
    perm_result <- adonis2(bc_dist ~ perm_meta[[group_col]], permutations=999)
    sink(file.path(outdir, "permanova_results.txt"))
    print(perm_result)
    sink()
    cat("PERMANOVA R2:", perm_result$R2[1], "p-value:", perm_result$`Pr(>F)`[1], "\n")
  }
} else {
  cat("NOTE: Need >= 3 samples for beta diversity. Skipping PCoA.\n")
}

# ── Top 20 species (mean across samples) ────────────────────────────────────
mean_abund <- data.frame(
  species = feature_names,
  mean_abundance = rowMeans(abundance %>% select(-feature))
)
top20 <- mean_abund %>% arrange(desc(mean_abundance)) %>% head(20)
write_csv(top20, file.path(outdir, "top20_species.csv"))

p_top20 <- ggplot(top20, aes(x=reorder(species, mean_abundance), y=mean_abundance)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Top 20 Species (Mean Across Samples)", x="Species", y="Mean Relative Abundance") +
  theme_minimal(base_size=12) +
  theme(axis.text.y = element_text(face="italic"))
ggsave(file.path(outdir, "top20_species_barplot.pdf"), p_top20, width=10, height=7)

# ── ANCOM-BC (differential abundance) ───────────────────────────────────────
if (nrow(abund_mat) >= 4 && !is.null(group_col)) {
  tryCatch({
    library(ANCOMBC)

    otu_mat <- round(t(abund_mat) * 1e6)
    mode(otu_mat) <- "integer"
    OTU <- otu_table(otu_mat, taxa_are_rows=TRUE)

    meta_df <- metadata %>%
      filter(sample_id %in% rownames(abund_mat)) %>%
      tibble::column_to_rownames("sample_id")
    SAMP <- sample_data(meta_df)

    ps <- phyloseq(OTU, SAMP)

    ancom_out <- ancombc2(
      data = ps,
      fix_formula = group_col,
      p_adj_method = "BH",
      alpha = opts$fdr_ancombc
    )
    write_csv(ancom_out$res, file.path(outdir, "ancombc_results.csv"))
    cat("ANCOM-BC complete\n")
  }, error = function(e) {
    cat("ANCOM-BC skipped:", conditionMessage(e), "\n")
  })
} else {
  cat("NOTE: Need >= 4 samples for ANCOM-BC.\n")
}

cat("Diversity analysis complete. Results in:", outdir, "\n")
