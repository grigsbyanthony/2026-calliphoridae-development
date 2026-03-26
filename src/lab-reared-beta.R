library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(umap)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  install.packages("pairwiseAdonis")
}
library(pairwiseAdonis)

theme_pub <- function(base_size=14, base_family="cmu sans serif") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

ps <- qiime2R::qza_to_phyloseq(
  features = "data/filtered-dada-table-nmnc.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
)

ps_insects <- subset_samples(ps, Class == "Insect")

ps_lab <- subset_samples(ps_insects, Group == "EGG")

cat("=== AVAILABLE STAGES IN LAB-REARED SAMPLES ===\n")
available_stages <- unique(sample_data(ps_lab)$Stage)
cat("Available stages:", paste(available_stages, collapse = ", "), "\n")

egg_stages <- c("Egg (Generation 0)", "Egg (Generation 1)")
non_egg_stages <- available_stages[!available_stages %in% egg_stages]
cat("Stages included in analysis (excluding eggs):", paste(non_egg_stages, collapse = ", "), "\n")

ps_lab_stages <- subset_samples(ps_lab, Stage %in% non_egg_stages)

ps_lab_stages <- prune_taxa(taxa_sums(ps_lab_stages) > 0, ps_lab_stages)

stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
actual_stages <- unique(sample_data(ps_lab_stages)$Stage)
stage_order <- stage_order[stage_order %in% actual_stages]

sample_data(ps_lab_stages)$Stage <- factor(sample_data(ps_lab_stages)$Stage, levels = stage_order)

stage_colors <- brewer.pal(length(stage_order), "Set1")
names(stage_colors) <- stage_order

dist_bray <- phyloseq::distance(ps_lab_stages, method = "bray")
dist_unifrac <- phyloseq::distance(ps_lab_stages, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps_lab_stages, method = "wunifrac")

perform_permanova <- function(dist_matrix, metadata) {
  metadata_df <- as.data.frame(metadata)
  permanova <- adonis2(dist_matrix ~ metadata_df$Stage, permutations = 999, by = "terms")
  return(permanova)
}

metadata <- as.data.frame(sample_data(ps_lab_stages))

cat("=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
permanova_bray <- perform_permanova(dist_bray, metadata)
print(permanova_bray)

cat("\n=== PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
permanova_unifrac <- perform_permanova(dist_unifrac, metadata)
print(permanova_unifrac)

cat("\n=== PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
permanova_wunifrac <- perform_permanova(dist_wunifrac, metadata)
print(permanova_wunifrac)

format_pairwise_adonis_results <- function(pairwise_results) {
  results <- data.frame(
    Stage1 = character(),
    Stage2 = character(),
    R2 = numeric(),
    p_value = numeric(),
    significance = character(),
    consecutive = logical(),
    stringsAsFactors = FALSE
  )

  stages <- levels(sample_data(ps_lab_stages)$Stage)

  for (i in 1:nrow(pairwise_results)) {
    pair <- as.character(pairwise_results$pairs[i])
    pair_parts <- strsplit(pair, "_vs_")[[1]]
    stage1 <- pair_parts[1]
    stage2 <- pair_parts[2]

    r2 <- pairwise_results$R2[i]
    p_value <- pairwise_results$p.value[i]

    sig <- ""
    if (p_value < 0.001) sig <- "***"
    else if (p_value < 0.01) sig <- "**"
    else if (p_value < 0.05) sig <- "*"
    else sig <- "ns"

    stage1_idx <- which(stages == stage1)
    stage2_idx <- which(stages == stage2)
    consecutive <- abs(stage1_idx - stage2_idx) == 1

    results <- rbind(results, data.frame(
      Stage1 = stage1,
      Stage2 = stage2,
      R2 = r2,
      p_value = p_value,
      significance = sig,
      consecutive = consecutive,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")

metadata_df <- data.frame(sample_data(ps_lab_stages))
metadata_df$Stage <- as.factor(metadata_df$Stage)
rownames(metadata_df) <- sample_names(ps_lab_stages)

tryCatch({
  pairwise_bray_results <- pairwiseAdonis::pairwise.adonis2(
    dist_bray ~ Stage,
    data = metadata_df,
    permutations = 999
  )
  cat("SUCCESS: pairwiseAdonis with formula worked!\n")
}, error = function(e) {
  cat("pairwiseAdonis formula failed:", e$message, "\n")

  tryCatch({
    pairwise_bray_results <- pairwiseAdonis::pairwise.adonis2(
      x = dist_bray,
      factors = metadata_df$Stage,
      permutations = 999
    )
    cat("SUCCESS: pairwiseAdonis with factors worked!\n")
  }, error = function(e2) {
    cat("pairwiseAdonis factors also failed:", e2$message, "\n")

    cat("Using manual pairwise approach...\n")

    stages <- levels(metadata_df$Stage)
    pairwise_results <- data.frame(
      pairs = character(),
      F.Model = numeric(),
      R2 = numeric(),
      p.value = numeric(),
      p.adjusted = numeric(),
      stringsAsFactors = FALSE
    )

    p_values <- c()
    for (i in 1:(length(stages) - 1)) {
      for (j in (i + 1):length(stages)) {
        stage1 <- stages[i]
        stage2 <- stages[j]

        subset_samples <- rownames(metadata_df)[metadata_df$Stage %in% c(stage1, stage2)]
        subset_metadata <- metadata_df[subset_samples, , drop = FALSE]
        subset_metadata$Stage <- droplevels(subset_metadata$Stage)

        subset_dist <- as.dist(as.matrix(dist_bray)[subset_samples, subset_samples])

        result <- vegan::adonis2(subset_dist ~ Stage, data = subset_metadata, permutations = 999)

        pair_name <- paste(stage1, "vs", stage2)
        pairwise_results <- rbind(pairwise_results, data.frame(
          pairs = pair_name,
          F.Model = result$F[1],
          R2 = result$R2[1],
          p.value = result$`Pr(>F)`[1],
          p.adjusted = NA,
          stringsAsFactors = FALSE
        ))

        p_values <- c(p_values, result$`Pr(>F)`[1])
      }
    }

    pairwise_results$p.adjusted <- p.adjust(p_values, method = "fdr")
    pairwise_bray_results <- pairwise_results
    cat("SUCCESS: Manual pairwise approach completed!\n")
  })
})

if (exists("format_pairwise_adonis_results") && is.function(format_pairwise_adonis_results)) {
  pairwise_bray <- format_pairwise_adonis_results(pairwise_bray_results)
} else {
  pairwise_bray <- pairwise_bray_results
}

print(pairwise_bray)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")

metadata_df <- data.frame(sample_data(ps_lab_stages))
metadata_df$Stage <- as.factor(metadata_df$Stage)
rownames(metadata_df) <- sample_names(ps_lab_stages)

tryCatch({
  pairwise_unifrac_results <- pairwiseAdonis::pairwise.adonis2(
    dist_unifrac ~ Stage,
    data = metadata_df,
    permutations = 999
  )
  cat("SUCCESS: pairwiseAdonis with formula worked!\n")
}, error = function(e) {
  cat("pairwiseAdonis formula failed:", e$message, "\n")

  tryCatch({
    pairwise_unifrac_results <- pairwiseAdonis::pairwise.adonis2(
      x = dist_unifrac,
      factors = metadata_df$Stage,
      permutations = 999
    )
    cat("SUCCESS: pairwiseAdonis with factors worked!\n")
  }, error = function(e2) {
    cat("pairwiseAdonis factors also failed:", e2$message, "\n")

    cat("Using manual pairwise approach...\n")

    stages <- levels(metadata_df$Stage)
    pairwise_results <- data.frame(
      pairs = character(),
      F.Model = numeric(),
      R2 = numeric(),
      p.value = numeric(),
      p.adjusted = numeric(),
      stringsAsFactors = FALSE
    )

    p_values <- c()
    for (i in 1:(length(stages) - 1)) {
      for (j in (i + 1):length(stages)) {
        stage1 <- stages[i]
        stage2 <- stages[j]

        subset_samples <- rownames(metadata_df)[metadata_df$Stage %in% c(stage1, stage2)]
        subset_metadata <- metadata_df[subset_samples, , drop = FALSE]
        subset_metadata$Stage <- droplevels(subset_metadata$Stage)

        subset_dist <- as.dist(as.matrix(dist_unifrac)[subset_samples, subset_samples])

        result <- vegan::adonis2(subset_dist ~ Stage, data = subset_metadata, permutations = 999)

        pair_name <- paste(stage1, "vs", stage2)
        pairwise_results <- rbind(pairwise_results, data.frame(
          pairs = pair_name,
          F.Model = result$F[1],
          R2 = result$R2[1],
          p.value = result$`Pr(>F)`[1],
          p.adjusted = NA,
          stringsAsFactors = FALSE
        ))

        p_values <- c(p_values, result$`Pr(>F)`[1])
      }
    }

    pairwise_results$p.adjusted <- p.adjust(p_values, method = "fdr")
    pairwise_unifrac_results <- pairwise_results
    cat("SUCCESS: Manual pairwise approach completed!\n")
  })
})

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")

metadata_df <- data.frame(sample_data(ps_lab_stages))
metadata_df$Stage <- as.factor(metadata_df$Stage)
rownames(metadata_df) <- sample_names(ps_lab_stages)

tryCatch({
  pairwise_wunifrac_results <- pairwiseAdonis::pairwise.adonis2(
    dist_wunifrac ~ Stage,
    data = metadata_df,
    permutations = 999
  )
  cat("SUCCESS: pairwiseAdonis with formula worked!\n")
}, error = function(e) {
  cat("pairwiseAdonis formula failed:", e$message, "\n")

  tryCatch({
    pairwise_wunifrac_results <- pairwiseAdonis::pairwise.adonis2(
      x = dist_wunifrac,
      factors = metadata_df$Stage,
      permutations = 999
    )
    cat("SUCCESS: pairwiseAdonis with factors worked!\n")
  }, error = function(e2) {
    cat("pairwiseAdonis factors also failed:", e2$message, "\n")

    cat("Using manual pairwise approach...\n")

    stages <- levels(metadata_df$Stage)
    pairwise_results <- data.frame(
      pairs = character(),
      F.Model = numeric(),
      R2 = numeric(),
      p.value = numeric(),
      p.adjusted = numeric(),
      stringsAsFactors = FALSE
    )

    p_values <- c()
    for (i in 1:(length(stages) - 1)) {
      for (j in (i + 1):length(stages)) {
        stage1 <- stages[i]
        stage2 <- stages[j]

        subset_samples <- rownames(metadata_df)[metadata_df$Stage %in% c(stage1, stage2)]
        subset_metadata <- metadata_df[subset_samples, , drop = FALSE]
        subset_metadata$Stage <- droplevels(subset_metadata$Stage)

        subset_dist <- as.dist(as.matrix(dist_wunifrac)[subset_samples, subset_samples])

        result <- vegan::adonis2(subset_dist ~ Stage, data = subset_metadata, permutations = 999)

        pair_name <- paste(stage1, "vs", stage2)
        pairwise_results <- rbind(pairwise_results, data.frame(
          pairs = pair_name,
          F.Model = result$F[1],
          R2 = result$R2[1],
          p.value = result$`Pr(>F)`[1],
          p.adjusted = NA,
          stringsAsFactors = FALSE
        ))

        p_values <- c(p_values, result$`Pr(>F)`[1])
      }
    }

    pairwise_results$p.adjusted <- p.adjust(p_values, method = "fdr")
    pairwise_wunifrac_results <- pairwise_results
    cat("SUCCESS: Manual pairwise approach completed!\n")
  })
})

perform_pcoa <- function(ps_obj, distance_method) {
  dist_matrix <- phyloseq::distance(ps_obj, method = distance_method)

  pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)

  variance_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)

  pcoa_df <- data.frame(
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    Sample = rownames(pcoa$points)
  )

  metadata <- as.data.frame(sample_data(ps_obj))
  pcoa_df <- merge(pcoa_df, metadata, by.x = "Sample", by.y = "row.names")

  return(list(
    ordination = pcoa_df,
    variance_explained = variance_explained
  ))
}

perform_umap <- function(ps_obj, distance_method) {
  dist_matrix <- phyloseq::distance(ps_obj, method = distance_method)

  dist_matrix_mat <- as.matrix(dist_matrix)

  set.seed(42)
  umap_result <- umap(dist_matrix_mat, n_neighbors = min(15, nrow(dist_matrix_mat) - 1), min_dist = 0.1)

  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    Sample = rownames(dist_matrix_mat)
  )

  metadata <- as.data.frame(sample_data(ps_obj))
  umap_df <- merge(umap_df, metadata, by.x = "Sample", by.y = "row.names")

  return(umap_df)
}

library(ggExtra)

create_pcoa_plot <- function(pcoa_result, title, permanova_result) {
  pcoa_df <- pcoa_result$ordination
  variance_explained <- pcoa_result$variance_explained

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Stage, shape = RearingGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = stage_colors) +
    labs(
      title = title,
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme_pub() +
    theme(
      legend.position = "right"
    )

  return(p)
}

create_umap_plot <- function(umap_df, title, permanova_result) {
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Stage, shape = RearingGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = stage_colors) +
    labs(
      title = title,
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_pub() +
    theme(
      legend.position = "right"
    )

  return(p)
}

pcoa_bray <- perform_pcoa(ps_lab_stages, "bray")
umap_bray <- perform_umap(ps_lab_stages, "bray")

pcoa_unifrac <- perform_pcoa(ps_lab_stages, "unifrac")
umap_unifrac <- perform_umap(ps_lab_stages, "unifrac")

pcoa_wunifrac <- perform_pcoa(ps_lab_stages, "wunifrac")
umap_wunifrac <- perform_umap(ps_lab_stages, "wunifrac")

pcoa_bray_plot <- create_pcoa_plot(pcoa_bray, "PCoA - Bray-Curtis", permanova_bray)
umap_bray_plot <- create_umap_plot(umap_bray, "UMAP - Bray-Curtis", permanova_bray)

pcoa_unifrac_plot <- create_pcoa_plot(pcoa_unifrac, "PCoA - Unweighted UniFrac", permanova_unifrac)
umap_unifrac_plot <- create_umap_plot(umap_unifrac, "UMAP - Unweighted UniFrac", permanova_unifrac)

pcoa_wunifrac_plot <- create_pcoa_plot(pcoa_wunifrac, "PCoA - Weighted UniFrac", permanova_wunifrac)
umap_wunifrac_plot <- create_umap_plot(umap_wunifrac, "UMAP - Weighted UniFrac", permanova_wunifrac)

combine_plots <- function(pcoa_plot, umap_plot, distance_name) {
  legend <- cowplot::get_legend(pcoa_plot + theme(legend.position = "bottom"))

  pcoa_with_marginals <- ggExtra::ggMarginal(
    pcoa_plot +
      theme(legend.position = "none") +
      labs(title = paste0("PCoA - ", distance_name)),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  umap_with_marginals <- ggExtra::ggMarginal(
    umap_plot +
      theme(legend.position = "none") +
      labs(title = paste0("UMAP - ", distance_name)),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  combined_plot <- cowplot::plot_grid(
    pcoa_with_marginals, umap_with_marginals,
    ncol = 2,
    align = "h",
    axis = "tb"
  )

  final_plot <- cowplot::plot_grid(
    combined_plot, legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )

  return(final_plot)
}

create_stacked_plot <- function(bray_combined, unifrac_combined, wunifrac_combined) {
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Stage, shape = RearingGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = stage_colors) +
    theme_pub() +
    theme(legend.position = "bottom")

  legend <- cowplot::get_legend(legend_plot)

  stacked_plot <- cowplot::plot_grid(
    bray_combined,
    unifrac_combined,
    wunifrac_combined,
    ncol = 1,
    align = "v",
    axis = "lr"
  )

  final_plot <- cowplot::plot_grid(
    stacked_plot,
    legend,
    ncol = 1,
    rel_heights = c(0.95, 0.05)
  )

  return(final_plot)
}

create_stacked_plot_alt <- function() {
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Stage, shape = RearingGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = stage_colors) +
    theme_pub() +
    theme(legend.position = "bottom")

  legend <- cowplot::get_legend(legend_plot)

  pcoa_bray_marginal <- ggExtra::ggMarginal(
    pcoa_bray_plot +
      theme(legend.position = "none") +
      labs(title = "PCoA - Bray-Curtis"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  umap_bray_marginal <- ggExtra::ggMarginal(
    umap_bray_plot +
      theme(legend.position = "none") +
      labs(title = "UMAP - Bray-Curtis"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  pcoa_unifrac_marginal <- ggExtra::ggMarginal(
    pcoa_unifrac_plot +
      theme(legend.position = "none") +
      labs(title = "PCoA - Unweighted UniFrac"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  umap_unifrac_marginal <- ggExtra::ggMarginal(
    umap_unifrac_plot +
      theme(legend.position = "none") +
      labs(title = "UMAP - Unweighted UniFrac"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  pcoa_wunifrac_marginal <- ggExtra::ggMarginal(
    pcoa_wunifrac_plot +
      theme(legend.position = "none") +
      labs(title = "PCoA - Weighted UniFrac"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  umap_wunifrac_marginal <- ggExtra::ggMarginal(
    umap_wunifrac_plot +
      theme(legend.position = "none") +
      labs(title = "UMAP - Weighted UniFrac"),
    type = "density",
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )

  bray_row <- cowplot::plot_grid(
    pcoa_bray_marginal, umap_bray_marginal,
    ncol = 2,
    align = "h",
    axis = "tb"
  )

  unifrac_row <- cowplot::plot_grid(
    pcoa_unifrac_marginal, umap_unifrac_marginal,
    ncol = 2,
    align = "h",
    axis = "tb"
  )

  wunifrac_row <- cowplot::plot_grid(
    pcoa_wunifrac_marginal, umap_wunifrac_marginal,
    ncol = 2,
    align = "h",
    axis = "tb"
  )

  stacked_plot <- cowplot::plot_grid(
    bray_row,
    unifrac_row,
    wunifrac_row,
    ncol = 1,
    align = "v",
    axis = "lr"
  )

  final_plot <- cowplot::plot_grid(
    stacked_plot,
    legend,
    ncol = 1,
    rel_heights = c(0.95, 0.05)
  )

  return(final_plot)
}

bray_combined <- combine_plots(pcoa_bray_plot, umap_bray_plot, "Bray-Curtis")
unifrac_combined <- combine_plots(pcoa_unifrac_plot, umap_unifrac_plot, "Unweighted UniFrac")
wunifrac_combined <- combine_plots(pcoa_wunifrac_plot, umap_wunifrac_plot, "Weighted UniFrac")

tryCatch({
  all_distances_combined <- create_stacked_plot(bray_combined, unifrac_combined, wunifrac_combined)
}, error = function(e) {
  message("Using alternative stacked plot method")
  all_distances_combined <<- create_stacked_plot_alt()
})

legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Stage, shape = RearingGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +
  theme_pub() +
  theme(legend.position = "bottom")

legend_only <- cowplot::get_legend(legend_plot)

legend_plot_only <- cowplot::plot_grid(legend_only)

ggsave("lab_reared_bray_curtis.png", bray_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("lab_reared_unweighted_unifrac.png", unifrac_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("lab_reared_weighted_unifrac.png", wunifrac_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("lab_reared_all_distances.png", all_distances_combined + theme(legend.position = "none"), width = 12, height = 18, dpi = 300, bg = "white")

ggsave("lab_reared_legend.png", legend_plot_only, width = 12, height = 2, dpi = 300, bg = "white")

legend_plot_alt <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Stage, shape = RearingGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +
  guides(color = guide_legend(title = "Stage", override.aes = list(size = 5)),
         shape = guide_legend(title = "Rearing Group", override.aes = list(size = 5))) +
  theme_pub(base_family = "cmu sans serif") +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", family = "cmu sans serif"),
        legend.text = element_text(size = 12, family = "cmu sans serif"),
        legend.key.size = unit(1, "cm"),
        text = element_text(family = "cmu sans serif"))

ggsave("lab_reared_legend_alt.png", legend_plot_alt, width = 12, height = 3, dpi = 300, bg = "white")

cat("=== LAB-REARED BETA DIVERSITY ANALYSIS (EXCLUDING EGGS) ===\n")
cat("Number of samples:", nsamples(ps_lab_stages), "\n")
cat("Number of taxa:", ntaxa(ps_lab_stages), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")

cat("\nVisualization files created:\n")
cat("- lab_reared_bray_curtis.png (PCoA and UMAP for Bray-Curtis distance)\n")
cat("- lab_reared_unweighted_unifrac.png (PCoA and UMAP for Unweighted UniFrac distance)\n")
cat("- lab_reared_weighted_unifrac.png (PCoA and UMAP for Weighted UniFrac distance)\n")
