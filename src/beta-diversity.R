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
            legend.key.size = unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin = unit(c(10,5,5,5),"mm"),
            strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

ps <- qiime2R::qza_to_phyloseq(
  features = "filtered-dada-table-nmnc.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

ps_insects <- subset_samples(ps, Class == "Insect")

target_stages <- c("3rd instar", "Pupal", "Adult")
ps_stages <- subset_samples(ps_insects, Stage %in% target_stages)

cat("=== AVAILABLE STAGES IN FILTERED DATA ===\n")
available_stages <- unique(sample_data(ps_stages)$Stage)
cat("Available stages:", paste(available_stages, collapse = ", "), "\n")

cat("\n=== AVAILABLE REARING CONDITIONS (STATUS) ===\n")
available_status <- unique(sample_data(ps_stages)$Status)
cat("Available rearing conditions:", paste(available_status, collapse = ", "), "\n")

ps_stages <- prune_taxa(taxa_sums(ps_stages) > 0, ps_stages)

stage_order <- c("3rd instar", "Pupal", "Adult")
sample_data(ps_stages)$Stage <- factor(sample_data(ps_stages)$Stage, levels = stage_order)

sample_data(ps_stages)$Status <- factor(sample_data(ps_stages)$Status)
available_status <- levels(sample_data(ps_stages)$Status)

n_colors <- max(3, length(available_status))
status_colors <- brewer.pal(n_colors, "Set1")[1:length(available_status)]
names(status_colors) <- available_status

cat("\nStatus levels (as factors):", paste(levels(sample_data(ps_stages)$Status), collapse = ", "), "\n")

dist_bray <- phyloseq::distance(ps_stages, method = "bray")
dist_unifrac <- phyloseq::distance(ps_stages, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps_stages, method = "wunifrac")

perform_permanova <- function(dist_matrix, metadata) {
  metadata_df <- as.data.frame(metadata)

  permanova <- adonis2(dist_matrix ~ metadata_df$Status, permutations = 999, by = "terms")

  return(permanova)
}

metadata <- as.data.frame(sample_data(ps_stages))

cat("=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
permanova_bray <- perform_permanova(dist_bray, metadata)
print(permanova_bray)

cat("\n=== PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
permanova_unifrac <- perform_permanova(dist_unifrac, metadata)
print(permanova_unifrac)

cat("\n=== PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
permanova_wunifrac <- perform_permanova(dist_wunifrac, metadata)
print(permanova_wunifrac)

metadata_df <- data.frame(sample_data(ps_stages))
metadata_df$Stage <- as.factor(metadata_df$Stage)
rownames(metadata_df) <- sample_names(ps_stages)

pairwise_dist_bray_results <- pairwiseAdonis::pairwise.adonis2(
  dist_bray ~ Stage,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_bray_results)

pairwise_dist_unifrac_results <- pairwiseAdonis::pairwise.adonis2(
  dist_unifrac ~ Stage,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_unifrac_results)

pairwise_dist_wunifrac_results <- pairwiseAdonis::pairwise.adonis2(
  dist_wunifrac ~ Stage,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_wunifrac_results)

pairwise_dist_bray_status_results <- pairwiseAdonis::pairwise.adonis2(
  dist_bray ~ Status,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_bray_status_results)

pairwise_dist_unifrac_status_results <- pairwiseAdonis::pairwise.adonis2(
  dist_unifrac ~ Status,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_unifrac_status_results)

pairwise_dist_wunifrac_status_results <- pairwiseAdonis::pairwise.adonis2(
  dist_wunifrac ~ Status,
  data = metadata_df,
  permutations = 999
)

print(pairwise_dist_wunifrac_status_results)

perform_pairwise_permanova <- function(ps_obj, distance_method) {
  status_values <- as.character(sample_data(ps_obj)$Status)
  conditions <- unique(status_values)

  cat("Conditions for pairwise PERMANOVA:", paste(conditions, collapse = ", "), "\n")

  results <- data.frame(
    Condition1 = character(),
    Condition2 = character(),
    R2 = numeric(),
    p_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )

  if (length(conditions) > 1) {
    for (i in 1:(length(conditions) - 1)) {
      for (j in (i+1):length(conditions)) {
        cond1 <- conditions[i]
        cond2 <- conditions[j]

        cat("Comparing conditions:", cond1, "vs", cond2, "\n")

        ps_subset <- subset_samples(ps_obj, Status %in% c(cond1, cond2))

        if (nsamples(ps_subset) < 3) {
          cat("Skipping due to insufficient samples\n")
          next
        }

        dist_matrix <- phyloseq::distance(ps_subset, method = distance_method)

        metadata_subset <- as.data.frame(sample_data(ps_subset))

        permanova <- adonis2(dist_matrix ~ metadata_subset$Status, permutations = 999, by = "terms")

        r2 <- permanova["Status", "R2"]
        p_value <- permanova["Status", "Pr(>F)"]

        sig <- ""
        if (p_value < 0.001) sig <- "***"
        else if (p_value < 0.01) sig <- "**"
        else if (p_value < 0.05) sig <- "*"
        else sig <- "ns"

        results <- rbind(results, data.frame(
          Condition1 = cond1,
          Condition2 = cond2,
          R2 = r2,
          p_value = p_value,
          significance = sig,
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    cat("Not enough conditions for pairwise comparison\n")
  }

  return(results)
}

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
pairwise_bray <- perform_pairwise_permanova(ps_stages, "bray")
print(pairwise_bray)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_unifrac <- perform_pairwise_permanova(ps_stages, "unifrac")
print(pairwise_unifrac)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_wunifrac <- perform_pairwise_permanova(ps_stages, "wunifrac")
print(pairwise_wunifrac)

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

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Status, shape = Stage)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = status_colors) +
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
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Status, shape = Stage)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = status_colors) +
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

pcoa_bray <- perform_pcoa(ps_stages, "bray")
umap_bray <- perform_umap(ps_stages, "bray")

pcoa_unifrac <- perform_pcoa(ps_stages, "unifrac")
umap_unifrac <- perform_umap(ps_stages, "unifrac")

pcoa_wunifrac <- perform_pcoa(ps_stages, "wunifrac")
umap_wunifrac <- perform_umap(ps_stages, "wunifrac")

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
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Status, shape = Stage)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = status_colors) +
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
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Status, shape = Stage)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = status_colors) +
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

legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Status, shape = Stage)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  theme_pub() +
  theme(legend.position = "bottom")

legend_only <- cowplot::get_legend(legend_plot)

legend_plot_only <- cowplot::plot_grid(legend_only)

ggsave("all_conditions_bray_curtis.png", bray_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_unweighted_unifrac.png", unifrac_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_weighted_unifrac.png", wunifrac_combined + theme(legend.position = "none"), width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_all_distances.png", all_distances_combined + theme(legend.position = "none"), width = 12, height = 18, dpi = 300, bg = "white")

ggsave("all_conditions_legend.png", legend_plot_only, width = 12, height = 2, dpi = 300, bg = "white")

legend_plot_alt <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Status, shape = Stage)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  guides(color = guide_legend(title = "Rearing Condition", override.aes = list(size = 5)),
         shape = guide_legend(title = "Stage", override.aes = list(size = 5))) +
  theme_pub(base_family = "cmu sans serif") +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", family = "cmu sans serif"),
        legend.text = element_text(size = 12, family = "cmu sans serif"),
        legend.key.size = unit(1, "cm"),
        text = element_text(family = "cmu sans serif"))

ggsave("all_conditions_legend_alt.png", legend_plot_alt, width = 12, height = 3, dpi = 300, bg = "white")

cat("=== ALL CONDITIONS BETA DIVERSITY ANALYSIS (3RD INSTAR, PUPAL, ADULT) ===\n")
cat("Number of samples:", nsamples(ps_stages), "\n")
cat("Number of taxa:", ntaxa(ps_stages), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")
cat("Rearing conditions included:", paste(available_status, collapse = ", "), "\n")

cat("\nVisualization files created:\n")
cat("- all_conditions_bray_curtis.png (PCoA and UMAP for Bray-Curtis distance)\n")
cat("- all_conditions_unweighted_unifrac.png (PCoA and UMAP for Unweighted UniFrac distance)\n")
cat("- all_conditions_weighted_unifrac.png (PCoA and UMAP for Weighted UniFrac distance)\n")
