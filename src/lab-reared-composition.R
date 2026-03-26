library(qiime2R)
library(phyloseq)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidyr)
library(ape)
library(RColorBrewer)
library(viridis)
library(ggalluvial)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(scales)


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

scale_fill_pub <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_pub <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
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

stage_counts <- table(sample_data(ps_lab)$Stage)
cat("Sample counts per stage:\n")
print(stage_counts)

egg_stages <- c("Egg (Generation 0)", "Egg (Generation 1)")
non_egg_stages <- available_stages[!available_stages %in% egg_stages]
cat("Stages included in analysis (excluding eggs):", paste(non_egg_stages, collapse = ", "), "\n")

ps_lab_stages <- subset_samples(ps_lab, Stage %in% non_egg_stages)

ps_lab_stages <- prune_taxa(taxa_sums(ps_lab_stages) > 0, ps_lab_stages)

ps_rel <- transform_sample_counts(ps_lab_stages, function(x) x / sum(x))

prepare_composition_data <- function(ps_obj, tax_level = "Genus", top_n = 15) {
  otu_table_df <- as.data.frame(otu_table(ps_obj))
  tax_table_df <- as.data.frame(tax_table(ps_obj))
  sample_data_df <- as.data.frame(sample_data(ps_obj))

  if(taxa_are_rows(ps_obj)) {
    otu_table_df <- t(otu_table_df)
  }

  tax_table_df <- tax_table_df %>%
    mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
    mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

  tax_col <- tax_table_df[[tax_level]]

  agg_data <- data.frame()

  for(sample in rownames(otu_table_df)) {
    sample_abundances <- otu_table_df[sample, ]

    tax_abundances <- tapply(sample_abundances, tax_col, sum, na.rm = TRUE)

    sample_df <- data.frame(
      Sample = sample,
      Taxon = names(tax_abundances),
      Abundance = as.numeric(tax_abundances),
      stringsAsFactors = FALSE
    )

    agg_data <- rbind(agg_data, sample_df)
  }

  agg_data <- merge(agg_data, sample_data_df, by.x = "Sample", by.y = "row.names")

  agg_data <- agg_data %>%
    group_by(Sample) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ungroup()

  top_taxa <- agg_data %>%
    filter(Taxon != "Unknown") %>%
    group_by(Taxon) %>%
    summarise(MeanAbundance = mean(RelativeAbundance)) %>%
    arrange(desc(MeanAbundance)) %>%
    head(top_n) %>%
    pull(Taxon)

  agg_data <- agg_data %>%
    mutate(Taxon_grouped = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
    group_by(Sample, Stage, Taxon_grouped) %>%
    summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")

  agg_data <- agg_data %>%
    mutate(Taxon_grouped = case_when(
      Taxon_grouped == "Lactococcus_A_343473" ~ "Lactococcus A.1",
      Taxon_grouped == "Vagococcus_B" ~ "Vagococcus B",
      Taxon_grouped == "Vagococcus_A" ~ "Vagococcus A",
      Taxon_grouped == "Lactococcus_A_346120" ~ "Lactococcus A.2",
      Taxon_grouped == "Jeotgalicoccus_A_310962" ~ "Jeotgalicoccus A",
      Taxon_grouped == "Leuconostoc_B" ~ "Leuconostoc B",
      Taxon_grouped == "Mammaliicoccus_319278" ~ "Mammaliicoccus",
      TRUE ~ Taxon_grouped
    ))

  return(agg_data)
}

genus_data <- prepare_composition_data(ps_rel, "Genus", 15)

create_stacked_bar_chart <- function(data, title = "Community Composition by Developmental Stage") {
  stage_means <- data %>%
    group_by(Stage, Taxon_grouped) %>%
    summarise(MeanRelativeAbundance = mean(RelativeAbundance), .groups = "drop")

  actual_stages <- unique(stage_means$Stage)

  all_stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")

  stage_order <- all_stage_order[all_stage_order %in% actual_stages]

  stage_means$Stage <- factor(stage_means$Stage, levels = stage_order)

  unique_taxa <- unique(stage_means$Taxon_grouped)
  n_taxa <- length(unique_taxa)

  if(n_taxa <= 11) {
    base_colors <- brewer.pal(max(3, n_taxa), "Spectral")
  } else {
    base_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n_taxa)
  }

  colors <- base_colors[1:n_taxa]
  names(colors) <- unique_taxa
  if("Other" %in% unique_taxa) {
    colors["Other"] <- "#808080"
  }

  p <- ggplot(stage_means, aes(x = Stage, y = MeanRelativeAbundance, fill = Taxon_grouped)) +
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Genus") +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title = title,
      x = "Developmental Stage",
      y = "Mean Relative Abundance"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(ncol = 1))

  return(p)
}

stacked_bar_plot <- create_stacked_bar_chart(genus_data)
print(stacked_bar_plot)

ggsave("lab_reared_stacked_bar_chart.png", stacked_bar_plot,
       width = 12, height = 8, dpi = 300, bg = "white")

prepare_alluvial_data <- function(data, min_abundance = 0.01) {
  stage_means <- data %>%
    group_by(Stage, Taxon_grouped) %>%
    summarise(MeanRelativeAbundance = mean(RelativeAbundance), .groups = "drop") %>%
    filter(MeanRelativeAbundance >= min_abundance)

  actual_stages <- unique(stage_means$Stage)
  all_stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
  stage_order <- all_stage_order[all_stage_order %in% actual_stages]

  stage_means$Stage <- factor(stage_means$Stage, levels = stage_order)

  alluvial_wide <- stage_means %>%
    select(Stage, Taxon_grouped, MeanRelativeAbundance) %>%
    pivot_wider(names_from = Stage, values_from = MeanRelativeAbundance, values_fill = 0)

  alluvial_long <- alluvial_wide %>%
    pivot_longer(cols = -Taxon_grouped, names_to = "Stage", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    mutate(Stage = factor(Stage, levels = stage_order))

  return(alluvial_long)
}

create_alluvial_plot <- function(data, title = "Microbial Community Flow Across Developmental Stages") {
  alluvial_data <- prepare_alluvial_data(data, min_abundance = 0.01)

  taxa_abundance <- alluvial_data %>%
    group_by(Taxon_grouped) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))

  if("Other" %in% taxa_abundance$Taxon_grouped) {
    taxa_abundance <- taxa_abundance %>%
      filter(Taxon_grouped != "Other") %>%
      bind_rows(taxa_abundance %>% filter(Taxon_grouped == "Other"))
  }

  ordered_taxa <- taxa_abundance$Taxon_grouped
  n_taxa <- length(ordered_taxa)

  if(n_taxa <= 11) {
    base_colors <- brewer.pal(max(3, n_taxa), "Set3")
  } else {
    base_colors <- colorRampPalette(brewer.pal(11, "Set3"))(n_taxa)
  }

  colors <- base_colors[1:n_taxa]
  names(colors) <- ordered_taxa
  if("Other" %in% ordered_taxa) {
    colors["Other"] <- "#808080"
  }

  alluvial_data$Taxon_grouped <- factor(alluvial_data$Taxon_grouped,
                                       levels = ordered_taxa)

  p <- ggplot(alluvial_data,
              aes(x = Stage, stratum = Taxon_grouped, alluvium = Taxon_grouped,
                  y = Abundance, fill = Taxon_grouped)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray", alpha = 0.7) +
    geom_stratum(alpha = 0.8) +
    scale_fill_manual(values = colors, name = "Genus", breaks = ordered_taxa) +
    labs(
      x = "Developmental Stage",
      y = "Relative Abundance"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.title = element_text(angle = 90, face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0), labels = percent_format())

  return(p)
}

alluvial_plot <- create_alluvial_plot(genus_data)
print(alluvial_plot)

ggsave("lab_reared_alluvial_plot.png", alluvial_plot,
       width = 14, height = 10, dpi = 300, bg = "white")

prepare_heatmap_data <- function(ps_obj, top_n = 25) {
  heatmap_data <- prepare_composition_data(ps_obj, "Genus", top_n)

  heatmap_wide <- heatmap_data %>%
    select(Sample, Stage, Taxon_grouped, RelativeAbundance) %>%
    pivot_wider(names_from = Taxon_grouped, values_from = RelativeAbundance, values_fill = 0)

  stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
  actual_stages <- unique(heatmap_wide$Stage)
  stage_order <- stage_order[stage_order %in% actual_stages]

  heatmap_wide <- heatmap_wide %>%
    mutate(Stage = factor(Stage, levels = stage_order)) %>%
    arrange(Stage)

  sample_names <- heatmap_wide$Sample
  stage_info <- heatmap_wide$Stage
  heatmap_matrix <- as.matrix(heatmap_wide[, -c(1, 2)])
  rownames(heatmap_matrix) <- sample_names

  sample_annotation <- data.frame(Stage = stage_info, row.names = sample_names)

  return(list(matrix = heatmap_matrix, annotation = sample_annotation))
}

create_pheatmap <- function(ps_obj, title = "Microbial Community Heatmap") {
  heatmap_data <- prepare_heatmap_data(ps_obj, 25)

  log_matrix <- log10(heatmap_data$matrix + 1e-6)

  actual_stages <- unique(heatmap_data$annotation$Stage)

  all_stage_colors <- c("1st instar" = "#e41a1c", "2nd instar" = "#377eb8",
                       "3rd instar" = "#4daf4a", "Pupal" = "#984ea3",
                       "Adult" = "#ff7f00")

  stage_colors <- all_stage_colors[names(all_stage_colors) %in% actual_stages]
  annotation_colors <- list(Stage = stage_colors)

  heatmap_colors <- colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(100)

  p <- pheatmap::pheatmap(
    t(log_matrix),
    annotation_col = heatmap_data$annotation,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    color = heatmap_colors,
    main = title,
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 6
  )

  return(p)
}

heatmap_plot <- create_pheatmap(ps_lab_stages)

png("lab_reared_heatmap.png", width = 12, height = 10, units = "in", res = 300)
print(heatmap_plot)
dev.off()

create_complex_heatmap <- function(ps_obj) {
  heatmap_data <- prepare_heatmap_data(ps_obj, 25)

  log_matrix <- log10(heatmap_data$matrix + 1e-6)

  scaled_matrix <- scale(t(log_matrix))

  actual_stages <- unique(heatmap_data$annotation$Stage)

  all_stage_colors <- c("1st instar" = "#e41a1c", "2nd instar" = "#377eb8",
                       "3rd instar" = "#4daf4a", "Pupal" = "#984ea3",
                       "Adult" = "#ff7f00")

  stage_colors <- all_stage_colors[names(all_stage_colors) %in% actual_stages]

  col_annotation <- HeatmapAnnotation(
    Stage = heatmap_data$annotation$Stage,
    col = list(Stage = stage_colors),
    annotation_name_gp = gpar(fontsize = 10)
  )

  col_fun <- colorRamp2(c(-2, 0, 2), c("#0d0887", "#f0f921", "#cc4778"))

  ht <- Heatmap(
    scaled_matrix,
    name = "Z-score",
    col = col_fun,
    top_annotation = col_annotation,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    column_title = "Lab-Reared Samples by Developmental Stage (Excluding Eggs)",
    row_title = "Bacterial Genera",
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )

  return(ht)
}

complex_heatmap <- create_complex_heatmap(ps_lab_stages)

png("lab_reared_complex_heatmap.png", width = 14, height = 10, units = "in", res = 300)
draw(complex_heatmap)
dev.off()

create_stage_averaged_heatmap <- function(ps_obj) {
  heatmap_data <- prepare_heatmap_data(ps_obj, 25)

  log_matrix <- log10(heatmap_data$matrix + 1e-6)

  scaled_matrix <- scale(t(log_matrix))

  scaled_df <- as.data.frame(scaled_matrix)
  scaled_df$Genus <- rownames(scaled_df)

  stage_info <- heatmap_data$annotation$Stage
  names(stage_info) <- rownames(heatmap_data$annotation)

  scaled_long <- scaled_df %>%
    pivot_longer(cols = -Genus, names_to = "Sample", values_to = "ZScore") %>%
    mutate(Stage = stage_info[Sample])

  stage_means <- scaled_long %>%
    group_by(Genus, Stage) %>%
    summarise(MeanZScore = mean(ZScore, na.rm = TRUE), .groups = "drop")

  stage_wide <- stage_means %>%
    pivot_wider(names_from = Stage, values_from = MeanZScore, values_fill = 0)

  genus_names <- stage_wide$Genus
  stage_matrix <- as.matrix(stage_wide[, -1])
  rownames(stage_matrix) <- genus_names

  stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
  actual_stages <- colnames(stage_matrix)
  stage_order <- stage_order[stage_order %in% actual_stages]
  stage_matrix <- stage_matrix[, stage_order, drop = FALSE]

  col_fun <- colorRamp2(c(-2, 0, 2), c("#0d0887", "#f0f921", "#cc4778"))

  ht <- Heatmap(
    stage_matrix,
    name = "Mean Z-score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10, fontfamily = "cmu sans serif"),
    column_names_gp = gpar(fontsize = 12, fontfamily = "cmu sans serif", rot = 45),
    column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "cmu sans serif"),
    row_title = "Bacterial Genera",
    row_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "cmu sans serif"),
    heatmap_legend_param = list(
      title = "Mean Z-score",
      title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "cmu sans serif"),
      labels_gp = gpar(fontsize = 10, fontfamily = "cmu sans serif"),
      legend_direction = "vertical",
      legend_height = unit(4, "cm"),
      title_position = "topcenter"
    ),
    border = TRUE,
    width = unit(6, "cm"),
    height = unit(12, "cm")
  )

  return(ht)
}

stage_averaged_heatmap <- create_stage_averaged_heatmap(ps_lab_stages)

png("lab_reared_stage_averaged_heatmap.png", width = 12, height = 14, units = "in", res = 300)
draw(stage_averaged_heatmap)
dev.off()

cat("=== LAB-REARED MICROBIOME COMPOSITION ANALYSIS (EXCLUDING EGGS) ===\n")
cat("Number of samples:", nsamples(ps_lab_stages), "\n")
cat("Number of taxa:", ntaxa(ps_lab_stages), "\n")
cat("Developmental stages included:", paste(unique(sample_data(ps_lab_stages)$Stage), collapse = ", "), "\n")

diversity_by_stage <- genus_data %>%
  group_by(Stage) %>%
  summarise(
    Shannon = -sum(RelativeAbundance * log(RelativeAbundance + 1e-10)),
    Simpson = 1 - sum(RelativeAbundance^2),
    Richness = sum(RelativeAbundance > 0),
    .groups = "drop"
  )

cat("\nDiversity metrics by developmental stage:\n")
print(diversity_by_stage)

cat("\nTop 15 most abundant genera (used in bar charts and alluvial plots):\n")
top_genera_summary <- genus_data %>%
  filter(Taxon_grouped != "Other") %>%
  group_by(Taxon_grouped) %>%
  summarise(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance))

print(top_genera_summary)

cat("\nVisualization files created:\n")
cat("- lab_reared_stacked_bar_chart.png (top 15 genera)\n")
cat("- lab_reared_alluvial_plot.png (top 15 genera)\n")
cat("- lab_reared_heatmap.png (top 25 genera, samples grouped by stage, viridis-inspired palette)\n")
cat("- lab_reared_complex_heatmap.png (top 25 genera, samples grouped by stage, plasma palette)\n")
cat("- lab_reared_stage_averaged_heatmap.png (top 25 genera, mean Z-scores by stage, plasma palette)\n")
