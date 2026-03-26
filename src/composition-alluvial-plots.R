library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)
library(RColorBrewer)
library(viridis)
library(ggalluvial)
library(scales)
library(patchwork)

theme_pub <- function(base_size = 14, base_family = "helvetica", legend_pos = "bottom") {
  library(ggthemes)
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.2, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}

ps <- qiime2R::qza_to_phyloseq(
  features = "data/filtered-dada-table-nmnc.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
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

status_order <- c("Lab-reared", "Carrion-reared")
status_order <- status_order[status_order %in% unique(sample_data(ps_stages)$Status)]
other_status <- setdiff(unique(sample_data(ps_stages)$Status), status_order)
status_order <- c(status_order, other_status)

sample_data(ps_stages)$Status <- factor(sample_data(ps_stages)$Status, levels = status_order)
available_status <- levels(sample_data(ps_stages)$Status)

n_colors <- max(3, length(available_status))
status_colors <- brewer.pal(n_colors, "Set1")[1:length(available_status)]
names(status_colors) <- available_status

cat("\nStatus levels (as factors):", paste(levels(sample_data(ps_stages)$Status), collapse = ", "), "\n")

ps_rel <- transform_sample_counts(ps_stages, function(x) x / sum(x))

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
    group_by(Sample, Status, Stage, Taxon_grouped) %>%
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

prepare_alluvial_data <- function(data, min_abundance = 0.01) {
  status_stage_means <- data %>%
    group_by(Status, Stage, Taxon_grouped) %>%
    summarise(MeanRelativeAbundance = mean(RelativeAbundance), .groups = "drop")

  filtered_sum <- status_stage_means %>%
    filter(MeanRelativeAbundance < min_abundance) %>%
    group_by(Status, Stage) %>%
    summarise(SumFilteredOut = sum(MeanRelativeAbundance), .groups = "drop")

  status_stage_means_filtered <- status_stage_means %>%
    filter(MeanRelativeAbundance >= min_abundance)

  other_rows <- filtered_sum %>%
    mutate(Taxon_grouped = "Other")

  for (i in 1:nrow(other_rows)) {
    current_status <- other_rows$Status[i]
    current_stage <- other_rows$Stage[i]
    current_sum <- other_rows$SumFilteredOut[i]

    existing_other <- status_stage_means_filtered %>%
      filter(Status == current_status, Stage == current_stage, Taxon_grouped == "Other")

    if (nrow(existing_other) > 0) {
      status_stage_means_filtered <- status_stage_means_filtered %>%
        mutate(MeanRelativeAbundance = ifelse(
          Status == current_status & Stage == current_stage & Taxon_grouped == "Other",
          MeanRelativeAbundance + current_sum,
          MeanRelativeAbundance
        ))
    } else if (current_sum > 0) {
      new_other <- data.frame(
        Status = current_status,
        Stage = current_stage,
        Taxon_grouped = "Other",
        MeanRelativeAbundance = current_sum
      )
      status_stage_means_filtered <- rbind(status_stage_means_filtered, new_other)
    }
  }

  actual_stages <- unique(status_stage_means_filtered$Stage)
  all_stage_order <- c("3rd instar", "Pupal", "Adult")
  stage_order <- all_stage_order[all_stage_order %in% actual_stages]

  status_stage_means_filtered$Stage <- factor(status_stage_means_filtered$Stage, levels = stage_order)

  return(status_stage_means_filtered)
}

create_faceted_alluvial_plot <- function(data, title = NULL) {
  alluvial_data <- prepare_alluvial_data(data, min_abundance = 0.01)

  taxa_abundance <- alluvial_data %>%
    group_by(Taxon_grouped) %>%
    summarise(MeanAbundance = mean(MeanRelativeAbundance)) %>%
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
                  y = MeanRelativeAbundance, fill = Taxon_grouped)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray", alpha = 0.7) +
    geom_stratum(alpha = 0.8) +
    scale_fill_manual(values = colors, name = "Genus", breaks = ordered_taxa) +
    facet_wrap(~ Status, scales = "free_y") +
    labs(
      x = "Developmental Stage",
      y = "Mean Relative Abundance"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(angle = 90, face = "bold", hjust = 0.5),
      legend.key.height = unit(0.5, "cm"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0), labels = percent_format()) +
    guides(fill = guide_legend(ncol = 4))

  return(p)
}

faceted_alluvial_plot <- create_faceted_alluvial_plot(genus_data)
print(faceted_alluvial_plot)

ggsave("all_conditions_faceted_alluvial_plot.png", faceted_alluvial_plot,
       width = 14, height = 10, dpi = 300, bg = "white")

create_individual_alluvial_plots <- function(data) {
  alluvial_data <- prepare_alluvial_data(data, min_abundance = 0.01)

  taxa_abundance <- alluvial_data %>%
    group_by(Taxon_grouped) %>%
    summarise(MeanAbundance = mean(MeanRelativeAbundance)) %>%
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

  status_values <- unique(alluvial_data$Status)

  plot_list <- list()

  for (status in status_values) {
    status_data <- alluvial_data %>%
      filter(Status == status)

    p <- ggplot(status_data,
                aes(x = Stage, stratum = Taxon_grouped, alluvium = Taxon_grouped,
                    y = MeanRelativeAbundance, fill = Taxon_grouped)) +
      geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray", alpha = 0.7) +
      geom_stratum(alpha = 0.8) +
      scale_fill_manual(values = colors, name = "Genus", breaks = ordered_taxa) +
      labs(
        x = "Developmental Stage",
        y = "Mean Relative Abundance"
      ) +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm")
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0), labels = percent_format()) +
      guides(fill = guide_legend(ncol = 3))

    plot_list[[status]] <- p

    ggsave(paste0("all_conditions_alluvial_plot_", gsub(" ", "_", tolower(status)), ".png"), p,
           width = 12, height = 8, dpi = 300, bg = "white")
  }

  return(plot_list)
}

individual_alluvial_plots <- create_individual_alluvial_plots(genus_data)

create_combined_alluvial_plot <- function(plot_list) {
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))

  plot_list_no_legend <- lapply(plot_list, function(p) {
    p + theme(legend.position = "none")
  })

  plot_names <- names(plot_list)

  if (length(plot_names) >= 2 && "Lab-reared" %in% plot_names && "Carrion-reared" %in% plot_names) {
    ordered_names <- c("Lab-reared", "Carrion-reared")
    other_names <- setdiff(plot_names, ordered_names)
    ordered_names <- c(ordered_names, other_names)
    ordered_names <- ordered_names[ordered_names %in% plot_names]
    plot_list_no_legend <- plot_list_no_legend[ordered_names]
  }

  if (length(plot_list_no_legend) == 2) {
    combined_plot <- cowplot::plot_grid(plotlist = plot_list_no_legend, ncol = 2)
  } else {
    combined_plot <- cowplot::plot_grid(plotlist = plot_list_no_legend, ncol = 2)
  }

  final_plot <- cowplot::plot_grid(
    combined_plot, legend,
    ncol = 1,
    rel_heights = c(1, 0.2)
  )

  return(final_plot)
}

combined_alluvial_plot <- create_combined_alluvial_plot(individual_alluvial_plots)
ggsave("all_conditions_combined_alluvial_plot.png", combined_alluvial_plot,
       width = 16, height = 12, dpi = 300, bg = "white")

cat("=== ALL CONDITIONS MICROBIOME COMPOSITION ANALYSIS ===\n")
cat("Number of samples:", nsamples(ps_stages), "\n")
cat("Number of taxa:", ntaxa(ps_stages), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")
cat("Rearing conditions included:", paste(available_status, collapse = ", "), "\n")

diversity_by_status_stage <- genus_data %>%
  group_by(Status, Stage) %>%
  summarise(
    Shannon = -sum(RelativeAbundance * log(RelativeAbundance + 1e-10)),
    Simpson = 1 - sum(RelativeAbundance^2),
    Richness = sum(RelativeAbundance > 0),
    .groups = "drop"
  )

cat("\nDiversity metrics by rearing condition and developmental stage:\n")
print(diversity_by_status_stage)

cat("\nTop 15 most abundant genera (used in alluvial plots):\n")
top_genera_summary <- genus_data %>%
  filter(Taxon_grouped != "Other") %>%
  group_by(Taxon_grouped) %>%
  summarise(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance))

print(top_genera_summary)

cat("\nVisualization files created:\n")
cat("- all_conditions_faceted_alluvial_plot.png (alluvial plots faceted by Status)\n")
for (status in unique(genus_data$Status)) {
  cat(paste0("- all_conditions_alluvial_plot_", gsub(" ", "_", tolower(status)), ".png (individual alluvial plot for ", status, ")\n"))
}
cat("- all_conditions_combined_alluvial_plot.png (combined individual alluvial plots with shared legend)\n")
