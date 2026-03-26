library(qiime2R)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(viridis)
library(patchwork)

theme_pub <- function(base_size=14, base_family="Polymath") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5,
                                      family = "Polymath"),
            text = element_text(family = "Polymath"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1), family = "Polymath"),
            axis.title.y = element_text(angle=90, vjust=2, family = "Polymath"),
            axis.title.x = element_text(vjust = -0.2, family = "Polymath"),
            axis.text = element_text(family = "Polymath"),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_text(face="italic", family = "Polymath"),
            legend.text = element_text(family = "Polymath"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", family = "Polymath")
    ))
}

flies_ps <- qza_to_phyloseq(
  features="data/filtered-dada-table-nmnc.qza",
  tree="data/rooted-tree.qza",
  "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
)

flies_metadata <- data.frame(sample_data(flies_ps))
flies_metadata$Sample <- rownames(flies_metadata)

stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")

group_colors <- c(
  "Lab-reared" = "#4E6851",
  "Carcass-reared" = "#B83A2D"
)

observed_richness <- estimate_richness(flies_ps, measures = "Observed")
observed_richness$Sample <- rownames(observed_richness)

shannon_diversity <- estimate_richness(flies_ps, measures = "Shannon")
shannon_diversity$Sample <- rownames(shannon_diversity)

faith_pd <- picante::pd(t(otu_table(flies_ps)), phy_tree(flies_ps), include.root = FALSE)
faith_pd$Sample <- rownames(faith_pd)

observed_richness$Sample <- gsub("\\.", "-", observed_richness$Sample)
shannon_diversity$Sample <- gsub("\\.", "-", shannon_diversity$Sample)
faith_pd$Sample <- gsub("\\.", "-", faith_pd$Sample)

alpha_diversity <- observed_richness %>%
  select(Sample, Observed) %>%
  left_join(shannon_diversity %>% select(Sample, Shannon), by = "Sample") %>%
  left_join(faith_pd %>% select(Sample, PD), by = "Sample") %>%
  left_join(flies_metadata %>% select(Sample, Stage, Group, Class), by = "Sample")

cat("Unique groups in metadata:", paste(unique(flies_metadata$Group), collapse = ", "), "\n")

filtered_alpha <- alpha_diversity %>%
  filter(
    Class == "Insect" &
    !is.na(Stage) &
    Stage %in% stage_order &
    Group %in% c("EGG", "LM")
  ) %>%
  mutate(Stage = factor(Stage, levels = stage_order),
         Group = case_when(
           Group == "EGG" ~ "Lab-reared",
           Group == "LM" ~ "Carcass-reared",
           TRUE ~ Group
         ))

cat("Number of samples by group:\n")
print(table(filtered_alpha$Group))

if (nrow(filtered_alpha) == 0) {
  stop("No samples found")
}

cat("Number of samples found:", nrow(filtered_alpha), "\n")
cat("Stages represented:", paste(unique(filtered_alpha$Stage), collapse = ", "), "\n")
cat("Groups represented:", paste(unique(filtered_alpha$Group), collapse = ", "), "\n")

alpha_long <- filtered_alpha %>%
  pivot_longer(
    cols = c(Observed, Shannon, PD),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed", "Shannon", "PD"),
                   labels = c("Observed ASVs", "Shannon Diversity", "Faith's PD"))
  )

alpha_summary <- alpha_long %>%
  group_by(Stage, Group, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

library(ggsignif)
library(rstatix)

get_significance_data <- function(data, metric_name) {
  test_data <- data %>%
    filter(Metric == metric_name & Stage %in% c("3rd instar", "Pupal", "Adult"))

  results <- data.frame()

  for (stage in c("3rd instar", "Pupal", "Adult")) {
    stage_data <- test_data %>% filter(Stage == stage)

    if (length(unique(stage_data$Group)) < 2) {
      next
    }

    t_test_result <- t.test(Value ~ Group, data = stage_data)

    if (t_test_result$p.value < 0.05) {
      signif <- ""
      if (t_test_result$p.value < 0.001) signif <- "***"
      else if (t_test_result$p.value < 0.01) signif <- "**"
      else if (t_test_result$p.value < 0.05) signif <- "*"

      group_means <- stage_data %>%
        group_by(Group) %>%
        summarise(Mean = mean(Value, na.rm = TRUE), .groups = "drop")

      y_pos <- max(group_means$Mean) * 1.1

      results <- rbind(results, data.frame(
        Stage = stage,
        p.value = t_test_result$p.value,
        annotation = signif,
        y.position = y_pos
      ))
    }
  }

  return(results)
}

create_bar_plot <- function(data, metric_name, raw_data) {
  plot_data <- data %>% filter(Metric == metric_name)

  sig_data <- get_significance_data(raw_data, metric_name)

  p <- ggplot(plot_data, aes(x = Stage, y = Mean, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SE, ymax = Mean + SE),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = metric_name,
      x = "Life Stage",
      y = paste("Mean", metric_name),
      fill = "Rearing Method"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    coord_cartesian(expand = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  if (nrow(sig_data) > 0) {
    for (i in 1:nrow(sig_data)) {
      stage <- sig_data$Stage[i]
      p <- p +
        geom_signif(
          xmin = which(levels(plot_data$Stage) == stage) - 0.2,
          xmax = which(levels(plot_data$Stage) == stage) + 0.2,
          y_position = sig_data$y.position[i],
          annotation = sig_data$annotation[i],
          tip_length = 0.01
        )
    }
  }

  return(p)
}

observed_plot <- create_bar_plot(alpha_summary, "Observed ASVs", alpha_long)
shannon_plot <- create_bar_plot(alpha_summary, "Shannon Diversity", alpha_long)
faith_pd_plot <- create_bar_plot(alpha_summary, "Faith's PD", alpha_long)

combined_plot <- observed_plot / shannon_plot / faith_pd_plot +
  plot_layout(heights = c(1, 1, 1), guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(family = "Polymath", face = "bold", size = 16, hjust = 0.5),
      legend.position = "bottom"
    )
  )

print(combined_plot)

ggsave("figures/alpha_diversity_mean_barplot.pdf", combined_plot, width = 10, height = 15)
ggsave("figures/alpha_diversity_mean_barplot.png", combined_plot, width = 10, height = 15, dpi = 300)

cat("Mean alpha diversity bar plots created and saved to the figures directory.\n")
