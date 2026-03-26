library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)
library(picante)
library(scales)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(patchwork)

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
            legend.title = element_text(face="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
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

calculate_alpha_diversity <- function(ps_obj) {
  shannon <- estimate_richness(ps_obj, measures = "Shannon")

  tree <- phy_tree(ps_obj)

  otu_mat <- as(otu_table(ps_obj), "matrix")

  if(taxa_are_rows(ps_obj)) {
    otu_mat <- t(otu_mat)
  }

  pa_mat <- 1 * (otu_mat > 0)

  tree_tips <- tree$tip.label
  otu_taxa <- colnames(pa_mat)

  common_taxa <- intersect(tree_tips, otu_taxa)

  tree_pruned <- ape::keep.tip(tree, common_taxa)
  pa_mat_pruned <- pa_mat[, common_taxa]

  faith_pd_result <- pd(pa_mat_pruned, tree_pruned, include.root = FALSE)

  faith_pd_df <- data.frame(
    Sample = rownames(faith_pd_result),
    FaithPD = faith_pd_result$PD,
    stringsAsFactors = FALSE
  )

  shannon_df <- data.frame(
    Sample = rownames(shannon),
    Shannon = shannon$Shannon,
    stringsAsFactors = FALSE
  )

  alpha_div <- merge(shannon_df, faith_pd_df, by = "Sample")

  sample_data_df <- as.data.frame(sample_data(ps_obj))
  alpha_div <- merge(alpha_div, sample_data_df, by.x = "Sample", by.y = "row.names")

  return(alpha_div)
}

alpha_diversity <- calculate_alpha_diversity(ps_stages)

alpha_by_status_stage <- alpha_diversity %>%
  group_by(Status, Stage) %>%
  summarise(
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    FaithPD_mean = mean(FaithPD, na.rm = TRUE),
    FaithPD_se = sd(FaithPD, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

alpha_by_status <- alpha_diversity %>%
  group_by(Status) %>%
  summarise(
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    FaithPD_mean = mean(FaithPD, na.rm = TRUE),
    FaithPD_se = sd(FaithPD, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

alpha_by_stage <- alpha_diversity %>%
  group_by(Stage) %>%
  summarise(
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    FaithPD_mean = mean(FaithPD, na.rm = TRUE),
    FaithPD_se = sd(FaithPD, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("=== ALPHA DIVERSITY BY STATUS AND DEVELOPMENTAL STAGE ===\n")
print(alpha_by_status_stage)

cat("\n=== ALPHA DIVERSITY BY STATUS ===\n")
print(alpha_by_status)

cat("\n=== ALPHA DIVERSITY BY DEVELOPMENTAL STAGE ===\n")
print(alpha_by_stage)

perform_wilcoxon_tests <- function(data, metric, group_var) {
  groups <- unique(data[[group_var]])

  results <- data.frame(
    Group1 = character(),
    Group2 = character(),
    p_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:(length(groups) - 1)) {
    for (j in (i+1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]

      data1 <- data[data[[group_var]] == group1, metric]
      data2 <- data[data[[group_var]] == group2, metric]

      wilcox_test <- wilcox.test(data1, data2)
      p_value <- wilcox_test$p.value

      sig <- ""
      if (p_value < 0.001) sig <- "***"
      else if (p_value < 0.01) sig <- "**"
      else if (p_value < 0.05) sig <- "*"
      else sig <- "ns"

      results <- rbind(results, data.frame(
        Group1 = group1,
        Group2 = group2,
        p_value = p_value,
        significance = sig,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}

perform_status_tests_by_stage <- function(data, metric) {
  stage_groups <- unique(data$Stage)

  results <- data.frame(
    Stage = character(),
    Status1 = character(),
    Status2 = character(),
    p_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )

  for (stage in stage_groups) {
    stage_data <- data[data$Stage == stage, ]

    status_groups <- unique(stage_data$Status)

    for (i in 1:(length(status_groups) - 1)) {
      for (j in (i+1):length(status_groups)) {
        status1 <- status_groups[i]
        status2 <- status_groups[j]

        data1 <- stage_data[stage_data$Status == status1, metric]
        data2 <- stage_data[stage_data$Status == status2, metric]

        if (length(data1) < 2 || length(data2) < 2) {
          next
        }

        wilcox_test <- wilcox.test(data1, data2)
        p_value <- wilcox_test$p.value

        sig <- ""
        if (p_value < 0.001) sig <- "***"
        else if (p_value < 0.01) sig <- "**"
        else if (p_value < 0.05) sig <- "*"
        else sig <- "ns"

        results <- rbind(results, data.frame(
          Stage = stage,
          Status1 = status1,
          Status2 = status2,
          p_value = p_value,
          significance = sig,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results)
}

shannon_status_tests <- perform_wilcoxon_tests(alpha_diversity, "Shannon", "Status")
cat("=== SHANNON DIVERSITY TESTS BY STATUS ===\n")
print(shannon_status_tests)

faithpd_status_tests <- perform_wilcoxon_tests(alpha_diversity, "FaithPD", "Status")
cat("\n=== FAITH'S PD TESTS BY STATUS ===\n")
print(faithpd_status_tests)

shannon_stage_tests <- perform_wilcoxon_tests(alpha_diversity, "Shannon", "Stage")
cat("\n=== SHANNON DIVERSITY TESTS BY STAGE ===\n")
print(shannon_stage_tests)

faithpd_stage_tests <- perform_wilcoxon_tests(alpha_diversity, "FaithPD", "Stage")
cat("\n=== FAITH'S PD TESTS BY STAGE ===\n")
print(faithpd_stage_tests)

shannon_status_by_stage_tests <- perform_status_tests_by_stage(alpha_diversity, "Shannon")
cat("\n=== SHANNON DIVERSITY TESTS BETWEEN STATUS GROUPS BY STAGE ===\n")
print(shannon_status_by_stage_tests)

faithpd_status_by_stage_tests <- perform_status_tests_by_stage(alpha_diversity, "FaithPD")
cat("\n=== FAITH'S PD TESTS BETWEEN STATUS GROUPS BY STAGE ===\n")
print(faithpd_status_by_stage_tests)

shannon_status_kruskal <- kruskal.test(Shannon ~ Status, data = alpha_diversity)
shannon_stage_kruskal <- kruskal.test(Shannon ~ Stage, data = alpha_diversity)
faithpd_status_kruskal <- kruskal.test(FaithPD ~ Status, data = alpha_diversity)
faithpd_stage_kruskal <- kruskal.test(FaithPD ~ Stage, data = alpha_diversity)

cat("\n=== KRUSKAL-WALLIS TESTS ===\n")
cat("Shannon diversity by Status: chi-squared =", shannon_status_kruskal$statistic,
    ", df =", shannon_status_kruskal$parameter,
    ", p-value =", shannon_status_kruskal$p.value, "\n")
cat("Shannon diversity by Stage: chi-squared =", shannon_stage_kruskal$statistic,
    ", df =", shannon_stage_kruskal$parameter,
    ", p-value =", shannon_stage_kruskal$p.value, "\n")
cat("Faith's PD by Status: chi-squared =", faithpd_status_kruskal$statistic,
    ", df =", faithpd_status_kruskal$parameter,
    ", p-value =", faithpd_status_kruskal$p.value, "\n")
cat("Faith's PD by Stage: chi-squared =", faithpd_stage_kruskal$statistic,
    ", df =", faithpd_stage_kruskal$parameter,
    ", p-value =", faithpd_stage_kruskal$p.value, "\n")

shannon_color <- "#1f77b4"
faith_pd_color <- "#ff7f0e"

create_boxplot <- function(data, metric, metric_name, color, sig_tests) {
  p <- ggplot(data, aes_string(x = "Stage", y = metric, fill = "Status")) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2,
                 position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = status_colors) +
    labs(
      x = "Developmental Stage",
      y = metric_name,
      fill = "Rearing Condition"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.margin = unit(c(40, 20, 20, 20), "pt")
    )

  sig_tests_filtered <- sig_tests[sig_tests$significance != "ns", ]

  if (nrow(sig_tests_filtered) > 0) {
    y_range <- layer_scales(p)$y$range$range
    if (is.null(y_range)) {
      y_range <- range(data[[metric]], na.rm = TRUE)
    }
    y_range_diff <- diff(y_range)

    stage_groups <- unique(sig_tests_filtered$Stage)

    for (i in seq_along(stage_groups)) {
      stage <- stage_groups[i]

      stage_tests <- sig_tests_filtered[sig_tests_filtered$Stage == stage, ]

      for (j in 1:nrow(stage_tests)) {
        status1 <- stage_tests$Status1[j]
        status2 <- stage_tests$Status2[j]

        stage_pos <- which(levels(data$Stage) == stage)

        status1_pos <- which(levels(data$Status) == status1)
        status2_pos <- which(levels(data$Status) == status2)

        group_width <- 0.8

        box_width <- group_width / length(unique(data$Status))

        dodge_width1 <- (status1_pos - 1) * box_width - (group_width - box_width) / 2
        dodge_width2 <- (status2_pos - 1) * box_width - (group_width - box_width) / 2

        x1 <- stage_pos + dodge_width1
        x2 <- stage_pos + dodge_width2

        y_values <- c(
          data[data$Status == status1 & data$Stage == stage, metric],
          data[data$Status == status2 & data$Stage == stage, metric]
        )

        bar_height <- max(y_values, na.rm = TRUE) + 0.1 * y_range_diff

        bar_offset <- (i - 1) * 0.05 * y_range_diff

        p <- p +
          annotate("segment",
                  x = x1,
                  xend = x2,
                  y = bar_height + bar_offset,
                  yend = bar_height + bar_offset,
                  color = "black",
                  size = 0.5) +
          annotate("segment",
                  x = x1,
                  xend = x1,
                  y = bar_height + bar_offset - 0.02 * y_range_diff,
                  yend = bar_height + bar_offset,
                  color = "black",
                  size = 0.5) +
          annotate("segment",
                  x = x2,
                  xend = x2,
                  y = bar_height + bar_offset - 0.02 * y_range_diff,
                  yend = bar_height + bar_offset,
                  color = "black",
                  size = 0.5) +
          annotate("text",
                  x = (x1 + x2) / 2,
                  y = bar_height + bar_offset + 0.02 * y_range_diff,
                  label = stage_tests$significance[j],
                  color = "black")
      }
    }
  }

  return(p)
}

create_line_plot <- function(data, metric_mean, metric_se, metric_name, color) {
  p <- ggplot(data, aes(x = Stage, y = .data[[metric_mean]], color = Status, group = Status)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = .data[[metric_mean]] - .data[[metric_se]],
                      ymax = .data[[metric_mean]] + .data[[metric_se]]),
                  width = 0.2) +
    scale_color_manual(values = status_colors) +
    labs(
      title = paste(metric_name, "by Rearing Condition and Developmental Stage"),
      x = "Developmental Stage",
      y = metric_name,
      color = "Rearing Condition"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  return(p)
}

create_status_bar_plot <- function(data, metric_mean, metric_se, metric_name, color) {
  p <- ggplot(data, aes(x = Status, y = .data[[metric_mean]], fill = Status)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_errorbar(aes(ymin = .data[[metric_mean]] - .data[[metric_se]],
                      ymax = .data[[metric_mean]] + .data[[metric_se]]),
                  width = 0.2) +
    scale_fill_manual(values = status_colors) +
    labs(
      title = paste(metric_name, "by Rearing Condition"),
      x = "Rearing Condition",
      y = metric_name
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(p)
}

create_stage_bar_plot <- function(data, metric_mean, metric_se, metric_name, color) {
  p <- ggplot(data, aes(x = Stage, y = .data[[metric_mean]], fill = Stage)) +
    geom_bar(stat = "identity", alpha = 0.7, fill = color) +
    geom_errorbar(aes(ymin = .data[[metric_mean]] - .data[[metric_se]],
                      ymax = .data[[metric_mean]] + .data[[metric_se]]),
                  width = 0.2) +
    labs(
      title = paste(metric_name, "by Developmental Stage"),
      x = "Developmental Stage",
      y = metric_name
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(p)
}

shannon_boxplot <- create_boxplot(alpha_diversity, "Shannon", "Shannon Diversity", shannon_color, shannon_status_by_stage_tests)
faithpd_boxplot <- create_boxplot(alpha_diversity, "FaithPD", "Faith's Phylogenetic Diversity", faith_pd_color, faithpd_status_by_stage_tests)

shannon_line_plot <- create_line_plot(alpha_by_status_stage, "Shannon_mean", "Shannon_se", "Shannon Diversity", shannon_color)
faithpd_line_plot <- create_line_plot(alpha_by_status_stage, "FaithPD_mean", "FaithPD_se", "Faith's Phylogenetic Diversity", faith_pd_color)

shannon_status_plot <- create_status_bar_plot(alpha_by_status, "Shannon_mean", "Shannon_se", "Shannon Diversity", shannon_color)
faithpd_status_plot <- create_status_bar_plot(alpha_by_status, "FaithPD_mean", "FaithPD_se", "Faith's Phylogenetic Diversity", faith_pd_color)

shannon_stage_plot <- create_stage_bar_plot(alpha_by_stage, "Shannon_mean", "Shannon_se", "Shannon Diversity", shannon_color)
faithpd_stage_plot <- create_stage_bar_plot(alpha_by_stage, "FaithPD_mean", "FaithPD_se", "Faith's Phylogenetic Diversity", faith_pd_color)

combined_boxplots <- shannon_boxplot / faithpd_boxplot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_line_plots <- shannon_line_plot / faithpd_line_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_status_plots <- shannon_status_plot / faithpd_status_plot

combined_stage_plots <- shannon_stage_plot / faithpd_stage_plot

ggsave("all_conditions_alpha_boxplots.png", combined_boxplots,
       width = 6, height = 18, dpi = 300, bg = "white")
ggsave("all_conditions_alpha_line_plots.png", combined_line_plots,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave("all_conditions_alpha_by_status.png", combined_status_plots,
       width = 8, height = 10, dpi = 300, bg = "white")
ggsave("all_conditions_alpha_by_stage.png", combined_stage_plots,
       width = 8, height = 10, dpi = 300, bg = "white")

ggsave("all_conditions_shannon_boxplot.png", shannon_boxplot,
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_faithpd_boxplot.png", faithpd_boxplot,
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_shannon_line_plot.png", shannon_line_plot,
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave("all_conditions_faithpd_line_plot.png", faithpd_line_plot,
       width = 12, height = 6, dpi = 300, bg = "white")

cat("=== ALL CONDITIONS ALPHA DIVERSITY ANALYSIS ===\n")
cat("Number of samples:", nrow(alpha_diversity), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")
cat("Rearing conditions included:", paste(available_status, collapse = ", "), "\n")

cat("\nAlpha diversity metrics by rearing condition and developmental stage:\n")
print(alpha_by_status_stage)

cat("\nVisualization files created:\n")
cat("- all_conditions_alpha_boxplots.png (boxplots of Shannon and Faith's PD by Status and Stage)\n")
cat("- all_conditions_alpha_line_plots.png (line plots of Shannon and Faith's PD by Status and Stage)\n")
cat("- all_conditions_alpha_by_status.png (bar plots of Shannon and Faith's PD by Status)\n")
cat("- all_conditions_alpha_by_stage.png (bar plots of Shannon and Faith's PD by Stage)\n")
cat("- all_conditions_shannon_boxplot.png (boxplot of Shannon diversity by Status and Stage)\n")
cat("- all_conditions_faithpd_boxplot.png (boxplot of Faith's PD by Status and Stage)\n")
cat("- all_conditions_shannon_line_plot.png (line plot of Shannon diversity by Status and Stage)\n")
cat("- all_conditions_faithpd_line_plot.png (line plot of Faith's PD by Status and Stage)\n")
