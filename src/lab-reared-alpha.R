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

alpha_diversity <- calculate_alpha_diversity(ps_lab_stages)

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

stage_order <- c("1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
actual_stages <- unique(alpha_by_stage$Stage)
stage_order <- stage_order[stage_order %in% actual_stages]

alpha_by_stage$Stage <- factor(alpha_by_stage$Stage, levels = stage_order)
alpha_diversity$Stage <- factor(alpha_diversity$Stage, levels = stage_order)

cat("=== ALPHA DIVERSITY BY DEVELOPMENTAL STAGE ===\n")
print(alpha_by_stage)

perform_pairwise_tests <- function(data, metric) {
  stages <- levels(data$Stage)

  results <- data.frame(
    Stage1 = character(),
    Stage2 = character(),
    p_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:(length(stages) - 1)) {
    stage1 <- stages[i]
    stage2 <- stages[i + 1]

    data1 <- data[data$Stage == stage1, metric]
    data2 <- data[data$Stage == stage2, metric]

    wilcox_test <- wilcox.test(data1, data2)
    p_value <- wilcox_test$p.value

    sig <- ""
    if (p_value < 0.001) sig <- "***"
    else if (p_value < 0.01) sig <- "**"
    else if (p_value < 0.05) sig <- "*"
    else sig <- "ns"

    results <- rbind(results, data.frame(
      Stage1 = stage1,
      Stage2 = stage2,
      p_value = p_value,
      significance = sig,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

shannon_tests <- perform_pairwise_tests(alpha_diversity, "Shannon")
cat("=== SHANNON DIVERSITY PAIRWISE TESTS ===\n")
print(shannon_tests)

faithpd_tests <- perform_pairwise_tests(alpha_diversity, "FaithPD")
cat("=== FAITH'S PD PAIRWISE TESTS ===\n")
print(faithpd_tests)

shannon_color <- "#1f77b4"
faith_pd_color <- "#ff7f0e"

create_dual_axis_plot <- function(data, shannon_tests, faithpd_tests) {
  shannon_range <- range(data$Shannon_mean)
  faith_pd_range <- range(data$FaithPD_mean)
  scale_factor <- diff(shannon_range) / diff(faith_pd_range)

  data$FaithPD_scaled <- (data$FaithPD_mean - faith_pd_range[1]) * scale_factor + shannon_range[1]
  data$FaithPD_se_scaled <- data$FaithPD_se * scale_factor

  p <- ggplot(data, aes(x = Stage)) +
    geom_ribbon(aes(ymin = Shannon_mean - Shannon_se,
                    ymax = Shannon_mean + Shannon_se,
                    fill = "Shannon"),
                alpha = 0.2) +
    geom_line(aes(y = Shannon_mean, group = 1, color = "Shannon"), size = 1) +
    geom_point(aes(y = Shannon_mean, color = "Shannon"), size = 3) +
    geom_errorbar(aes(ymin = Shannon_mean - Shannon_se,
                      ymax = Shannon_mean + Shannon_se,
                      color = "Shannon"),
                  width = 0.2) +
    geom_ribbon(aes(ymin = FaithPD_scaled - FaithPD_se_scaled,
                    ymax = FaithPD_scaled + FaithPD_se_scaled,
                    fill = "Faith's PD"),
                alpha = 0.2) +
    geom_line(aes(y = FaithPD_scaled, group = 1, color = "Faith's PD"), size = 1) +
    geom_point(aes(y = FaithPD_scaled, color = "Faith's PD"), size = 3) +
    geom_errorbar(aes(ymin = FaithPD_scaled - FaithPD_se_scaled,
                      ymax = FaithPD_scaled + FaithPD_se_scaled,
                      color = "Faith's PD"),
                  width = 0.2) +
    scale_y_continuous(
      name = "Shannon Diversity",
      sec.axis = sec_axis(
        ~ (. - shannon_range[1]) / scale_factor + faith_pd_range[1],
        name = "Faith's Phylogenetic Diversity"
      )
    ) +
    scale_color_manual(
      name = "Alpha Diversity Metric",
      values = c("Shannon" = shannon_color, "Faith's PD" = faith_pd_color)
    ) +
    scale_fill_manual(
      name = "Alpha Diversity Metric",
      values = c("Shannon" = shannon_color, "Faith's PD" = faith_pd_color)
    ) +
    labs(
      x = "Developmental Stage"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y.left = element_text(color = shannon_color),
      axis.title.y.right = element_text(color = faith_pd_color),
      legend.position = "bottom",
      plot.margin = margin(t = 40, r = 20, b = 20, l = 20, unit = "pt")
    )

  all_shannon_tests <- data.frame()
  all_faithpd_tests <- data.frame()

  stages <- levels(data$Stage)
  for (i in 1:(length(stages) - 1)) {
    for (j in (i+1):length(stages)) {
      stage1 <- stages[i]
      stage2 <- stages[j]

      data1 <- alpha_diversity[alpha_diversity$Stage == stage1, "Shannon"]
      data2 <- alpha_diversity[alpha_diversity$Stage == stage2, "Shannon"]
      wilcox_test <- wilcox.test(data1, data2)
      p_value <- wilcox_test$p.value
      sig <- ""
      if (p_value < 0.001) sig <- "***"
      else if (p_value < 0.01) sig <- "**"
      else if (p_value < 0.05) sig <- "*"
      else sig <- "ns"

      all_shannon_tests <- rbind(all_shannon_tests, data.frame(
        Stage1 = stage1,
        Stage2 = stage2,
        p_value = p_value,
        significance = sig,
        consecutive = (j == i+1),
        stringsAsFactors = FALSE
      ))

      data1 <- alpha_diversity[alpha_diversity$Stage == stage1, "FaithPD"]
      data2 <- alpha_diversity[alpha_diversity$Stage == stage2, "FaithPD"]
      wilcox_test <- wilcox.test(data1, data2)
      p_value <- wilcox_test$p.value
      sig <- ""
      if (p_value < 0.001) sig <- "***"
      else if (p_value < 0.01) sig <- "**"
      else if (p_value < 0.05) sig <- "*"
      else sig <- "ns"

      all_faithpd_tests <- rbind(all_faithpd_tests, data.frame(
        Stage1 = stage1,
        Stage2 = stage2,
        p_value = p_value,
        significance = sig,
        consecutive = (j == i+1),
        stringsAsFactors = FALSE
      ))
    }
  }

  shannon_tests_to_plot <- all_shannon_tests[all_shannon_tests$consecutive & all_shannon_tests$significance != "ns", ]
  faithpd_tests_to_plot <- all_faithpd_tests[all_faithpd_tests$consecutive & all_faithpd_tests$significance != "ns", ]

  cat("=== ALL SHANNON DIVERSITY PAIRWISE TESTS ===\n")
  print(all_shannon_tests)
  cat("\n=== SHANNON TESTS SHOWN IN PLOT (CONSECUTIVE STAGES ONLY) ===\n")
  print(shannon_tests_to_plot)

  cat("\n=== ALL FAITH'S PD PAIRWISE TESTS ===\n")
  print(all_faithpd_tests)
  cat("\n=== FAITH'S PD TESTS SHOWN IN PLOT (CONSECUTIVE STAGES ONLY) ===\n")
  print(faithpd_tests_to_plot)

  for (i in 1:nrow(shannon_tests_to_plot)) {
    stage1_pos <- which(levels(data$Stage) == shannon_tests_to_plot$Stage1[i])
    stage2_pos <- which(levels(data$Stage) == shannon_tests_to_plot$Stage2[i])

    y1 <- data$Shannon_mean[data$Stage == shannon_tests_to_plot$Stage1[i]]
    y2 <- data$Shannon_mean[data$Stage == shannon_tests_to_plot$Stage2[i]]
    bar_height <- max(y1, y2) + 0.1 * diff(shannon_range)

    p <- p +
      annotate("segment",
              x = stage1_pos,
              xend = stage2_pos,
              y = bar_height,
              yend = bar_height,
              color = shannon_color,
              size = 0.5) +
      annotate("segment",
              x = stage1_pos,
              xend = stage1_pos,
              y = bar_height - 0.02 * diff(shannon_range),
              yend = bar_height,
              color = shannon_color,
              size = 0.5) +
      annotate("segment",
              x = stage2_pos,
              xend = stage2_pos,
              y = bar_height - 0.02 * diff(shannon_range),
              yend = bar_height,
              color = shannon_color,
              size = 0.5) +
      annotate("text",
              x = (stage1_pos + stage2_pos) / 2,
              y = bar_height + 0.02 * diff(shannon_range),
              label = shannon_tests_to_plot$significance[i],
              color = shannon_color)
  }

  for (i in 1:nrow(faithpd_tests_to_plot)) {
    stage1_pos <- which(levels(data$Stage) == faithpd_tests_to_plot$Stage1[i])
    stage2_pos <- which(levels(data$Stage) == faithpd_tests_to_plot$Stage2[i])

    y1 <- data$FaithPD_scaled[data$Stage == faithpd_tests_to_plot$Stage1[i]]
    y2 <- data$FaithPD_scaled[data$Stage == faithpd_tests_to_plot$Stage2[i]]
    bar_height <- max(y1, y2) + 0.2 * diff(shannon_range)

    p <- p +
      annotate("segment",
              x = stage1_pos,
              xend = stage2_pos,
              y = bar_height,
              yend = bar_height,
              color = faith_pd_color,
              size = 0.5) +
      annotate("segment",
              x = stage1_pos,
              xend = stage1_pos,
              y = bar_height - 0.02 * diff(shannon_range),
              yend = bar_height,
              color = faith_pd_color,
              size = 0.5) +
      annotate("segment",
              x = stage2_pos,
              xend = stage2_pos,
              y = bar_height - 0.02 * diff(shannon_range),
              yend = bar_height,
              color = faith_pd_color,
              size = 0.5) +
      annotate("text",
              x = (stage1_pos + stage2_pos) / 2,
              y = bar_height + 0.02 * diff(shannon_range),
              label = faithpd_tests_to_plot$significance[i],
              color = faith_pd_color)
  }

  return(p)
}

alpha_plot <- create_dual_axis_plot(alpha_by_stage, shannon_tests, faithpd_tests)
print(alpha_plot)

ggsave("lab_reared_alpha_diversity.png", alpha_plot,
       width = 12, height = 8, dpi = 300, bg = "white")

create_shannon_plot <- function(data) {
  ribbon_data <- data.frame(
    Stage = data$Stage,
    y = data$Shannon_mean,
    ymin = data$Shannon_mean - data$Shannon_se,
    ymax = data$Shannon_mean + data$Shannon_se
  )

  p <- ggplot() +
    geom_ribbon(data = ribbon_data,
                aes(x = Stage, ymin = ymin, ymax = ymax),
                fill = shannon_color, alpha = 0.5) +
    geom_line(data = data, aes(x = Stage, y = Shannon_mean, group = 1),
              color = shannon_color, size = 1) +
    geom_point(data = data, aes(x = Stage, y = Shannon_mean),
               color = shannon_color, size = 3) +
    labs(
      title = "Shannon Diversity Across Developmental Stages",
      x = "Developmental Stage",
      y = "Shannon Diversity"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

create_faith_pd_plot <- function(data) {
  ribbon_data <- data.frame(
    Stage = data$Stage,
    y = data$FaithPD_mean,
    ymin = data$FaithPD_mean - data$FaithPD_se,
    ymax = data$FaithPD_mean + data$FaithPD_se
  )

  p <- ggplot() +
    geom_ribbon(data = ribbon_data,
                aes(x = Stage, ymin = ymin, ymax = ymax),
                fill = faith_pd_color, alpha = 0.5) +
    geom_line(data = data, aes(x = Stage, y = FaithPD_mean, group = 1),
              color = faith_pd_color, size = 1) +
    geom_point(data = data, aes(x = Stage, y = FaithPD_mean),
               color = faith_pd_color, size = 3) +
    labs(
      title = "Faith's Phylogenetic Diversity Across Developmental Stages",
      x = "Developmental Stage",
      y = "Faith's PD"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

shannon_plot <- create_shannon_plot(alpha_by_stage)
faith_pd_plot <- create_faith_pd_plot(alpha_by_stage)

ggsave("lab_reared_shannon_diversity.png", shannon_plot,
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("lab_reared_faith_pd.png", faith_pd_plot,
       width = 8, height = 6, dpi = 300, bg = "white")

shannon_plot <- shannon_plot +
  scale_fill_manual(values = c(shannon_color))

faith_pd_plot <- faith_pd_plot +
  scale_fill_manual(values = c(faith_pd_color))

library(patchwork)
combined_plot <- shannon_plot + faith_pd_plot + plot_layout(ncol = 2)

ggsave("lab_reared_alpha_diversity_combined.png", combined_plot,
       width = 14, height = 6, dpi = 300, bg = "white")

cat("=== LAB-REARED ALPHA DIVERSITY ANALYSIS (EXCLUDING EGGS) ===\n")
cat("Number of samples:", nrow(alpha_diversity), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")

cat("\nAlpha diversity metrics by developmental stage:\n")
print(alpha_by_stage)

cat("\nVisualization files created:\n")
cat("- lab_reared_alpha_diversity.png (dual y-axis plot with Shannon and Faith's PD)\n")
cat("- lab_reared_shannon_diversity.png (Shannon diversity only)\n")
cat("- lab_reared_faith_pd.png (Faith's PD only)\n")
cat("- lab_reared_alpha_diversity_combined.png (Shannon and Faith's PD side by side)\n")
