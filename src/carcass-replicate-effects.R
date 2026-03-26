library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(picante)
library(scales)
library(gridExtra)
library(grid)
library(patchwork)
library(RColorBrewer)
library(car)
library(lme4)
library(lmerTest)

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

cat("============================================================================\n")
cat("REARING GROUP (CARCASS) EFFECTS ANALYSIS\n")
cat("============================================================================\n\n")

ps <- qiime2R::qza_to_phyloseq(
  features = "filtered-dada-table-nmnc.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
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

cat("\n=== REARING GROUPS (CARCASS REPLICATES) ===\n")
rearing_groups <- unique(sample_data(ps_lab_stages)$RearingGroup)
cat("Rearing groups:", paste(rearing_groups, collapse = ", "), "\n")
sample_data(ps_lab_stages)$RearingGroup <- factor(sample_data(ps_lab_stages)$RearingGroup)

cat("\n=== SAMPLE COUNTS BY STAGE AND REARING GROUP ===\n")
metadata_df <- as.data.frame(sample_data(ps_lab_stages))
sample_counts <- table(metadata_df$Stage, metadata_df$RearingGroup)
print(sample_counts)

cat("\n============================================================================\n")
cat("PART 1: ALPHA DIVERSITY - REARING GROUP EFFECTS\n")
cat("============================================================================\n\n")

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

alpha_diversity$Stage <- factor(alpha_diversity$Stage, levels = stage_order)
alpha_diversity$RearingGroup <- factor(alpha_diversity$RearingGroup)

cat("=== TWO-WAY ANOVA: SHANNON DIVERSITY ===\n")
cat("Testing: Stage, RearingGroup, and Stage × RearingGroup interaction\n\n")

shannon_aov <- aov(Shannon ~ Stage * RearingGroup, data = alpha_diversity)
shannon_anova_results <- Anova(shannon_aov, type = "II")

cat("Type II ANOVA Results:\n")
print(shannon_anova_results)

shannon_ss <- summary(shannon_aov)[[1]]
shannon_ss_total <- sum(shannon_ss$`Sum Sq`)
shannon_eta_sq <- data.frame(
  Factor = rownames(shannon_ss),
  Sum_Sq = shannon_ss$`Sum Sq`,
  Eta_Squared = shannon_ss$`Sum Sq` / shannon_ss_total,
  Partial_Eta_Sq = shannon_ss$`Sum Sq` / (shannon_ss$`Sum Sq` + shannon_ss$`Sum Sq`[length(shannon_ss$`Sum Sq`)])
)
shannon_eta_sq$Partial_Eta_Sq[length(shannon_eta_sq$Partial_Eta_Sq)] <- NA

cat("\nEffect Sizes (Eta-squared):\n")
print(shannon_eta_sq)

cat("\n=== TWO-WAY ANOVA: FAITH'S PHYLOGENETIC DIVERSITY ===\n")
cat("Testing: Stage, RearingGroup, and Stage × RearingGroup interaction\n\n")

faithpd_aov <- aov(FaithPD ~ Stage * RearingGroup, data = alpha_diversity)
faithpd_anova_results <- Anova(faithpd_aov, type = "II")

cat("Type II ANOVA Results:\n")
print(faithpd_anova_results)

faithpd_ss <- summary(faithpd_aov)[[1]]
faithpd_ss_total <- sum(faithpd_ss$`Sum Sq`)
faithpd_eta_sq <- data.frame(
  Factor = rownames(faithpd_ss),
  Sum_Sq = faithpd_ss$`Sum Sq`,
  Eta_Squared = faithpd_ss$`Sum Sq` / faithpd_ss_total,
  Partial_Eta_Sq = faithpd_ss$`Sum Sq` / (faithpd_ss$`Sum Sq` + faithpd_ss$`Sum Sq`[length(faithpd_ss$`Sum Sq`)])
)
faithpd_eta_sq$Partial_Eta_Sq[length(faithpd_eta_sq$Partial_Eta_Sq)] <- NA

cat("\nEffect Sizes (Eta-squared):\n")
print(faithpd_eta_sq)

cat("\n=== MIXED EFFECTS MODELS (RearingGroup as random effect) ===\n")
cat("This treats RearingGroup as a random blocking factor\n\n")

cat("Shannon Diversity Mixed Model:\n")
shannon_mixed <- lmer(Shannon ~ Stage + (1|RearingGroup), data = alpha_diversity)
print(summary(shannon_mixed))
cat("\nANOVA table for fixed effects:\n")
print(anova(shannon_mixed))

cat("\nFaith's PD Mixed Model:\n")
faithpd_mixed <- lmer(FaithPD ~ Stage + (1|RearingGroup), data = alpha_diversity)
print(summary(faithpd_mixed))
cat("\nANOVA table for fixed effects:\n")
print(anova(faithpd_mixed))

cat("\n=== VARIANCE COMPONENTS ===\n")
shannon_var <- as.data.frame(VarCorr(shannon_mixed))
cat("Shannon Diversity:\n")
cat("  RearingGroup variance:", shannon_var$vcov[1], "\n")
cat("  Residual variance:", shannon_var$vcov[2], "\n")
cat("  ICC (proportion of variance due to RearingGroup):",
    shannon_var$vcov[1] / (shannon_var$vcov[1] + shannon_var$vcov[2]), "\n")

faithpd_var <- as.data.frame(VarCorr(faithpd_mixed))
cat("\nFaith's PD:\n")
cat("  RearingGroup variance:", faithpd_var$vcov[1], "\n")
cat("  Residual variance:", faithpd_var$vcov[2], "\n")
cat("  ICC (proportion of variance due to RearingGroup):",
    faithpd_var$vcov[1] / (faithpd_var$vcov[1] + faithpd_var$vcov[2]), "\n")

rearing_colors <- brewer.pal(max(3, length(rearing_groups)), "Set2")[1:length(rearing_groups)]
names(rearing_colors) <- levels(alpha_diversity$RearingGroup)

shannon_rearing_plot <- ggplot(alpha_diversity, aes(x = Stage, y = Shannon, fill = RearingGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = RearingGroup), position = position_jitterdodge(jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  scale_fill_manual(values = rearing_colors) +
  scale_color_manual(values = rearing_colors) +
  labs(
    title = "Shannon Diversity by Developmental Stage and Rearing Group",
    x = "Developmental Stage",
    y = "Shannon Diversity",
    fill = "Rearing Group",
    color = "Rearing Group"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

faithpd_rearing_plot <- ggplot(alpha_diversity, aes(x = Stage, y = FaithPD, fill = RearingGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = RearingGroup), position = position_jitterdodge(jitter.width = 0.2),
             alpha = 0.6, size = 2) +
  scale_fill_manual(values = rearing_colors) +
  scale_color_manual(values = rearing_colors) +
  labs(
    title = "Faith's PD by Developmental Stage and Rearing Group",
    x = "Developmental Stage",
    y = "Faith's Phylogenetic Diversity",
    fill = "Rearing Group",
    color = "Rearing Group"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

alpha_rearing_combined <- shannon_rearing_plot + faithpd_rearing_plot +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave("rearing_group_alpha_diversity.png", alpha_rearing_combined,
       width = 14, height = 7, dpi = 300, bg = "white")

cat("\n============================================================================\n")
cat("PART 2: BETA DIVERSITY - REARING GROUP EFFECTS (PERMANOVA)\n")
cat("============================================================================\n\n")

metadata_df <- data.frame(sample_data(ps_lab_stages))
metadata_df$Stage <- factor(metadata_df$Stage, levels = stage_order)
metadata_df$RearingGroup <- factor(metadata_df$RearingGroup)
rownames(metadata_df) <- sample_names(ps_lab_stages)

dist_bray <- phyloseq::distance(ps_lab_stages, method = "bray")
dist_unifrac <- phyloseq::distance(ps_lab_stages, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps_lab_stages, method = "wunifrac")

perform_permanova_full <- function(dist_matrix, metadata, formula_str) {
  formula <- as.formula(paste("dist_matrix ~", formula_str))

  permanova <- adonis2(formula, data = metadata, permutations = 999, by = "terms")

  return(permanova)
}

cat("=== PERMANOVA: BRAY-CURTIS DISTANCE ===\n")
cat("Model: Distance ~ Stage * RearingGroup\n\n")

permanova_bray_full <- perform_permanova_full(dist_bray, metadata_df, "Stage * RearingGroup")
print(permanova_bray_full)

bray_r2_total <- sum(permanova_bray_full$R2[1:3], na.rm = TRUE)
cat("\nVariance explained:\n")
cat("  Stage:", round(permanova_bray_full$R2[1] * 100, 2), "%\n")
cat("  RearingGroup:", round(permanova_bray_full$R2[2] * 100, 2), "%\n")
cat("  Stage × RearingGroup:", round(permanova_bray_full$R2[3] * 100, 2), "%\n")
cat("  Residual:", round(permanova_bray_full$R2[4] * 100, 2), "%\n")

cat("\n=== PERMANOVA: UNWEIGHTED UNIFRAC DISTANCE ===\n")
cat("Model: Distance ~ Stage * RearingGroup\n\n")

permanova_unifrac_full <- perform_permanova_full(dist_unifrac, metadata_df, "Stage * RearingGroup")
print(permanova_unifrac_full)

cat("\nVariance explained:\n")
cat("  Stage:", round(permanova_unifrac_full$R2[1] * 100, 2), "%\n")
cat("  RearingGroup:", round(permanova_unifrac_full$R2[2] * 100, 2), "%\n")
cat("  Stage × RearingGroup:", round(permanova_unifrac_full$R2[3] * 100, 2), "%\n")
cat("  Residual:", round(permanova_unifrac_full$R2[4] * 100, 2), "%\n")

cat("\n=== PERMANOVA: WEIGHTED UNIFRAC DISTANCE ===\n")
cat("Model: Distance ~ Stage * RearingGroup\n\n")

permanova_wunifrac_full <- perform_permanova_full(dist_wunifrac, metadata_df, "Stage * RearingGroup")
print(permanova_wunifrac_full)

cat("\nVariance explained:\n")
cat("  Stage:", round(permanova_wunifrac_full$R2[1] * 100, 2), "%\n")
cat("  RearingGroup:", round(permanova_wunifrac_full$R2[2] * 100, 2), "%\n")
cat("  Stage × RearingGroup:", round(permanova_wunifrac_full$R2[3] * 100, 2), "%\n")
cat("  Residual:", round(permanova_wunifrac_full$R2[4] * 100, 2), "%\n")

cat("\n=== PERMANOVA WITH REARING GROUP AS STRATA (BLOCKING FACTOR) ===\n")
cat("This approach tests Stage effect while controlling for RearingGroup\n\n")

cat("Bray-Curtis (stratified by RearingGroup):\n")
permanova_bray_strata <- adonis2(dist_bray ~ Stage,
                                 data = metadata_df,
                                 permutations = 999,
                                 strata = metadata_df$RearingGroup,
                                 by = "terms")
print(permanova_bray_strata)

cat("\nUnweighted UniFrac (stratified by RearingGroup):\n")
permanova_unifrac_strata <- adonis2(dist_unifrac ~ Stage,
                                    data = metadata_df,
                                    permutations = 999,
                                    strata = metadata_df$RearingGroup,
                                    by = "terms")
print(permanova_unifrac_strata)

cat("\nWeighted UniFrac (stratified by RearingGroup):\n")
permanova_wunifrac_strata <- adonis2(dist_wunifrac ~ Stage,
                                     data = metadata_df,
                                     permutations = 999,
                                     strata = metadata_df$RearingGroup,
                                     by = "terms")
print(permanova_wunifrac_strata)

cat("\n=== BETADISPER: HOMOGENEITY OF DISPERSIONS ===\n")
cat("Testing whether dispersion differs between RearingGroups\n\n")

cat("Bray-Curtis - Dispersion by RearingGroup:\n")
betadisp_bray_rearing <- betadisper(dist_bray, metadata_df$RearingGroup)
print(anova(betadisp_bray_rearing))
print(permutest(betadisp_bray_rearing, permutations = 999))

cat("\nUnweighted UniFrac - Dispersion by RearingGroup:\n")
betadisp_unifrac_rearing <- betadisper(dist_unifrac, metadata_df$RearingGroup)
print(anova(betadisp_unifrac_rearing))
print(permutest(betadisp_unifrac_rearing, permutations = 999))

cat("\nWeighted UniFrac - Dispersion by RearingGroup:\n")
betadisp_wunifrac_rearing <- betadisper(dist_wunifrac, metadata_df$RearingGroup)
print(anova(betadisp_wunifrac_rearing))
print(permutest(betadisp_wunifrac_rearing, permutations = 999))

cat("\n============================================================================\n")
cat("SUMMARY: REARING GROUP EFFECTS\n")
cat("============================================================================\n\n")

summary_table <- data.frame(
  Metric = c("Bray-Curtis", "Bray-Curtis", "Bray-Curtis",
             "Unweighted UniFrac", "Unweighted UniFrac", "Unweighted UniFrac",
             "Weighted UniFrac", "Weighted UniFrac", "Weighted UniFrac"),
  Factor = rep(c("Stage", "RearingGroup", "Stage × RearingGroup"), 3),
  R2 = c(permanova_bray_full$R2[1:3],
         permanova_unifrac_full$R2[1:3],
         permanova_wunifrac_full$R2[1:3]),
  F_value = c(permanova_bray_full$F[1:3],
              permanova_unifrac_full$F[1:3],
              permanova_wunifrac_full$F[1:3]),
  p_value = c(permanova_bray_full$`Pr(>F)`[1:3],
              permanova_unifrac_full$`Pr(>F)`[1:3],
              permanova_wunifrac_full$`Pr(>F)`[1:3])
)

summary_table$R2_percent <- round(summary_table$R2 * 100, 2)
summary_table$Significance <- ifelse(summary_table$p_value < 0.001, "***",
                                     ifelse(summary_table$p_value < 0.01, "**",
                                            ifelse(summary_table$p_value < 0.05, "*", "ns")))

cat("=== PERMANOVA SUMMARY TABLE ===\n")
print(summary_table[, c("Metric", "Factor", "R2_percent", "F_value", "p_value", "Significance")])

alpha_summary <- data.frame(
  Metric = c("Shannon", "Shannon", "Shannon",
             "Faith's PD", "Faith's PD", "Faith's PD"),
  Factor = rep(c("Stage", "RearingGroup", "Stage × RearingGroup"), 2),
  F_value = c(shannon_anova_results$`F value`[1:3],
              faithpd_anova_results$`F value`[1:3]),
  p_value = c(shannon_anova_results$`Pr(>F)`[1:3],
              faithpd_anova_results$`Pr(>F)`[1:3])
)

alpha_summary$Significance <- ifelse(alpha_summary$p_value < 0.001, "***",
                                     ifelse(alpha_summary$p_value < 0.01, "**",
                                            ifelse(alpha_summary$p_value < 0.05, "*", "ns")))

cat("\n=== ALPHA DIVERSITY ANOVA SUMMARY ===\n")
print(alpha_summary)
