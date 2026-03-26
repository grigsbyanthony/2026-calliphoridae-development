library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(caret)
library(ROCR)
library(RColorBrewer)
library(gridExtra)

theme_pub <- function(base_size=14, base_family="helvetica") {
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
  features = "data/filtered-dada-table-nmnc.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
)

ps_insects <- subset_samples(ps, Class == "Insect")

cat("=== DATASET SUMMARY ===\n")
cat("Total samples:", nsamples(ps_insects), "\n")
cat("Total ASVs:", ntaxa(ps_insects), "\n")

target_stages <- c("3rd instar", "Pupal", "Adult")
ps_stages <- subset_samples(ps_insects, Stage %in% target_stages)

cat("\n=== STATUS VALUES IN THE DATASET ===\n")
status_values <- unique(sample_data(ps_stages)$Status)
cat("Status values:", paste(status_values, collapse = ", "), "\n")

prepare_rf_data <- function(ps_obj) {
  otu_table <- as.data.frame(t(otu_table(ps_obj)))
  sample_data <- as.data.frame(sample_data(ps_obj))

  if(!all(rownames(otu_table) == rownames(sample_data))) {
    stop("Sample IDs in OTU table and sample data do not match")
  }

  sample_data$Status <- factor(sample_data$Status,
                              levels = c("Lab Reared", "Carrion Reared"))

  rf_data <- cbind(otu_table, Status = sample_data$Status)

  return(rf_data)
}

run_random_forest <- function(data, title, ntree = 500, seed = 123) {
  set.seed(seed)

  train_indices <- createDataPartition(data$Status, p = 0.7, list = FALSE)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]

  rf_model <- randomForest(
    Status ~ .,
    data = train_data,
    ntree = ntree,
    importance = TRUE
  )

  predictions <- predict(rf_model, test_data)
  conf_matrix <- confusionMatrix(predictions, test_data$Status)

  rf_probs <- predict(rf_model, test_data, type = "prob")

  if("Carrion Reared" %in% colnames(rf_probs)) {
    pred_obj <- prediction(rf_probs[, "Carrion Reared"], test_data$Status)
    perf_obj <- performance(pred_obj, "tpr", "fpr")
    auc_obj <- performance(pred_obj, "auc")
    auc_value <- auc_obj@y.values[[1]]

    roc_data <- data.frame(
      FPR = perf_obj@x.values[[1]],
      TPR = perf_obj@y.values[[1]]
    )
  } else {
    cat("Warning: Could not calculate proper ROC curve (possibly only one class in predictions)\n")
    roc_data <- NULL
    auc_value <- NA
  }

  var_importance <- importance(rf_model)
  var_importance_df <- as.data.frame(var_importance)
  var_importance_df$ASV <- rownames(var_importance_df)

  var_importance_df <- var_importance_df %>%
    arrange(desc(MeanDecreaseGini))

  top_vars <- head(var_importance_df, 20)

  tax_table_df <- as.data.frame(tax_table(ps_stages))
  tax_table_df <- tax_table_df %>%
    mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
    mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

  top_vars$Phylum <- tax_table_df[top_vars$ASV, "Phylum"]
  top_vars$Family <- tax_table_df[top_vars$ASV, "Family"]
  top_vars$Genus <- tax_table_df[top_vars$ASV, "Genus"]

  importance_plot <- ggplot(top_vars, aes(x = reorder(ASV, MeanDecreaseGini), y = MeanDecreaseGini, fill = Phylum)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste("Top 20 Important ASVs -", title),
      x = "ASV",
      y = "Mean Decrease in Gini Index"
    ) +
    theme_pub() +
    theme(legend.position = "right")

  if(!is.null(roc_data)) {
    roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      labs(
        title = paste("ROC Curve -", title),
        subtitle = paste("AUC =", round(auc_value, 3)),
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      theme_pub()
  } else {
    roc_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "ROC curve could not be calculated\n(possibly only one class in predictions)") +
      labs(
        title = paste("ROC Curve -", title),
        subtitle = "AUC = NA"
      ) +
      theme_pub() +
      theme(panel.border = element_rect(color = "gray", fill = NA))
  }

  return(list(
    model = rf_model,
    confusion_matrix = conf_matrix,
    auc = auc_value,
    importance = top_vars,
    importance_plot = importance_plot,
    roc_plot = roc_plot
  ))
}

all_data <- prepare_rf_data(ps_stages)
cat("\n=== PREPARED DATA FOR ALL SAMPLES ===\n")
cat("Number of samples:", nrow(all_data), "\n")
cat("Number of features:", ncol(all_data) - 1, "\n")
cat("Class distribution:\n")
print(table(all_data$Status))

cat("\n=== RUNNING RANDOM FOREST FOR ALL SAMPLES ===\n")
all_results <- run_random_forest(all_data, "All Stages")

cat("\n=== RANDOM FOREST RESULTS FOR ALL SAMPLES ===\n")
print(all_results$confusion_matrix)
cat("AUC:", all_results$auc, "\n")

stages <- target_stages
stage_results <- list()

for(stage in stages) {
  cat(paste("\n=== RUNNING RANDOM FOREST FOR", stage, "===\n"))

  ps_stage <- subset_samples(ps_stages, Stage == stage)

  status_counts <- table(sample_data(ps_stage)$Status)
  if(length(status_counts) < 2 || any(status_counts < 5)) {
    cat("Not enough samples for both conditions in stage:", stage, "\n")
    cat("Status counts:", paste(names(status_counts), status_counts, sep = ": ", collapse = ", "), "\n")
    next
  }

  stage_data <- prepare_rf_data(ps_stage)
  cat("Number of samples:", nrow(stage_data), "\n")
  cat("Class distribution:\n")
  print(table(stage_data$Status))

  stage_results[[stage]] <- run_random_forest(stage_data, stage)

  cat(paste("\n=== RANDOM FOREST RESULTS FOR", stage, "===\n"))
  print(stage_results[[stage]]$confusion_matrix)
  cat("AUC:", stage_results[[stage]]$auc, "\n")
}

roc_plots <- list()
if(!is.null(all_results$roc_plot)) {
  roc_plots[[length(roc_plots) + 1]] <- all_results$roc_plot
}

for(stage in names(stage_results)) {
  if(!is.null(stage_results[[stage]]$roc_plot)) {
    roc_plots[[length(roc_plots) + 1]] <- stage_results[[stage]]$roc_plot
  }
}

if(length(roc_plots) > 0) {
  combined_roc_plot <- grid.arrange(
    grobs = roc_plots,
    ncol = min(2, length(roc_plots)),
    top = "ROC Curves for Predicting Rearing Condition"
  )

  ggsave("random_forest_roc_curves.png", combined_roc_plot, width = 12, height = 10, dpi = 300)
} else {
  cat("Warning: No valid ROC curves to plot\n")
}

importance_plots <- list(all_results$importance_plot)
for(stage in names(stage_results)) {
  importance_plots[[length(importance_plots) + 1]] <- stage_results[[stage]]$importance_plot
}

combined_importance_plot <- grid.arrange(
  grobs = importance_plots,
  ncol = 2,
  top = "Important ASVs for Predicting Rearing Condition"
)

ggsave("random_forest_importance.png", combined_importance_plot, width = 16, height = 14, dpi = 300)

summary_results <- data.frame(
  Analysis = "All Stages",
  Accuracy = all_results$confusion_matrix$overall["Accuracy"],
  Sensitivity = all_results$confusion_matrix$byClass["Sensitivity"],
  Specificity = all_results$confusion_matrix$byClass["Specificity"],
  AUC = all_results$auc,
  stringsAsFactors = FALSE
)

for(stage in names(stage_results)) {
  summary_results <- rbind(summary_results, data.frame(
    Analysis = stage,
    Accuracy = stage_results[[stage]]$confusion_matrix$overall["Accuracy"],
    Sensitivity = stage_results[[stage]]$confusion_matrix$byClass["Sensitivity"],
    Specificity = stage_results[[stage]]$confusion_matrix$byClass["Specificity"],
    AUC = stage_results[[stage]]$auc,
    stringsAsFactors = FALSE
  ))
}

cat("\n=== SUMMARY OF RANDOM FOREST RESULTS ===\n")
print(summary_results)

accuracy_plot <- ggplot(summary_results, aes(x = Analysis, y = Accuracy, fill = Analysis)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Accuracy, 3)), vjust = -0.5) +
  labs(
    title = "Random Forest Accuracy by Life Stage",
    x = "Analysis",
    y = "Accuracy"
  ) +
  ylim(0, 1.1) +
  theme_pub() +
  theme(legend.position = "none")

ggsave("random_forest_accuracy.png", accuracy_plot, width = 10, height = 6, dpi = 300)

top_asvs_all <- all_results$importance %>% head(10) %>% select(ASV, Genus, MeanDecreaseGini)
colnames(top_asvs_all)[3] <- "All_Stages"

combined_importance <- top_asvs_all %>% select(ASV, Genus)
combined_importance$All_Stages <- top_asvs_all$All_Stages

for(stage in names(stage_results)) {
  top_asvs_stage <- stage_results[[stage]]$importance %>% head(10) %>% select(ASV, MeanDecreaseGini)
  colnames(top_asvs_stage)[2] <- stage

  combined_importance <- left_join(combined_importance, top_asvs_stage, by = "ASV")
}

combined_importance[is.na(combined_importance)] <- 0

importance_long <- combined_importance %>%
  pivot_longer(cols = -c(ASV, Genus), names_to = "Analysis", values_to = "Importance")

heatmap_plot <- ggplot(importance_long, aes(x = Analysis, y = paste(ASV, Genus, sep = " - "), fill = Importance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Importance of Top ASVs Across Analyses",
    x = "Analysis",
    y = "ASV - Genus"
  ) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

ggsave("random_forest_importance_heatmap.png", heatmap_plot, width = 12, height = 10, dpi = 300)

cat("\n=== RANDOM FOREST ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. random_forest_roc_curves.png - ROC curves for all analyses\n")
cat("2. random_forest_importance.png - Important ASVs for all analyses\n")
cat("3. random_forest_accuracy.png - Accuracy by life stage\n")
cat("4. random_forest_importance_heatmap.png - Heatmap of important ASVs across analyses\n")
