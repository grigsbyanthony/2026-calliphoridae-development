cat("Step 1: Loading required libraries...\n")
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

cat("Step 2: Setting up plotting theme...\n")
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

cat("Step 3: Importing QIIME2 artifacts...\n")
ps <- qiime2R::qza_to_phyloseq(
  features = "data/filtered-dada-table-nmnc.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
)

cat("Step 4: Examining the phyloseq object...\n")
cat("Number of samples:", nsamples(ps), "\n")
cat("Number of ASVs:", ntaxa(ps), "\n")

cat("Step 5: Filtering to include only insect samples...\n")
ps_insects <- subset_samples(ps, Class == "Insect")
cat("After filtering for insects - Number of samples:", nsamples(ps_insects), "\n")

cat("Step 6: Filtering to include only specific life stages (3rd instar, Pupal, Adult)...\n")
target_stages <- c("3rd instar", "Pupal", "Adult")
ps_stages <- subset_samples(ps_insects, Stage %in% target_stages)
cat("After filtering for target stages - Number of samples:", nsamples(ps_stages), "\n")

cat("Step 7: Checking rearing conditions in the dataset...\n")
status_values <- unique(sample_data(ps_stages)$Status)
cat("Status values (rearing conditions):", paste(status_values, collapse = ", "), "\n")

cat("Step 8: Defining function to prepare data for random forest...\n")
prepare_rf_data <- function(ps_obj) {
  otu_table <- as.data.frame(t(otu_table(ps_obj)))

  sample_data <- as.data.frame(sample_data(ps_obj))

  if(!all(rownames(otu_table) == rownames(sample_data))) {
    stop("Sample IDs in OTU table and sample data do not match")
  }

  sample_data$Status <- factor(sample_data$Status,
                              levels = c("Lab Reared", "Carrion Reared"))

  original_asv_ids <- colnames(otu_table)

  new_asv_names <- paste0("ASV_", 1:ncol(otu_table))
  colnames(otu_table) <- new_asv_names

  asv_mapping <- data.frame(
    ASV_Name = new_asv_names,
    Original_ASV_ID = original_asv_ids,
    stringsAsFactors = FALSE
  )

  tax_table_df <- as.data.frame(tax_table(ps_obj))
  tax_table_df <- tax_table_df %>%
    mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
    mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

  asv_mapping$Kingdom <- tax_table_df[asv_mapping$Original_ASV_ID, "Kingdom"]
  asv_mapping$Phylum <- tax_table_df[asv_mapping$Original_ASV_ID, "Phylum"]
  asv_mapping$Class <- tax_table_df[asv_mapping$Original_ASV_ID, "Class"]
  asv_mapping$Order <- tax_table_df[asv_mapping$Original_ASV_ID, "Order"]
  asv_mapping$Family <- tax_table_df[asv_mapping$Original_ASV_ID, "Family"]
  asv_mapping$Genus <- tax_table_df[asv_mapping$Original_ASV_ID, "Genus"]
  asv_mapping$Species <- tax_table_df[asv_mapping$Original_ASV_ID, "Species"]

  cat("First few ASVs with their taxonomic information:\n")
  print(head(asv_mapping[, c("ASV_Name", "Original_ASV_ID", "Phylum", "Family", "Genus")], 10))

  rf_data <- cbind(otu_table, Status = sample_data$Status)

  cat("First few column names after cleaning:\n")
  print(head(colnames(rf_data), 10))

  return(list(
    data = rf_data,
    asv_mapping = asv_mapping
  ))
}

cat("Step 9: Preparing data for random forest analysis...\n")
rf_result <- prepare_rf_data(ps_stages)
all_data <- rf_result$data
asv_mapping <- rf_result$asv_mapping

cat("Number of samples:", nrow(all_data), "\n")
cat("Number of features (ASVs):", ncol(all_data) - 1, "\n")
cat("Class distribution:\n")
print(table(all_data$Status))

write.csv(asv_mapping, "asv_taxonomy_mapping.csv", row.names = FALSE)
cat("ASV taxonomy mapping saved to asv_taxonomy_mapping.csv\n")

cat("Step 10: Defining function to run random forest and evaluate performance...\n")
run_random_forest <- function(data, title, ntree = 500, seed = 123) {
  set.seed(seed)

  train_indices <- createDataPartition(data$Status, p = 0.7, list = FALSE)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]

  cat("Training set size:", nrow(train_data), "\n")
  cat("Testing set size:", nrow(test_data), "\n")

  cat("Checking for problematic column names...\n")
  problematic_cols <- grep("[^A-Za-z0-9_.]", colnames(train_data))
  if(length(problematic_cols) > 0) {
    cat("Found problematic column names. Cleaning them...\n")
    colnames(train_data)[problematic_cols] <- paste0("ASV_", problematic_cols)
    colnames(test_data)[problematic_cols] <- paste0("ASV_", problematic_cols)
  }

  cat("Training random forest model with", ntree, "trees...\n")
  formula_str <- paste("Status ~", paste(colnames(train_data)[colnames(train_data) != "Status"], collapse = " + "))
  cat("Using formula:", formula_str, "\n")

  rf_model <- randomForest(
    as.formula(formula_str),
    data = train_data,
    ntree = ntree,
    importance = TRUE
  )

  cat("Random Forest Model Summary:\n")
  print(rf_model)

  cat("Making predictions on test data...\n")
  predictions <- predict(rf_model, test_data)

  cat("Calculating confusion matrix...\n")
  conf_matrix <- confusionMatrix(predictions, test_data$Status)

  cat("Calculating ROC curve and AUC...\n")
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

  cat("Extracting variable importance...\n")
  var_importance <- importance(rf_model)
  var_importance_df <- as.data.frame(var_importance)
  var_importance_df$ASV <- rownames(var_importance_df)

  var_importance_df <- var_importance_df %>%
    arrange(desc(MeanDecreaseGini))

  top_vars <- head(var_importance_df, 20)

  top_vars <- top_vars %>%
    left_join(asv_mapping, by = c("ASV" = "ASV_Name")) %>%
    rename(ASV_Name = ASV)

  cat("Creating importance plot...\n")
  importance_plot <- ggplot(top_vars, aes(x = reorder(ASV_Name, MeanDecreaseGini), y = MeanDecreaseGini, fill = Phylum)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste("Top 20 Important ASVs -", title),
      x = "ASV",
      y = "Mean Decrease in Gini Index"
    ) +
    theme_pub() +
    theme(legend.position = "right")

  cat("Creating ROC curve plot...\n")
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

  cat("Returning results...\n")
  return(list(
    model = rf_model,
    confusion_matrix = conf_matrix,
    auc = auc_value,
    importance = top_vars,
    importance_plot = importance_plot,
    roc_plot = roc_plot
  ))
}

cat("\nStep 11: Running random forest for all samples...\n")
all_results <- run_random_forest(all_data, "All Stages")

cat("\nStep 12: Printing results for all samples...\n")
print(all_results$confusion_matrix)
cat("AUC:", all_results$auc, "\n")

cat("\nStep 13: Saving importance plot for all samples...\n")
ggsave("random_forest_importance_all.png", all_results$importance_plot, width = 10, height = 8, dpi = 300)

cat("\nStep 14: Saving ROC curve plot for all samples...\n")
ggsave("random_forest_roc_all.png", all_results$roc_plot, width = 8, height = 6, dpi = 300)

cat("\nStep 15: Creating a detailed table of top important ASVs...\n")
top_asvs_table <- all_results$importance %>%
  select(ASV_Name, MeanDecreaseGini, Original_ASV_ID, Phylum, Family, Genus, Species)

write.csv(top_asvs_table, "top_important_asvs_all.csv", row.names = FALSE)

cat("\nStep 16: Performing random forest analysis for each life stage...\n")

life_stages <- c("3rd instar", "Pupal", "Adult")
stage_results <- list()
stage_data_list <- list()

dir.create("life_stage_results", showWarnings = FALSE)

for (stage in life_stages) {
  cat("\n=== Analyzing", stage, "===\n")

  ps_stage <- subset_samples(ps_stages, Stage == stage)

  status_counts <- table(sample_data(ps_stage)$Status)
  cat("Sample counts by rearing condition:\n")
  print(status_counts)

  if (length(status_counts) < 2 || any(status_counts < 3)) {
    cat("Not enough samples for both conditions in stage:", stage, "\n")
    cat("Skipping random forest analysis for this stage.\n")
    next
  }

  cat("Preparing data for random forest analysis...\n")
  stage_rf_result <- prepare_rf_data(ps_stage)
  stage_data <- stage_rf_result$data
  stage_asv_mapping <- stage_rf_result$asv_mapping

  stage_data_list[[stage]] <- list(
    data = stage_data,
    asv_mapping = stage_asv_mapping
  )

  cat("Running random forest for", stage, "...\n")
  stage_results[[stage]] <- run_random_forest(stage_data, stage)

  cat("\nRandom forest results for", stage, ":\n")
  print(stage_results[[stage]]$confusion_matrix)
  cat("AUC:", stage_results[[stage]]$auc, "\n")

  importance_file <- paste0("life_stage_results/random_forest_importance_", gsub(" ", "_", stage), ".png")
  cat("Saving importance plot to", importance_file, "...\n")
  ggsave(importance_file, stage_results[[stage]]$importance_plot, width = 10, height = 8, dpi = 300)

  roc_file <- paste0("life_stage_results/random_forest_roc_", gsub(" ", "_", stage), ".png")
  cat("Saving ROC curve plot to", roc_file, "...\n")
  ggsave(roc_file, stage_results[[stage]]$roc_plot, width = 8, height = 6, dpi = 300)

  top_asvs_file <- paste0("life_stage_results/top_important_asvs_", gsub(" ", "_", stage), ".csv")
  cat("Saving top ASVs table to", top_asvs_file, "...\n")
  top_stage_asvs <- stage_results[[stage]]$importance %>%
    select(ASV_Name, MeanDecreaseGini, Original_ASV_ID, Phylum, Family, Genus, Species)
  write.csv(top_stage_asvs, top_asvs_file, row.names = FALSE)
}

cat("\nStep 17: Creating visualizations of confusion matrices...\n")

create_confusion_matrix_plot <- function(conf_matrix, title) {
  cm <- as.matrix(conf_matrix$table)

  cm_df <- as.data.frame(cm)
  colnames(cm_df) <- c("Reference", "Prediction", "Freq")

  accuracy <- conf_matrix$overall["Accuracy"]
  sensitivity <- conf_matrix$byClass["Sensitivity"]
  specificity <- conf_matrix$byClass["Specificity"]

  p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "white", size = 8, family = "CMU Sans Serif") +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = title,
      subtitle = paste0(
        "Accuracy: ", round(accuracy, 3),
        " | Sensitivity: ", round(sensitivity, 3),
        " | Specificity: ", round(specificity, 3)
      ),
      x = "Reference (True)",
      y = "Prediction"
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      text = element_text(family = "CMU Sans Serif"),
      plot.title = element_text(family = "CMU Sans Serif", face = "bold", hjust = 0.5),
      plot.subtitle = element_text(family = "CMU Sans Serif", hjust = 0.5),
      axis.title = element_text(family = "CMU Sans Serif"),
      axis.text = element_text(family = "CMU Sans Serif"),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    ) +
    coord_fixed()

  return(p)
}

all_cm_plot <- create_confusion_matrix_plot(all_results$confusion_matrix, "All Stages")

stage_cm_plots <- list()
for (stage in names(stage_results)) {
  stage_cm_plots[[stage]] <- create_confusion_matrix_plot(
    stage_results[[stage]]$confusion_matrix,
    stage
  )
}

if (!requireNamespace("extrafont", quietly = TRUE)) {
  install.packages("extrafont")
  library(extrafont)
} else {
  library(extrafont)
}

all_cm_plots <- c(list(all_cm_plot), stage_cm_plots)
combined_cm_plot <- gridExtra::grid.arrange(
  grobs = all_cm_plots,
  ncol = 2,
  nrow = 2
)

ggsave("random_forest_confusion_matrices.png", combined_cm_plot, width = 13, height = 13, dpi = 300)

cat("\nStep 18: Creating a summary comparison of results across life stages...\n")

summary_results <- data.frame(
  Analysis = "All Stages",
  Accuracy = all_results$confusion_matrix$overall["Accuracy"],
  Sensitivity = all_results$confusion_matrix$byClass["Sensitivity"],
  Specificity = all_results$confusion_matrix$byClass["Specificity"],
  AUC = all_results$auc,
  stringsAsFactors = FALSE
)

for (stage in names(stage_results)) {
  summary_results <- rbind(summary_results, data.frame(
    Analysis = stage,
    Accuracy = stage_results[[stage]]$confusion_matrix$overall["Accuracy"],
    Sensitivity = stage_results[[stage]]$confusion_matrix$byClass["Sensitivity"],
    Specificity = stage_results[[stage]]$confusion_matrix$byClass["Specificity"],
    AUC = stage_results[[stage]]$auc,
    stringsAsFactors = FALSE
  ))
}

cat("\nSummary of random forest results across life stages:\n")
print(summary_results)

write.csv(summary_results, "random_forest_summary_by_stage.csv", row.names = FALSE)

accuracy_plot <- ggplot(summary_results, aes(x = Analysis, y = Accuracy, fill = Analysis)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Accuracy, 3)), vjust = -0.5) +
  labs(
    title = "Random Forest Accuracy by Life Stage",
    x = "Life Stage",
    y = "Accuracy"
  ) +
  ylim(0, 1.1) +
  theme_pub() +
  theme(legend.position = "none")

ggsave("random_forest_accuracy_by_stage.png", accuracy_plot, width = 10, height = 6, dpi = 300)

cat("\nStep 19: Comparing important ASVs across life stages...\n")

top_asvs_all <- all_results$importance %>% head(10) %>% select(ASV_Name, MeanDecreaseGini)
colnames(top_asvs_all)[2] <- "All_Stages"

combined_importance <- top_asvs_all %>% select(ASV_Name)

combined_importance$All_Stages <- top_asvs_all$All_Stages

for (stage in names(stage_results)) {
  top_asvs_stage <- stage_results[[stage]]$importance %>%
    head(10) %>%
    select(ASV_Name, MeanDecreaseGini)
  colnames(top_asvs_stage)[2] <- gsub(" ", "_", stage)

  combined_importance <- full_join(combined_importance, top_asvs_stage, by = "ASV_Name")
}

combined_importance[is.na(combined_importance)] <- 0

combined_importance <- left_join(combined_importance, asv_mapping, by = c("ASV_Name" = "ASV_Name"))

importance_long <- combined_importance %>%
  select(-c(Original_ASV_ID, Kingdom, Class, Order, Species)) %>%
  pivot_longer(cols = -c(ASV_Name, Phylum, Family, Genus),
               names_to = "Analysis",
               values_to = "Importance")

heatmap_plot <- ggplot(importance_long,
                      aes(x = Analysis,
                          y = paste0(ASV_Name, " (", Genus, ")"),
                          fill = Importance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Importance of Top ASVs Across Life Stages",
    x = "Analysis",
    y = "ASV - Genus"
  ) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

ggsave("random_forest_importance_heatmap_by_stage.png", heatmap_plot, width = 12, height = 10, dpi = 300)

cat("\nRandom forest analysis completed successfully!\n")
cat("Output files:\n")
cat("1. random_forest_importance_all.png - Importance plot for all samples\n")
cat("2. random_forest_roc_all.png - ROC curve plot for all samples\n")
cat("3. asv_taxonomy_mapping.csv - Complete mapping of ASVs to their taxonomic classification\n")
cat("4. top_important_asvs_all.csv - Detailed table of top important ASVs with their taxonomic information\n")
cat("5. life_stage_results/ - Directory containing results for each life stage\n")
cat("6. random_forest_confusion_matrices.png - Combined visualization of confusion matrices\n")
cat("7. random_forest_summary_by_stage.csv - Summary of random forest results across life stages\n")
cat("8. random_forest_accuracy_by_stage.png - Bar plot of accuracy by life stage\n")
cat("9. random_forest_importance_heatmap_by_stage.png - Heatmap of important ASVs across life stages\n")
