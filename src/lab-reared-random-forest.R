library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(caret)
library(pROC)
library(RColorBrewer)

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

ps_insects <- subset_samples(ps_insects, Group == "EGG")

prepare_data_for_rf <- function(ps_obj) {
  otu_table <- as.data.frame(t(otu_table(ps_obj)))

  otu_table[is.na(otu_table)] <- 0

  metadata <- as.data.frame(sample_data(ps_obj))

  data <- cbind(metadata, otu_table)

  if (any(is.na(data))) {
    cat("Warning: Data still contains NA values after processing.\n")
    cat("Number of NA values:", sum(is.na(data)), "\n")

    na_cols <- colnames(data)[colSums(is.na(data)) > 0]
    cat("Columns with NA values:", paste(na_cols, collapse=", "), "\n")

    metadata_cols <- colnames(metadata)
    data <- data[complete.cases(data[, metadata_cols]), ]
    cat("Removed rows with NA values in metadata. Remaining rows:", nrow(data), "\n")
  }

  return(data)
}

normalize_counts <- function(otu_table) {
  rel_abund <- sweep(otu_table, 1, rowSums(otu_table), "/")
  return(rel_abund)
}

split_data <- function(data, target_var, train_prop = 0.7, seed = 42) {
  set.seed(seed)

  data[[target_var]] <- as.factor(data[[target_var]])

  train_indices <- createDataPartition(data[[target_var]], p = train_prop, list = FALSE)

  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]

  return(list(train = train_data, test = test_data))
}

train_rf_model <- function(train_data, target_var, ntree = 500, mtry = NULL, seed = 42) {
  set.seed(seed)

  features <- train_data[, !colnames(train_data) %in% c("Sample", "Class", "Group", "Stage", "RearingGroup", "Generation")]
  target <- train_data[[target_var]]

  if (any(is.na(features))) {
    cat("Warning: Features contain NA values. Replacing with 0...\n")
    features[is.na(features)] <- 0
  }

  if (any(is.na(target))) {
    cat("Warning: Target variable contains NA values. Removing corresponding samples...\n")
    valid_indices <- !is.na(target)
    features <- features[valid_indices, ]
    target <- target[valid_indices]
  }

  if (is.null(mtry)) {
    mtry <- floor(sqrt(ncol(features)))
  }

  rf_model <- randomForest(
    x = features,
    y = target,
    ntree = ntree,
    mtry = mtry,
    importance = TRUE,
    na.action = na.omit
  )

  return(rf_model)
}

evaluate_model <- function(model, test_data, target_var) {
  features <- test_data[, !colnames(test_data) %in% c("Sample", "Class", "Group", "Stage", "RearingGroup", "Generation")]
  actual <- test_data[[target_var]]

  predictions <- predict(model, features)

  accuracy <- sum(predictions == actual) / length(actual)

  all_levels <- unique(c(levels(actual), levels(predictions)))
  actual <- factor(actual, levels = all_levels)
  predictions <- factor(predictions, levels = all_levels)

  conf_matrix <- tryCatch({
    confusionMatrix(predictions, actual)
  }, error = function(e) {
    cat("Error creating confusion matrix:", e$message, "\n")
    cat("Continuing without confusion matrix...\n")
    return(NULL)
  })

  return(list(
    accuracy = accuracy,
    confusion_matrix = conf_matrix
  ))
}

plot_variable_importance <- function(model, n_top = 20, taxonomy = NULL) {
  var_imp <- importance(model)

  var_imp_sorted <- var_imp[order(var_imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]

  top_vars <- rownames(var_imp_sorted)[1:min(n_top, nrow(var_imp_sorted))]
  top_imp <- var_imp_sorted[top_vars, "MeanDecreaseGini", drop = FALSE]

  imp_df <- data.frame(
    Feature = rownames(top_imp),
    Importance = top_imp[, "MeanDecreaseGini"]
  )

  if (!is.null(taxonomy)) {
    otu_ids <- imp_df$Feature

    tax_df <- taxonomy[otu_ids, ]

    tax_labels <- paste0(
      tax_df$Genus, " (",
      tax_df$Phylum, "; ",
      tax_df$Class, "; ",
      tax_df$Order, "; ",
      tax_df$Family, ")"
    )

    tax_labels <- gsub("NA", "Unknown", tax_labels)

    imp_df$Taxonomy <- tax_labels
  } else {
    imp_df$Taxonomy <- imp_df$Feature
  }

  p <- ggplot(imp_df, aes(x = reorder(Taxonomy, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Top Features by Importance",
      x = "",
      y = "Mean Decrease in Gini"
    ) +
    theme_pub()

  return(p)
}

available_stages <- unique(sample_data(ps_insects)$Stage)
cat("Available stages:", paste(available_stages, collapse = ", "), "\n")

stage_order <- c("Egg (Generation 0)", "Egg (Generation 1)", "1st instar", "2nd instar", "3rd instar", "Pupal", "Adult")
actual_stages <- stage_order[stage_order %in% available_stages]

stage_results <- list()

for (stage in actual_stages) {
  cat("\n=== RANDOM FOREST MODEL FOR STAGE:", stage, "===\n")

  ps_stage <- subset_samples(ps_insects, Stage == stage)

  if (nsamples(ps_stage) < 10) {
    cat("Not enough samples for stage:", stage, "\n")
    next
  }

  rearing_groups <- unique(sample_data(ps_stage)$RearingGroup)
  if (length(rearing_groups) < 2) {
    cat("Only one rearing group for stage:", stage, "\n")
    next
  }

  data <- prepare_data_for_rf(ps_stage)

  data_split <- split_data(data, "RearingGroup", train_prop = 0.7)

  rf_model <- train_rf_model(data_split$train, "RearingGroup")

  evaluation <- evaluate_model(rf_model, data_split$test, "RearingGroup")

  cat("Accuracy:", evaluation$accuracy, "\n")
  if (!is.null(evaluation$confusion_matrix)) {
    print(evaluation$confusion_matrix)
  } else {
    cat("Confusion matrix not available due to error.\n")
  }

  var_imp_plot <- plot_variable_importance(rf_model, n_top = 15)

  ggsave(paste0("rf_rearing_group_", gsub(" ", "_", stage), ".png"), var_imp_plot, width = 10, height = 8, dpi = 300, bg = "white")

  stage_results[[stage]] <- list(
    model = rf_model,
    evaluation = evaluation,
    var_imp_plot = var_imp_plot
  )
}

cat("\n=== RANDOM FOREST MODEL FOR PREDICTING STAGE ===\n")

data <- prepare_data_for_rf(ps_insects)

data_split <- split_data(data, "Stage", train_prop = 0.7)

rf_model_stage <- train_rf_model(data_split$train, "Stage")

evaluation_stage <- evaluate_model(rf_model_stage, data_split$test, "Stage")

cat("Accuracy:", evaluation_stage$accuracy, "\n")
if (!is.null(evaluation_stage$confusion_matrix)) {
  print(evaluation_stage$confusion_matrix)
} else {
  cat("Confusion matrix not available due to error.\n")
}

var_imp_plot_stage <- plot_variable_importance(rf_model_stage, n_top = 20)

ggsave("rf_stage_prediction.png", var_imp_plot_stage, width = 10, height = 8, dpi = 300, bg = "white")

perform_cv <- function(data, target_var, k = 5, ntree = 500, seed = 42) {
  set.seed(seed)

  data[[target_var]] <- as.factor(data[[target_var]])

  folds <- createFolds(data[[target_var]], k = k)

  accuracy <- numeric(k)

  for (i in 1:k) {
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]

    rf_model <- train_rf_model(train_data, target_var, ntree = ntree)

    evaluation <- evaluate_model(rf_model, test_data, target_var)

    accuracy[i] <- evaluation$accuracy
  }

  mean_accuracy <- mean(accuracy)
  sd_accuracy <- sd(accuracy)

  return(list(
    accuracy = accuracy,
    mean_accuracy = mean_accuracy,
    sd_accuracy = sd_accuracy
  ))
}

cat("\n=== CROSS-VALIDATION FOR STAGE PREDICTION ===\n")
cv_stage <- perform_cv(data, "Stage", k = 5)
cat("Mean accuracy:", cv_stage$mean_accuracy, "\n")
cat("Standard deviation:", cv_stage$sd_accuracy, "\n")

cat("\n=== SUMMARY OF RANDOM FOREST MODELS ===\n")

cat("Stage-specific models for predicting rearing group:\n")
for (stage in names(stage_results)) {
  cat("- Stage:", stage, "- Accuracy:", stage_results[[stage]]$evaluation$accuracy, "\n")
}

cat("\nModel for predicting stage:\n")
cat("- Accuracy:", evaluation_stage$accuracy, "\n")
cat("- Cross-validation mean accuracy:", cv_stage$mean_accuracy, "\n")

for (stage in names(stage_results)) {
  saveRDS(stage_results[[stage]]$model, paste0("rf_model_rearing_group_", gsub(" ", "_", stage), ".rds"))
}

saveRDS(rf_model_stage, "rf_model_stage.rds")

cat("\nModels saved to disk.\n")
