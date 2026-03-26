library(qiime2R)
library(phyloseq)
library(randomForest)
library(rpart)
library(rpart.plot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
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

  rf_data <- cbind(otu_table, Status = sample_data$Status)

  return(list(
    data = rf_data,
    asv_mapping = asv_mapping
  ))
}

visualize_tree <- function(rf_model, tree_num = 1, asv_mapping = NULL, max_depth = 5) {
  tree <- getTree(rf_model, k = tree_num, labelVar = TRUE)

  tree_df <- as.data.frame(tree)

  if (!is.null(asv_mapping)) {
    var_names <- unique(tree_df$`split var`)
    var_names <- var_names[!is.na(var_names)]

    var_mapping <- data.frame(
      original = var_names,
      new = var_names,
      stringsAsFactors = FALSE
    )

    for (i in 1:nrow(var_mapping)) {
      var_name <- var_mapping$original[i]
      if (var_name %in% asv_mapping$ASV_Name) {
        asv_info <- asv_mapping[asv_mapping$ASV_Name == var_name, ]
        new_name <- paste0(var_name, " (", asv_info$Genus, ")")
        var_mapping$new[i] <- new_name
      }
    }

    for (i in 1:nrow(var_mapping)) {
      tree_df$`split var`[tree_df$`split var` == var_mapping$original[i]] <- var_mapping$new[i]
    }
  }

  data <- rf_model$call$data
  if (is.null(data)) {
    cat("Warning: Original data not available in the model. Using a simplified approach.\n")

    plot_tree_structure(tree_df, max_depth)
    return(invisible(NULL))
  }

  formula_str <- paste("Status ~", paste(colnames(data)[colnames(data) != "Status"], collapse = " + "))
  rpart_model <- rpart(as.formula(formula_str), data = data, method = "class", maxdepth = max_depth)

  rpart.plot(rpart_model,
             box.palette = "RdBu",
             shadow.col = "gray",
             nn = TRUE,
             main = paste("Decision Tree", tree_num, "from Random Forest"))

  return(invisible(NULL))
}

plot_tree_structure <- function(tree_df, max_depth = 5) {
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "",
       main = "Simplified Tree Structure",
       axes = FALSE)

  plot_node(tree_df, node = 1, x = 0.5, y = 0.9, width = 0.4, depth = 0, max_depth = max_depth)
}

plot_node <- function(tree_df, node, x, y, width, depth, max_depth) {
  if (depth >= max_depth) {
    return()
  }

  node_info <- tree_df[node, ]

  if (!is.na(node_info$`split var`) && node_info$`split var` == "<leaf>") {
    rect(x - width/4, y - 0.03, x + width/4, y + 0.03, col = "lightgreen", border = "black")
    text(x, y, paste("Class", node_info$prediction), cex = 0.8)
    return()
  }

  if (is.na(node_info$`split var`)) {
    return()
  }

  rect(x - width/3, y - 0.04, x + width/3, y + 0.04, col = "lightblue", border = "black")
  text(x, y, paste0(substr(node_info$`split var`, 1, 10), " < ", round(node_info$`split point`, 2)), cex = 0.8)

  left_x <- x - width/2
  right_x <- x + width/2
  child_y <- y - 0.15

  lines(c(x, left_x), c(y - 0.04, child_y + 0.04))
  lines(c(x, right_x), c(y - 0.04, child_y + 0.04))

  text(x - width/4, y - 0.08, "Yes", cex = 0.7)
  text(x + width/4, y - 0.08, "No", cex = 0.7)

  plot_node(tree_df, node = node_info$`left daughter`, x = left_x, y = child_y,
            width = width * 0.8, depth = depth + 1, max_depth = max_depth)
  plot_node(tree_df, node = node_info$`right daughter`, x = right_x, y = child_y,
            width = width * 0.8, depth = depth + 1, max_depth = max_depth)
}

visualize_multiple_trees <- function(rf_model, n_trees = 3, asv_mapping = NULL, max_depth = 4) {
  pdf("random_forest_trees.pdf", width = 12, height = 10)

  for (i in 1:min(n_trees, rf_model$ntree)) {
    cat("Visualizing tree", i, "of", min(n_trees, rf_model$ntree), "\n")
    visualize_tree(rf_model, tree_num = i, asv_mapping = asv_mapping, max_depth = max_depth)
  }

  dev.off()

  cat("Trees visualized and saved to random_forest_trees.pdf\n")
}

visualize_tree_ggplot <- function(rf_model, tree_num = 1, asv_mapping = NULL, max_depth = 5) {
  tree <- getTree(rf_model, k = tree_num, labelVar = TRUE)

  tree_df <- as.data.frame(tree)

  nodes <- data.frame(
    id = 1:nrow(tree_df),
    x = NA_real_,
    y = NA_real_,
    label = NA_character_,
    type = NA_character_,
    stringsAsFactors = FALSE
  )

  nodes$x[1] <- 0.5
  nodes$y[1] <- 1

  if (!is.null(asv_mapping)) {
    for (i in 1:nrow(tree_df)) {
      var_name <- tree_df$`split var`[i]
      if (!is.na(var_name) && var_name != "<leaf>" && var_name %in% asv_mapping$ASV_Name) {
        asv_info <- asv_mapping[asv_mapping$ASV_Name == var_name, ]
        tree_df$`split var`[i] <- paste0(var_name, " (", asv_info$Genus, ")")
      }
    }
  }

  for (i in 1:nrow(tree_df)) {
    if (!is.na(tree_df$`split var`[i]) && tree_df$`split var`[i] == "<leaf>") {
      nodes$label[i] <- paste("Class:", tree_df$prediction[i])
      nodes$type[i] <- "leaf"
    } else if (!is.na(tree_df$`split var`[i])) {
      nodes$label[i] <- paste0(tree_df$`split var`[i], " < ", round(tree_df$`split point`[i], 2))
      nodes$type[i] <- "split"
    } else {
      nodes$label[i] <- "NA"
      nodes$type[i] <- "unknown"
    }
  }

  queue <- 1
  visited <- rep(FALSE, nrow(tree_df))
  visited[1] <- TRUE

  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]

    if (is.na(tree_df$`split var`[current]) || tree_df$`split var`[current] == "<leaf>") {
      next
    }

    left_child <- tree_df$`left daughter`[current]
    right_child <- tree_df$`right daughter`[current]

    nodes$x[left_child] <- nodes$x[current] - 0.2 * (max_depth - floor(log2(current)))
    nodes$x[right_child] <- nodes$x[current] + 0.2 * (max_depth - floor(log2(current)))

    depth <- floor(log2(current)) + 1
    nodes$y[left_child] <- nodes$y[current] - 0.15
    nodes$y[right_child] <- nodes$y[current] - 0.15

    if (!visited[left_child]) {
      visited[left_child] <- TRUE
      queue <- c(queue, left_child)
    }
    if (!visited[right_child]) {
      visited[right_child] <- TRUE
      queue <- c(queue, right_child)
    }
  }

  edges <- data.frame(
    from = rep(NA_integer_, nrow(tree_df) * 2),
    to = rep(NA_integer_, nrow(tree_df) * 2),
    label = rep(NA_character_, nrow(tree_df) * 2),
    stringsAsFactors = FALSE
  )

  edge_count <- 0
  for (i in 1:nrow(tree_df)) {
    if (tree_df$`split var`[i] != "<leaf>") {
      edge_count <- edge_count + 1
      edges$from[edge_count] <- i
      edges$to[edge_count] <- tree_df$`left daughter`[i]
      edges$label[edge_count] <- "Yes"

      edge_count <- edge_count + 1
      edges$from[edge_count] <- i
      edges$to[edge_count] <- tree_df$`right daughter`[i]
      edges$label[edge_count] <- "No"
    }
  }

  edges <- edges[1:edge_count, ]

  p <- ggplot() +
    geom_segment(data = edges,
                aes(x = nodes$x[from], y = nodes$y[from],
                    xend = nodes$x[to], yend = nodes$y[to]),
                arrow = arrow(length = unit(0.2, "cm"))) +
    geom_text(data = edges,
              aes(x = (nodes$x[from] + nodes$x[to])/2,
                  y = (nodes$y[from] + nodes$y[to])/2,
                  label = label),
              size = 3, vjust = -0.5) +
    geom_point(data = nodes,
               aes(x = x, y = y, color = type),
               size = 15, alpha = 0.7) +
    geom_text(data = nodes,
              aes(x = x, y = y, label = label),
              size = 2.5) +
    scale_color_manual(values = c("leaf" = "lightgreen", "split" = "lightblue")) +
    labs(title = paste("Decision Tree", tree_num, "from Random Forest"),
         subtitle = "Simplified Visualization") +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12))

  return(p)
}

visualize_multiple_trees_ggplot <- function(rf_model, n_trees = 3, asv_mapping = NULL, max_depth = 4) {
  plot_list <- list()

  for (i in 1:min(n_trees, rf_model$ntree)) {
    cat("Visualizing tree", i, "of", min(n_trees, rf_model$ntree), "\n")
    plot_list[[i]] <- visualize_tree_ggplot(rf_model, tree_num = i, asv_mapping = asv_mapping, max_depth = max_depth)
  }

  for (i in 1:length(plot_list)) {
    ggsave(paste0("random_forest_tree_", i, ".png"), plot_list[[i]], width = 10, height = 8, dpi = 300)
  }

  if (length(plot_list) > 1) {
    combined_plot <- gridExtra::grid.arrange(
      grobs = plot_list,
      ncol = min(2, length(plot_list)),
      top = "Sample Trees from Random Forest"
    )

    ggsave("random_forest_trees_combined.png", combined_plot, width = 16, height = 12, dpi = 300)
  }

  cat("Trees visualized and saved as PNG files\n")
}

cat("Loading data and models...\n")

ps <- qiime2R::qza_to_phyloseq(
  features = "data/filtered-dada-table-nmnc.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata.tsv"
)

ps_insects <- subset_samples(ps, Class == "Insect")

target_stages <- c("3rd instar", "Pupal", "Adult")
ps_stages <- subset_samples(ps_insects, Stage %in% target_stages)

cat("Preparing data for random forest...\n")
rf_result <- prepare_rf_data(ps_stages)
all_data <- rf_result$data
asv_mapping <- rf_result$asv_mapping

cat("Training random forest model...\n")
set.seed(123)
rf_model <- randomForest(
  Status ~ .,
  data = all_data,
  ntree = 50,
  importance = TRUE
)

cat("\nRandom Forest Model Summary:\n")
print(rf_model)

dir.create("tree_visualizations", showWarnings = FALSE)

cat("\nVisualizing trees using different methods...\n")

cat("\nMethod 1: Using getTree and custom plotting...\n")
pdf("tree_visualizations/random_forest_trees_method1.pdf", width = 12, height = 10)
for (i in 1:3) {
  cat("Visualizing tree", i, "of 3\n")
  tree <- getTree(rf_model, k = i, labelVar = TRUE)
  plot_tree_structure(as.data.frame(tree), max_depth = 5)
  title(paste("Decision Tree", i, "from Random Forest"))
}
dev.off()

cat("\nMethod 2: Using ggplot2 for visualization...\n")
visualize_multiple_trees_ggplot(rf_model, n_trees = 3, asv_mapping = asv_mapping, max_depth = 5)

cat("\nTree visualizations completed!\n")
cat("Output files:\n")
cat("1. tree_visualizations/random_forest_trees_method1.pdf - Basic tree visualizations\n")
cat("2. random_forest_tree_1.png, random_forest_tree_2.png, random_forest_tree_3.png - Individual tree visualizations using ggplot2\n")
cat("3. random_forest_trees_combined.png - Combined visualization of multiple trees\n")
