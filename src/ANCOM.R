library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ANCOMBC)
library(RColorBrewer)
library(patchwork)
library(scales)

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

perform_ancombc <- function(ps_obj, stage) {
  ps_stage <- subset_samples(ps_obj, Stage == stage)

  status_counts <- table(sample_data(ps_stage)$Status)
  if (any(status_counts < 3)) {
    cat("Warning: Stage", stage, "has fewer than 3 samples for some Status groups.\n")
    cat("Status counts:", paste(names(status_counts), status_counts, sep = ": ", collapse = ", "), "\n")
    cat("Skipping ANCOM-BC analysis for this stage.\n")
    return(NULL)
  }

  tryCatch({
    ancombc_result <- ancombc(
      phyloseq = ps_stage,
      formula = "Status",
      p_adj_method = "BH",
      group = "Status",
      struc_zero = TRUE,
      neg_lb = TRUE,
      tol = 1e-5,
      max_iter = 100,
      conserve = TRUE,
      alpha = 0.05,
      global = TRUE
    )
  }, error = function(e) {
    cat("First attempt failed with error:", e$message, "\n")
    cat("Trying alternative parameters...\n")

    ancombc_result <- ancombc(
      data = ps_stage,
      formula = "Status",
      p_adj_method = "BH"
    )

    return(ancombc_result)
  })

  return(ancombc_result)
}

ancombc_results <- list()
for (stage in stage_order) {
  cat("\n=== PERFORMING ANCOM-BC ANALYSIS FOR STAGE:", stage, "===\n")
  ancombc_results[[stage]] <- perform_ancombc(ps_stages, stage)
}

cat("\n=== SUMMARY OF ANCOM-BC RESULTS ===\n")
for (stage in stage_order) {
  cat("\nStage:", stage, "\n")

  if (is.null(ancombc_results[[stage]])) {
    cat("No ANCOM-BC results available for this stage.\n")
    next
  }

  cat("Structure of ANCOM-BC results:\n")
  print(names(ancombc_results[[stage]]))

  if (!is.null(ancombc_results[[stage]]$feature_table)) {
    cat("Number of taxa analyzed:", nrow(ancombc_results[[stage]]$feature_table), "\n")
  }

  if (!is.null(ancombc_results[[stage]]$res)) {
    res <- ancombc_results[[stage]]$res
    cat("Components in results:", paste(names(res), collapse = ", "), "\n")

    if (!is.null(res$p_val) && ncol(res$p_val) >= 2) {
      sig_count <- sum(res$p_val[, 2] < 0.05, na.rm = TRUE)
      cat("Number of taxa with p < 0.05:", sig_count, "\n")
    }

    if (!is.null(res$q_val) && ncol(res$q_val) >= 2) {
      sig_count <- sum(res$q_val[, 2] < 0.05, na.rm = TRUE)
      cat("Number of taxa with q < 0.05 (FDR-corrected):", sig_count, "\n")
    }

    if (!is.null(res$diff_abn) && ncol(res$diff_abn) >= 2) {
      diff_count <- sum(res$diff_abn[, 2], na.rm = TRUE)
      cat("Number of differentially abundant taxa:", diff_count, "\n")
    }

    if (!is.null(res$p_val) && ncol(res$p_val) >= 2) {
      cat("\nTop 5 taxa by p-value:\n")

      taxa_df <- NULL

      tryCatch({
        taxa_df <- data.frame(
          taxon = rownames(res$p_val),
          p_val = res$p_val[, 2],
          q_val = res$q_val[, 2],
          lfc = res$beta[, 2],
          Genus = rownames(res$p_val),
          stringsAsFactors = FALSE
        )

        if (!is.null(ancombc_results[[stage]]$tax_info) &&
            nrow(ancombc_results[[stage]]$tax_info) > 0 &&
            nrow(ancombc_results[[stage]]$tax_info) == nrow(taxa_df)) {

          tax_info <- ancombc_results[[stage]]$tax_info %>%
            mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
            mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

          if ("Genus" %in% colnames(tax_info)) {
            taxa_df$Genus <- tax_info$Genus
          }
        }
      }, error = function(e) {
        cat("Error creating data frame:", e$message, "\n")
        cat("Will try a simpler approach...\n")
        taxa_df <- NULL
      })

      if (is.null(taxa_df)) {
        tryCatch({
          taxa_df <- data.frame(
            taxon = rownames(res$p_val),
            p_val = res$p_val[, 2],
            stringsAsFactors = FALSE
          )

          if (!is.null(res$q_val) && ncol(res$q_val) >= 2) {
            taxa_df$q_val <- res$q_val[, 2]
          } else {
            taxa_df$q_val <- NA
          }

          if (!is.null(res$beta) && ncol(res$beta) >= 2) {
            taxa_df$lfc <- res$beta[, 2]
          } else {
            taxa_df$lfc <- NA
          }

          taxa_df$Genus <- rownames(res$p_val)

        }, error = function(e) {
          cat("Error creating simplified data frame:", e$message, "\n")
          cat("Unable to create a data frame for this stage.\n")
          taxa_df <- NULL
        })
      }

      if (is.null(taxa_df)) {
        cat("Skipping top taxa display for this stage due to data frame creation errors.\n")
        next
      }

      top_taxa <- taxa_df %>%
        arrange(p_val) %>%
        head(5)

      for (i in 1:nrow(top_taxa)) {
        cat("  ", top_taxa$Genus[i], "(LFC =", round(top_taxa$lfc[i], 2),
            ", p =", format(top_taxa$p_val[i], digits = 3),
            ", q =", format(top_taxa$q_val[i], digits = 3), ")\n")
      }
    }
  }
}

cat("\n=== OVERALL SUMMARY ===\n")
for (stage in stage_order) {
  if (!is.null(ancombc_results[[stage]]) && !is.null(ancombc_results[[stage]]$res) &&
      !is.null(ancombc_results[[stage]]$res$q_val) && ncol(ancombc_results[[stage]]$res$q_val) >= 2) {

    sig_count <- sum(ancombc_results[[stage]]$res$q_val[, 2] < 0.05, na.rm = TRUE)
    cat("Stage", stage, ":", sig_count, "differentially abundant taxa (q < 0.05)\n")
  } else {
    cat("Stage", stage, ": No results available\n")
  }
}

cat("\n=== ALL CONDITIONS ANCOM ANALYSIS SUMMARY ===\n")
cat("Number of samples:", nsamples(ps_stages), "\n")
cat("Number of taxa:", ntaxa(ps_stages), "\n")
cat("Developmental stages included:", paste(stage_order, collapse = ", "), "\n")
cat("Rearing conditions included:", paste(available_status, collapse = ", "), "\n")

cat("\nDifferentially abundant taxa by stage (q < 0.05):\n")
for (stage in stage_order) {
  if (!is.null(ancombc_results[[stage]]) && !is.null(ancombc_results[[stage]]$res) &&
      !is.null(ancombc_results[[stage]]$res$q_val) && ncol(ancombc_results[[stage]]$res$q_val) >= 2) {

    sig_count <- sum(ancombc_results[[stage]]$res$q_val[, 2] < 0.05, na.rm = TRUE)
    cat("Stage", stage, ":", sig_count, "differentially abundant taxa\n")

    if (sig_count > 0) {
      taxa_df <- NULL

      tryCatch({
        taxa_df <- data.frame(
          taxon = rownames(ancombc_results[[stage]]$res$q_val),
          p_val = ancombc_results[[stage]]$res$p_val[, 2],
          q_val = ancombc_results[[stage]]$res$q_val[, 2],
          lfc = ancombc_results[[stage]]$res$beta[, 2],
          Genus = rownames(ancombc_results[[stage]]$res$q_val),
          stringsAsFactors = FALSE
        )

        if (!is.null(ancombc_results[[stage]]$tax_info) &&
            nrow(ancombc_results[[stage]]$tax_info) > 0 &&
            nrow(ancombc_results[[stage]]$tax_info) == nrow(taxa_df)) {

          tax_info <- ancombc_results[[stage]]$tax_info %>%
            mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
            mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

          if ("Genus" %in% colnames(tax_info)) {
            taxa_df$Genus <- tax_info$Genus
          }
        }
      }, error = function(e) {
        cat("Error creating data frame:", e$message, "\n")
        cat("Will try a simpler approach...\n")
        taxa_df <- NULL
      })

      if (is.null(taxa_df)) {
        tryCatch({
          taxa_df <- data.frame(
            taxon = rownames(ancombc_results[[stage]]$res$q_val),
            p_val = ancombc_results[[stage]]$res$p_val[, 2],
            stringsAsFactors = FALSE
          )

          if (!is.null(ancombc_results[[stage]]$res$q_val) && ncol(ancombc_results[[stage]]$res$q_val) >= 2) {
            taxa_df$q_val <- ancombc_results[[stage]]$res$q_val[, 2]
          } else {
            taxa_df$q_val <- NA
          }

          if (!is.null(ancombc_results[[stage]]$res$beta) && ncol(ancombc_results[[stage]]$res$beta) >= 2) {
            taxa_df$lfc <- ancombc_results[[stage]]$res$beta[, 2]
          } else {
            taxa_df$lfc <- NA
          }

          taxa_df$Genus <- rownames(ancombc_results[[stage]]$res$q_val)

        }, error = function(e) {
          cat("Error creating simplified data frame:", e$message, "\n")
          cat("Unable to create a data frame for this stage.\n")
          taxa_df <- NULL
        })
      }

      if (is.null(taxa_df)) {
        cat("Skipping top taxa display for this stage due to data frame creation errors.\n")
        next
      }

      tryCatch({
        sig_taxa <- taxa_df %>%
          filter(q_val < 0.05) %>%
          arrange(q_val) %>%
          head(5)

        cat("  Top 5 taxa by significance:\n")
        for (i in 1:nrow(sig_taxa)) {
          cat("    ", sig_taxa$Genus[i], "(LFC =", round(sig_taxa$lfc[i], 2),
              ", p =", format(sig_taxa$p_val[i], digits = 3),
              ", q =", format(sig_taxa$q_val[i], digits = 3), ")\n")
        }
      }, error = function(e) {
        cat("Error filtering significant taxa:", e$message, "\n")
        cat("Unable to display top taxa for this stage.\n")
      })
    }
  } else {
    cat("Stage", stage, ": No results available\n")
  }
}
