library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(VennDiagram)
library(gridExtra)
library(grid)
library(circlize)

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

cat("\n=== AVAILABLE REARING CONDITIONS (STATUS) ===\n")
available_status <- unique(sample_data(ps_insects)$Status)
cat("Available rearing conditions:", paste(available_status, collapse = ", "), "\n")

target_stages <- c("3rd instar", "Pupal", "Adult")
ps_stages <- subset_samples(ps_insects, Stage %in% target_stages)

cat("\n=== STATUS VALUES IN THE DATASET ===\n")
status_values <- unique(sample_data(ps_stages)$Status)
cat("Status values:", paste(status_values, collapse = ", "), "\n")

ps_lab <- subset_samples(ps_stages, Status == "Lab Reared")
ps_carrion <- subset_samples(ps_stages, Status == "Carrion Reared")

cat("\n=== FILTERED DATASET SUMMARY ===\n")
cat("Lab-reared samples:", nsamples(ps_lab), "\n")
cat("Carrion-reared samples:", nsamples(ps_carrion), "\n")

ps_lab <- prune_taxa(taxa_sums(ps_lab) > 0, ps_lab)
ps_carrion <- prune_taxa(taxa_sums(ps_carrion) > 0, ps_carrion)

lab_asvs <- taxa_names(ps_lab)
carrion_asvs <- taxa_names(ps_carrion)

shared_asvs <- intersect(lab_asvs, carrion_asvs)
lab_only_asvs <- setdiff(lab_asvs, carrion_asvs)
carrion_only_asvs <- setdiff(carrion_asvs, lab_asvs)

cat("\n=== ASV OVERLAP SUMMARY ===\n")
cat("Lab-reared ASVs:", length(lab_asvs), "\n")
cat("Carrion-reared ASVs:", length(carrion_asvs), "\n")
cat("Shared ASVs:", length(shared_asvs), "\n")
cat("Lab-reared only ASVs:", length(lab_only_asvs), "\n")
cat("Carrion-reared only ASVs:", length(carrion_only_asvs), "\n")

lab_shared_percent <- round(length(shared_asvs) / length(lab_asvs) * 100, 1)
carrion_shared_percent <- round(length(shared_asvs) / length(carrion_asvs) * 100, 1)

cat("\n=== PERCENTAGE OVERLAP ===\n")
cat("Percentage of Lab-reared ASVs shared with Carrion-reared:", lab_shared_percent, "%\n")
cat("Percentage of Carrion-reared ASVs shared with Lab-reared:", carrion_shared_percent, "%\n")

all_tax_table <- as.data.frame(tax_table(ps))

all_tax_table <- all_tax_table %>%
  mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
  mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

lab_abundance_matrix <- as.matrix(otu_table(ps_lab))
carrion_abundance_matrix <- as.matrix(otu_table(ps_carrion))

lab_total_abundance <- rowSums(lab_abundance_matrix)
carrion_total_abundance <- rowSums(carrion_abundance_matrix)

asv_data <- data.frame(
  ASV = taxa_names(ps),
  LabAbundance = 0,
  CarrionAbundance = 0,
  Phylum = all_tax_table$Phylum,
  Family = all_tax_table$Family,
  Genus = all_tax_table$Genus,
  stringsAsFactors = FALSE
)

asv_data$LabAbundance[match(names(lab_total_abundance), asv_data$ASV)] <- lab_total_abundance
asv_data$CarrionAbundance[match(names(carrion_total_abundance), asv_data$ASV)] <- carrion_total_abundance
asv_data$TotalAbundance <- asv_data$LabAbundance + asv_data$CarrionAbundance

simplify_genus_name <- function(genus_name) {
  if (grepl("_", genus_name)) {
    parts <- strsplit(genus_name, "_")[[1]]
    if (length(parts) > 1 && nchar(parts[2]) > 0) {
      simplified <- paste0(parts[1], "_", substr(parts[2], 1, 1))
    } else {
      simplified <- parts[1]
    }
    return(simplified)
  } else {
    return(genus_name)
  }
}

asv_data$SimplifiedGenus <- sapply(asv_data$Genus, simplify_genus_name)

genus_abundance_all <- asv_data %>%
  group_by(SimplifiedGenus) %>%
  summarize(
    LabAbundance = sum(LabAbundance),
    CarrionAbundance = sum(CarrionAbundance),
    TotalAbundance = sum(TotalAbundance)
  ) %>%
  filter(SimplifiedGenus != "Unknown")

overlap_genera <- genus_abundance_all %>%
  filter(LabAbundance > 0 & CarrionAbundance > 0) %>%
  arrange(desc(TotalAbundance)) %>%
  head(10) %>%
  pull(SimplifiedGenus)

lab_exclusive_genera <- genus_abundance_all %>%
  filter(LabAbundance > 0 & CarrionAbundance == 0) %>%
  arrange(desc(LabAbundance)) %>%
  head(5) %>%
  pull(SimplifiedGenus)

carrion_exclusive_genera <- genus_abundance_all %>%
  filter(LabAbundance == 0 & CarrionAbundance > 0) %>%
  arrange(desc(CarrionAbundance)) %>%
  head(5) %>%
  pull(SimplifiedGenus)

top_genera_combined <- c(overlap_genera, lab_exclusive_genera, carrion_exclusive_genera)

genus_abundance <- genus_abundance_all %>%
  filter(SimplifiedGenus %in% top_genera_combined) %>%
  arrange(desc(TotalAbundance))

cat("\n=== SELECTED GENERA FOR CHORD DIAGRAM ===\n")
cat("Number of overlapping genera (present in both conditions):", length(overlap_genera), "\n")
cat("Number of lab-reared exclusive genera:", length(lab_exclusive_genera), "\n")
cat("Number of carrion-reared exclusive genera:", length(carrion_exclusive_genera), "\n")
cat("Total number of genera in combined list:", length(top_genera_combined), "\n")

cat("\n=== TOP 20 GENERA BY TOTAL ABUNDANCE ===\n")
print(genus_abundance)

chord_data <- data.frame(
  source = c(rep("Lab-reared", nrow(genus_abundance)), rep("Carrion-reared", nrow(genus_abundance))),
  target = c(genus_abundance$SimplifiedGenus, genus_abundance$SimplifiedGenus),
  value = c(genus_abundance$LabAbundance, genus_abundance$CarrionAbundance),
  stringsAsFactors = FALSE
)

chord_data <- chord_data[chord_data$value > 0, ]

cat("\n=== CHORD DIAGRAM DATA ===\n")
print(head(chord_data))

n_genera <- nrow(genus_abundance)
genera_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n_genera)
names(genera_colors) <- genus_abundance$SimplifiedGenus

condition_colors <- c("Lab-reared" = "#E41A1C", "Carrion-reared" = "#377EB8")

all_colors <- c(genera_colors, condition_colors)

png("asv_genera_chord_diagram.png", width = 16, height = 16, units = "in", res = 300, bg = "white")

circos.clear()

circos.par(
  gap.degree = c(
    rep(2, nrow(genus_abundance)),
    10, 10
  ),
  start.degree = 90
)

chordDiagram(
  chord_data,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = TRUE,
  grid.col = all_colors,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.lwd = 0.8,
  link.lty = 1,
  reduce = 0,
  group = c(
    setNames(rep("Genera", nrow(genus_abundance)), genus_abundance$SimplifiedGenus),
    c("Lab-reared" = "Conditions", "Carrion-reared" = "Conditions")
  ),
  order = c(
    genus_abundance$SimplifiedGenus,
    "Lab-reared", "Carrion-reared"
  )
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")

    if(sector.name %in% c("Lab-reared", "Carrion-reared")) {
      circos.text(
        mean(xlim),
        ylim[1] + 0.5,
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.5,
        font = 2,
        family = "cmu sans serif"
      )
    } else {
      is_exclusive <- sector.name %in% c(lab_exclusive_genera, carrion_exclusive_genera)
      font_size <- if(is_exclusive) 0.9 else 1.2

      circos.text(
        mean(xlim),
        ylim[1] + 0.1,
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = font_size,
        font = 3,
        family = "cmu sans serif"
      )
    }
  },
  bg.border = NA
)

dev.off()

tax_table_df <- as.data.frame(tax_table(ps)[shared_asvs, ])

tax_table_df <- tax_table_df %>%
  mutate(across(everything(), ~gsub("^[kpcofgs]__", "", .))) %>%
  mutate(across(everything(), ~ifelse(. == "" | is.na(.), "Unknown", .)))

phylum_counts <- table(tax_table_df$Phylum)
phylum_counts <- sort(phylum_counts, decreasing = TRUE)

genus_counts <- table(tax_table_df$Genus)
genus_counts <- sort(genus_counts, decreasing = TRUE)

cat("\n=== SHARED ASVs BY TAXONOMY ===\n")
cat("Top phyla in shared ASVs:\n")
print(head(phylum_counts, 10))

cat("\nTop genera in shared ASVs:\n")
print(head(genus_counts, 10))

top_genera <- names(head(genus_counts, 10))
top_genera_counts <- as.data.frame(genus_counts[top_genera])
colnames(top_genera_counts) <- c("Genus", "Count")

genus_plot <- ggplot(top_genera_counts, aes(x = reorder(Genus, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Top 10 Genera in Shared ASVs",
    x = "Genus",
    y = "Number of ASVs"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("shared_asvs_top_genera.png", genus_plot, width = 10, height = 6, dpi = 300)

lab_abundance <- otu_table(ps_lab)
carrion_abundance <- otu_table(ps_carrion)

lab_rel_abundance <- apply(lab_abundance, 1, function(x) mean(x / sum(x)))
carrion_rel_abundance <- apply(carrion_abundance, 1, function(x) mean(x / sum(x)))

shared_abundance <- data.frame(
  ASV = shared_asvs,
  Lab = lab_rel_abundance[shared_asvs],
  Carrion = carrion_rel_abundance[shared_asvs]
)

shared_abundance$Phylum <- tax_table_df$Phylum
shared_abundance$Genus <- tax_table_df$Genus

top_shared <- shared_abundance %>%
  mutate(TotalAbundance = Lab + Carrion) %>%
  arrange(desc(TotalAbundance)) %>%
  head(20)

abundance_plot <- ggplot(top_shared, aes(x = Lab, y = Carrion, color = Phylum, label = Genus)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(hjust = -0.2, vjust = 0.5, size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Relative Abundance of Top 20 Shared ASVs",
    x = "Lab-reared Mean Relative Abundance",
    y = "Carrion-reared Mean Relative Abundance"
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave("shared_asvs_abundance.png", abundance_plot, width = 12, height = 8, dpi = 300)

if (!requireNamespace("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
}
library(eulerr)

euler_counts <- c(
  "Lab-reared" = length(lab_only_asvs),
  "Carrion-reared" = length(carrion_only_asvs),
  "Lab-reared&Carrion-reared" = length(shared_asvs)
)

euler_plot <- euler(euler_counts)

euler_diagram <- plot(
  euler_plot,
  quantities = TRUE,
  fills = c("#E41A1C", "#377EB8"),
  alpha = 0.7,
  labels = c("Lab-reared", "Carrion-reared"),
  main = "ASV Overlap Between Lab-reared and Carrion-reared Samples"
)

png("asv_overlap_proportional.png", width = 10, height = 8, units = "in", res = 300)
plot(euler_diagram)
dev.off()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. asv_genera_chord_diagram.png - Chord diagram showing connections between top 20 genera and rearing conditions\n")
cat("2. shared_asvs_top_genera.png - Bar plot of top genera in shared ASVs\n")
cat("3. shared_asvs_abundance.png - Scatter plot of relative abundances of top shared ASVs\n")
cat("4. asv_overlap_proportional.png - Proportional Venn diagram of ASV overlap\n")
