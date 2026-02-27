### Figure 2
### Author: Yunjeong So

# Fig2 a. Phylum area ==========================
# Load libraries
library(tidyverse)
library(ggplot2)

# Load data
metadata <- "path to file"
taxa_abundance <- "path to file"
gene_abundance <- "path to file"
mapping <- "path to file"
grouped_cazyme <- "path to file"

# Reshape data
taxa_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  filter( ID %in% metadata$ID) %>%
  separate(taxa,
    into = c("kingdom","phylum","class","order","family","genus","species"),
    sep = ";",
    fill = "right",
    remove = TRUE) %>%
  mutate(across(kingdom:species, ~ stringr::str_remove(.x, "^[dpcofgsk]__"))
  ) 

# Join metadata + summarise at phylum level
p_levels <- c("Actinobacteriota", "Firmicutes", "Bacteroidota", "Proteobacteria", "Others", "Unclassified")

filtered_data <- metadata %>%
  inner_join(taxa_long, by = "ID") %>%
  filter(infant == "I", substrate == "F") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    phylum = fct_other(phylum, keep = c("Actinobacteriota", "Firmicutes", "Bacteroidota", "Proteobacteria", "Others", "Unclassified"), other_level = "Others"),
    phylum = factor(phylum, levels = p_levels)
  ) %>%
  group_by(participant, timepoint, phylum) %>%
  summarise(taxa_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(timepoint_numeric = as.integer(timepoint))

# Colours 
p_colours <- c(
  "Actinobacteriota" = "#171b60",
  "Bacteroidota"     = "#065e06",
  "Firmicutes"       = "#a54846",
  "Proteobacteria"   = "wheat3",
  "Others"           = "#8a9eb4",
  "Unclassified"     = "grey80"
)

# Plot 
one_tp <- filtered_data %>%
  distinct(participant, timepoint_numeric) %>%
  count(participant, name = "n_tp") %>%
  filter(n_tp == 1) %>%
  pull(participant)

p <- ggplot(filtered_data, aes(x = timepoint_numeric, y = taxa_abundance, fill = phylum)) +
  geom_area(alpha = 0.85, linewidth = 0, colour = NA) +
  geom_col(
    data = filtered_data %>% filter(participant %in% one_tp),
    width = 0.9, alpha = 0.85
  ) +
  scale_fill_manual(values = p_colours, drop = FALSE) +
  facet_wrap(~participant, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = 1:3, labels = c("T1", "T2", "T3"),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank()
  ) +
  coord_cartesian(xlim = c(0.9, 3.1))

print(p)

# Print and save plot
filename <- paste0("figures/area-Infant-faeces_others.svg")
#ggsave(plot = p, filename = filename, width = 8, height = 5, bg = "transparent")

# Fig2 b. Bifidobacterium species heatmap ==========================
# Load libraries
library(pheatmap)

# Reshape data
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", ""))) %>%
  replace_na(list(genus = "unassigned_genus", species = "unassigned_species")) %>% 
  filter(taxa_abundance >= 0.01)

# Merge metadata and taxa abundance data
merged_data <- merge(metadata, taxa_abundance_long, by = "ID")
merged_data$timepoint <- factor(merged_data$timepoint, levels = c('T1', 'T2', 'T3'))

# Define a function
class_func <- function(classes, species) {
  if (!(classes %in% c("Actinomycetia", "Bacilli", "Bacteroidia", "Clostridia", "Gammaproteobacteria", "Negativicutes", "Verrucomicrobiae", "Unclassified"))) 
  { return("Others") 
  } else if (classes == "Actinomycetia" & startsWith(as.character(species), "Bifidobacterium" ) & !(species %in% c("Bifidobacterium animalis", "Bifidobacterium angulatum")) ) 
  { return(species)
  } else 
  { return(classes)}}

phylum_func <- function(phyla, fill_taxa) {
  if (fill_taxa%in% c("Others", "Unclassified")) {
    return("Others") 
  } else {
    return(phyla)}}

# Define the order of taxa
taxa_levels <- c("Bifidobacterium longum.infantis", "Bifidobacterium longum.longum", "Bifidobacterium bifidum", "Bifidobacterium breve",  
                 "Bifidobacterium adolescentis", "Bifidobacterium angulatum", "Bifidobacterium animalis", "Bifidobacterium catenulatum", "Bifidobacterium dentium", "Bifidobacterium pseudocatenulatum",
                 "Actinomycetia", "Bacilli", "Clostridia", "Negativicutes", "Bacteroidia", "Gammaproteobacteria", "Verrucomicrobiae", "Others", "Unclassified")

# format the abundance data
bif_phyla <- merged_data %>%
  filter(infant == "I", substrate %in% c("F")) %>% 
  mutate(
    fill_taxa = mapply(class_func, classes, species),
    newID = paste(participant, timepoint, sep = "_")) %>%
  mutate(anno_phyla = mapply(phylum_func, phyla, fill_taxa)) %>%
  select(newID, fill_taxa, anno_phyla, taxa_abundance) 

bif_phyla_log <- bif_phyla %>% 
  select(newID, fill_taxa, anno_phyla, taxa_abundance) %>%
  group_by(newID, fill_taxa) %>%
  summarise(fill_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = 'drop')  %>% 
  mutate(fill_abundance = ((fill_abundance))) %>% 
  pivot_wider(names_from = fill_taxa, values_from = fill_abundance) %>%
  column_to_rownames(var = "newID") 

# Create a phyla annotation
phyla_anno <- bif_phyla %>%
  select(fill_taxa, anno_phyla) %>%               
  distinct(fill_taxa, .keep_all = TRUE) %>%  
  mutate(fill_taxa = factor(fill_taxa, levels = taxa_levels)) %>% 
  arrange(fill_taxa) %>% 
  column_to_rownames(var = "fill_taxa") 

# Set colours
my_colours <- colorRampPalette(c("snow2", "darkmagenta"))(100)
p_colours <- c(
  "Actinobacteriota" = "#171b60", 
  "Bacteroidota" = "#065e06",  
  "Firmicutes" = "#a54846", 
  "Proteobacteria" = "wheat3", 
  "Verrucomicrobiota" = "#4e193e",
  "Others" = "grey80")

# Generate heatmap with the correct phyla annotation
matrix <- as.matrix(t(bif_phyla_log))
species_names <- rownames(matrix)
matrix <- matrix[match(rownames(phyla_anno), rownames(matrix)), ]
matrix[is.na(matrix)] <- 0

# Define breaks and labels manually
breaks <- seq(0, max(matrix, na.rm = TRUE), length.out = 100)
label_matrix <- ifelse(matrix == 0, "0", "")

# Plot the heatmap
p <- pheatmap(matrix, 
              border_color = FALSE, 
              col = my_colours,
              breaks = breaks,
              cluster_rows = F, 
              cluster_cols = F, 
              annotation_row = phyla_anno, 
              annotation_colors = list(anno_phyla = p_colours), 
              gaps_col = c(3,6,9,11,14,15),
              gaps_row = c(9,12,13,14,15),
              fontsize_row = 6,         
              fontsize_col = 7,          
              angle_col = 45, 
              show_colnames = T, 
              show_rownames = T,
              annotation_legend = T,
              treeheight_row = 10,  
              treeheight_col = 20,
              display_numbers = label_matrix, 
              number_color = "grey50")

print(p)
filename <- paste0("figures/heatmap_bif_class_by_subject.svg")
#ggsave(plot = p, filename = filename, width = 5, height = 2.5, bg = "transparent")

# Fig2 c, e. Violin plot - Counts of cazyme over weaning ==========================
# Merge metadata
md_infant_F <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, sep = "_")) %>%
  left_join(mapping, by = "ID") %>%
  filter(substrate == "F" & infant == "I")

order_substrates <- c(
  "HMOs", "HMOs/host glycans", "O-glycans", "Glycosaminoglycans",
  "Starch", "Xylan-backbone", "Pectin-backbone", "β-glucans", "Food additives")

# Filter gene abundance data based on CAZyme_ID 
filtered_abundance <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>% 
  filter(!is.na(prefix)) %>%
  select(mappedID, prefix, timepoint, participant, all_of(grouped_cazyme$cazyme)) %>% 
  pivot_longer(cols = -c(mappedID, prefix, timepoint, participant), names_to = "cazyme", values_to = "abundance") %>% 
  left_join(grouped_cazyme, by = "cazyme") %>% 
  filter(cazyme_sub_category %in% order_substrates) %>% 
  group_by(prefix, timepoint, participant, mappedID, cazyme_sub_category)  %>% 
  summarise(total_abundance = (sum(abundance)), .groups = "drop") %>% 
  mutate(cazyme_sub_category = factor(cazyme_sub_category, levels = order_substrates))

participant_colours <- c("A" = "#4c89cd",
                         "B" = "#e99215",
                         "C" = "#478c2c",
                         "D" = "#d64a28",
                         "E" = "#713d91",
                         "F" = "#8c564b",
                         "G" = "#ce9ab9")

log10p1_trans <- trans_new(
  name      = "log10p1",
  transform = function(x) log10(x + 1),
  inverse   = function(x) 10^x - 1
)

p <- ggplot(filtered_abundance, aes(x = timepoint, y = total_abundance)) +
  geom_violin(fill = "#998675", alpha = 0.1) +
  stat_summary(fun = median, geom = "crossbar", colour = "black", fatten = 1) +
  geom_point(aes(colour = participant),
             position = position_dodge2(width = 0.8), size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = participant_colours) +
  scale_y_continuous(
    trans  = log10p1_trans,
    breaks = c(0, 5, 10, 100, 500, 1000, 5000),
    labels = c("0", "5", "10", "100", "500", "1000", "5000")
  ) +
  labs(y = "GPM") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    strip.background = element_blank(),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black")
  ) +
  facet_wrap(~ cazyme_sub_category, nrow = 1) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.format",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))
  )

print(p)

# Save the plot
filename <- "path to file"
#ggsave(plot = p, filename = filename, width = 13, height = 4, bg = "transparent")

# Fig2 d. Heatmap of cazyme ==========================
# Map mappedID -> prefix 
id2prefix <- md_infant_F %>%
  distinct(mappedID, prefix) %>%
  tibble::deframe()

# Set group colours
cazy_colours <- c(
  "HMOs/host glycans" = "yellow",
  "Host glycans" = "skyblue3",
  "HMOs/host glycans/plant" = "red",
  "Starch" = "#ebc4c1",
  "Dietary fibres" = "darkgreen",
  "Food additives" = "grey"
)

# Order of cazyme_sub_category
target_cazyme_sub_category <- c(
  "HMOs", "HMOs/host glycans", "HMOs/plant",
  "O-glycans", "Glycosaminoglycans", "Starch",
  "Inulin", "Xylan-backbone",
  "Pectin-backbone",
  "β-glucans", "Mannan", "Gum Arabic"
)

# Set substrate colours
substrate_colours <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(target_cazyme_sub_category)),
  target_cazyme_sub_category
)

# abundance_matrix stays RAW (no log transform)
abundance_raw <- gene_abundance %>%
  filter(mappedID %in% md_infant_F$mappedID) %>%
  column_to_rownames("mappedID") %>%
  select(all_of(grouped_cazyme$cazyme)) %>%
  mutate(across(everything(), ~ as.numeric(.)))

# Rename sample IDs 
rownames(abundance_raw) <- id2prefix[rownames(abundance_raw)]

# TOP 4 genes per subcategory
gene_totals <- abundance_raw %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "cazyme", values_to = "total_abundance")

top_genes <- grouped_cazyme %>%
  filter(cazyme_sub_category %in% target_cazyme_sub_category) %>%
  inner_join(gene_totals, by = "cazyme") %>%
  group_by(cazyme_category, cazyme_sub_category) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  pull(cazyme)

abundance_raw <- abundance_raw %>%
  select(all_of(top_genes))

abundance_matrix <- as.matrix(t(abundance_raw))
abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix <- abundance_matrix[, order(colnames(abundance_matrix)), drop = FALSE]

# Create row annotation
cazy_anno <- grouped_cazyme %>%
  as.data.frame() %>%
  filter(cazyme %in% colnames(abundance_raw)) %>%
  column_to_rownames("cazyme") %>%
  mutate(
    category = factor(cazyme_category, levels = names(cazy_colours)),
    cazyme_sub_category = factor(cazyme_sub_category, levels = target_cazyme_sub_category)
  ) %>%
  select(category, cazyme_sub_category) %>%
  arrange(category, cazyme_sub_category)

# Reorder matrix rows to match annotation order, and drop any unmatched
keep_rows <- intersect(rownames(cazy_anno), rownames(abundance_matrix))
cazy_anno <- cazy_anno[keep_rows, , drop = FALSE]
abundance_matrix <- abundance_matrix[keep_rows, , drop = FALSE]

# LOG-SCALED COLOR MAPPING
n_breaks <- 100
max_val <- max(abundance_matrix, na.rm = TRUE)
my_colours <- colorRampPalette(c("snow2", "#b22222"))(n_breaks)
breaks <- 10^(seq(log10(1), log10(max_val + 1), length.out = n_breaks)) - 1

# Legend in RAW units
legend_values <- c(0, 1, 10, 100, 500, 1000)
legend_values <- legend_values[legend_values <= max_val]
legend_breaks <- legend_values

# Annotation color mapping
annotation_colors <- list(
  category = cazy_colours,
  Figure_1e = substrate_colours
)
label_matrix <- ifelse(abundance_matrix == 0, "0", "")

# Generate the heatmap
p <- pheatmap(
  abundance_matrix,
  color = my_colours,
  breaks = breaks,
  border_color = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(3, 6, 9, 11, 14, 15),
  gaps_row = c(2, 5, 11, 12, 15),
  annotation_row = cazy_anno,
  annotation_colors = annotation_colors,
  annotation_names_row = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_legend = TRUE
  
)

# Save 
print(p)

filename <- "path to file"
#ggsave(plot = p, filename = filename, width = 6.65, height = 7, bg = "transparent")

# Fig2 f. Shannon diversity of microbiota and cazyme ==========================
all_ids <- md_infant_F %>% select(ID) %>% distinct()

# taxa abundance
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"), sep = ";", fill = "right") %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", ""))) %>%
  mutate(
    genus = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = ifelse(is.na(species) | species == "", "unassigned_species", species)
  ) %>%
  filter(taxa_abundance >= 0.01) %>%
  filter(ID %in% md_infant_F$ID)

filtered_taxa <- taxa_abundance_long %>%
  select(ID, species, taxa_abundance) %>%
  filter(species != "Unclassified") %>%
  group_by(ID, species) %>%
  summarise(total_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop")

shannon_taxa <- filtered_taxa %>%
  group_by(ID) %>%
  summarise(shannon = diversity(total_abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "taxa abundance")

# PLANT gene abundance (Starch + Dietary fibres)
plant_cazyme <- grouped_cazyme %>%
  filter(cazyme_category %in%c("Starch", "Dietary fibres"))

filtered_abundance_plant <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(plant_cazyme$cazyme)) %>%
  pivot_longer(cols = -c(ID, prefix, timepoint, participant), 
               names_to = "cazyme", values_to = "abundance") %>%
  left_join(grouped_cazyme, by = "cazyme")

shannon_gene <- filtered_abundance_plant %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "Plant gene abundance")


HG_cazyme <- grouped_cazyme %>%
  filter(cazyme_category %in%c("Host glycans"))

filtered_abundance_HG <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(HG_cazyme$cazyme)) %>%
  pivot_longer(cols = -c(ID, prefix, timepoint, participant), 
               names_to = "cazyme", values_to = "abundance") %>%
  left_join(grouped_cazyme, by = "cazyme")

shannon_HG <- filtered_abundance_HG %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "HG abundance")

# FOOD ADDITIVES gene abundance
additives_cazyme <- grouped_cazyme %>%
  filter(cazyme_category %in%c("Food additives"))

filtered_abundance_add <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(additives_cazyme$cazyme)) %>%
  pivot_longer(cols = -c(ID, prefix, timepoint, participant), 
               names_to = "cazyme", values_to = "abundance") %>%
  left_join(grouped_cazyme, by = "cazyme")

shannon_additives <- filtered_abundance_add %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(shannon = replace_na(shannon, 0),
         source = "Food additives abundance")

HMO_HG_cazyme <- grouped_cazyme %>%
  filter(cazyme_category %in%c("HMOs", "HMOs/host glycans"))

filtered_abundance_HMO_HG <- gene_abundance %>%
  left_join(md_infant_F, by = "mappedID") %>%
  filter(!is.na(prefix)) %>%
  select(ID, prefix, timepoint, participant, all_of(HMO_HG_cazyme$cazyme)) %>%
  pivot_longer(
    cols = -c(ID, prefix, timepoint, participant),
    names_to = "cazyme",
    values_to = "abundance"
  ) %>%
  left_join(grouped_cazyme, by = "cazyme")

shannon_HMO_HG <- filtered_abundance_HMO_HG %>%
  group_by(ID) %>%
  summarise(shannon = diversity(abundance, index = "shannon"), .groups = "drop") %>%
  right_join(all_ids, by = "ID") %>%
  mutate(
    shannon = replace_na(shannon, 0),
    source = "HMO + Host glycans"
  )
# COMBINE and plot
combined_shannon <- bind_rows(
  shannon_taxa,
  shannon_gene,
  shannon_additives,
  shannon_HG,
  shannon_HMO_HG
)

# Add metadata
combined_shannon <- combined_shannon %>%
  left_join(md_infant_F %>% select(ID, timepoint, participant), by = "ID")

# Plot
p <- ggplot(combined_shannon, aes(x = timepoint, y = shannon)) +
  geom_violin(alpha = 0.1, lwd = 0.5, fill = "#998675") +
  stat_summary(fun = median,
               geom = "crossbar",
               colour = "black",
               fatten = 1) +  
  geom_point(aes(colour = participant),
             position = position_jitter(width = 0.25),
             size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = participant_colours) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),
    plot.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(),
    axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(y = "Shannon Index") +
  facet_wrap(~ source, nrow = 1) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    label = "p.format",
    label.y = c(3.47, 3.59, 3.71))

print(p)

# save
#ggsave("figures/shannon_diversity_taxa_gene-with-additives.svg", plot = p, width = 6, height = 3)

# Fig2 g. NMDS with cazyme ==========================
# Taxa -> long -> species matrix
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  separate(
    taxa,
    into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"),
    sep = ";", fill = "right"
  ) %>%
  mutate(across(c(kingdom, phyla, classes, orders, family, genus, species),
                ~ str_replace(., "^[a-z]__", ""))) %>%
  mutate(
    genus   = ifelse(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = ifelse(is.na(species) | species == "", "unassigned_species", species)
  ) %>%
  filter(ID %in% md_infant_F$ID) %>%
  filter(taxa_abundance >= 0.01) %>%
  select(ID, species, taxa_abundance) %>%
  filter(species != "Unclassified")

df_abundance <- taxa_abundance_long %>%
  pivot_wider(names_from = species, values_from = taxa_abundance, values_fill = 0) %>%
  column_to_rownames("ID")

abundance_matrix <- as.matrix(df_abundance)

# NMDS
set.seed(123)
nmds <- metaMDS(abundance_matrix, distance = "bray")

data.scores <- as.data.frame(scores(nmds)$sites)
data.scores2 <- merge(data.scores, md_infant_F, by.x = 0, by.y = "ID")

# grouped_cazyme: keep only desired rows
grouped_cazyme <- grouped_cazyme %>%
  mutate(
    cazyme = as.character(cazyme),
    cazyme_category  = as.character(cazyme_category)
  )

# Filter gene abundance data to cazyme in grouped_cazyme
gene_filtered_abundance <- gene_abundance %>%
  filter(mappedID %in% md_infant_F$mappedID) %>%
  left_join(mapping, by = "mappedID") %>%
  select(ID, any_of(grouped_cazyme$cazyme)) %>%
  mutate(across(-ID, ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  mutate(across(-ID, ~ replace_na(., 0))) %>%
  column_to_rownames("ID") %>%
  select(where(~ sum(.) > 0))  # drop all-zero genes

# Envfit + BH(FDR) adjustment + join category
en <- envfit(nmds, gene_filtered_abundance, permutations = 999, na.rm = TRUE)

en_coord_cont <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cont$pvals <- en$vectors$pvals
en_coord_cont$r     <- en$vectors$r
en_coord_cont$r2    <- en_coord_cont$r^2
en_coord_cont$padj  <- p.adjust(en_coord_cont$pvals, method = "BH")

# thresholds
fdr_threshold <- 0.05
r2_threshold  <- 0.00

en_coord_cont_plot <- en_coord_cont %>%
  rownames_to_column("cazyme") %>%
  left_join(grouped_cazyme %>% select(cazyme, cazyme_category), by = "cazyme") %>%
  filter(!is.na(padj), padj < 0.05, r2 >= r2_threshold)

en_coord_cont_plot2 <- en_coord_cont %>%
  rownames_to_column("cazyme") %>%
  left_join(grouped_cazyme %>% select(cazyme, cazyme_category), by = "cazyme") %>%
  filter(!is.na(padj), padj < 0.015, r2 >= r2_threshold)

# Colours
participant_colours <- c(
  "A" = "#4c89cd", "B" = "#e99215", "C" = "#478c2c", "D" = "#d64a28",
  "E" = "#713d91", "F" = "#8c564b", "G" = "#ce9ab9"
)

# envfit colours based on grouped_cazyme$category
cazyme_colors <- c(
  "HMOs/host glycans"       = "red",
  "HMOs/host glycans/plant" = "blue",
  "Host glycans"            = "#82addc",
  "Dietary fibres"          = "#4c9141",
  "Food additives"          = "grey50",
  "Starch"                  = "orange"
)

# Plot
p <- ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(group = timepoint),
               type = "t", level = 0.68, alpha = 0.2) +
  geom_point(aes(colour = as.factor(participant), shape = timepoint),
             size = 4, alpha = 0.8, stroke = 0.5) +
  scale_colour_manual(values = participant_colours, name = "Participant") +
  
  ggnewscale::new_scale_colour() +
  
  geom_segment(
    data = en_coord_cont_plot,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, colour = cazyme_category),
    linewidth = 0.5, alpha = 0.5,
    arrow = arrow(length = unit(0.1, "cm"))
  ) +
  geom_text_repel(
    data = en_coord_cont_plot2,
    aes(x = NMDS1, y = NMDS2, colour = cazyme_category, label = cazyme),
    size = 3,
    show.legend = TRUE,
    max.overlaps = 50
  ) +
  scale_colour_manual(values = cazyme_colors, name = "CAZyme group") +
  
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

print(p)

cat("NMDS stress:", nmds$stress, "\n")
cat("Envfit vectors plotted:", nrow(en_coord_cont_plot), "\n")

# PERMANOVA
set.seed(2010)
model <- adonis2(abundance_matrix ~ timepoint,
                 data = md_infant_F,
                 method = "bray",
                 permutations = 9999,
                 by = "terms")
print(model)

md2 <- md_infant_F %>% filter(participant != "F")
abundance2 <- abundance_matrix[md2$ID, ]
adonis2(abundance2 ~ participant, data = md2)

model3 <- adonis2(abundance_matrix ~ participant + timepoint,
                  data = md_infant_F,
                  method = "bray",
                  permutations = 9999,
                  by = "terms")
print(model3)
