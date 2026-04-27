### Figure 7
### Author: Yunjeong So

# Figure7 b. Bubble plot------------------------------------------------------
# Load necessary libraries
library(tidyverse)
library(RColorBrewer)

# Load the metadata and taxa abundance data
metadata <- "path to file"
taxa_abundance <- "path to file"

# Define substrate list
sub_list <- c("F", "HMOs", "MUC")

# Prepare metadata
metadata <- metadata %>%
  filter(substrate %in% sub_list, infant == "M", timepoint == "T2") %>% 
  mutate(
    participant = as.factor(participant),
    prefix = paste(participant, timepoint, substrate, sep = "_"), 
    newID = paste(participant, timepoint, sep = "_")
  ) %>%
  arrange(timepoint, match(substrate, sub_list), participant) %>%  
  mutate(prefix = factor(prefix, levels = unique(paste(participant, timepoint, substrate, sep = "_")))) 

# Define class order and key genus
key_classes <- c(
  "Actinomycetia", "Coriobacteriia", 
  "Bacilli", "Clostridia", "Negativicutes", 
  "Bacteroidia", "Verrucomicrobiae", "Gammaproteobacteria", "Others", "Unclassified")

# Reshape taxa abundance data to long format
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(-ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  filter(taxa_abundance >= 0.01, ID %in% metadata$ID) %>% 
  separate(taxa, into = c("kingdom", "phyla", "classes", "orders", "family", "genus", "species"),
           sep = ";", fill = "right") %>%
  mutate(across(kingdom:species, ~str_remove(.x, "^[dpcofgs]__"))) %>%
  mutate(
    genus = str_remove(genus, "_.*"),
    species = str_remove(species, "_[A-Za-z]+"),
    genus = if_else(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = if_else(is.na(species) | species == "", "unassigned_species", species)
  )

# Merge metadata and taxa abundance data
merged_data <- metadata %>%
  left_join(taxa_abundance_long, by = "ID") %>%
  mutate(
    timepoint = factor(timepoint, levels = c('T1', 'T2', 'T3'))
  )

# Filter faeces
F_data <- merged_data %>%
  filter(substrate == "F") %>%
  group_by(newID, classes, genus, species) %>%  
  summarise(F_abundance = sum(taxa_abundance, na.rm = TRUE)) %>%
  ungroup()

# Filter HMOs and mucin group
enrichment_data <- merged_data %>%
  filter(substrate %in% sub_list) %>%
  group_by(prefix, newID, participant, timepoint, classes, genus, species, substrate) %>%
  summarise(enrichment_abundance = sum(taxa_abundance, na.rm = TRUE)) %>%
  ungroup()

# Merge faeces and non-faeces data for fold change calculation
merged_subset <- enrichment_data %>% 
  left_join(F_data, by = c("newID", "classes", "genus", "species")) %>% 
  mutate(
    F_abundance = coalesce(F_abundance, 0),
    enrichment_abundance = coalesce(enrichment_abundance, 0))

# Calculate fold change (log2 transformation)
filtered_subset <- merged_subset %>%
  mutate(
    enrichment_abundance = round(as.numeric(enrichment_abundance), 2),  
    F_abundance = round(as.numeric(F_abundance), 2),  
    fold_change = round(log2(enrichment_abundance + 1) - log2(F_abundance + 1), 2),
    substrate = factor(substrate, levels = sub_list),
    classes = factor(classes, levels = key_classes)
  ) 

# Species list
filter_list <- filtered_subset %>%
  filter(fold_change >= 1, enrichment_abundance >= 1) %>%
  distinct(species) %>%
  pull(species)

# Filter data
filtered_subset2 <- merged_subset %>%
  mutate(
    classes = if_else(classes %in% key_classes, classes, 'Others'),
    grouping = case_when(
      classes == "Bacilli" & genus == "Clostridium" ~ paste0('z-', classes),
      genus %in% c("Enterococcus", "Streptococcus") ~ paste0(genus),
      !species %in% filter_list ~ paste0('z-', classes),
      classes == "Gammaproteobacteria" ~ paste0('z-', classes),
      TRUE ~ species)) %>%
  group_by(prefix, classes, substrate, grouping) %>%
  summarise(
    enrichment_abundance = sum(enrichment_abundance),
    F_abundance = sum(F_abundance),
    .groups = "drop") %>%
  mutate(
    fold_change = log2(enrichment_abundance + 1) - log2(F_abundance + 1),
    classes = factor(classes, levels = key_classes),
    substrate = factor(substrate, levels = sub_list)) %>%
  arrange(classes, grouping) %>% 
  mutate(grouping = factor(grouping, levels = rev(unique(grouping))))


# Define the color scale with grey80 at the midpoint
my_colours <- c("#003171", "white", "#ae3127")
fold_change_min <- min(filtered_subset$fold_change, na.rm = TRUE)
fold_change_max <- max(filtered_subset$fold_change, na.rm = TRUE)
fold_change_range <- c(fold_change_min, 0, fold_change_max)

# Plot with substrate and participant facets
p <- ggplot(filtered_subset2, aes(x = prefix, y = grouping, size = enrichment_abundance, fill = fold_change)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(0.5, 6), limits = c(0.5, max(filtered_subset2$enrichment_abundance))) +
  scale_fill_gradientn(
    colours = my_colours,
    limits = c(fold_change_min, fold_change_max),
    values = scales::rescale(fold_change_range, to = c(0, 1))
  ) +
  theme_minimal() +
  labs(size = "Relative abundance (%)", fill = "Log2 fold change") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "white", color = "black", size = 1)
  ) 

# Display the plot
print(p)

# Save the bubble plot
#ggsave(plot = p, filename = "figures/2025.mother-bubbles.0.01.svg", width = 6, height = 8, bg = "transparent")


# Figure7 c. CAZymes gene enrichment ------------------------------------------------------
library(tidyverse)
library(ggpubr)

# Read data
metadata        <- "path to file"
gene_abundance  <- "path to file"
mapping         <- "path to file"
taxa_abundance  <- "path to file"
out_dir         <- "path to file"

# Settings
substrate_levels <- c("F", "HMOs", "MUC")
timepoint_levels <- c("T1", "T2", "T3")

# Selected HMO utilisation genes
hmos_genes <- c("GH112","GH136","GH20","GH42","GH29","GH95","GH33")

# Participant colours
participant_colours <- c(
  "A" = "#4c89cd",
  "B" = "#e99215",
  "C" = "#478c2c",
  "D" = "#d64a28",
  "E" = "#713d91",
  "F" = "#8c564b",
  "G" = "#ce9ab9"
)

# CAZy colours
hmos_colors <- c(
  "GH112" = "#ffa500",
  "GH136" = "dodgerblue",
  "GH20"  = "pink",
  "GH42"  = "#aa98a9",
  "GH29"  = "#a12232",
  "GH95"  = "#af593e",
  "GH33"  = "purple4"
)

# Process metadata
metadata <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "M", substrate %in% substrate_levels, timepoint =="T2") %>%
  mutate(
    timepoint = factor(timepoint, levels = timepoint_levels),
    substrate = factor(substrate, levels = substrate_levels),
    participant = factor(participant, levels = names(participant_colours))
  )

# Gene abundance long table
abundance_long <- gene_abundance %>%
  filter(mappedID %in% metadata$mappedID) %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(
    metadata %>% select(mappedID, timepoint, substrate, participant),
    by = "mappedID"
  ) %>%
  filter(cazymes %in% hmos_genes) %>%
  mutate(
    cazymes   = factor(cazymes, levels = hmos_genes),
    timepoint = factor(timepoint, levels = timepoint_levels),
    substrate = factor(substrate, levels = substrate_levels)
  )

# function 
extract_target_genes <- function(target_substrate, target_genes) {
  abundance_long %>%
    filter(substrate %in% c("F", target_substrate)) %>%
    group_by(participant, timepoint, substrate, cazymes) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(
      target = target_substrate,
      substrate_cmp = if_else(substrate == "F", "F", "enrichment"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "enrichment"))
    )
}

# run
all_targets <- bind_rows(
  extract_target_genes("HMOs", hmos_genes),
  extract_target_genes("MUC", hmos_genes)
)

# filter out non paired samples
paired <- all_targets %>%
  group_by(target, participant, cazymes) %>%
  filter(n_distinct(substrate_cmp) == 2) %>% 
  ungroup() %>%
  mutate(prefix = paste(target, participant, sep = "_")) 

# Aggregate 
all_plot <- paired %>%
  group_by(target, participant, prefix, substrate_cmp, cazymes) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Build wide table per CAZyme
wide_cazyme <- all_plot %>%
  select(target, participant, prefix, cazymes, substrate_cmp, abundance) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = abundance,
    values_fill = 0
  ) %>%
  mutate(
    both_zero = (F == 0 & enrichment == 0),
    log2FC = log2(enrichment + 1) - log2(F + 1)
  )

# Filter out both-zero CAZymes
keep_keys <- wide_cazyme %>%
  filter(!both_zero) %>%
  select(target, prefix, cazymes)

all_plot2 <- all_plot %>%
  semi_join(keep_keys, by = c("target", "prefix", "cazymes")) %>% 
  mutate(cazymes = factor(cazymes, levels = hmos_genes))

# Keys to star (log2FC >= 1)
star_keys <- wide_cazyme %>%
  filter(!both_zero, log2FC >= 1) %>%
  select(target, prefix, cazymes) %>%
  mutate(star = TRUE)

star_layer <- all_plot2 %>%
  filter(substrate_cmp == "enrichment") %>%
  left_join(star_keys, by = c("target", "prefix", "cazymes")) %>%
  mutate(label = if_else(!is.na(star) & star, "*", NA_character_)) %>%
  mutate(cazymes = factor(cazymes, levels = hmos_genes)) %>%
  arrange(target, prefix, cazymes)

p <- ggplot(all_plot2, aes(x = prefix, y = abundance, fill = cazymes)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", liparticipantth = 0.25) +
  scale_fill_manual(values = hmos_colors, drop = FALSE) +
  facet_wrap(target ~ substrate_cmp, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  geom_text(
    data = star_layer,
    aes(label = label),
    position = position_stack(vjust = 0.4),
    na.rm = TRUE,
    size = 5,
    colour = "white"
  )

print(p)

out_file="path to file"
#ggsave(filename = out_file, plot = p_bar,  width = 7, height = 3, bg = "transparent")


# Merge selected genes into one HMO utilisation metric per sample
hmo_utilisation <- abundance_long %>%
  group_by(mappedID, timepoint, substrate, participant) %>%
  filter(substrate != "MUC") %>% 
  summarise(
    hmo_utilisation_total = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

w <- hmo_utilisation %>%
  filter(substrate %in% c("F","HMOs")) %>%
  select(participant, timepoint, substrate, hmo_utilisation_total) %>%
  pivot_wider(names_from = substrate,
              values_from = hmo_utilisation_total) %>%
  filter(!is.na(F) & !is.na(HMOs))

# Paired Wilcoxon
wilcox.test(w$F, w$HMOs, paired = TRUE, exact = FALSE)

# Medians (paired samples only)
cat("Median F:", median(w$F, na.rm = TRUE), "\n")
cat("Median HMOs:", median(w$HMOs, na.rm = TRUE), "\n")
