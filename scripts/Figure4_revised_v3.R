### Figure 4
### Author: Yunjeong So

#====== Figure 4a Bubble plot (F, HMOs, Mucin)-=======
library(tidyverse)
library(RColorBrewer)

# Load data
metadata <- "path to file"
taxa_abundance <- "path to file"

# Substrates
sub_list <- c("F", "HMOs", "MUC")

# Prepare metadata
metadata <- metadata %>%
  filter(substrate %in% sub_list, infant == "I") %>%
  mutate(
    participant = as.factor(participant),
    prefix = paste(participant, timepoint, substrate, sep = "_"),
    newID  = paste(participant, timepoint, sep = "_"),
    substrate = factor(substrate, levels = sub_list),
    timepoint = factor(timepoint, levels = c("T1","T2","T3"))
  ) %>%
  arrange(timepoint, match(substrate, sub_list), participant) %>%
  mutate(prefix = factor(prefix, levels = unique(prefix)))

# Class order
key_classes <- c(
  "Actinomycetia", "Coriobacteriia",
  "Bacilli", "Clostridia", "Negativicutes",
  "Bacteroidia",
  "Gammaproteobacteria",
  "Others", "Unclassified"
)

# Long-format genome abundance
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(-ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  filter(taxa_abundance >= 0.01, ID %in% metadata$ID) %>%
  separate(
    taxa,
    into = c("kingdom","phyla","classes","orders","family","genus","species"),
    sep = ";",
    fill = "right") %>%
  mutate(across(kingdom:species, ~ str_remove(.x, "^[dpcofgs]__"))) %>%
  mutate(
    genus   = str_remove(genus, "_.*"),
    species = str_remove(species, "_[A-Za-z]+"),
    genus   = if_else(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = if_else(is.na(species) | species == "", "unassigned_species", species),
    classes = if_else(classes %in% key_classes, classes, "Others")
  )

# Merge metadata + abundances
merged_data <- metadata %>%
  left_join(taxa_abundance_long, by = "ID")

# Faeces baseline for fold-change
F_data <- merged_data %>%
  filter(substrate == "F") %>%
  group_by(newID, classes, genus, species) %>%
  summarise(F_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop")

# Include ALL substrates in plotted data
all_data <- merged_data %>%
  group_by(prefix, newID, participant, substrate, timepoint, classes, genus, species) %>%
  summarise(enrichment_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop")

# Merge + fold change
merged_subset <- all_data %>%
  left_join(F_data, by = c("newID","classes","genus","species")) %>%
  mutate(
    F_abundance     = coalesce(F_abundance, 0),
    enrichment_abundance = coalesce(enrichment_abundance, 0),
    fold_change     = log2(enrichment_abundance + 1) - log2(F_abundance + 1)
  )

# Name infantis and longum 
target_infantis <-"Bifidobacterium longum.infantis"
target_longum <- "Bifidobacterium longum.longum"
append_species <- c(target_longum, target_infantis)

# Species list
filter_list <- merged_subset %>%
  filter((fold_change >= 1 & enrichment_abundance >= 10 | (species %in% c("Bifidobacterium breve", target_longum)))) %>%
  distinct(species) %>%
  pull(species)

# Regroup
filtered_subset <- merged_subset %>%
  mutate(
    grouping = case_when(
      species %in% append_species ~ species,
      genus == "Enterococcus" ~ genus,
      genus == "Streptococcus" ~ genus,
      classes == "Bacilli" & genus == "Clostridium" ~ paste0("z-", classes),
      classes == "Negativicutes" ~ classes,
      !species %in% filter_list ~ paste0("z-", classes),
      TRUE ~ species)
  ) %>%
  group_by(prefix, newID, classes, substrate, grouping) %>%
  summarise(
    enrichment_abundance = sum(enrichment_abundance),
    F_abundance     = sum(F_abundance),
    .groups = "drop"
  ) %>%
  mutate(
    fold_change = log2(enrichment_abundance + 1) - log2(F_abundance + 1),
    classes = factor(classes, levels = key_classes),
    substrate = factor(substrate, levels = sub_list)
  ) %>%
  arrange(classes, grouping) %>%
  mutate(grouping = factor(grouping, levels = rev(unique(grouping))))

# Keep original order, then append longum + infantis AFTER it
base_levels <- levels(filtered_subset$grouping)
base_levels_no_append <- base_levels[!base_levels %in% append_species]
final_levels <- c(base_levels_no_append, append_species)

# Facet group: with vs without infantis
infantis_prefix <- merged_subset %>%
  filter(species == target_infantis, enrichment_abundance > 0) %>%
  distinct(prefix) %>%
  pull(prefix)

infantis_prefix_FHMO <- merged_subset %>%
  filter(substrate %in% c("F","HMOs"),
         species == target_infantis,
         enrichment_abundance > 0) %>%
  distinct(prefix) %>%
  pull(prefix)


filtered_subset <- filtered_subset %>%
  mutate(
    grouping = factor(as.character(grouping), levels = final_levels),
    substrate_box = if_else(substrate == "MUC", "MUC", "F + HMOs") %>%
      factor(levels = c("F + HMOs", "MUC")),
    infantis_group = case_when(
      substrate == "MUC" ~ "NA",
      prefix %in% infantis_prefix_FHMO ~ "+ infantis",
      TRUE ~ "- infantis"
    ) %>% factor(levels = c("+ infantis", "- infantis", "NA")),
    
    facet_panel = paste(substrate_box, infantis_group, sep = " | ") %>%
      factor(levels = c(
        "F + HMOs | + infantis",
        "F + HMOs | - infantis",
        "MUC | NA"
      ))
  )

# Colour scale
my_colours <- c("#003171", "white", "#ae3127")
fc_min <- min(filtered_subset$fold_change, na.rm = TRUE)
fc_max <- max(filtered_subset$fold_change, na.rm = TRUE)

# Plot
p <- ggplot(filtered_subset,
            aes(x = prefix, y = grouping, size = enrichment_abundance, fill = fold_change)) +
  geom_point(shape = 21, colour = "black", stroke = 0.3) +
  scale_size(range = c(0.5, 6)) +
  scale_fill_gradientn(
    colours = my_colours,
    limits = c(fc_min, fc_max),
    values = scales::rescale(c(fc_min, 0, fc_max))
  ) +
  facet_grid(. ~ facet_panel, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(size = "Relative abundance (%)", fill = "Log2 fold change", x = NULL) +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7, face = "italic"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "white", colour = "black", size = 1)
  )

print(p)

file_out="path to file"
#ggsave(filename = file_out, plot = p,  device = "svg",  width = 12, height = 4,  units = "in",  bg = "transparent")

#===== Figure 4c HMO gene abundance in stacked bar with Wilcoxon (infantis vs non-infantis) ===========
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

# Build B. infantis presence/absence by sample
taxa_abundance_long <- taxa_abundance %>%
  pivot_longer(-ID, names_to = "taxa", values_to = "taxa_abundance") %>%
  filter(taxa_abundance >= 0.01, ID %in% metadata$ID) %>%
  separate(
    taxa,
    into = c("kingdom","phyla","classes","orders","family","genus","species"),
    sep = ";",
    fill = "right") %>%
  mutate(across(kingdom:species, ~ str_remove(.x, "^[dpcofgs]__"))) %>%
  mutate(
    genus   = str_remove(genus, "_.*"),
    species = str_remove(species, "_[A-Za-z]+"),
    genus   = if_else(is.na(genus) | genus == "", "unassigned_genus", genus),
    species = if_else(is.na(species) | species == "", "unassigned_species", species)
  )

# merge metadata
metadata <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I", substrate %in% substrate_levels) %>%
  mutate(
    participant = as.factor(participant),
    substrate = factor(substrate, levels = substrate_levels),
    timepoint = factor(timepoint, levels = timepoint_levels),
    participant = factor(participant, levels = names(participant_colours)),
    newID = paste(timepoint, participant, sep = "_")
  ) %>%
  distinct(mappedID, substrate, newID, participant, timepoint, .keep_all = TRUE)

# Merge metadata + abundances
merged_data <- metadata %>%
  left_join(taxa_abundance_long %>%
              filter(ID %in% metadata$ID), by = "ID")

# create infantis groups
infantis_ids <- merged_data %>%
  filter(substrate %in% c("F","HMOs"),
         species == "Bifidobacterium longum.infantis",
         taxa_abundance > 0) %>%
  pull(ID) 

infantis_group_df <- metadata %>%
  distinct(ID, substrate) %>%
  left_join(mapping, by = "ID") %>%
  transmute(
    mappedID,
    infantis_group = factor(
      case_when(
        substrate == "MUC" ~ "NA",
        ID %in% infantis_ids ~ "+ infantis",
        TRUE ~ "- infantis"),
      levels = c("+ infantis", "- infantis", "NA")
    )
  )

metadata2 <- metadata %>%
  left_join(infantis_group_df, by = "mappedID")

# Gene abundance long table
abundance_long <- gene_abundance %>%
  filter(mappedID %in% metadata2$mappedID) %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(metadata2,by = "mappedID") %>%
  filter(cazymes %in% hmos_genes) 

# function 
extract_target_genes <- function(target_substrate, target_genes) {
  abundance_long %>%
    filter(substrate %in% c("F", target_substrate)) %>%
    group_by(newID, participant, timepoint, substrate, cazymes, infantis_group) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(
      target = target_substrate,
      substrate_cmp = if_else(substrate == "F", "F", "enrichment"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "enrichment"))
    )
}

# run for HMOs ######
all_targets <- bind_rows(
  extract_target_genes("HMOs", hmos_genes)
)

# filter out non paired samples
paired <- all_targets %>%
  group_by(target, newID, cazymes, infantis_group) %>%
  filter(n_distinct(substrate_cmp) == 2) %>% 
  ungroup() %>%
  mutate(prefix = paste(target, newID, sep = "_")) 

# Aggregate 
all_plot <- paired %>%
  group_by(target, newID, prefix, substrate_cmp, cazymes, infantis_group) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Build wide table per CAZyme
wide_cazyme <- all_plot %>%
  select(target, newID, prefix, cazymes, substrate_cmp, abundance, infantis_group) %>%
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

p_HMOs <- ggplot(all_plot2, aes(x = prefix, y = abundance, fill = cazymes)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = hmos_colors, drop = FALSE) +
  facet_wrap(target ~ infantis_group + substrate_cmp, scales = "free_x", space = "free_x") +
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

print(p_HMOs)

# Run for mucin #####
all_targets <- bind_rows(
  extract_target_genes("MUC", hmos_genes)
)

# filter out non paired samples
paired <- all_targets %>%
  group_by(target, newID, cazymes) %>%
  filter(n_distinct(substrate_cmp) == 2) %>% 
  ungroup() %>%
  mutate(prefix = paste(target, newID, sep = "_")) 

# Aggregate 
all_plot <- paired %>%
  group_by(target, newID, prefix, substrate_cmp, cazymes) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Build wide table per CAZyme
wide_cazyme <- all_plot %>%
  select(target, newID, prefix, cazymes, substrate_cmp, abundance) %>%
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

p_MUC <- ggplot(all_plot2, aes(x = prefix, y = abundance, fill = cazymes)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = hmos_colors, drop = FALSE) +
  facet_wrap( ~ substrate_cmp, scales = "free_x", space = "free_x") +
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

print(p_MUC)

out_file="path to file"
#ggsave(filename = out_file, plot = p,  width = 7, height = 3, bg = "transparent")


# violin plot for comparisoin between infantis groups
hmo_utilisation <- abundance_long %>%
  group_by(mappedID, timepoint, substrate, participant, infantis_group) %>%
  summarise(
    hmo_utilisation_total = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

infantis_violin <- ggplot(
  hmo_utilisation,
  aes(x = infantis_group, y = hmo_utilisation_total)
) +
  geom_violin(fill = "#998675", alpha = 0.1) +
  stat_summary(fun = median, geom = "crossbar", colour = "black", fatten = 1) +
  geom_point(
    aes(shape = timepoint, fill = participant),
    colour = "black",
    position = position_jitter(width = 0.12, height = 0),
    size = 1.5, alpha = 0.8
  ) +
  scale_fill_manual(values = participant_colours, drop = FALSE) +
  scale_shape_manual(values = c(T1 = 21, T2 = 22, T3 = 24)) +
  scale_y_log10(
    breaks = c(1, 6, 11, 101, 501, 1001, 5001),
    labels = c("0", "5", "10", "100", "500", "1000", "5000")
  ) +
  labs(x = NULL, y = "HMO utilisation genes (GPM)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    strip.background = element_blank(),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black")
  ) +
  facet_wrap(~ substrate, nrow = 1) +
  stat_compare_means(
    comparisons = list(c("+ infantis", "- infantis")),
    method = "wilcox.test",
    method.args = list(paired = FALSE, exact = FALSE),
    label = "p"
  )

print(infantis_violin)

median_table <- hmo_utilisation %>%
  group_by(substrate, infantis_group) %>%
  summarise(
    median_hmo = median(hmo_utilisation_total, na.rm = TRUE),
    .groups = "drop"
  )

print(median_table)

#===== Figure 4d Differences in HMO gene abundance among substrates Wilcoxon (F vs HMOs vs Mucin) ===========
library(tidyverse)
library(ggpubr)

# Read data
metadata        <- "path to file"
gene_abundance  <- "path to file"
mapping         <- "path to file"
taxa_abundance  <- "path to file"

# Settings
substrate_levels <- c("F", "HMOs", "MUC")

# Selected HMO utilisation genes
hmos_genes <- c("GH112","GH136","GH20","GH42","GH29","GH95","GH33")

participant_colours <- c("A" = "#4c89cd",
                         "B" = "#e99215",
                         "C" = "#478c2c",
                         "D" = "#d64a28",
                         "E" = "#713d91",
                         "F" = "#8c564b",
                         "G" = "#ce9ab9")

# Metadata processing
metadata2 <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I", substrate %in% substrate_levels) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("T1","T2","T3")),
    substrate = factor(substrate, levels = substrate_levels)
  )

# Gene level long table
abundance_long <- gene_abundance %>%
  filter(mappedID %in% metadata2$mappedID) %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(
    metadata2 %>% select(mappedID, timepoint, substrate, participant),
    by = "mappedID"
  ) %>%
  filter(cazymes %in% hmos_genes) %>%
  mutate(
    cazymes   = factor(cazymes, levels = hmos_genes),
    timepoint = factor(timepoint, levels = c("T1","T2","T3")),
    substrate = factor(substrate, levels = substrate_levels)
  )

# Merge selected genes into one HMO utilisation metric per sample
hmo_utilisation <- abundance_long %>%
  group_by(mappedID, timepoint, substrate, participant) %>%
  summarise(
    hmo_utilisation_total = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Plot 
sub_comparisons <- list(c("F","HMOs"), c("HMOs","MUC"), c("F","MUC"))

p_substrate <- ggplot(hmo_utilisation, aes(x = substrate, y = hmo_utilisation_total)) +
  geom_violin(fill = "#998675", alpha = 0.1) +
  stat_summary(fun = median, geom = "crossbar", colour = "black", fatten = 1) +
  geom_point(
    aes(shape = timepoint, fill = participant),
    colour = "black",
    position = position_jitter(width = 0.12, height = 0),
    size = 1.5, alpha = 0.8
  ) +
  scale_fill_manual(values = participant_colours) +
  scale_shape_manual(values = c(T1 = 21, T2 = 22, T3 = 24)) +  # <- key line
  scale_y_log10(
    breaks = c(1, 6, 11, 101, 501, 1001, 5001),
    labels = c("0", "5", "10", "100", "500", "1000", "5000")
  ) +
  labs(x = NULL, y = "HMO utilisation genes (GPM)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    strip.background = element_blank(),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black")
  ) +
  stat_compare_means(
    comparisons = sub_comparisons,
    method = "wilcox.test",
    method.args = list(paired = FALSE, exact = FALSE),
    label = "p.format",
    label.y.npc = "top",
    step.increase = 0.10,
    accuracy = 1e-5
  )

print(p_substrate)

hmo_utilisation %>%
  group_by(substrate) %>%
  summarise(median = median(hmo_utilisation_total, na.rm = TRUE))

out_file="path to file"
#ggsave(filename = out_file, plot = p_substrate,  width = 3, height = 1.5, bg = "transparent")

