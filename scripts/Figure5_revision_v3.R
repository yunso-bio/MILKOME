### Figure 5
### Author: Yunjeong So

# ======= Figure5 a. Bubble plot  =======
# Load libraries
library(tidyverse)
library(RColorBrewer)

# Load tdata
metadata <- "path to the file"
taxa_abundance <- "path to the file"

# Define substrate list
sub_list <- c("F", "AX", "PEC", "XG")

# Prepare metadata
metadata <- metadata %>%
  filter(infant == "I", substrate %in% sub_list) %>% 
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
  "Bacteroidia", "Gammaproteobacteria", "Others", "Unclassified")

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

# Species list
filter_list <- merged_subset %>%
  filter((fold_change >= 1 & enrichment_abundance >= 10)) %>%
  distinct(species) %>%
  pull(species)

# Regroup
filtered_subset <- merged_subset %>%
  mutate(
    grouping = case_when(
      genus == "Enterococcus" ~ genus,
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
    fold_change = log2(enrichment_abundance + ps) - log2(F_abundance + ps),
    classes = factor(classes, levels = key_classes),
    substrate = factor(substrate, levels = sub_list)
  ) %>%
  arrange(classes, grouping) %>%
  mutate(grouping = factor(grouping, levels = rev(unique(grouping))))

# Define the color scale with grey80 at the midpoint
my_colours <- c("#003171", "white", "#ae3127")
fold_change_min <- min(filtered_subset$fold_change, na.rm = TRUE)
fold_change_max <- max(filtered_subset$fold_change, na.rm = TRUE)
fold_change_range <- c(fold_change_min, 0, fold_change_max)

# Plot with substrate and participant facets
p <- ggplot(filtered_subset, aes(x = prefix, y = grouping, size = enrichment_abundance, fill = fold_change)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(0.5, 6), limits = c(0.5, max(filtered_subset$enrichment_abundance))) +
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

file_out="path to the file"
# Save the bubble plot
#ggsave(plot = p, filename = file_out, width = 12.25, height = 4.65, bg = "transparent")

# ======= Figure5 b. CAZymes gene enrichment =======
library(data.table)
library(tidyverse)
library(ggpubr)

# Read data
metadata <- "path to the file"
gene_abundance <-"path to the file"
grouped_cazymes <- "path to the file"
mapping <- "path to the file"

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

# Define colors
AX_colors <- c(
  "GH10" = "red",
  "GH11" = "orchid",
  "GH5_21" = "#87789c",
  "GH5_35" = "#994D78",
  "GH8" = "#ed9897",
  "GH9" = "#8b5d2f")

PEC_colors <- c(
  "GH28" = "#6f8520",
  "GH105" = "#989279",
  "PL9_1" = "black",
  "PL" = "#006400")

XG_colors <- c(
  "GH74" = "#1C0CA1",
  "GH5_4" = "#5599cc",
  "GH12" = "#6FE8FD",
  "GH8" = "#ed9897",
  "GH9" = "#8b5d2f")

# combine target genes and colors
target_colors <- c(AX_colors, PEC_colors, XG_colors)
targets <- list(AX = names(AX_colors), PEC = names(PEC_colors), XG = names(XG_colors))

# merge metadata
metadata2 <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I") %>%
  mutate(
    participant = as.factor(participant),
    substrate = factor(substrate, levels = c("F", "AX", "PEC", "XG")),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    newID = paste(timepoint, participant, sep = "_")
  ) %>%
  distinct(mappedID, substrate, newID, participant, timepoint, .keep_all = TRUE)

# gene abundance to long format
abundance_long <- gene_abundance %>%
  semi_join(metadata2, by = "mappedID") %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(metadata2 %>% select(mappedID, substrate, newID, participant, timepoint), by = "mappedID") %>%
  left_join(grouped_cazymes %>% select(cazymes, Figure_5), by = "cazymes") %>%
  filter(!is.na(Figure_5), substrate %in% c("F", names(targets)))

# function 
extract_target_genes <- function(target_substrate, target_genes) {
  abundance_long %>%
    filter(Figure_5 %in% target_genes) %>%
    filter(substrate %in% c("F", target_substrate)) %>%
    group_by(newID, participant, timepoint, substrate, Figure_5) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(
      target = target_substrate,
      substrate_cmp = if_else(substrate == "F", "F", "enrichment"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "enrichment"))
    )
}

# run the function for each substrate
all_targets <- bind_rows(
  extract_target_genes("AX", targets$AX),
  extract_target_genes("PEC", targets$PEC),
  extract_target_genes("XG", targets$XG)
)

# filter out non paired samples
paired <- all_targets %>%
  group_by(target, newID, Figure_5) %>%
  filter(all(c("F", "enrichment") %in% substrate_cmp)) %>%
  ungroup() %>%
  mutate(
    prefix = paste(target, newID, sep = "_"))

# Aggregate 
all_plot <- paired %>%
  group_by(target, newID, prefix, substrate_cmp, Figure_5) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Build wide table per CAZyme
wide_cazyme <- all_plot %>%
  select(target, newID, prefix, Figure_5, substrate_cmp, abundance) %>%
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
  select(target, prefix, Figure_5)

all_plot2 <- all_plot %>%
  semi_join(keep_keys, by = c("target", "prefix", "Figure_5"))

# Keys to star (log2FC >= 1)
star_keys <- wide_cazyme %>%
  filter(!both_zero, log2FC >= 1) %>%
  select(target, prefix, Figure_5) %>%
  mutate(star = TRUE)

star_layer <- all_plot2 %>%
  filter(substrate_cmp == "enrichment") %>%
  left_join(star_keys, by = c("target", "prefix", "Figure_5")) %>%
  mutate(label = if_else(!is.na(star) & star, "*", NA_character_)) %>%
  arrange(target, prefix, Figure_5)

p <- ggplot(all_plot2, aes(x = prefix, y = abundance, fill = Figure_5)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = target_colors, drop = FALSE) +
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
    aes(label = label, group = Figure_5),
    position = position_stack(vjust = 0.4),
    na.rm = TRUE,
    size = 5,
    colour = "white"
  )

print(p)

file_out="path to the file"
# Save the bubble plot
#ggsave(plot = p, filename = file_out, width = 12, height = 4.5, bg = "transparent")


### ============ Statistics ==============
# function for summing target gene abundance
target_genes_sum <- function(target_substrate, target_genes) {
  abundance_long %>%
    filter(Figure_5 %in% target_genes) %>%
    filter(substrate %in% c("F", target_substrate)) %>%
    group_by(newID, participant, timepoint, substrate) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      target = target_substrate,
      substrate_cmp = if_else(substrate == "F", "F", "enrichment"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "enrichment"))
    )
}

# run the function for each substrate
totals_all_targets <- bind_rows(
  target_genes_sum("AX", targets$AX),
  target_genes_sum("PEC", targets$PEC),
  target_genes_sum("XG", targets$XG)
)

# filter out non paired samples
totals_paired <- totals_all_targets %>%
  group_by(target, newID) %>%
  filter(all(c("F", "enrichment") %in% substrate_cmp)) %>%
  ungroup()

# wide 
wide_for_tests <- totals_paired %>%
  select(target, newID, substrate_cmp, total) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = total,
    values_fill = list(total = NA_real_)
  ) %>%
  filter(!is.na(F) & !is.na(enrichment))

# paired wilcoxon 
pvals <- wide_for_tests %>%
  group_by(target) %>%
  group_modify(~{
    df <- .x
    tibble(
      n_pairs = nrow(df),
      p = if (nrow(df) >= 2) wilcox.test(df$F, df$enrichment, paired = TRUE, exact = FALSE)$p.value else NA_real_
    )
  }) %>%
  ungroup() %>%
  mutate(label = if_else(is.na(p), NA_character_, rstatix::p_format(p)))

# plot setting 
y_top <- 3400
ann <- pvals %>%
  filter(!is.na(p)) %>%
  mutate(
    x1 = 1,
    x2 = 2,
    y = y_top * 0.85,
    tick = y_top * 0.82,
    xmid = (x1 + x2) / 2,
    ytext = y_top * 0.90
  )

# plot 
p_F_vs_target <- ggplot(totals_paired, aes(x = substrate_cmp, y = total)) +
  geom_violin(alpha = 0.15) +
  stat_summary(fun = median, geom = "crossbar", colour = "black", width = 0.6) +
  geom_point(
    aes(colour = participant, shape = timepoint),
    position = position_jitter(width = 0.12, height = 0),
    size = 1.5, alpha = 0.85
  ) +
  scale_colour_manual(values = participant_colours, drop = FALSE) +
  facet_wrap(~ target, nrow = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right"
  ) +
  scale_y_log10(
    limits = c(10, 3500),
    breaks = c(10, 200, 500, 1000, 2000, 3500),
    labels = c("10", "200", "500", "1000", "2000", "3500")
  ) +
  labs(x = NULL, y = "Gene abundance") +
  geom_segment(data = ann, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE) +
  geom_segment(data = ann, aes(x = x1, xend = x1, y = tick, yend = y), inherit.aes = FALSE) +
  geom_segment(data = ann, aes(x = x2, xend = x2, y = tick, yend = y), inherit.aes = FALSE) +
  geom_text(
    data = ann,
    aes(x = xmid, y = ytext, label = label),
    inherit.aes = FALSE,
    vjust = 0
  )

print(p_F_vs_target)

# print the medians
medians_F_vs_target <- totals_paired %>%
  group_by(target, substrate_cmp) %>%
  summarise(
    n = n(),
    median_total = median(total, na.rm = TRUE),
    .groups = "drop"
  )

print(medians_F_vs_target, n = Inf)

