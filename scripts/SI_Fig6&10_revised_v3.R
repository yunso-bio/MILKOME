### SI Figure 6 & 10
### Author: Yunjeong So

# SI Figure6 c. Bubble plot------------------------------------------------------
# Load  libraries
library(tidyverse)
library(RColorBrewer)

# Load the metadata and genome abundance data
metadata <- "path to file"
taxa_abundance <- "path to file"
gene_abundance <- "path to file"
mapping <- "path to file"
grouped_cazymes <- "path to file"
out_dir <- "path to file"

# Define substrate list
sub_list <- c("F", "GM", "GA")

# Prepare metadata
metadata <- metadata %>%
  filter(substrate %in% sub_list) %>% 
  mutate(
    participant = as.factor(participant),
    prefix = paste(substrate, participant, timepoint,  sep = "_"), 
    newID = paste(participant, timepoint, sep = "_")
  ) %>%
  arrange(timepoint, match(substrate, sub_list), participant) %>%  
  mutate(prefix = factor(prefix, levels = unique(paste(substrate, participant, timepoint, sep = "_")))) 

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
  summarise(non_F_abundance = sum(taxa_abundance, na.rm = TRUE), .groups = "drop")

# Merge + fold change
ps = 1
merged_subset <- all_data %>%
  left_join(F_data, by = c("newID","classes","genus","species")) %>%
  mutate(
    F_abundance     = coalesce(F_abundance, 0),
    non_F_abundance = coalesce(non_F_abundance, 0),
    fold_change     = log2(non_F_abundance + ps) - log2(F_abundance + ps)
  )

# species list
filter_list <- merged_subset %>%
  filter((fold_change >= 1 & non_F_abundance >= 5)) %>%
  distinct(species) %>%
  pull(species)

# Regroup
filtered_subset <- merged_subset %>%
  filter(substrate != "F") %>% 
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
    non_F_abundance = sum(non_F_abundance),
    F_abundance     = sum(F_abundance),
    .groups = "drop"
  ) %>%
  mutate(
    fold_change = log2(non_F_abundance + ps) - log2(F_abundance + ps),
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
p <- ggplot(filtered_subset, aes(x = prefix, y = grouping, size = non_F_abundance, fill = fold_change)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(0.5, 6), limits = c(0.5, max(filtered_subset$non_F_abundance))) +
  scale_fill_gradientn(
    colours = my_colours,
    limits = c(fold_change_min, fold_change_max),
    values = scales::rescale(fold_change_range, to = c(0, 1))
  ) +
  theme_minimal() +
  labs(size = "Relative abundance (%)", fill = "Log2 fold change") +
  theme(
    GMis.text.x = element_text(angle = 45, hjust = 1),
    GMis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "white", color = "black", size = 1)
  ) 

# Display the plot
print(p)

# SI Figure6 f. CAZymes gene enrichment ------------------------------------------------------
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
GM_colors <- c(
  "GH5_7" = "#363700",
  "GH5_8" = "#4A8A4C",
  "GH5_41" = "#E5EAE7",
  "GH26" = "#9AC3AD",
  "GH113" = "#9C6769",
  "GH130_1" = "#2B4C40",
  "GH130_2" = "#D8C978",
  "GH130_3" = "#515B50",
  "GH130_4" = "#FDD4D3",
  "GH130_5" = "#1F2020")

GA_colors <- c(
  "GH39" = "#D5D254",
  "GH43_24" = "#A46F53",
  "GH154" = "#7AA1C3",
  "PL27" = "#6B7566",
  "PL42" = "#4656A2")

# combine target genes and colors
target_colors <- c(GM_colors, GA_colors)
targets <- list(GM = names(GM_colors), GA = names(GA_colors))

# merge metadata
metadata2 <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I") %>%
  mutate(
    participant = as.factor(participant),
    substrate = factor(substrate, levels = c("F", "GM", "GA")),
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
      substrate_cmp = if_else(substrate == "F", "F", "Target"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "Target"))
    )
}

# run the function for each substrate
totals_all_targets <- bind_rows(
  extract_target_genes("GM", targets$GM),
  extract_target_genes("GA", targets$GA)
)

# Pairing filter + prefix + factor levels
totals_paired <- totals_all_targets %>%
  group_by(target, newID, Figure_5) %>%
  filter(all(c("F", "Target") %in% substrate_cmp)) %>%
  ungroup() %>%
  mutate(
    prefix = paste(target, newID, sep = "_"),
    Figure_5 = factor(Figure_5, levels = names(target_colors))
  )

# Aggregate to one row per segment
totals_plot <- totals_paired %>%
  group_by(target, newID, prefix, substrate_cmp, Figure_5) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Build wide table per CAZyme to:
wide_cazyme <- totals_plot %>%
  select(target, newID, prefix, Figure_5, substrate_cmp, abundance) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = abundance,
    values_fill = 0
  ) %>%
  mutate(
    both_zero = (F == 0 & Target == 0),
    log2FC = log2(Target + 1) - log2(F + 1)
  )

# Filter out both-zero CAZymes FROM THE PLOT too
keep_keys <- wide_cazyme %>%
  filter(!both_zero) %>%
  select(target, prefix, Figure_5)

totals_plot2 <- totals_plot %>%
  semi_join(keep_keys, by = c("target", "prefix", "Figure_5"))

# Keys to star (log2FC >= 1), after removing both-zero
star_keys <- wide_cazyme %>%
  filter(!both_zero, log2FC >= 1) %>%
  select(target, prefix, Figure_5) %>%
  mutate(star = TRUE)

star_layer <- totals_plot2 %>%
  filter(substrate_cmp == "Target") %>%
  left_join(star_keys, by = c("target", "prefix", "Figure_5")) %>%
  mutate(label = if_else(!is.na(star) & star, "*", NA_character_)) %>%
  arrange(target, prefix, Figure_5)

p <- ggplot(totals_plot2, aes(x = prefix, y = abundance, fill = Figure_5)) +
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

### Statistics ###
# function for summing target gene abundance
target_genes_sum <- function(target_substrate, target_genes) {
  abundance_long %>%
    filter(Figure_5 %in% target_genes) %>%
    filter(substrate %in% c("F", target_substrate)) %>%
    group_by(newID, participant, timepoint, substrate) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      target = target_substrate,
      substrate_cmp = if_else(substrate == "F", "F", "Target"),
      substrate_cmp = factor(substrate_cmp, levels = c("F", "Target"))
    )
}

# run the function for each substrate
totals_all_targets <- bind_rows(
  target_genes_sum("GM", targets$GM),
  target_genes_sum("GA", targets$GA)
)

# filter out non paired samples
totals_paired <- totals_all_targets %>%
  group_by(target, newID) %>%
  filter(all(c("F", "Target") %in% substrate_cmp)) %>%
  ungroup()

# wide 
wide_for_tests <- totals_paired %>%
  select(target, newID, substrate_cmp, total) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = total,
    values_fill = list(total = NA_real_)
  ) %>%
  filter(!is.na(F) & !is.na(Target))

# paired wilcoxon 
pvals <- wide_for_tests %>%
  group_by(target) %>%
  group_modify(~{
    df <- .x
    tibble(
      n_pairs = nrow(df),
      p = if (nrow(df) >= 2) wilcox.test(df$F, df$Target, paired = TRUE, exact = FALSE)$p.value else NA_real_
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
    GMis.line = element_line(color = "black"),
    GMis.ticks = element_line(color = "black"),
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

# SI Figrue 10 ------------------------------------------------------
GM_colors <- c(
  "GH5_7" = "#363700",
  "GH5_8" = "#4A8A4C",
  "GH26" = "#9AC3AD",
  "GH113" = "#9C6769",
  "GH130_1" = "#2B4C40",
  "GH130_2" = "#D8C978",
  "GH130_3" = "#515B50",
  "GH130_4" = "#FDD4D3"
)

GA_colors <- c(
  "GH39" = "#D5D254",
  "GH43_24" = "#A46F53",
  "GH154" = "#7AA1C3"
)

fill_colors <- c(
  setNames(GM_colors, paste("GM", names(GM_colors), sep = "::")),
  setNames(GA_colors, paste("GA", names(GA_colors), sep = "::"))
)

target_groups <- list(
  GM = names(GM_colors),
  GA = names(GA_colors)
)

metadata2 <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I") %>%
  mutate(
    participant = factor(participant),
    substrate = factor(substrate, levels = c("F", "GM", "GA")),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    newID = paste(timepoint, participant, sep = "_")
  ) %>%
  distinct(mappedID, substrate, newID, participant, timepoint, .keep_all = TRUE)

abundance_long <- gene_abundance %>%
  semi_join(metadata2, by = "mappedID") %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(metadata2 %>% select(mappedID, substrate, newID, participant, timepoint), by = "mappedID") %>%
  left_join(grouped_cazymes %>% select(cazymes, Figure_5), by = "cazymes") %>%
  filter(
    substrate %in% c("F", names(target_groups)),
    Figure_5 %in% c(target_groups$GM, target_groups$GA)
  )

extract_target_genes <- function(target_substrate, target_group_levels) {
  abundance_long %>%
    filter(
      Figure_5 %in% target_group_levels,
      substrate %in% c("F", target_substrate)
    ) %>%
    group_by(newID, participant, timepoint, substrate, cazymes, Figure_5) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    mutate(
      target = target_substrate,
      substrate_cmp = factor(if_else(substrate == "F", "F", "enrichment"), levels = c("F", "enrichment"))
    )
}

all_target_genes <- bind_rows(
  extract_target_genes("GM", target_groups$GM),
  extract_target_genes("GA", target_groups$GA)
)

all_pairs <- all_target_genes %>%
  group_by(target, newID, cazymes) %>%
  filter(all(c("F", "enrichment") %in% substrate_cmp)) %>%
  ungroup()

wide_abundance <- all_pairs %>%
  select(target, participant, timepoint, newID, substrate_cmp, cazymes, abundance, Figure_5) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = abundance,
    values_fill = list(abundance = 0)
  ) %>%
  mutate(
    fold_change = log2((enrichment + 1) / (F + 1)),
    target = factor(target, levels = c("GM", "GA")),
    key = paste(target, Figure_5, sep = "::")
  ) %>%
  group_by(target, cazymes) %>%
  mutate(n_pairs = n()) %>%
  ungroup()

wide_for_plot <- wide_abundance %>%
  group_by(target, cazymes) %>%
  filter(any(enrichment > 0, na.rm = TRUE) & any(F > 0, na.rm = TRUE)) %>%
  ungroup()

thresh <- 1

samples_gt_thresh <- wide_for_plot %>%
  filter(fold_change >= thresh) %>%
  distinct(target, newID, participant, timepoint) %>%
  arrange(target, newID)

print(samples_gt_thresh, n = Inf)

paired_tests <- wide_abundance %>%
  filter(cazymes %in% wide_for_plot$cazymes) %>%
  group_by(target, cazymes, Figure_5) %>%
  summarise(
    n_pairs = n(),
    median_fc = median(fold_change, na.rm = TRUE),
    mean_fc = mean(fold_change, na.rm = TRUE),
    p_value = if (n() < 3) NA_real_ else {
      suppressWarnings(wilcox.test(fold_change, mu = 0, alternative = "two.sided", exact = FALSE)$p.value)
    },
    .groups = "drop"
  ) %>%
  group_by(target) %>%
  mutate(
    P.adj = p.adjust(p_value, method = "BH"),
    is_sig = !is.na(P.adj) & P.adj < 0.05
  ) %>%
  ungroup()

ann_tbl <- paired_tests %>%
  mutate(
    label = paste0(
      "Med. = ", if_else(is.finite(median_fc), sprintf("%.2f", median_fc), "NA"),
      "   P = ", if_else(is.finite(p_value), sprintf("%.3g", p_value), "NA"),
      "   P.adj=", if_else(is.finite(P.adj), sprintf("%.3g", P.adj), "NA")
    ),
    label_col = if_else(is_sig & median_fc > 0, "red", "black")
  )

make_target_plot <- function(df_plot, target_name, group_levels, ann_tbl_target, ymax, thresh = 1, show_per_sample = TRUE) {
  df_t <- df_plot %>%
    filter(target == target_name) %>%
    mutate(group_rank = match(as.character(Figure_5), group_levels)) %>%
    arrange(group_rank, cazymes, newID, fold_change) %>%
    mutate(
      newID = factor(newID, levels = unique(newID)),
      cazymes_in_target = factor(
        paste(target, cazymes, sep = " :: "),
        levels = unique(paste(target, cazymes, sep = " :: "))
      )
    )
  
  gene_counts <- df_t %>%
    group_by(cazymes_in_target) %>%
    summarise(
      n_gt_gene = n_distinct(newID[fold_change >= thresh]),
      n_total_gene = n_distinct(newID),
      .groups = "drop"
    ) %>%
    mutate(label_gene = paste0("n (log2FC≥", thresh, ") = ", n_gt_gene, " / ", n_total_gene))
  
  sample_counts_target <- df_t %>%
    group_by(newID) %>%
    summarise(
      n_gt_sample = sum(fold_change >= thresh, na.rm = TRUE),
      n_total_sample = sum(is.finite(fold_change)),
      .groups = "drop"
    ) %>%
    mutate(label_sample = paste0(n_gt_sample, "/", n_total_sample))
  
  ann_facet <- ann_tbl_target %>%
    mutate(
      cazymes_in_target = paste(target_name, cazymes, sep = " :: "),
      cazymes_in_target = factor(cazymes_in_target, levels = levels(df_t$cazymes_in_target))
    )
  
  p <- ggplot(df_t, aes(x = newID, y = fold_change, fill = key)) +
    geom_col(alpha = 0.85, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_text(
      data = ann_facet,
      aes(x = 1, y = ymax, label = label, color = label_col),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.2,
      size = 3
    ) +
    geom_text(
      data = gene_counts,
      aes(x = Inf, y = ymax, label = label_gene),
      inherit.aes = FALSE,
      hjust = 1.02, vjust = 1.2,
      size = 3
    ) +
    coord_cartesian(ylim = c(-ymax, ymax), clip = "off") +
    scale_color_identity() +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    labs(title = target_name, x = NULL, y = "log2 fold change (enrichment vs F)") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      plot.margin = margin(5.5, 25, 5.5, 5.5),
      legend.position = "none"
    ) +
    facet_wrap(~cazymes_in_target, nrow = 2)
  
  if (isTRUE(show_per_sample)) {
    p <- p +
      geom_text(
        data = sample_counts_target,
        aes(x = newID, y = -ymax, label = label_sample),
        inherit.aes = FALSE,
        vjust = -0.2,
        size = 2.3
      )
  }
  
  p
}

ann_GM <- ann_tbl %>% filter(target == "GM")
ann_GA <- ann_tbl %>% filter(target == "GA")

ymax <- ceiling(max(abs(wide_for_plot$fold_change), na.rm = TRUE) * 1.05)

p_GM <- make_target_plot(wide_for_plot, "GM", names(GM_colors), ann_GM, ymax, thresh = thresh, show_per_sample = TRUE)
p_GA <- make_target_plot(wide_for_plot, "GA", names(GA_colors), ann_GA, ymax, thresh = thresh, show_per_sample = TRUE)

p_all <- (p_GM | p_GA)
print(p_all)
