### SI Figrue 9
### Author: Yunjeong So

# Library
library(tidyverse)
library(data.table)
library(stringr)
library(forcats)
library(ggpubr)
library(patchwork)

# Read data
metadata <- "path to file"
gene_abundance <- "path to file"
mapping <- "path to file"
grouped_cazymes <- "path to file"
out_dir <- "path to file"

# target genes and colors
AX_colors <- c(
  "GH10" = "red",
  "GH11" = "orchid",
  "GH5_21" = "#87789c",
  "GH5_35" = "#994D78",
  "GH8" = "#ed9897",
  "xylan/cellulose/b-glucans" = "#8b5d2f",
  "Xylan-side-chains" = "#e4ddd4"
)

PEC_colors <- c(
  "GH28" = "#6f8520",
  "GH105" = "#989279",
  "PL9_1" = "black",
  "PL" = "#006400",
  "Pectin-side-chains" = "#E7E178"
)

XG_colors <- c(
  "GH74" = "#1C0CA1",
  "GH5_4" = "#5599cc",
  "GH12" = "#6FE8FD",
  "GH8" = "#ed9897",
  "xylan/cellulose/b-glucans" = "#8b5d2f",
  "Xylan-side-chains" = "#e4ddd4"
)

fill_colors <- c(
  setNames(AX_colors, paste("AX", names(AX_colors), sep = "::")),
  setNames(PEC_colors, paste("PEC", names(PEC_colors), sep = "::")),
  setNames(XG_colors, paste("XG", names(XG_colors), sep = "::"))
)

target_groups <- list(
  AX = names(AX_colors),
  PEC = names(PEC_colors),
  XG = names(XG_colors)
)

# filter metadata
metadata2 <- metadata %>%
  left_join(mapping, by = "ID") %>%
  filter(infant == "I") %>%
  mutate(
    participant = factor(participant),
    substrate = factor(substrate, levels = c("F", "AX", "PEC", "XG")),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    newID = paste(timepoint, participant, sep = "_")
  ) %>%
  distinct(mappedID, substrate, newID, participant, timepoint, .keep_all = TRUE)

# gene abundance in a long format
abundance_long <- gene_abundance %>%
  semi_join(metadata2, by = "mappedID") %>%
  pivot_longer(cols = -mappedID, names_to = "cazymes", values_to = "abundance") %>%
  left_join(
    metadata2 %>% select(mappedID, substrate, newID, participant, timepoint),
    by = "mappedID"
  ) %>%
  left_join(grouped_cazymes %>% select(cazymes, Figure_5), by = "cazymes") %>%
  filter(
    substrate %in% c("F", names(target_groups)),
      Figure_5 %in% target_groups$AX |
      Figure_5 %in% target_groups$PEC |
      Figure_5 %in% target_groups$XG
  )

# extract target gene abundance
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
      substrate_cmp = factor(if_else(substrate == "F", "F", "enrichment"),
                             levels = c("F", "enrichment")))
}

# run extraction
all_target_genes <- bind_rows(
  extract_target_genes("AX", target_groups$AX),
  extract_target_genes("PEC", target_groups$PEC),
  extract_target_genes("XG", target_groups$XG)
)

# extract all pairs
all_pairs <- all_target_genes %>%
  group_by(target, newID, cazymes) %>%
  filter(all(c("F", "enrichment") %in% substrate_cmp)) %>%
  ungroup()

# wide
wide_abundance <- all_pairs %>%
  select(target, participant, timepoint, newID, substrate_cmp, cazymes, abundance, Figure_5) %>%
  pivot_wider(
    names_from = substrate_cmp,
    values_from = abundance,
    values_fill = list(abundance = 0)) %>%
  mutate(
    fold_change = log2((enrichment + 1) / (F + 1)),
    target = factor(target, levels = c("AX", "PEC", "XG")),
    key = paste(target, Figure_5, sep = "::") ) %>%
  group_by(target, cazymes) %>%
  mutate(n_pairs = n()) %>%
  ungroup()

# Filter out genes that never appear in both conditions within each target
wide_for_plot <- wide_abundance %>%
  group_by(target, cazymes) %>%
  filter(any(enrichment > 0, na.rm = TRUE) & any(F > 0, na.rm = TRUE)) %>%
  ungroup()

# One-sample Wilcoxon test on fold_change vs 0 (median fold_change > 0)
paired_tests <- wide_abundance %>%
  filter(cazymes %in% wide_for_plot$cazymes) %>% 
  group_by(target, cazymes, Figure_5) %>%
  summarise(
    n_pairs = n(),
    median_fc = median(fold_change, na.rm = TRUE),
    mean_fc = mean(fold_change, na.rm = TRUE),
    p_value = if (n() < 3) NA_real_ else {
      suppressWarnings(
        wilcox.test(fold_change, mu = 0, alternative = "two.sided", exact = FALSE)$p.value)},
    .groups = "drop") %>%
  group_by(target) %>%
  mutate(
    P.adj = p.adjust(p_value, method = "BH"),
    is_sig = !is.na(P.adj) & P.adj < 0.05
  ) %>%
  ungroup()

# stat format
ann_tbl <- paired_tests %>%
  mutate(
    label = paste0(
      "Med. = ",
      if_else(is.finite(median_fc), sprintf("%.2f", median_fc), "NA"),
      "  P = ",
      if_else(is.finite(p_value), sprintf("%.3g", p_value), "NA"),
      "\n.                    P.adj=",
      if_else(is.finite(P.adj), sprintf("%.3g", P.adj), "NA")
    ),
    label_col = if_else(is_sig & median_fc > 0, "red", "black")
  )


# make plots
make_target_plot <- function(df_plot, target_name, group_levels, ann_tbl_target, ymax) {
  df_t <- df_plot %>%
    filter(target == target_name) %>%
    mutate(group_rank = match(as.character(Figure_5), group_levels)) %>%
    arrange(group_rank, cazymes, newID) %>%
    mutate(
      newID = factor(newID, levels = unique(newID)),
      cazymes_in_target = factor(
        paste(target, cazymes, sep = " :: "),
        levels = unique(paste(target, cazymes, sep = " :: "))
      )
    )
  
  ann_facet <- ann_tbl_target %>%
    mutate(
      cazymes_in_target = paste(target_name, cazymes, sep = " :: "),
      cazymes_in_target = factor(cazymes_in_target, levels = levels(df_t$cazymes_in_target))
    )
  
  ggplot(df_t, aes(x = newID, y = fold_change, fill = key)) +
    geom_col(alpha = 0.85, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_text(
      data = ann_facet,
      aes(x = 1, y = ymax, label = label, color = label_col),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.2,
      size = 5.5
    ) +
    coord_cartesian(ylim = c(-ymax, ymax), clip = "off") +
    scale_color_identity() +
    scale_fill_manual(values = fill_colors) +
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
    facet_wrap(~cazymes_in_target, nrow = 10)
}

# Extract data 
ann_AX <- ann_tbl %>% filter(target == "AX")
ann_PEC <- ann_tbl %>% filter(target == "PEC")
ann_XG <- ann_tbl %>% filter(target == "XG")

# plot
ymax <- ceiling(max(abs(wide_for_plot$fold_change), na.rm = TRUE) * 1.05)
p_AX <- make_target_plot(wide_for_plot, "AX", names(AX_colors), ann_AX, ymax)
p_PEC <- make_target_plot(wide_for_plot, "PEC", names(PEC_colors), ann_PEC, ymax)
p_XG <- make_target_plot(wide_for_plot, "XG", names(XG_colors), ann_XG, ymax)

p_all <- (p_AX | p_PEC | p_XG)
print(p_all)
