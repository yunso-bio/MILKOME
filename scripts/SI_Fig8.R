### SI Figure 8 Backhed
### Author: Yunjeong So

# Extended Figure 5. Violin plot------------------------------------------------------
# Load libraries 
library(DT)
library(curatedMetagenomicData)
library(data.table)
library(tidyverse)
library(paletteer)

# Load data
metadata <- "path to file"
gene_abundance <- "path to file"
grouped_cazymes <- "path to file"
out_dir <- "path to file"

exc_bm_backhed <- sampleMetadata |> 
  filter(study_name == "BackhedF_2015" &
           feeding_practice == "exclusively_breastfeeding" & 
           body_site =="stool" & disease == "healthy" &  
           antibiotics_current_use == "no" &
           born_method != "c_section")

order_substrates <- c(
  "HMOs" = "#252e9a", 
  "HMOs/host glycans" = "#70a6d9",
  "Host glycans" = "#a0cac7",
  "Starch" = "#ffff99", 
  "Inulin" = "darkgreen", 
  "Xylan-backbone" = "#e99215",
  "Pectin-backbone" = "#6f9940",
  "β-glucans" = "#b19cd9", 
  "β-glucans, Xylan-backbone, Chitosanase" = "#b19cd9",
  "β-glucans, Xylan-backbone, Food additives" = "#b19cd9",
  "Food additives" = "#ddc5c1")

# order months
month_order <- c("B", "4M")
new_id_levels <- c()

for (substrate in names(order_substrates)) {
  for (month in month_order) {
    new_id_levels <- c(new_id_levels, paste(substrate, month, sep = "_"))
  }
}

# Filter gene abundance data based on CAZyme_ID 
filtered_abundance <- gene_abundance %>%
  filter(ID %in% exc_bm_backhed$NCBI_accession) %>%
  left_join(metadata, by = c("ID")) %>%
  pivot_longer(cols = -c(ID, Month, Sample_name, personalID),
               names_to = "cazyme", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazyme") %>%
  mutate(
    cazyme_sub_category = recode(
      cazyme_sub_category,
      "β-glucans, Xylan-backbone, Chitosanase" = "β-glucans",
      "β-glucans, Xylan-backbone, Food additives" = "β-glucans"
    )
  ) %>%
  filter(cazyme_sub_category %in% names(order_substrates)) %>%
  mutate(new_ID = paste(cazyme_sub_category, Month, sep = "_")) %>%
  group_by(new_ID, ID, cazyme_sub_category, Month) %>%
  summarise(total_abundance = log10(sum(abundance, na.rm = TRUE) + 1), .groups = "drop") %>%
  mutate(new_ID = factor(new_ID, levels = new_id_levels))

# Calculate Wilcoxon p-values per group (new_ID)
wilcox_results <- filtered_abundance %>%
  group_by(new_ID) %>%
  summarise(
    p_value = wilcox.test(total_abundance, mu = 0, alternative = "greater")$p.value,
    .groupd = "drop") %>%
  mutate(p_label = ifelse(p_value < 0.01, 
                          formatC(p_value, format = "e", digits = 2),  
                          round(p_value, 4)))

# Plot with p-value labels
violin_p <- ggplot(filtered_abundance, aes(x = new_ID, y = total_abundance, fill = cazyme_sub_category)) +
  geom_violin(alpha = 0.5) +
  stat_summary(fun = median,
               geom = "crossbar",
               colour = "black",
               fatten = 1,
               width = 0.5) +
  geom_point(colour = "black", size = 1, alpha = 0.2, position = position_jitter(width = 0.1)) +
  scale_y_continuous(
    breaks = log10(c(0, 10, 100, 1000, 10000) + 1),
    labels = c("0", "10", "100", "1,000", "10,000")
  ) +
  scale_fill_manual(values = order_substrates) +
  theme_minimal() +
  labs(y = "Genes per million") +
  geom_text(data = wilcox_results,
            aes(x = new_ID, y = log10(10000), label = p_label),
            inherit.aes = FALSE,
            size = 4,
            color = "black")

print(violin_p)

wide <- filtered_abundance %>%
  select(ID, new_ID, total_abundance) %>%
  pivot_wider(
    names_from = new_ID,
    values_from = total_abundance,
    values_fill = 0
  )

#write_tsv(wide,file.path(out_dir, "backhed_filtered_abundance_wide.tsv"))
