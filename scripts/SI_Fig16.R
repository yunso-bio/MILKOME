### SI Figure 16
### Author: Yunjeong So

# SI Figure16 a. NMDS with infants mothers------------------------------------------------------
# Load libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(paletteer)
library(data.table)

# Read data using data.table for better memory management
metadata <- "path to file"
mp_abundance <- "path to file"
mapping <- "path to file"
gene_abundance <- "path to file"
genome_abundance <- "path to file"

# Process metadata
new_metadata <- metadata %>%
  mutate(prefix = paste(infant, participant, timepoint, substrate, sep = "_")) %>%
  left_join(mapping, by = "ID")  %>% 
  filter(substrate %in% c("F", "HMOs", "MUC", "PEC", "AX", "XG", "GM", "GA")) %>% 
  mutate(age_substrate = case_when(
    substrate %in% c("GM", "GA") ~ paste(infant, "Food additives", sep = "_"),
    substrate %in% c("PEC", "AX", "XG") ~ paste(infant, "Dietary fibres", sep = "_"),
    TRUE ~ paste(infant, substrate, sep = "_")))

list <- new_metadata %>% 
  filter(!(new_metadata$prefix %in% new_metadata$prefix))

# Filter the abundance data efficiently
filtered_abundance <- mp_abundance %>%
  filter(ID %in% new_metadata$ID) %>%
  select(-Unclassified) %>%
  mutate(across(-ID, ~ round(as.numeric(as.character(.)) * 10000, digits = 0))) %>% 
  select(ID, where(~ sum(. != 0, na.rm = TRUE) > 0)) 
filtered_abundance <- as.data.frame(filtered_abundance)
rownames(filtered_abundance) <- filtered_abundance[,1]
filtered_abundance <- filtered_abundance[,-1]

# NMDS analysis using parallel processing
ga_matrix <- as.matrix(filtered_abundance)
set.seed(123)
nmds <- metaMDS(ga_matrix, distance = "bray")

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column(var = "ID") %>%
  left_join(new_metadata, by = "ID") %>%
  mutate(
    shape_group = case_when(
      infant == "I" & timepoint == "T1" ~ "T1",
      infant == "I" & timepoint == "T2" ~ "T2",
      infant == "I" & timepoint == "T3" ~ "T3",
      infant == "M" ~ "Mother",
      TRUE ~ NA_character_))

custom_shapes <- c("T1" = 21, "T2" = 24, "T3" = 22, "Mother" = 8)

# Define specific colors for the substrates
substrate_colors <- c("F" = "#5a3d35",
                      "HMOs" = "#ffcd41",  
                      "MUC" = "#82addc",
                      "AX" = "#4c9141",
                      "PEC" = "#4c9141",
                      "XG" = "#4c9141",
                      "GM" = "pink",
                      "GA" = "pink")  

# Plot with ggplot2 using custom colors
p <- ggplot(data = nmds_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(fill = as.factor(substrate), colour = as.factor(substrate), shape = shape_group), size = 3, alpha = 0.8, stroke = 0.5) +
  scale_fill_manual(values = substrate_colors) +
  scale_colour_manual(values = substrate_colors) +
  scale_shape_manual(values = custom_shapes) +
  stat_ellipse(aes(group = age_substrate), 
               type = "t", level = 0.68, alpha = 0.2) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

print(p)


# SI Figure16 b. Phyla bar plot------------------------------------------------------
# Reshape genome abundance data to long format
genome_abundance_long <- genome_abundance %>%
  pivot_longer(cols = -ID, names_to = "taxa", values_to = "genome_abundance") %>%
  separate(taxa, into = c("phylum", "species"), sep = ";") %>%
  mutate(
    phylum  = str_replace(phylum, "^p__", ""),
    species = str_replace(species, "^s__", "")
  ) 

# Function to group phyla
group_others <- function(phylum) {
  keep <- c("Actinobacteriota", "Firmicutes", "Bacteroidota",
            "Proteobacteria", "Verrucomicrobiota", "Others", "Unclassified")
  if (phylum %in% keep) return(phylum)
  "Others"
}

merged_data <- genome_abundance_long %>% 
  left_join(metadata, by = "ID")

# Filter data
filtered_data <- merged_data %>%
  filter(infant == "M", timepoint == "T2") %>%
  mutate(
    phylum = sapply(phylum, group_others),
    prefix = paste(participant, substrate, sep = "_")
  ) %>%
  group_by(phylum, prefix) %>%
  summarise(
    total_abundance = sum(genome_abundance, na.rm = TRUE),
    .groups = "drop"
  )

wide_data <- filtered_data %>%
  pivot_wider(
    names_from = phylum,
    values_from = total_abundance,
    values_fill = 0
  )

# Phylum order and colours
p_levels <- c("Actinobacteriota", "Firmicutes", "Bacteroidota",
              "Proteobacteria", "Verrucomicrobiota", "Others", "Unclassified")

p_colours <- c(
  "Actinobacteriota" = "#171b60",
  "Bacteroidota"      = "#065e06",
  "Firmicutes"        = "#a54846",
  "Proteobacteria"    = "wheat3",
  "Verrucomicrobiota" = "#5b264b",
  "Others"            = "#8a9eb4",
  "Unclassified"      = "grey80"
)

filtered_data$phylum <- factor(filtered_data$phylum, levels = p_levels)

# Plot
p <- ggplot(filtered_data, aes(x = prefix, y = total_abundance, fill = phylum)) +
  geom_col(alpha = 0.85, linewidth = 0.0, colour = NA) +
  scale_fill_manual(values = p_colours) +
  labs(y = "Relative Abundance") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text()
  ) 

print(p)
