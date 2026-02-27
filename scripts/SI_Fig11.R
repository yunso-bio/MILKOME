### SI Figure 11 — UpSet plot: Backhed vs Milkome ===========
### Author: Yunjeong So

# Load library
library(tidyverse)
library(data.table)
library(ComplexUpset)
library(curatedMetagenomicData)

# Load data
milkome_metadata <- "path to file"
milkome_abundance <- "path to file"
mapping <- "path to file"
backhed_metadata <- "path to file"
backhed_abundance <- "path to file"
grouped_cazymes <- "path to file"

# Filter Backhed samples
exc_bm_backhed <- sampleMetadata |> 
  filter(study_name == "BackhedF_2015",
         feeding_practice == "exclusively_breastfeeding",
         body_site == "stool",
         disease == "healthy",
         antibiotics_current_use == "no",
         born_method != "c_section")

# Substrate color order
order_substrates <- c(
  "Inulin" = "darkgreen", 
  "Xylan-backbone" = "#e99215",
  "Pectin-backbone" = "#6f9940",
  "β-glucans" = "#b19cd9", 
  "Food additives" = "#ddc5c1"
)

order_substrates <- c(
  "HMOs" = "darkgreen", 
  "GH20" = "pink",
  "HMOs/host glycans" = "red",
  "Host glycans" = "#e99215",
  "HMOs/plant" = "blue"
)

# Backhed
filtered_backhed <- backhed_abundance %>%
  filter(ID %in% exc_bm_backhed$NCBI_accession) %>%
  left_join(backhed_metadata, by = c("ID")) %>%
  pivot_longer(
    cols = -c(ID, Month, Sample_name, personalID),
    names_to = "cazymes", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazymes") %>% 
  mutate(
    Figure_1c = dplyr::case_when(
      cazymes %in% c("GH112", "GH136") ~ "HMOs",
      TRUE ~ Figure_1c
    )
  ) %>%
  filter(Figure_1c %in% names(order_substrates), abundance > 0) 

backhed_B_long <- filtered_backhed %>%
  filter(Month == "B") %>% 
  group_by(ID, Figure_1c) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>% 
  mutate(cohort = "backhed_B")

backhed_4M_long <- filtered_backhed %>%
  filter(Month == "4M") %>% 
  group_by(ID, Figure_1c) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>% 
  mutate(cohort = "backhed_4M")

# Milkome
milkome_metadata2 <- left_join(milkome_metadata, mapping, by = "ID")
milkome_long <- milkome_abundance %>%
  left_join(milkome_metadata2, by = c("mappedID")) %>%
  filter(infant == "I" & timepoint == "T1" & substrate == "F") %>% 
  pivot_longer(
    cols = -c(mappedID, ID, infant, participant, timepoint, substrate),
    names_to = "cazymes", values_to = "abundance") %>%
  left_join(grouped_cazymes, by = "cazymes") %>% 
  mutate(
    Figure_1c = dplyr::case_when(
      cazymes %in% c("GH112", "GH136") ~ "HMOs",
      TRUE ~ Figure_1c
    )
  ) %>%
  filter(Figure_1c %in% names(order_substrates), abundance > 0) %>% 
  group_by(ID, Figure_1c) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>% 
  mutate(cohort = "milkome")

# Combine
combined <- bind_rows(backhed_B_long, backhed_4M_long, milkome_long)
combined_wide <- combined %>%
  distinct(ID, cohort, Figure_1c) %>%
  pivot_wider(
    names_from = Figure_1c,
    values_from = Figure_1c,
    values_fill = FALSE,
    values_fn = function(x) TRUE)

# UpsetPlot
for (study in unique(combined_wide$cohort)) {
  combined_wide2 <- combined_wide %>% 
    filter(cohort == study)
  
  p <- upset(
    combined_wide2,
    names(order_substrates),
    width_ratio = 0.1,
    base_annotations=list(
      'Intersection size' = intersection_size(
        text_mapping = aes(label = paste0(!!upset_text_percentage()))
      )
    )
  ) +
    ggtitle(study) +
    theme_minimal() +                  
    theme(
      axis.title.x = element_text(),
      axis.title.y = element_text()
    )
  print(p)
}

