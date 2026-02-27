### Figure 3b, Figure 6a, Extended Figure 4b.
### Author: Yunjeong So

# Second Enrichment - infants------------------------------------------------------
library(tidyverse);library(reshape2);library(vegan);library(ggpubr);library(plyr);library(readr)

data_ph <- "path to file"

data_ph$sample_type <- gsub(".*_", "",data_ph$Substrate)
data_ph$sample_type <- ifelse(data_ph$sample_type == "NC", "NC", ifelse(data_ph$sample_type == "plantHMOs", "plant-HMOs", "condition"))
data_ph$Substrate <- gsub("_.*", "",data_ph$Substrate)

data_ph_m <- melt(data_ph, value.vars = c("Substrate", "sample_type"))
data_ph_m <- separate(data_ph_m, col = variable, into = c("subject", "replicate"), sep = "\\.")
data_ph_m <- data_ph_m[grep("^M", data_ph_m$subject, invert = T),]

stats <- compare_means(value ~ sample_type, group.by = c("Substrate", "subject"), data = data_ph_m, method = "t.test")

stats$p.adj <- NA
for (i in unique(stats$subject)) {
  stats[stats$subject == i,]$p.adj <- p.adjust(stats[stats$subject == i,]$p, method = "BH")}
stats$prefix <-  ifelse(stats$group1 == "NC", paste0(stats$Substrate, "_", stats$group2), "")

stats <- stats %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE             ~ "ns")) %>%
  ungroup()

summary_data <- ddply(data_ph_m, .(sample_type, Substrate, subject), summarize, mean = mean(value), SD = sd(value))

summary_data <- summary_data %>% 
  mutate(prefix = paste0(Substrate, "_", sample_type)) %>% 
  left_join(stats %>% select(prefix, subject, p, p.adj, p.adj.signif, method), 
            by = c("prefix", "subject"))

summary_data2 <- summary_data %>%
  mutate(
    label = paste0(
      ifelse(!is.na(p.adj.signif) & p.adj.signif != "ns", p.adj.signif,  "ns")), 
    prefix = factor(prefix, levels = c("HMOs_NC", "HMOs_condition",
                                       "Mucin_NC", "Mucin_condition",
                                       "AX_NC", "AX_condition", "AX_plant-HMOs",
                                       "PEC_NC", "PEC_condition", "PEC_plant-HMOs",
                                       "XG_NC", "XG_condition", "XG_plant-HMOs",
                                       "GM_NC", "GM_condition", "GM_plant-HMOs",
                                       "GA_NC", "GA_condition", "GA_plant-HMOs")))

# Create a matrix for the heatmap
mean_matrix <- reshape2::acast(summary_data2, prefix ~ subject, value.var = "mean")

# Optional: create a matrix of labels
label_matrix <- reshape2::acast(summary_data2, prefix ~ subject, value.var = "label")

# Plot with gaps after rows 2 and 4 (example)
limit <- max(abs(mean_matrix), na.rm = TRUE)
heatmap_p <- pheatmap(
  mean_matrix,
  color = colorRampPalette(c("#d1837b", "snow2", "#629fbd"))(100),
  breaks = seq(-limit, limit, length.out = 101), 
  display_numbers = label_matrix,
  gaps_row = c(2, 4, 7, 10, 13, 16),
  gaps_col = c(3, 6, 9, 11, 14, 15),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_number = 12
)

print(heatmap_p)
outdir <-  "path to file"
#ggsave(paste0(outdir, "2025.10.13_combined_heatmap_infant-pH-only.svg"), plot = heatmap_p, width = 5.45, height = 4.65)

# Second Enrichment - mothers------------------------------------------------------
data_ph <-  "path to file"

data_ph$sample_type <- gsub(".*_", "",data_ph$Substrate)
data_ph$sample_type <- ifelse(data_ph$sample_type == "NC", "NC", "condition")
data_ph$Substrate <- gsub("_.*", "",data_ph$Substrate)

data_ph_m <- melt(data_ph, value.vars = c("Substrate", "sample_type"))
data_ph_m <- separate(data_ph_m, col = variable, into = c("subject", "replicate"), sep = "\\.")
data_ph_m <- data_ph_m[grep("^I", data_ph_m$subject, invert = T),]
data_ph_m <- data_ph_m %>% 
  filter(Substrate %in% c("HMOs", "Mucin"))

stats <- compare_means(value ~ sample_type, group.by = c("Substrate", "subject"), data = data_ph_m, method = "t.test")

stats$p.adj <- NA
for (i in unique(stats$subject)) {
  stats[stats$subject == i,]$p.adj <- p.adjust(stats[stats$subject == i,]$p, method = "BH")}
stats$prefix <-  ifelse(stats$group1 == "NC", paste0(stats$Substrate, "_", stats$group2), "")

stats <- stats %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE             ~ "ns")) %>%
  ungroup()

summary_data <- ddply(data_ph_m, .(sample_type, Substrate, subject), summarize, mean = mean(value), SD = sd(value))

summary_data <- summary_data %>% 
  mutate(prefix = paste0(Substrate, "_", sample_type)) %>% 
  left_join(stats %>% select(prefix, subject, p, p.adj, p.adj.signif, method), 
            by = c("prefix", "subject"))

summary_data2 <- summary_data %>%
  mutate(
    label = paste0(
      ifelse(!is.na(p.adj.signif) & p.adj.signif != "ns", p.adj.signif,  "ns")), 
    prefix = factor(prefix, levels = c("HMOs_NC", "HMOs_condition",
                                       "Mucin_NC", "Mucin_condition")))

# Create a matrix for the heatmap
mean_matrix <- reshape2::acast(summary_data2, prefix ~ subject, value.var = "mean")

library(pheatmap);

# Plot with gaps after rows 2 and 4 (example)
limit <- max(abs(mean_matrix), na.rm = TRUE)
heatmap_p <- pheatmap(
  mean_matrix,
  color = colorRampPalette(c("#d1837b", "snow2", "#629fbd"))(100),
  breaks = seq(-limit, limit, length.out = 101), 
  display_numbers = label_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_number = 12
)

print(heatmap_p)
outdir <-  "path to file"
#ggsave(paste0(outdir, "2025.10.13_combined_heatmap_infant-pH-only.svg"), plot = heatmap_p, width = 5.45, height = 4.65)


