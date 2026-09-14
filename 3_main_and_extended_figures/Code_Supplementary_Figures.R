#Code for Supplementary Figures

library(ggplot2)
library(dplyr)

# Get list of all 40 proteins identified in Model 1 and 2
protein_list <- unique(id)
#edit and change ACP3 and KLK3 to be males only
replace_protein = c("SeqId_8468_19","SeqId_19317_114")
protein_list <- setdiff(protein_list,replace_protein)
selected_ids <- v2$SampleId
# Create long format data for ALL proteins
create_long_data <- function(protein_list, selected_ids) {
  all_data <- data.frame()
  
  for(protein in protein_list) {
    protein_data <- rbind(
      data.frame(SampleId = v2$SampleId, 
                 Visit = "V2", 
                 Protein_Value = v2[[protein]],
                 Protein = protein),
      data.frame(SampleId = v3$SampleId, 
                 Visit = "V3", 
                 Protein_Value = v3[[protein]],
                 Protein = protein),
      data.frame(SampleId = v5$SampleId, 
                 Visit = "V5", 
                 Protein_Value = v5[[protein]],
                 Protein = protein)
    )
    all_data <- rbind(all_data, protein_data)
  }
  
  return(all_data %>%
           filter(SampleId %in% selected_ids) %>%
           mutate(Visit = factor(Visit, levels = c("V2", "V3", "V5"))))
}

# Create the combined dataset
all_long_data <- create_long_data(protein_list, selected_ids)

all_long_data$Protein_Target = paste0(anno$`Entrez Gene Name`[match(all_long_data$Protein,anno$id)]," (",anno$`UniProt ID`[match(all_long_data$Protein,anno$id)],")")


create_long_data_male <- function(protein_list, selected_ids) {
  all_data <- data.frame()
  
  for(protein in protein_list) {
    protein_data <- rbind(
      data.frame(SampleId = v2_male$SampleId, 
                 Visit = "V2", 
                 Protein_Value = v2_male[[protein]],
                 Protein = protein),
      data.frame(SampleId = v3_male$SampleId, 
                 Visit = "V3", 
                 Protein_Value = v3_male[[protein]],
                 Protein = protein),
      data.frame(SampleId = v5_male$SampleId, 
                 Visit = "V5", 
                 Protein_Value = v5_male[[protein]],
                 Protein = protein)
    )
    all_data <- rbind(all_data, protein_data)
  }
  
  return(all_data %>%
           filter(SampleId %in% selected_ids) %>%
           mutate(Visit = factor(Visit, levels = c("V2", "V3", "V5"))))
}

# Create the combined dataset
male_long_data <- create_long_data_male(replace_protein, male_ids)

male_long_data$Protein_Target = paste0(anno$`Entrez Gene Name`[match(male_long_data$Protein,anno$id)]," (",anno$`UniProt ID`[match(male_long_data$Protein,anno$id)],")")

all_long_data = rbind(all_long_data,male_long_data)


# Calculate summary statistics for each protein
summary_stats_all <- all_long_data %>%
  group_by(Protein_Target, Visit) %>%
  summarise(
    Mean = mean(Protein_Value, na.rm = TRUE),
    SD = sd(Protein_Value, na.rm = TRUE),
    .groups = "drop"
  )


# Create the vertical panel plot
panel_plot <- ggplot(all_long_data, aes(x = Visit, y = Protein_Value)) +
  # Individual lines
  geom_line(aes(group = SampleId), alpha = 0.3, color = "lightblue", size = 0.3) +
  # Box plots at each visit
  geom_boxplot(aes(group = Visit), alpha = 0.6, width = 0.4, 
               outlier.shape = NA, fill = "darkblue", size = 0.3) +
  # Mean line
  geom_line(data = summary_stats_all, aes(x = Visit, y = Mean, group = 1), 
            color = "red", size = 0.8) +
  geom_point(data = summary_stats_all, aes(x = Visit, y = Mean), 
             color = "red", size = 1.5) +
  # Facet by protein - 5 columns, 8 rows for vertical page
  facet_wrap(~Protein_Target, scales = "free_y", ncol = 5, nrow = 8) +
  labs(title = "Longitudinal Profiles Across Cancer-Risk Associated Proteins",
       subtitle = paste("n =", length(selected_ids), "individuals per protein"),
       x = "Visit",
       y = "Protein ANML standardized log2(RFU) level") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9),           # Protein names
    axis.text = element_text(size = 7),            # Axis labels
    axis.title = element_text(size = 11),          # Axis titles
    plot.title = element_text(size = 14),          # Main title
    plot.subtitle = element_text(size = 12),       # Subtitle
    panel.spacing = unit(0.3, "lines")             # Space between panels
  )

# Display the plot
print(panel_plot)

# Save as high-resolution image optimized for vertical layout
ggsave("./collaboration/ARIC_Cancer/manuscript/edit/top_proteins_longitudinal_final.png", plot = panel_plot, 
       width = 12, height = 20, dpi = 320)  # Taller than wide for vertical page

# Alternative: Save as PDF for better quality in documents
ggsave("./collaboration/ARIC_Cancer/manuscript/edit/top_proteins_longitudinal_final.pdf", plot = panel_plot, 
       width = 12, height = 20, dpi = 320)





#----------------------------------------------------------------------------------

#heatmap
#only keep males for ACP3 and KLK3
library(reshape2)
library(viridis)

# Extract just the correlation values for visualization
correlation_matrix_data <- combined_wide_anno %>%
  select(Protein,Correlation.cor_V2V3, Correlation.cor_V2V5, Correlation.cor_V3V5) %>%
  melt(id.vars = "Protein", variable.name = "Comparison", value.name = "Correlation") %>%
  mutate(Comparison = case_when(
    Comparison == "Correlation.cor_V2V3" ~ "V2 vs V3",
    Comparison == "Correlation.cor_V2V5" ~ "V2 vs V5", 
    Comparison == "Correlation.cor_V3V5" ~ "V3 vs V5"
  ))

correlation_matrix_data_top = correlation_matrix_data[which(correlation_matrix_data$Protein %in% protein_list),]
correlation_matrix_data_top$Protein_Target = paste0(anno$`Entrez Gene Name`[match(correlation_matrix_data_top$Protein,anno$id)]," (",anno$`UniProt ID`[match(correlation_matrix_data_top$Protein,anno$id)],")")

load("./collaboration/ARIC_Cancer/manuscript/edit/protein_cor_v2v3v5_male.RData")
correlation_matrix_data_male <- combined_wide_male_anno %>%
  select(Protein,Correlation.cor_V2V3, Correlation.cor_V2V5, Correlation.cor_V3V5) %>%
  melt(id.vars = "Protein", variable.name = "Comparison", value.name = "Correlation") %>%
  mutate(Comparison = case_when(
    Comparison == "Correlation.cor_V2V3" ~ "V2 vs V3",
    Comparison == "Correlation.cor_V2V5" ~ "V2 vs V5", 
    Comparison == "Correlation.cor_V3V5" ~ "V3 vs V5"
  ))

correlation_matrix_data_top_male = correlation_matrix_data_male[which(correlation_matrix_data_male$Protein %in% replace_protein),]
correlation_matrix_data_top_male$Protein_Target = paste0(anno$`Entrez Gene Name`[match(correlation_matrix_data_top_male$Protein,anno$id)]," (",anno$`UniProt ID`[match(correlation_matrix_data_top_male$Protein,anno$id)],")")

correlation_matrix_data_top = rbind(correlation_matrix_data_top,correlation_matrix_data_top_male)

# Order proteins by V2 vs V3 correlation (largest to smallest)
protein_order = correlation_matrix_data_top$Protein_Target[which(correlation_matrix_data_top$Comparison == "V2 vs V3")][order(correlation_matrix_data_top$Correlation[which(correlation_matrix_data_top$Comparison == "V2 vs V3")], decreasing = F)]

# Update correlation matrix data with ordered factor
correlation_matrix_data_top <- correlation_matrix_data_top %>%
  mutate(Protein_Target = factor(Protein_Target, levels = protein_order))


# Create heatmap
correlation_heatmap <- ggplot(correlation_matrix_data_top, 
                              aes(x = Comparison, y = Protein_Target, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  geom_text(aes(label = round(Correlation, 2)), 
            size = 2.5, 
            color = "black",#ifelse(abs(correlation_matrix_data_top$Correlation) > 0.5, "white", "black"), 
            fontface = "bold") +
  scale_fill_gradient2(name = "Correlation\nCoefficient", 
                       low = "white", 
                       mid = "red", 
                       high = "darkred",
                       midpoint = 0.5,  
                       limits = c(0, 1),  
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       labels = c("0.0", "0.25", "0.50", "0.75", "1.0")) +
  labs(title = "Correlations of Cancer-Risk Associated Proteins Across Visits",
       x = "Visit Comparison",
       y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 14),
    legend.title = element_text(size = 10)
  )

print(correlation_heatmap)

ggsave("correlation_heatmap_final.png", 
       path="./collaboration/ARIC_Cancer/manuscript/edit/",
       width = 8, height = 10, dpi = 320)
ggsave("correlation_heatmap_final.pdf", 
       path="./collaboration/ARIC_Cancer/manuscript/edit/",
       width = 8, height = 10, dpi = 320)


