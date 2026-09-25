#Per reviewer's request, we add a forest plot in the supplementary figure
library(ggplot2)
library(dplyr)

load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
results = res_final[which(res_final$fdr<0.05),]
results <- results %>%
  mutate(
   # HR = exp(coef),  # References the 'coef' column
   # lower_CI = exp(coef - qnorm(0.975) * `se(coef)`),  # Use backticks for special characters
   # upper_CI = exp(coef + qnorm(0.975) * `se(coef)`)
    lower_CI = coef - qnorm(0.975) * `se(coef)`,  # Use backticks for special characters
    upper_CI = coef + qnorm(0.975) * `se(coef)`
  )

# Define cancer colors
cancer_colors <- c(
  "Liver" = "#E41A1C",
  "Lung" = "#377EB8", 
  "Kidney" = "#4DAF4A",
  "Prostate" = "#984EA3"
)


# Create ordering and add cancer header rows
results_structured <- results %>%
  mutate(
    cancer = factor(cancer, levels = c("Lung", "Kidney", "Prostate", "Liver"))
  ) %>%
  arrange(cancer, desc(coef)) %>%
  mutate(
    protein_label = paste0(target, " (", `Entrez Gene Name`, ")")
  )

# Create header rows for each cancer type with n_cases
cancer_headers <- results_structured %>%
  group_by(cancer) %>%
  slice(1) %>%  # Get first row of each cancer to get n_cases
  ungroup() %>%
  distinct(cancer, n_cases) %>%
  mutate(
    is_header = TRUE,
    protein_label = paste0(as.character(cancer), " (n=", n_cases, ")"),
    coef = NA,
    lower_CI = NA,
    upper_CI = NA,
    target = NA,
    `Entrez Gene Name` = NA
  )

# Add header indicator to protein rows
results_structured <- results_structured %>%
  mutate(is_header = FALSE)

# Combine headers and data
plot_data <- bind_rows(cancer_headers, results_structured) %>%
  arrange(cancer, desc(is_header), desc(coef)) %>%
  mutate(
    y_pos = rev(row_number()),
    # Indent protein names
    display_label = if_else(is_header, 
                            protein_label,  # Use the label with n_cases for headers
                            paste0("     ", protein_label))
  )

# Check the data
print(plot_data %>% select(cancer, is_header, display_label, y_pos, coef, n_cases) %>% head(10))

# Create the forest plot
p <- ggplot(plot_data, aes(x = coef, y = y_pos)) +
  # Add vertical line at log HR = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
  # Add confidence intervals (only for non-header rows)
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = lower_CI, xmax = upper_CI, color = cancer), 
                 height = 0.3, linewidth = 0.6) +
  
  # Add point estimates (only for non-header rows)
  geom_point(data = filter(plot_data, !is_header),
             aes(color = cancer), size = 2.5, shape = 18) +
  
  # Scales and labels
  scale_color_manual(values = cancer_colors, name = "Cancer Type") +
  scale_x_continuous(
    breaks = c(-2,-1.5, -1, -0.5,0, 0.5,1,1.5, 2)
   # labels = c("-2","-1","0", "1", "2","3")
  )+
  scale_y_continuous(
    breaks = plot_data$y_pos,
    labels = plot_data$display_label,
    expand = expansion(add = c(1, 1))  # Changed from 0.5 to 1 for more space at top/bottom
  ) +
  
  # Labels
  labs(
    x = "log Hazard Ratio (95% CI)",
    y = NULL,
    title = NULL#"Significant Protein-Cancer Associations After Common and Site-Specific Risk Factor Adjustments"
  ) +
  
  # Theme
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7, hjust = 0),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 7),
    plot.margin = margin(10, 5, 5, 10)  # Add more margin (top, right, bottom, left)
  )

print(p)


ggsave("forest_plot_combined_model2.pdf", p, dpi=450,path = "./collaboration/ARIC_Cancer/manuscript/edit/Nature Comm/revision/plots/",
       width = 170, height = 170,units = "mm")
ggsave("forest_plot_combined_model2.png", p, dpi=450,path = "./collaboration/ARIC_Cancer/manuscript/edit/Nature Comm/revision/plots/",
       width = 170, height = 170,units = "mm")

