library(ggplot2)
library(ggrepel)
library(ggpubr)
tab = read.csv("./collaboration/ARIC_Cancer/manuscript/UKB/table2.csv") #data downloaded from Papier et al., Nat Comm
colnames(tab)[c(4,5,6,7)] = c("HR","lower_CI","upper_CI","pvalue")
tab$Cancer.site =  gsub(" Cancer", "", tab$Cancer.site)


load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
df = res_final[which(res_final$fdr < 0.05),]

#merge two datasets
df_plot = merge(df,tab, by.x = c("Entrez Gene Name","cancer"),by.y = c("Protein","Cancer.site"))

# Define significance: for example, based on study 1 p-values
df_plot$significant <- df_plot$pvalue < 0.05
colnames(df_plot)[1]="protein"

# # Compute correlation and concordance
# correlation <- cor(df_plot$HR.x, df_plot$HR.y)
# concordant <- sum(sign(df_plot$HR.x) == sign(df_plot$HR.y))
# concordance_rate <- concordant / nrow(df_plot)


# Step 1: Compute correlation and concordance per cancer type
stat_labels <- df_plot %>%
  mutate(
    sign_x = sign(HR.x - 1),
    sign_y = sign(HR.y - 1),
    concordant = sign_x == sign_y & sign_x != 0  # exclude HR=1 (neutral)
  ) %>%
  group_by(cancer) %>%
  summarise(
    correlation = cor(HR.x, HR.y, use = "complete.obs", method = "pearson"),
    concordance = mean(concordant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r: %.2f\nConcordance: %.2f", correlation, concordance),
    x = 2.2,
    y = 0.1
  )

df_plot$significant = factor(df_plot$significant, levels = c("TRUE","FALSE"))

df_plot$cancer <- factor(df_plot$cancer, 
                      levels = c("Liver", "Kidney", "Lung"))

remove_facets <- function(plot, layout) {
  layout <- strsplit(layout, split = '\n')[[1]]
  layout <- lapply(layout, trimws)
  layout <- matrix(unlist(sapply(layout, strsplit, "")),
                   nrow = length(layout), byrow = T)
  layout <- which(layout == "#", arr.ind = TRUE)
  prm <- apply(layout,1,\(x) {
    c(glue::glue("panel-{x[1]}-{x[2]}"),
      glue::glue("strip-t-{x[2]}-{x[1]}"))
  })
  g <- ggplot2::ggplotGrob(plot)
  rm_grobs <- g$layout$name %in% prm
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  ggpubr::as_ggplot(g)
}

# Step 2: Plot with facets and in-panel stats
p1 = ggplot(df_plot, aes(x = HR.x, y = HR.y, color = significant, label = protein)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(0,3.2)) +
  scale_y_continuous(breaks = c(0,1,2,3), limits = c(0,3.2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_text_repel(show.legend = FALSE, max.overlaps = 10) + 
  scale_color_manual(values = c("dark blue", "dark red")) +
  theme_minimal(base_size = 12) +
  labs(
    x = "HR (ARIC)",
    y = "HR (UKB)",
    color = "UKB P < 0.05",
    title = NULL
  ) +
  facet_wrap(~ cancer, scales = "fixed") +
  geom_text(
    data = stat_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "gray90")
    
  )

#with EPIC
#results comparison with EPIC
setwd("./collaboration/ARIC_Cancer/manuscript/EPIC_validation") #results received from EPIC collaborators
library(data.table)
prostate = read.csv("prostate_validation.csv")
lung = read.csv("lung_validation.csv")
kidney = read.csv("kidney_validation.csv")
liver = read.csv("liver_res.csv")


liver <- liver %>% mutate(
  HR = exp(`logHR`),
  lower_CI = `HRmin`,
  upper_CI = `HRmax`,
  pvalue = pval,
  cancer = "liver",
  seqid_in_sample = ID)


kidney$z = kidney$log.HR./kidney$robust.se.log.HR..
kidney$pvalue = 2 * (1 - pnorm(abs(kidney$z)))
format(kidney$pvalue, scientific = TRUE)
kidney$pvalue[1] = c(3.54e-19)



kidney <- kidney %>% mutate(
  HR = exp(`log.HR.`),
  lower_CI = exp(`log.HR.` - 1.96 * `robust.se.log.HR..`),
  upper_CI = exp(`log.HR.` + 1.96 * `robust.se.log.HR..`))

lung <- lung %>% mutate(
  HR = exp(`log.HR.`),
  lower_CI = exp(`log.HR.` - 1.96 * `robust.se`),
  upper_CI = exp(`log.HR.` + 1.96 * `robust.se`))

prostate = prostate[which(prostate$Model == "Fully Adjusted"),]
prostate = prostate[which(prostate$EntrezGeneSymbol %in% c("ACP3","KLK3")),]
prostate$cancer = "prostate"

final_out = rbind(kidney[,c("seqid_in_sample","HR","lower_CI","upper_CI","pvalue","cancer")],
                  liver[,c("seqid_in_sample","HR","lower_CI","upper_CI","pvalue","cancer")],
                  lung[,c("seqid_in_sample","HR","lower_CI","upper_CI","pvalue","cancer")],
                  prostate[,c("seqid_in_sample","HR","lower_CI","upper_CI","pvalue","cancer")])

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

final_out$cancer = capitalize_first(final_out$cancer)

load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
df = res_final[which(res_final$fdr < 0.05),]

#merge two datasets
df_plot = merge(df,final_out, by= c("seqid_in_sample","cancer"))

# Define significance: for example, based on study 1 p-values
df_plot$significant <- df_plot$pvalue < 0.05
colnames(df_plot)[12]="protein"
#df_plot$protein = paste0(df_plot$protein,"\n",df_plot$cancer)


stat_labels <- df_plot %>%
  mutate(
    sign_x = sign(HR.x - 1),
    sign_y = sign(HR.y - 1),
    concordant = sign_x == sign_y & sign_x != 0  # exclude HR=1 (neutral)
  ) %>%
  group_by(cancer) %>%
  summarise(
    correlation = cor(HR.x, HR.y, use = "complete.obs", method = "pearson"),
    concordance = mean(concordant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r: %.2f\nConcordance: %.2f", correlation, concordance),
    x = 2.2,
    y = 0.1
  )

df_plot$significant = factor(df_plot$significant, levels = c("TRUE","FALSE"))

# Step 2: Plot with facets and in-panel stats
p2 = ggplot(df_plot, aes(x = HR.x, y = HR.y, color = significant, label = protein)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(0,3.2)) +
  scale_y_continuous(breaks = c(0,1,2,3), limits = c(0,3.2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_text_repel(show.legend = FALSE, max.overlaps = 10) + 
  scale_color_manual(values = c("dark blue", "dark red")) +
  theme_minimal(base_size = 12) +
  labs(
    x = "HR (ARIC)",
    y = "HR (EPIC)",
    color = "EPIC P < 0.05",
    title = NULL
  ) +
  facet_wrap(~ cancer, scales = "fixed",nrow=1) +
  geom_text(
    data = stat_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "gray90")
    
  )


# Create completely empty plot
empty_plot <- ggplot() + theme_void()

# Remove legend from p1
p1_no_legend <- p1 + theme(legend.position = "none")

legend_plot <- as_ggplot(get_legend(p1))
white_space <- ggplot() + theme_void()

# Stack white space above legend
right_panel <- ggarrange(white_space, legend_plot, 
                         ncol = 2, widths = c(0.5, 0.5))

# Combine using ggarrange

p <-ggarrange(p2, ggarrange(p1, empty_plot,ncol=2,labels=NULL,widths = c(0.75,0.25),common.legend = T,legend="right"),
                         nrow = 2, labels = c("a", "b"),
                         heights = c(0.5,0.5))

in.dir = "./collaboration/ARIC_Cancer/manuscript/EPIC_validation"
ggsave(filename=paste0("Figure4.pdf"),
       plot=p, device="pdf",
       path=in.dir,
       width=450, height=300, units="mm", dpi=320)

ggsave(filename=paste0("Figure4.png"),
       plot=p, device="png",
       path=in.dir,
       width=450, height=300, units="mm", dpi=320)

