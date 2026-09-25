#Figure 3

# Prepare data for plots
load("./collaboration/ARIC_Cancer/all/original/pathway/pathway_liver.RData")
go_df <- as.data.frame(go_enrich)
go_enrich_bp_df <- as.data.frame(go_enrich_bp)
go_enrich_mf_df <- as.data.frame(go_enrich_mf)
reactome_df <-as.data.frame(reactome)

reactome_df$pathway = "Reactome"
go_enrich_bp_df$pathway = "GO Biological Process"
go_enrich_mf_df$pathway = "GO Molecular Function"
go_df$pathway = "GO Cellular Component"

full_results_print = rbind(reactome_df,go_enrich_bp_df,go_enrich_mf_df,go_df)
full_results_print$GeneRatio = as.character(full_results_print$GeneRatio)
full_results_print$GeneRatio <- ifelse(grepl("^\\d+/\\d+$", full_results_print$GeneRatio),
                                       paste0("'", full_results_print$GeneRatio),
                                       full_results_print$GeneRatio)


library(ggplot2)
library(forcats)
library(enrichplot)

library(stringr)

reactome@result$Description <- str_wrap(reactome@result$Description, width = 20)
edo <- pairwise_termsim(reactome)
#p1 <- emapplot(edo)
#p2 <- emapplot(edo, cex_category=1.5)
p1 <- emapplot(edo, layout="nicely",  node_label_size  = 1.7,
               size_category  = 1.4, size_edge = 0.7,color_edge = "grey60",
               # Instead of cex_label_category
               showCategory = 30)+  scale_color_viridis_c(
                 option = "plasma",
                 limits = c(1e-5, 0.05),
                 trans = "log10",
                 breaks = c(1e-5, 1e-3, 0.05),
                 labels = c("1e-5", "1e-3", "0.05"),
                 guide = "none"
               ) +
  theme(legend.position = "none")


# Wrap long descriptions to multiple lines
go_enrich_bp_df$Description <- str_wrap(go_enrich_bp_df$Description, width = 30)
reactome_df$Description <- str_wrap(reactome_df$Description, width = 30)
go_enrich_mf_df$Description <- str_wrap(go_enrich_mf_df$Description, width = 30)
go_df$Description <- str_wrap(go_df$Description, width = 30)

top_n <- 5
reactome_df <- reactome_df %>% arrange(p.adjust) %>% head(top_n)
p2 <- ggplot(reactome_df, aes(x = -log10(p.adjust), 
                              y = fct_reorder(Description, p.adjust, .desc = T))) + 
  geom_segment(aes(xend = 0, yend = Description), color = "gray80") +
  geom_point(aes(color = p.adjust), size = 6) +  # Larger points
  #scale_color_viridis_c(option = "plasma", guide = guide_colorbar(reverse = T)) +
  scale_color_viridis_c(
    option = "plasma",
    limits = c(1e-5, 0.05),
    trans = "log",
    breaks = c(1e-5, 1e-3, 0.05),  # Custom legend ticks
    labels = c("1e-5",  "1e-3", "0.05"),
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme_minimal(base_size = 12) +  # Increase base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  xlab(expression(-log[10]("Adjusted P value"))) +
  ylab(NULL) +ggtitle(NULL) 





top_n <- 5
go_enrich_bp_df <- go_enrich_bp_df %>% arrange(p.adjust) %>% head(top_n)
p3 <- ggplot(go_enrich_bp_df, aes(x = -log10(p.adjust), 
                                  y = fct_reorder(Description, p.adjust, .desc = T))) + 
  geom_segment(aes(xend = 0, yend = Description), color = "gray80") +
  geom_point(aes(color = p.adjust), size = 9) +  # Larger points
  # scale_color_viridis_c(option = "plasma", guide = guide_colorbar(reverse = T)) +
  scale_color_viridis_c(
    option = "plasma",
    limits = c(1e-5, 0.05),
    trans = "log",
    breaks = c(1e-5, 1e-3, 0.05),  # Custom legend ticks
    labels = c("1e-5",  "1e-3", "0.05"),
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  xlab(expression(-log[10]("Adjusted P value"))) +
  ylab(NULL) +
  ggtitle("Top 5 GO BP Pathways")



top_n <- 5
go_enrich_mf_df <- go_enrich_mf_df %>% arrange(p.adjust) %>% head(top_n)
p4 <- ggplot(go_enrich_mf_df, aes(x = -log10(p.adjust), 
                                  y = fct_reorder(Description, p.adjust, .desc = T))) + 
  geom_segment(aes(xend = 0, yend = Description), color = "gray80") +
  geom_point(aes(color = p.adjust), size = 9) +  # Larger points
  # scale_color_viridis_c(option = "plasma", guide = guide_colorbar(reverse = T)) +
  scale_color_viridis_c(
    option = "plasma",
    limits = c(1e-5, 0.05),
    trans = "log",
    breaks = c(1e-5, 1e-3, 0.05),  # Custom legend ticks
    labels = c("1e-5",  "1e-3", "0.05"),
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  xlab(expression(-log[10]("Adjusted P value"))) +
  ylab(NULL) +
  ggtitle("Top 5 GO MF Pathways")


go_df <- go_df %>% arrange(p.adjust) %>% head(top_n)
p5 <- ggplot(go_df, aes(x = -log10(p.adjust), 
                        y = fct_reorder(Description, p.adjust, .desc = T))) + 
  geom_segment(aes(xend = 0, yend = Description), color = "gray80") +
  geom_point(aes(color = p.adjust), size = 9) +  # Larger points
  # scale_color_viridis_c(option = "plasma", guide = guide_colorbar(reverse = T)) +
  scale_color_viridis_c(
    option = "plasma",
    limits = c(1e-5, 0.05),
    trans = "log",
    breaks = c(1e-5, 1e-3, 0.05),  # Custom legend ticks
    labels = c("1e-5",  "1e-3", "0.05"),
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  xlab(expression(-log[10]("Adjusted P value"))) +
  ylab(NULL) +
  ggtitle("Top 5 GO CC Pathways")


#put volcano plot together
#make the volcano plot
library(EnhancedVolcano)
library(latex2exp)
library(ggrepel)
in.dir = "./collaboration/ARIC_Cancer/all/original/"

#only plot liver cancer
load("./collaboration/ARIC_Cancer/all/original/results/cox_fit_fullmodel_liver_anno.RData")
library(ggrepel)
# plot adding up all layers we have seen so far
colnames(res_final)[5] = "P"
res_final$log2HR = log2(res_final$`exp(coef)`)
library(ggrepel)
# plot adding up all layers we have seen so far
res_final$diffexpressed = "NO"
res_final$diffexpressed[which(res_final$fdr < 0.05 & res_final$coef < 0 )] = "Down"
res_final$diffexpressed[which(res_final$fdr < 0.05 & res_final$coef > 0 )] = "Up"

res_final$delabel <- NA
res_final$delabel[res_final$diffexpressed != "NO"] <- res_final$`Entrez Gene Name`[res_final$diffexpressed != "NO"]
#p_thre = res_final[which.max(res_final$fdr[which(res_final$fdr < 0.05)]),"P"]#0.001793778
p_thre = 0.05
p6 <- ggplot(data=res_final, aes(x=log2HR, y=-log10(fdr), col=diffexpressed, label=delabel),max.overlaps =Inf) +
  #scale_y_continuous(breaks =c(0,0.5,1,1.5,2,40,60,80))  +
  scale_y_continuous(breaks =c(0,0.5,1,1.5,2,2.5,3),limits = c(0,3))  +
  scale_x_continuous(breaks =c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),limits = c(-2,2))  +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(p_thre), col="black",linetype="dashed") + theme(legend.position="none") +xlab("Log2 Hazard Ratio, Model 2 (adjusted for common and site-specific risk factors)") +ylab(TeX("$-log_{10}$ Adjusted P Value"))


#minimal model
load("./collaboration/ARIC_Cancer/all/original/results/cox_fit_minimal_liver_smk_anno.RData")
# plot adding up all layers we have seen so far
colnames(res_final)[5] = "P"
res_final$log2HR = log2(res_final$`exp(coef)`)
# plot adding up all layers we have seen so far
res_final$diffexpressed = "NO"
res_final$diffexpressed[which(res_final$fdr < 0.05 & res_final$coef < 0 )] = "Down"
res_final$diffexpressed[which(res_final$fdr < 0.05 & res_final$coef > 0 )] = "Up"

res_final$delabel <- NA
res_final$delabel[res_final$diffexpressed != "NO"] <- res_final$`Entrez Gene Name`[res_final$diffexpressed != "NO"]
#p_thre = res_final[which.max(res_final$fdr[which(res_final$fdr < 0.05)]),"P"]#0.001793778
p_thre = 0.05
p7 <- ggplot(data=res_final, aes(x=log2HR, y=-log10(fdr), col=diffexpressed, label=delabel),max.overlaps =Inf) +
  #scale_y_continuous(breaks =c(0,0.5,1,1.5,2,40,60,80))  +
  scale_y_continuous(breaks =c(0,0.5,1,1.5,2,2.5,3),limits = c(0,3))  +
  scale_x_continuous(breaks =c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),limits = c(-2,2))  +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(p_thre), col="black",linetype="dashed") + theme(legend.position="none") +xlab("Log2 Hazard Ratio, Model 1 (adjusted for common risk factors only)") +ylab(TeX("$-log_{10}$ Adjusted P Value"))
library(ggpubr)


p <- ggarrange(ggarrange(p6, p7,
                         ncol = 2, labels = c("a", "b"),
                         widths  = c(0.5,0.5)),
               ggarrange(p2,p1,ncol=2,labels = c("c","d"),widths = c(0.5,0.5)),
               # ggarrange(p2,p3,p4,p5,
               #      ncol = 4, labels = c("c", "d","e","f","g"),
               #      widths  = rep(0.2,5)),
               nrow = 2, labels = c(NULL,NULL),heights = c(0.6,0.4))

ggsave(filename=paste0("Figure3_liver_02232026.pdf"),
       plot=p, device="pdf",
       path=in.dir,
       width=350, height=250, units="mm", dpi=320)


ggsave(filename=paste0("Figure3_liver_02232026.png"),
       plot=p, device="png",
       path=in.dir,
       width=350, height=250, units="mm", dpi=320)



