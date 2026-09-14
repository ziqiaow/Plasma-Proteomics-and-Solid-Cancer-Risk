#timepoints

library(dplyr)
library(ggplot2)
library(ggpubr)
folder <- paste0("./collaboration/ARIC_Cancer/all/sensitivity/results/")
in.dir = "./collaboration/ARIC_Cancer/all/sensitivity/"

# Get list of all .RData files
cancer = c("prostate","kidney","liver","lung")
cancer_names <-  c("Prostate","Kidney","Liver","Lung")
tmp = NULL
h=0
for(i in cancer){
  h=h+1
  rdata_files <- list.files(paste0(folder,i), pattern = "_lag5yr_anno\\.RData$", full.names = TRUE)
  #rdata_files = rdata_files[c(1,3:8)]
  # Extract cancer name (after second-to-last underscore, remove trailing _anno)
 # cancer_names <-  c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")
  #sub("^.*?_.*?_.*?_.*?_(.*?)_.*$", "\\1", rdata_files)
  
  
  load(rdata_files[1])
  
  res_final <- res_final %>% mutate(
    chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
    
  ) %>%
    # Drop rows where CHR contains "chr" or is not numeric
    filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
    mutate(
      # Convert remaining CHR to numeric
      chr = as.numeric(chr),
      transcription_start_site  = as.numeric(transcription_start_site)
    )
  
  res_final$cancer = cancer_names[h]
  tmp = rbind(tmp,res_final)
  
  
  
}
res_final = tmp
save(res_final,file=paste0(in.dir,"summary_stats_combined_lag5yr.RData"))




# Get list of all .RData files
cancer = c("prostate","kidney","liver","lung")
cancer_names <-  c("Prostate","Kidney","Liver","Lung")
tmp = NULL
h=0
for(i in cancer){
  h=h+1
  rdata_files <- list.files(paste0(folder,i), pattern = "_lag10yr_anno\\.RData$", full.names = TRUE)
  #rdata_files = rdata_files[c(1,3:8)]
  # Extract cancer name (after second-to-last underscore, remove trailing _anno)
  # cancer_names <-  c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")
  #sub("^.*?_.*?_.*?_.*?_(.*?)_.*$", "\\1", rdata_files)
  
  
  load(rdata_files[1])
  
  res_final <- res_final %>% mutate(
    chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
    
  ) %>%
    # Drop rows where CHR contains "chr" or is not numeric
    filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
    mutate(
      # Convert remaining CHR to numeric
      chr = as.numeric(chr),
      transcription_start_site  = as.numeric(transcription_start_site)
    )
  
  res_final$cancer = cancer_names[h]
  tmp = rbind(tmp,res_final)
  
  
  
}
res_final = tmp
save(res_final,file=paste0(in.dir,"summary_stats_combined_lag10yr.RData"))


# Get list of all .RData files
cancer = c("prostate","kidney","lung")
cancer_names <-  c("Prostate","Kidney","Lung")
tmp = NULL
h=0
for(i in cancer){
  h=h+1
  
  rdata_files <- list.files(paste0(folder,i), pattern = "_trunc10yr_anno\\.RData$", full.names = TRUE)
  #rdata_files = rdata_files[c(1,3:8)]
  # Extract cancer name (after second-to-last underscore, remove trailing _anno)
  # cancer_names <-  c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")
  #sub("^.*?_.*?_.*?_.*?_(.*?)_.*$", "\\1", rdata_files)
  
  
  load(rdata_files[1])
  
  res_final <- res_final %>% mutate(
    chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
    
  ) %>%
    # Drop rows where CHR contains "chr" or is not numeric
    filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
    mutate(
      # Convert remaining CHR to numeric
      chr = as.numeric(chr),
      transcription_start_site  = as.numeric(transcription_start_site)
    )
  
  res_final$cancer = cancer_names[h]
  tmp = rbind(tmp,res_final)
  
  
  
}
res_final = tmp
save(res_final,file=paste0(in.dir,"summary_stats_combined_tunc10yr.RData"))



# Final plot, Add SLCO5A1 in kidney cancer
load("./collaboration/ARIC_Cancer/all/sensitivity/n_cases_lag5yr.RData")
load("./collaboration/ARIC_Cancer/all/sensitivity/n_cases_trunc10yr.RData")
load("./collaboration/ARIC_Cancer/all/sensitivity/summary_stats_combined_lag5yr.RData")
names(n_events_lag5yr) = paste0(toupper(substr(names(n_events_lag5yr), 1, 1)), substr(names(n_events_lag5yr), 2, nchar(names(n_events_lag5yr))))
names(n_events_trunc10yr) = paste0(toupper(substr(names(n_events_trunc10yr), 1, 1)), substr(names(n_events_trunc10yr), 2, nchar(names(n_events_trunc10yr))))
res_lag5yr = res_final
res_lag5yr$n_cases = sapply(res_lag5yr$cancer, function(x) n_events_lag5yr[[x]][2])
res_lag5yr$time = "> 5"
res_lag5yr$time =paste0(">= 5\n(n = ",res_lag5yr$n_cases,")")
load("./collaboration/ARIC_Cancer/all/sensitivity/summary_stats_combined_tunc10yr.RData")
res_trunc10yr = res_final
res_trunc10yr$n_cases = sapply(res_trunc10yr$cancer, function(x) n_events_trunc10yr[[x]][2])
res_trunc10yr$time = paste0("< 10\n(n = ",res_trunc10yr$n_cases,")")

load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_manhattan.RData")
res_final$time = "Full follow-up"
n_cases = c(136,	271,		96,	22,	416,	88,	588)
names(n_cases) = c("Bladder", "Colorectal", 	"Kidney", 	"Liver", 	"Lung", 	"Pancreatic", 	"Prostate") 
res_final$n_cases = n_cases[match(res_final$cancer,names(n_cases))]
res_final$time = paste0("Full follow-up\n(n = ",res_final$n_cases,")")


sig = res_final[which(res_final$fdr < 0.05),]
sig = rbind(sig,res_final[which(res_final$cancer == "Kidney" & res_final$`Entrez Gene Name` == "SLCO5A1"),])
sig$ID = paste0(sig$seqid_in_sample,"_",sig$cancer)
res_lag5yr$ID = paste0(res_lag5yr$seqid_in_sample,"_",res_lag5yr$cancer)
res_trunc10yr$ID = paste0(res_trunc10yr$seqid_in_sample,"_",res_trunc10yr$cancer)
df = rbind(sig[-which(sig$cancer == "Liver"),],res_lag5yr[match(sig$ID[-which(sig$cancer == "Liver")],res_lag5yr$ID),],res_trunc10yr[match(sig$ID[-which(sig$cancer == "Liver")],res_trunc10yr$ID),])
# Add HR and 95% CI
df <- df %>%
  mutate(
    HR = `exp(coef)`,
    lower_CI = exp(coef - 1.96 * `se(coef)`),
    upper_CI = exp(coef + 1.96 * `se(coef)`)
  )

# Faceted plot
df$protein = paste0(df$`Entrez Gene Name`," (",df$cancer,")")

#make the plot for each cancer separately
#kidney
df_plot = df[which(df$cancer == "Kidney"),]

p1 = ggplot(df_plot, aes(x = time, y = `exp(coef)`, group = protein, color = protein)) +
  # geom_line(size = 1.5) +                          # line for each protein
  # geom_point(size = 4) +                         # points at each timepoint
  geom_point(aes(shape = fdr < 0.05), size = 4, stroke = 2,position=position_dodge(width=0.5)) +  # thicker star for significant
  scale_shape_manual(values = c("FALSE" = 16,  #  non-significant
                                "TRUE"  = 8),guide="none") +  # star for significant
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.04,size = 1.5,position=position_dodge(width=0.5)) +  # CI
  geom_hline(yintercept = 1, linetype = "dashed", color = "black",size=1) +   # reference line
  # facet_wrap(~ cancer, ncol = 3)+#, #scales = "free_y") +
  scale_y_continuous(breaks =c(0,0.5,1,2,3,4,5),limits = c(0,5))  +
  scale_color_brewer(palette = "Set2") +   
  #coord_cartesian(ylim = c(0.5, 5)) +   # facet by cancer
  theme_minimal(base_size = 14) +
  labs(y = "Hazard Ratio (HR)", x = "Time to kidney cancer diagnosis (years)", color = NULL) +
  theme(legend.position = "top",
        axis.text.x = element_text(size=12,face="bold"), 
        legend.title = element_text(size = 12,face="bold"),
        legend.text = element_text(size = 12,face="bold"))


#prostate
df_plot = df[which(df$cancer == "Prostate"),]

p2 = ggplot(df_plot, aes(x = time, y = `exp(coef)`, group = protein, color = protein)) +
  # geom_line(size = 1.5) +                          # line for each protein
  geom_point(aes(shape = fdr < 0.05), size = 4, stroke = 2,position=position_dodge(width=0.5)) +  # thicker star for significant
  scale_shape_manual(values = c("FALSE" = 16,  #  non-significant
                                "TRUE"  = 8),guide="none") +  # star for significant
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.04,size = 1.5,position=position_dodge(width=0.5)) +  # CI
  geom_hline(yintercept = 1, linetype = "dashed", color = "black",size=1) +   # reference line
  # facet_wrap(~ cancer, ncol = 3)+#, #scales = "free_y") +
  scale_y_continuous(breaks =c(0.5,0,1,2,3,4,5),limits = c(0.5,5))  +
  scale_color_brewer(palette = "Set2") +   
  #coord_cartesian(ylim = c(0.5, 5)) +   # facet by cancer
  theme_minimal(base_size = 14) +
  labs(y = "Hazard Ratio (HR)", x = "Time to prostate cancer diagnosis (years)", color = NULL) +
  theme(legend.position = "top",
        axis.text.x = element_text(size=12,face="bold"), 
        legend.title = element_text(size = 12,face="bold"),
        legend.text = element_text(size = 12,face="bold"))


#lung
df_plot = df[which(df$cancer == "Lung"),]

p3 = ggplot(df_plot, aes(x = time, y = `exp(coef)`, group = protein, color = protein)) +
  #geom_line(size = 1.5) +                          # line for each protein
  geom_point(aes(shape = fdr < 0.05), size = 4, stroke = 2,position=position_dodge(width=0.5)) +  # thicker star for significant
  scale_shape_manual(values = c("FALSE" = 16,  #  non-significant
                                "TRUE"  = 8),guide="none") +  # star for significant
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.04,size = 1.5,position=position_dodge(width=0.5)) +  # CI
  geom_hline(yintercept = 1, linetype = "dashed", color = "black",size=1) +   # reference line
  # facet_wrap(~ cancer, ncol = 3)+#, #scales = "free_y") +
  scale_y_continuous(breaks =c(0.5,0,1,2,3),limits = c(0.5,2.5))  +
  scale_color_brewer(palette = "Set2") +   
  #coord_cartesian(ylim = c(0.5, 5)) +   # facet by cancer
  theme_minimal(base_size = 14) +
  labs(y = "Hazard Ratio (HR)", x = "Time to lung cancer diagnosis (years)", color = NULL) +
  theme(legend.position = "top",
        axis.text.x = element_text(size=12,face="bold"), 
        legend.title = element_text(size = 12,face="bold"),
        legend.text = element_text(size = 12,face="bold")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))  # wrap legend into 2 rows




p <- ggarrange(p1, p2,p3,
               ncol = 3, labels = c("a", "b","c"),
               widths  = c(0.33,0.33,0.33))


ggsave(filename=paste0("Figure5.pdf"),
       plot=p, device="pdf",
       path=in.dir,
       width=500, height=150, units="mm", dpi=320)



ggsave(filename=paste0("Figure5.png"),
       plot=p, device="png",
       path=in.dir,
       width=500, height=150, units="mm", dpi=320)

