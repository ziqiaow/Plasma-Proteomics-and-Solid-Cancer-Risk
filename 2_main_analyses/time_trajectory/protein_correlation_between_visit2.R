#Use somalogic ARIC proteomics data for v2, v3, v5
#Location: /dcs04/legacy-dcs01-arking/ARIC_static/ARIC_Data/Proteomics/ANML_Oct2021

#load the all cancer file
library(readstata13)
cancer <- read.dta13("./overall_cancer_2015_np.dta")
cancer = cancer[which(cancer$prvcancr == 0),] #14743 
table(cancer$anycancer)
# 0     1 
# 10130  4613 
id = cancer$id[which(cancer$anycancer == 0)]

#load the raw ANML normalized v2,v3,v5 protein data
library(data.table)
v5 = fread("./protein/Visit5/soma_visit_5_log2_ANML_SMP_updated.txt")
v5 = v5[,-c(1,3:33)]
v5 = v5[,-c(5286:5297)]

v3 = fread("./protein/Visit3/soma_visit_3_log2_ANML_SMP_updated.txt")
v3 = v3[,-c(1,2,4:34)]
v3 = v3[,-c(5286:5298)]

#for v2
#flag2=0 is applied, 4955 aptamers are used
load("./protein_flag2_all.RData")
v2 = protein_clean_final
id_intersect = intersect(v2$SampleId,id) #7517
id_v3 = intersect(v3$SampleId,id) #7336
id_v5 = intersect(v5$SampleId,id) #3710
id_v3v5 = intersect(id_v3, id_v5) #3324
id_all = intersect(id_v3v5,id_intersect) #3014

v2 = v2[match(id_all,v2$SampleId),]
v2= data.frame(v2)
col_id = match(colnames(v2),colnames(v3))
v3 = data.frame(v3)[match(id_all,v3$SampleId),col_id]
col_id = match(colnames(v2),colnames(v5))
v5 = data.frame(v5)[match(id_all,v5$SampleId),col_id]

#check pairwise correlations between v2, v3, and v5 for the top proteins
protein_cols <- setdiff(colnames(v2), "SampleId")

results <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    test <- cor.test(v2[[p]], v3[[p]])
    c(Correlation = test$estimate, 
      P_value = test$p.value,
      CI_lower = test$conf.int[1],
      CI_upper = test$conf.int[2])
  }))
)

print(results)


results_v2v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    test <- cor.test(v2[[p]], v5[[p]])
    c(Correlation = test$estimate, 
      P_value = test$p.value,
      CI_lower = test$conf.int[1],
      CI_upper = test$conf.int[2])
  }))
)
head(results_v2v5)


results_v3v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    test <- cor.test(v3[[p]], v5[[p]])
    c(Correlation = test$estimate, 
      P_value = test$p.value,
      CI_lower = test$conf.int[1],
      CI_upper = test$conf.int[2])
  }))
)
head(results_v3v5)

# Merge all results by Protein
combined_wide <- results %>%
  rename_with(~paste0(., "_V2V3"), -Protein) %>%
  full_join(
    results_v2v5 %>% rename_with(~paste0(., "_V2V5"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v3v5 %>% rename_with(~paste0(., "_V3V5"), -Protein),
    by = "Protein"
  )

print(combined_wide)


#annotation
load("./protein/anno_official_clean.RData")
combined_wide_anno = cbind(combined_wide,anno[match(combined_wide$Protein,anno$id),])
save(combined_wide_anno,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_cor_v2v3v5.RData")
write.csv(combined_wide_anno,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_cor_v2v3v5.csv")
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
id = res_final$seqid_in_sample[which(res_final$fdr < 0.05)]
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all_minimal_smk.RData")
id2= res_final$seqid_in_sample[which(res_final$fdr < 0.05)]
id_all = unique(c(id,id2))
combined_wide_anno_top = combined_wide_anno[match(id_all,combined_wide_anno$Protein),]
write.csv(combined_wide_anno_top,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_cor_v2v3v5_top.csv")



#summary statistics for 40 proteins
summary = correlation_matrix_data_top
summary(summary$Correlation[which(summary$Comparison == "V2 vs V3")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2389  0.5048  0.6297  0.5902  0.6874  0.7711 
summary(summary$Correlation[which(summary$Comparison == "V2 vs V5")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1569  0.4483  0.5499  0.5288  0.6335  0.7554 
sd(summary$Correlation[which(summary$Comparison == "V2 vs V3")])
#[1]  0.1296033
sd(summary$Correlation[which(summary$Comparison == "V2 vs V5")])
#[1] 0.1289984

#check correlations between two aptamers that map to IGF1 
#SeqId_2952_75 and SeqId_8406_17
cor.test(v2[,which(colnames(v2) == "SeqId_8406_17")],v2[,which(colnames(v2) == "SeqId_2952_75")])
#       cor 
#-0.06317022 
