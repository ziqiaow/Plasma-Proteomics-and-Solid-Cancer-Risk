library(irr)
#load the all cancer file
library(readstata13)
cancer <- read.dta13("./data/aric/overall_cancer_2015_np.dta")
cancer = cancer[which(cancer$prvcancr == 0),] #14743 
table(cancer$anycancer)
# 0     1 
# 10130  4613 
id = cancer$id[which(cancer$anycancer == 0)]

#load the raw ANML normalized v2,v3,v5 protein data
library(data.table)
v5 = fread("./data/aric/protein/Visit5/soma_visit_5_log2_ANML_SMP_updated.txt")
v5 = v5[,-c(1,3:33)]
v5 = v5[,-c(5286:5297)]

v3 = fread("./data/aric/protein/Visit3/soma_visit_3_log2_ANML_SMP_updated.txt")
v3 = v3[,-c(1,2,4:34)]
v3 = v3[,-c(5286:5298)]

#for v2
#flag2=0 is applied, 4955 proteins are used
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
#protein_cols <- setdiff(colnames(v2), "SampleId")
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
protein_cols = unique(res_final$seqid_in_sample[which(res_final$fdr < 0.05)])

results <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v3[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

print(results)



results_v2v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v3v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v3[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v2v3 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v3[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

# Merge all results by Protein
combined_wide <- results %>%
  rename_with(~paste0(., "_V2V3v5"), -Protein) %>%
  full_join(
    results_v2v5 %>% rename_with(~paste0(., "_V2V5"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v2v3 %>% rename_with(~paste0(., "_V2V3"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v3v5 %>% rename_with(~paste0(., "_V3V5"), -Protein),
    by = "Protein"
  )

print(combined_wide)


#edit and change ACP3 and KLK3 to be males only
replace_protein = c("SeqId_8468_19","SeqId_19317_114")
sex_info = cancer[match(id_all,cancer$id),]
# Separate by sex
male_ids <- v2$SampleId[sex_info$v1gender == "M"]  # Adjust column name if different
female_ids <- v2$SampleId[sex_info$v1gender == "F"]

v2_male <- v2[v2$SampleId %in% male_ids, ]
v3_male <- v3[v3$SampleId %in% male_ids, ]
v5_male <- v5[v5$SampleId %in% male_ids, ]

results <- data.frame(
  Protein = replace_protein,
  t(sapply(replace_protein, function(p) {
    tmp = cbind(v2_male[[p]],v3_male[[p]],v5_male[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)


results_v2v5 <- data.frame(
  Protein = replace_protein,
  t(sapply(replace_protein, function(p) {
    tmp = cbind(v2_male[[p]],v5_male[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v3v5 <- data.frame(
  Protein = replace_protein,
  t(sapply(replace_protein, function(p) {
    tmp = cbind(v3_male[[p]],v5_male[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v2v3 <- data.frame(
  Protein = replace_protein,
  t(sapply(replace_protein, function(p) {
    tmp = cbind(v2_male[[p]],v3_male[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

# Merge all results by Protein
combined_wide_male <- results %>%
  rename_with(~paste0(., "_V2V3v5"), -Protein) %>%
  full_join(
    results_v2v5 %>% rename_with(~paste0(., "_V2V5"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v2v3 %>% rename_with(~paste0(., "_V2V3"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v3v5 %>% rename_with(~paste0(., "_V3V5"), -Protein),
    by = "Protein"
  )

#replace these 2 proteins in the overall file
combined_wide[match(combined_wide_male$Protein,combined_wide$Protein),] = combined_wide_male
save(combined_wide,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_ICC_top.RData")
write.csv(combined_wide,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_ICC_top.csv")



#now look at all proteins in model 1 and 2
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
id = res_final$seqid_in_sample[which(res_final$fdr < 0.05)]
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all_minimal_smk.RData")
id2= res_final$seqid_in_sample[which(res_final$fdr < 0.05)]
protein_cols = unique(c(id,id2))


results <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v3[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

print(results)



results_v2v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v3v5 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v3[[p]],v5[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

results_v2v3 <- data.frame(
  Protein = protein_cols,
  t(sapply(protein_cols, function(p) {
    tmp = cbind(v2[[p]],v3[[p]])
    test <- icc(tmp,"twoway","consistency","single")
    c(ICC = test$value, 
      P_value = test$p.value,
      CI_lower = test$lbound,
      CI_upper = test$ubound)
  }))
)

# Merge all results by Protein
combined_wide <- results %>%
  rename_with(~paste0(., "_V2V3v5"), -Protein) %>%
  full_join(
    results_v2v5 %>% rename_with(~paste0(., "_V2V5"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v2v3 %>% rename_with(~paste0(., "_V2V3"), -Protein),
    by = "Protein"
  ) %>%
  full_join(
    results_v3v5 %>% rename_with(~paste0(., "_V3V5"), -Protein),
    by = "Protein"
  )
#replace these 2 proteins in the overall file
combined_wide[match(combined_wide_male$Protein,combined_wide$Protein),] = combined_wide_male
save(combined_wide,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_ICC_top_bothmodel.RData")
write.csv(combined_wide,file="./collaboration/ARIC_Cancer/manuscript/edit/protein_ICC_top_bothmodel.csv")


