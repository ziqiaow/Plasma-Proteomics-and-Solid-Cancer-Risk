library(data.table)
library(RNOmni, lib.loc = '/users/zwang4/R/4.3')
library(tibble)
#redo the proteomics data cleaning using the new PCs
load("./protein_flag2_all.RData")


pc_all <- fread("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/PCA/redo/ARIC_projection.sscore")
colnames(pc_all)[2] = "ID"
pc_all=pc_all[,-c(1,3,4,5)]
head(pc_all)

#load PEER factors 
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/ARIC_protein_peer10_v2.RData")

#load study center
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/clinical_keep_v2.RData")
clinical_keep = clinical_keep[,c("ID","v2center")]
covariates <- merge(factors,pc_all,by = "ID")
covariates <- merge(covariates,clinical_keep,by="ID")

aric.center <- covariates$v2center; aric.center <- model.matrix(~aric.center)
aric.center <- as_tibble(aric.center)
covariates <- cbind(covariates[,-ncol(covariates)], as.matrix(aric.center[,-1]))


prot_final = protein_clean_final[match(covariates$ID,protein_clean_final$SampleId),]
colnames(covariates)

identical(covariates$ID,prot_final$SampleId)

#################################################################

# inv-rank phenotype


prot_tmp = data.frame(prot_final[,-1])
prot_norm = data.frame(array(0,c(dim(prot_tmp)[1],dim(prot_tmp)[2])))

for(i in 1:dim(prot_tmp)[2]){
  
  dat <- data.frame(y=prot_tmp[,i], covariates[,-1])
  fit <- lm(y~., dat)
  prot_norm[,i] <- RankNorm(residuals(fit))
  
  print(i)
}
colnames(prot_norm) = colnames(prot_final)[-1]
rownames(prot_norm) = prot_final$SampleId


save(prot_norm,file="/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/data_clean/aric_protein_peer10_pc_redo.RData")
