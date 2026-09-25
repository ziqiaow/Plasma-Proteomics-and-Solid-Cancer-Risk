#Cox proportional hazards model
library(survival)
library(data.table)

cancer = c("prostate","bladder","colon","pancreatic","liver","kidney","lung")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
cancer = cancer[task_id]

in.dir <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/original/",cancer,"/")
in.dir0 <- "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/"



if(cancer %in% c("prostate","bladder","colon","pancreatic")){
e <- new.env()
obj_names <- load(paste0(in.dir0,cancer,"_incident_v2.RData"), envir = e)
obj_name <- obj_names[1]
cancer_incidence <- e[[obj_name]]

if(cancer == "colon"){
  colnames(cancer_incidence)[1] = "ID"
}
#load covariates
covariates = readRDS("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/covariates_clean_final.rds")
cancer_incidence = merge(cancer_incidence, covariates, by="ID",suffixes = c("", "_dup"))

} else {

  

  e <- new.env()
  obj_names <- load(paste0(in.dir0,cancer,"_discovery.RData"), envir = e)
  obj_name <- obj_names[1]
  cancer_incidence <- e[[obj_name]]
  
  cancer_incidence = cancer_incidence[,-which(colnames(cancer_incidence) == "race_center")]
  #load covariates
  covariates = readRDS("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/covariates_clean_final.rds")
  cancer_incidence = merge(cancer_incidence, covariates, by="ID",suffixes = c("", "_dup"))
  
  
}
#make alcohol and smoke as factors
cancer_incidence$cigt21 = factor(cancer_incidence$cigt21,levels = c(3,2,1))
#cigt21: 1.current ; 2.former ; 3. never(reference group)
cancer_incidence$drnkr21 = factor(cancer_incidence$drnkr21,levels = c(3,2,1))

#only use black individuals
cancer_incidence = cancer_incidence[which(cancer_incidence$racegrp == "B"),]


#load proteomics data
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/data_clean/aric_protein_peer10_pc_redo.RData")

prot_norm = prot_norm[match(cancer_incidence$ID,rownames(prot_norm)),]
identical(rownames(prot_norm),cancer_incidence$ID)

########################################################
########################################################
########################################################
#Fit cox proportional hazards model


#minimally adjusted model
if(cancer =="prostate"){
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x  + v2age22 +bmi21+cigt21 +packyr2 +packyr_missing_indicator+
                            as.factor(center) +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_smk_black.RData"))
  
  
} else if (cancer == "colon"){
  
  
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x  +factor(gender) +bmi21+ cigt21 +packyr2 +packyr_missing_indicator+ v2age22 +
                            as.factor(center) +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_colorectal_smk_black.RData"))
  
  
  #rectal only
  #check NAs
  #id = which(complete.cases(covariates) == T)
  #rectal alone
  covariates = cancer_incidence
  covariates$rectal_incidence_v2 = 0
  covariates$rectal_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$rectal_incidence==1))] = 1
  covariates$rectal_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$rectal_incidence==1))] = NA #prevalent cases or missing data at V2
  # #because race-center there is no B_F here, so we create a new race-cetner for rectal cancer, combine all black into one group
  # covariates$race_center_rectal = as.character(covariates$race_center)
  # covariates$race_center_rectal[which(covariates$race_center_rectal == "B_F")] = "B_J"
  # 
  # 
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(incidence_t2,rectal_incidence_v2)~ x  +factor(gender)+bmi21+ cigt21 +packyr2 +packyr_missing_indicator + v2age22 +
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_rectal_smk_black.RData"))
  
  
  #colon alone
  covariates = cancer_incidence
  covariates$colon_incidence_v2 = 0
  covariates$colon_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$colon_incidence==1))] = 1
  covariates$colon_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$colon_incidence==1))] = NA #prevalent cases or missing data at V2
  
  
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(incidence_t2,colon_incidence_v2)~ x +factor(gender)+bmi21+ cigt21 +packyr2 +packyr_missing_indicator + v2age22 +
                            as.factor(center) +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_colon_smk_black.RData"))
  
} else if (cancer == "bladder"){#for black, there is only 1 female with bladder cancer among black, rest 10 incidence were males, remove gender from the model
  
  cox.res = list()
  #zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x  +bmi21+ cigt21 +packyr2 +packyr_missing_indicator+ v2age22 +
                            as.factor(center) +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
   # zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res)  = colnames(prot_norm)
  save(cox.res, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_smk_black.RData"))
  
}  else {
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x  +factor(gender) +bmi21+ cigt21 +packyr2 +packyr_missing_indicator+ v2age22 +
                            as.factor(center) +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_smk_black.RData"))
  
}


