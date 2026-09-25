#Cox proportional hazards model
library(survival)
library(data.table)

cancer = c("prostate","bladder","colon","pancreatic")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
cancer = cancer[task_id]

in.dir <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/sensitivity/",cancer,"/")
in.dir0 <- "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/"

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


#make alcohol and smoke as factors
cancer_incidence$cigt21 = factor(cancer_incidence$cigt21,levels = c(3,2,1))
#cigt21: 1.current ; 2.former ; 3. never(reference group)
cancer_incidence$drnkr21 = factor(cancer_incidence$drnkr21,levels = c(3,2,1))



#load proteomics data
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/data_clean/aric_protein_peer10_pc_redo.RData")

prot_norm = prot_norm[match(cancer_incidence$ID,rownames(prot_norm)),]
identical(rownames(prot_norm),cancer_incidence$ID)

########################################################
########################################################
########################################################
#Fit cox proportional hazards model



#Subset of follow up >= 5 yrs
cancer_incidence2 <- subset(cancer_incidence, incidence_t2>=5)
cancer_incidence2$incidence_v2_after5yr <- 0
cancer_incidence2$incidence_v2_after5yr[cancer_incidence2$incidence_after_v2==1] <- 1
cancer_incidence2$follow_up_after5yr <- cancer_incidence2$incidence_t2
table(cancer_incidence2$incidence_v2_after5yr)
summary(cancer_incidence2$follow_up_after5yr)




########################################################################################

#subset follow up > 5 yrs


if(cancer == "pancreatic"){
  
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
   
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator  + bmi21  +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + drnkr21  + race_center +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
    
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_lag5yr.RData"))

  
} else if (cancer == "prostate"){
  
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x +cigt21+ packyr2 +packyr_missing_indicator  + bmi21 + v2age22 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_lag5yr.RData"))

  
  
} else if (cancer == "bladder"){
  
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x + v2age22 +factor(gender) +cigt21+ packyr2 +packyr_missing_indicator  + bmi21 
                          + race_center +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_lag5yr.RData"))

  
} else {
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colorectal_lag5yr.RData"))

  
  
  
  #rectal only
  #check NAs
  #id = which(complete.cases(covariates) == T)
  #rectal alone
  covariates = cancer_incidence
  covariates$rectal_incidence_v2 = 0
  covariates$rectal_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$rectal_incidence==1))] = 1
  covariates$rectal_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$rectal_incidence==1))] = NA #prevalent cases or missing data at V2
  
  #because race-center there is no B_F here, so we create a new race-cetner for rectal cancer, combine all black into one group
  covariates$race_center_rectal = as.character(covariates$race_center)
  covariates$race_center_rectal[which(covariates$race_center_rectal == "B_F")] = "B_J"
  
  #Subset of follow up >= 5 yrs
  covariates <- subset(covariates, incidence_t2>=5)
  covariates$incidence_v2_after5yr <- 0
  covariates$incidence_v2_after5yr[covariates$rectal_incidence_v2==1] <- 1
  covariates$follow_up_after5yr <- covariates$incidence_t2
  table(covariates$incidence_v2_after5yr)
  summary(covariates$follow_up_after5yr)
  
  
  
  
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center_rectal +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_rectal_lag5yr.RData"))
  table(dat$rectal_incidence_v2)
  
  
  #colon alone
  covariates = cancer_incidence
  covariates$colon_incidence_v2 = 0
  covariates$colon_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$colon_incidence==1))] = 1
  covariates$colon_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$colon_incidence==1))] = NA #prevalent cases or missing data at V2
  
  #Subset of follow up >= 5 yrs
  covariates <- subset(covariates, incidence_t2>=5)
  covariates$incidence_v2_after5yr <- 0
  covariates$incidence_v2_after5yr[covariates$colon_incidence_v2==1] <- 1
  covariates$follow_up_after5yr <- covariates$incidence_t2
  table(covariates$incidence_v2_after5yr)
  summary(covariates$follow_up_after5yr)
  
  
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(follow_up_after5yr, incidence_v2_after5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colon_only_lag5yr.RData"))
  table(dat$colon_incidence_v2)
  
}




