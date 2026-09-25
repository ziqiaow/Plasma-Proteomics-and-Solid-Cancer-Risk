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


#Truncated follow up < 5 yrs
cancer_incidence$incidence_v2_5yr <- 0                                                                                         
cancer_incidence$incidence_v2_5yr[which(cancer_incidence$incidence_after_v2==1 & cancer_incidence$incidence_t2<5)] <- 1                
cancer_incidence$follow_up_5yr <- cancer_incidence$incidence_t2                                         
cancer_incidence$follow_up_5yr[cancer_incidence$incidence_t2>=5] <- 5
table(cancer_incidence$incidence_v2_5yr)



#Truncated follow up < 10 yrs
#data is the original dataset
cancer_incidence$incidence_v2_10yr <- 0                                                                                          
cancer_incidence$incidence_v2_10yr[which(cancer_incidence$incidence_after_v2==1 & cancer_incidence$incidence_t2<10)] <- 1                
cancer_incidence$follow_up_10yr <- cancer_incidence$incidence_t2                                          
cancer_incidence$follow_up_10yr[cancer_incidence$incidence_t2>=10] <- 10
table(cancer_incidence$incidence_v2_10yr)






########################################################################################

#Truncated follow up < 5 yrs

                          
if(cancer == "pancreatic"){
  #because race-center there is no B_F here, so we create a new race-cetner for rectal cancer, combine all black into one group
  cancer_incidence$race_center_5yr = as.character(cancer_incidence$race_center)
  cancer_incidence$race_center_5yr[which(cancer_incidence$race_center_5yr == "B_F")] = "B_J"
  
cox.res = list()
zp = list()
for(i in 1:dim(prot_norm)[2]){
  
  dat <- data.frame(x=prot_norm[,i], cancer_incidence)
  cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) + race_center_5yr +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
  zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
  print(i)
 
}
names(cox.res) = names(zp) = colnames(prot_norm)
save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_trunc5yr.RData"))

} else if (cancer == "prostate"){
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x +cigt21+ packyr2 +packyr_missing_indicator  + bmi21 + v2age22 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc5yr.RData"))

  
} else if (cancer == "bladder"){
  #because race-center there is no B_F here, so we create a new race-cetner for rectal cancer, combine all black into one group
  cancer_incidence$race_center_5yr = as.character(cancer_incidence$race_center)
  cancer_incidence$race_center_5yr[which(cancer_incidence$race_center_5yr == "B_F")] = "B_J"
  cancer_incidence2 = cancer_incidence[-which(cancer_incidence$packyr_missing_indicator == 1),] #all cases have pack years, so we need to remove people who have missing data (all healthy individuals)
  table(cancer_incidence2$incidence_v2_5yr)
  
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) +cigt21+ packyr2  + bmi21 
                            + race_center_5yr +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc5yr.RData"))

} else {
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                           Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colorectal_trunc5yr.RData"))

  

        
        #colon alone
        covariates = cancer_incidence
        covariates$colon_incidence_v2 = 0
        covariates$colon_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$colon_incidence==1))] = 1
        covariates$colon_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$colon_incidence==1))] = NA #prevalent cases or missing data at V2
        
        
        #Truncated follow up < 5 yrs
        #data is the original dataset
        covariates$incidence_v2_5yr <- 0                                                                                         
        covariates$incidence_v2_5yr[which(covariates$colon_incidence_v2==1 & covariates$incidence_t2<5)] <- 1                
        covariates$follow_up_5yr <- covariates$incidence_t2                                         
        covariates$follow_up_5yr[covariates$incidence_t2>=5] <- 5
        table(covariates$incidence_v2_5yr)
        
        prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
        
        cox.res = list()
        zp = list()
        for(i in 1:dim(prot_norm)[2]){
          
          dat <- data.frame(x=prot_fit[,i], covariates)
          cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                                  DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                                  Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
          zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
          print(i)
        }
        names(cox.res) = names(zp) = colnames(prot_norm)
        save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colon_only_trunc5yr.RData"))
    
}







########################################################################################

#Truncated follow up < 10 yrs


if(cancer == "pancreatic"){
  
  
  #because race-center there is no B_F here, so we create a new race-cetner for rectal cancer, combine all black into one group
  cancer_incidence$race_center_10yr = as.character(cancer_incidence$race_center)
  cancer_incidence$race_center_10yr[which(cancer_incidence$race_center_10yr == "B_F")] = "B_J"
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) + race_center_10yr +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
    
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_trunc10yr.RData"))

} else if (cancer == "prostate"){
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x +cigt21+ packyr2 +packyr_missing_indicator  + bmi21 + v2age22 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc10yr.RData"))

  
} else if (cancer == "bladder"){
  
  cancer_incidence$race_center_10yr = as.character(cancer_incidence$race_center)
  cancer_incidence$race_center_10yr[which(cancer_incidence$race_center_10yr == "B_F")] = "B_J"
  cancer_incidence2 = cancer_incidence[-which(cancer_incidence$packyr_missing_indicator == 1),] #all cases but 1 have pack years, so we need to remove people who have missing data (all healthy individuals)
  table(cancer_incidence2$incidence_v2_10yr)
  
  prot_fit = prot_norm[match(cancer_incidence2$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) +cigt21+ packyr2  + bmi21 
                          + race_center_10yr +Factor1 + Factor2 + Factor3 + Factor4+
                            Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc10yr.RData"))

  
} else {
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colorectal_trunc10yr.RData"))

  
  
  
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
  # too many missing categories, directly use race
  
  #Truncated follow up < 10 yrs
  #data is the original dataset
  covariates$incidence_v2_10yr <- 0                                                                                         
  covariates$incidence_v2_10yr[which(covariates$rectal_incidence_v2==1 & covariates$incidence_t2<10)] <- 1                
  covariates$follow_up_10yr <- covariates$incidence_t2                                         
  covariates$follow_up_10yr[covariates$incidence_t2>=10] <- 10
  table(covariates$incidence_v2_10yr)
  
  
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + racegrp +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_rectal_trunc10yr.RData"))

  
  #colon alone
  covariates = cancer_incidence
  covariates$colon_incidence_v2 = 0
  covariates$colon_incidence_v2[which(!is.na(covariates$incidence_t2) & (covariates$colon_incidence==1))] = 1
  covariates$colon_incidence_v2[which(is.na(covariates$incidence_t2) & (covariates$prvcancr == 1 | covariates$colon_incidence==1))] = NA #prevalent cases or missing data at V2
  
  
  #Truncated follow up < 10 yrs
  #data is the original dataset
  covariates$incidence_v2_10yr <- 0                                                                                         
  covariates$incidence_v2_10yr[which(covariates$colon_incidence_v2==1 & covariates$incidence_t2<10)] <- 1                
  covariates$follow_up_10yr <- covariates$incidence_t2                                         
  covariates$follow_up_10yr[covariates$incidence_t2>=10] <- 10
  table(covariates$incidence_v2_10yr)
  
  prot_fit = prot_norm[match(covariates$ID, rownames(prot_norm)),]
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_fit[,i], covariates)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_colon_only_trunc10yr.RData"))

}


