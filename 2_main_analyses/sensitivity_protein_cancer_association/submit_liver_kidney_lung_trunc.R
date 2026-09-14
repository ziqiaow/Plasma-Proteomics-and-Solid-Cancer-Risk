#Cox proportional hazards model
library(survival)
library(data.table)

cancer = c("kidney","lung")
#liver only has 1 case in 10 yrs, none in 5yrs

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
cancer = cancer[task_id]

in.dir <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/sensitivity/",cancer,"/")
in.dir0 <- "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/"

e <- new.env()
obj_names <- load(paste0(in.dir0,cancer,"_discovery.RData"), envir = e)
obj_name <- obj_names[1]
cancer_incidence <- e[[obj_name]]

cancer_incidence = cancer_incidence[,-which(colnames(cancer_incidence) == "race_center")]
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
#data is the original dataset
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



  if(cancer == "lung"){
    cox.res = list()
    zp = list()
    cancer_incidence2 = cancer_incidence[-which(cancer_incidence$packyr_missing_indicator == 1),]
    for(i in 1:dim(prot_norm)[2]){
      
      dat <- data.frame(x=prot_norm[match(cancer_incidence2$ID,rownames(prot_norm)),i], cancer_incidence2)
      cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) +sum + cigt21 +packyr2 +drnkr21 + bmi21  +
                              DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center    +
                              Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
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
    cox.res[[i]] <- coxph(Surv(follow_up_5yr, incidence_v2_5yr)~ x + v2age22 +factor(gender) +  racegrp +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,"_trunc5yr.RData"))

  }
  





########################################################################################

#Truncated follow up < 10 yrs



if(cancer == "lung"){
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) +sum + cigt21 +packyr2 +packyr_missing_indicator +drnkr21 + bmi21  +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center    +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc10yr.RData"))
  
  
} else {
  cox.res = list()
  zp = list()
  
  cancer_incidence2 = cancer_incidence[-which(cancer_incidence$packyr_missing_indicator == 1),]
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[match(cancer_incidence2$ID,rownames(prot_norm)),i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(follow_up_10yr, incidence_v2_10yr)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_trunc10yr.RData"))
  
}

