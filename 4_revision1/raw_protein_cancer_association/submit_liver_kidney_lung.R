#Cox proportional hazards model
library(survival)
library(data.table)

cancer = c("liver","kidney","lung")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
cancer = cancer[task_id]

in.dir <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/raw_protein/",cancer,"/")
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
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/protein_flag2_all.RData")
prot_norm = data.frame(protein_clean_final[match(as.character(cancer_incidence$ID),as.character(protein_clean_final$SampleId)),])
rownames(prot_norm) = prot_norm$SampleId
identical(rownames(prot_norm),as.character(cancer_incidence$ID))
prot_norm = prot_norm[,-1]

########################################################
########################################################
########################################################
#Fit cox proportional hazards model

  if(cancer == "lung"){
    cox.res = list()
    zp = list()
    for(i in 1:dim(prot_norm)[2]){
      
      dat <- data.frame(x=prot_norm[,i], cancer_incidence)
      cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x + v2age22 +factor(gender) +sum + cigt21 +packyr2 +packyr_missing_indicator +drnkr21 + bmi21  +
                              DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center    +
                              Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
      zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
      print(i)
    }
    names(cox.res) = names(zp) = colnames(prot_norm)
    save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,".RData"))
    table(dat$incidence_after_v2)
    
    
  } else {
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                           Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,".RData"))
  table(dat$incidence_after_v2)
  }
  
  
#minimally adjusted model

  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_norm)[2]){
    
    dat <- data.frame(x=prot_norm[,i], cancer_incidence)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x  +factor(gender) + v2age22 +bmi21+ cigt21 +packyr2 +packyr_missing_indicator+
                            race_center +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_norm)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_minimal_",cancer,".RData"))
  


