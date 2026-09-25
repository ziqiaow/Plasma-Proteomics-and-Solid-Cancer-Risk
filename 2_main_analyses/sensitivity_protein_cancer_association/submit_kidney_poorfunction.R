#Cox proportional hazards model
library(survival)
library(data.table)

cancer = c("kidney")


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
covariates_sensitivity = readRDS("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/covariates_clean_final_sensitivity.rds")
cancer_incidence = merge(cancer_incidence, covariates_sensitivity , by.x = "ID",by.y = "subjectid")

#make alcohol and smoke as factors
cancer_incidence$cigt21 = factor(cancer_incidence$cigt21,levels = c(3,2,1))
#cigt21: 1.current ; 2.former ; 3. never(reference group)
cancer_incidence$drnkr21 = factor(cancer_incidence$drnkr21,levels = c(3,2,1))



#remove individuals with poor kidney function, meansured by eGFR2, use 60 as a threshold
cancer_incidence2 <- subset(cancer_incidence, value_egfr2 >= 60)

#load proteomics data
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/data_clean/aric_protein_peer10_pc_redo.RData")

prot_fit = prot_norm[match(cancer_incidence2$ID,rownames(prot_norm)),]
identical(rownames(prot_fit),cancer_incidence2$ID)

########################################################
########################################################
########################################################
#Fit cox proportional hazards model


  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_fit)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_poorfunction_60.RData"))
  table(dat$incidence_after_v2)
  cox.res[[1]]

  #N = 8801, cases = 91
  
  
  
  
#remove individuals with poor kidney function, meansured by eGFR2, use 90 as a threshold
 cancer_incidence2 <- subset(cancer_incidence, value_egfr2 >= 90)
  

  prot_fit = prot_norm[match(cancer_incidence2$ID,rownames(prot_norm)),]
  identical(rownames(prot_fit),cancer_incidence2$ID)
  
  ########################################################
  ########################################################
  ########################################################
  #Fit cox proportional hazards model
  
  
  cox.res = list()
  zp = list()
  for(i in 1:dim(prot_fit)[2]){
    
    dat <- data.frame(x=prot_fit[,i], cancer_incidence2)
    cox.res[[i]] <- coxph(Surv(incidence_t2, incidence_after_v2)~ x + v2age22 +factor(gender) + cigt21 +packyr2 +packyr_missing_indicator+anta01 +drnkr21 + bmi21 + wsthpr21 +
                            DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +work_i02 + sprt_i02 +lisr_i01   +
                            Factor1 + Factor2 + Factor3 + Factor4+ Factor5+Factor6+Factor7+Factor8+Factor9+Factor10+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, dat)
    zp[[i]] <- cox.zph(cox.res[[i]]) #For model diagnostics
    print(i)
  }
  names(cox.res) = names(zp) = colnames(prot_fit)
  save(cox.res,zp, file = paste0(in.dir,"cox_fit_fullmodel_",cancer,"_poorfunction_90.RData"))
  table(dat$incidence_after_v2)
  cox.res[[1]]
  
  #N = 6632, cases = 73
  
  
  