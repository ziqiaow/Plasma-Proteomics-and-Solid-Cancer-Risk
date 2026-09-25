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

folder <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/raw_protein/",cancer,"/")

# Get list of all .RData files
rdata_files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
#load proteomics data
load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/protein_flag2_all.RData")
prot_norm = protein_clean_final[,-1]

# Loop through files
for (file in rdata_files) {
  # Load the .RData file
  load(file)

  #Summarize the results
  res = array(0,c(4955,5))
  for(i in 1:length(cox.res)){
    res[i,] = summary(cox.res[[i]])$coef[1,]
    
    
  }
  res = data.frame(res)
  rownames(res) = colnames(prot_norm)
  colnames(res) = colnames(summary(cox.res[[1]])$coef)
  
  
  #FDR
  res$fdr = p.adjust(res$`Pr(>|z|)`,method="BH")
  length(which(res$fdr < 0.05))
  #[1] 0
  which(res$fdr < 0.05)
  
  #bonferonni
  res$bonferonni=0
  res$bonferonni[which(res$`Pr(>|z|)` < 0.05/length(cox.res))] = 1
  
  #ENT
  res$ENT=0
  M_eff = 2865.267 #check out /ARIC_Cancer/data_clean/Effective_number_tests.R
  res$ENT[which(res$`Pr(>|z|)` <  0.05/M_eff)] = 1
  
  
  res = res[order(res$`Pr(>|z|)`),]
  head(res)
  
  #add protein annotation
  load("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/anno_official_clean.RData")
  
  res_final = cbind(res , anno[match(rownames(res),anno$id),])
  out_file <- file.path(
    folder,
    paste0(tools::file_path_sans_ext(basename(file)), "_anno.RData")
  )
  
  save(res_final, file=out_file)
  
  out_file <- file.path(
    folder,
    paste0(tools::file_path_sans_ext(basename(file)), "_anno.csv")
  )
  
  write.csv(res_final, out_file, row.names = T)
  message("Exported: ", out_file)
  
  rm(list = c("cox.res"))
}


