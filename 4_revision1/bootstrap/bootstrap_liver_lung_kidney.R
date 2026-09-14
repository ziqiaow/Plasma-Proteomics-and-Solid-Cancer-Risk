#Cox proportional hazards model with Bootstrap
library(survival)
library(data.table)


cancer <-  c("kidney","liver","lung")
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)
print(task_id)
cancer = cancer[task_id]

in.dir <- paste0("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/all/original/",cancer,"/")
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



#load significant proteins
rdata_file =  list.files(in.dir, pattern = paste0("cox_fit_fullmodel_", cancer, "_anno.RData"), full.names = TRUE)
load(rdata_file)

sig_proteins <- res_final$seqid_in_sample[which(res_final$fdr< 0.05)]



########################################################
# BOOTSTRAP ANALYSIS FOR SIGNIFICANT PROTEINS
########################################################
bootstrap_cox_monitored <- function(data, prot_col, formula_vars, n_boot = 1000, 
                                    min_events = NULL, n_cores = 1) {
  
  n <- nrow(data)
  event_var <- "incidence_after_v2"
  
  # Set minimum events (e.g., 70% of original)
  if(is.null(min_events)) {
    min_events <- floor(0.7 * sum(data[[event_var]]))
  }
  
  n_events_orig <- sum(data[[event_var]])
  cat(sprintf("  Original events: %d, Minimum required: %d\n", n_events_orig, min_events))
  
  boot_iteration <- function(b) {
    set.seed(b)
    
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- data[boot_indices, ]
    n_events_boot <- sum(boot_data[[event_var]])
    
    # Skip if too few events
    if(n_events_boot < min_events) {
      return(c(coef = NA, HR = NA, 
               n_events = n_events_boot, skipped = 1))
    }
    
    tryCatch({
      fit <- coxph(formula_vars, data = boot_data)
      return(c(coef = coef(fit)[1], 
               HR = exp(coef(fit)[1]),
               n_events = n_events_boot,
               skipped = 0))
    }, error = function(e) {
      return(c(coef = NA, HR = NA,
               n_events = n_events_boot, skipped = 0))
    })
  }
  
  if(n_cores > 1) {
    boot_results <- mclapply(1:n_boot, boot_iteration, mc.cores = n_cores)
  } else {
    boot_results <- lapply(1:n_boot, boot_iteration)
  }
  
  boot_matrix <- do.call(rbind, boot_results)
  
  if(is.null(colnames(boot_matrix))) {
    colnames(boot_matrix) <- c("coef", "HR", "n_events", "skipped")
  }
  
  # Report statistics
  n_skipped <- sum(boot_matrix[, "skipped"], na.rm = TRUE)
  cat(sprintf("  Skipped %d/%d samples due to insufficient events\n", n_skipped, n_boot))
  cat(sprintf("  Event range: %d-%d (mean: %.1f)\n", 
              min(boot_matrix[, "n_events"], na.rm = TRUE),
              max(boot_matrix[, "n_events"], na.rm = TRUE),
              mean(boot_matrix[, "n_events"], na.rm = TRUE)))
  
  return(boot_matrix)
}



# Perform bootstrap for each significant protein
boot_results <- list()

for(prot in sig_proteins) {
  
  cat(sprintf("\nBootstrapping %s...\n", prot))
  
  # Prepare data
  prot_idx <- which(colnames(prot_norm) == prot)
  dat <- data.frame(x = prot_norm[, prot_idx], cancer_incidence)
  
  # Define formula based on cancer type
  if(cancer == "lung") {
    formula_boot <- Surv(incidence_t2, incidence_after_v2) ~ x + v2age22 + factor(gender) + sum + cigt21 + packyr2 + packyr_missing_indicator + drnkr21 + bmi21 +
      DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center +
      Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7 + Factor8 + Factor9 + Factor10 +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
  } else {
    formula_boot <- Surv(incidence_t2, incidence_after_v2) ~ x + v2age22 + factor(gender) + cigt21 + packyr2 + packyr_missing_indicator + anta01 + drnkr21 + bmi21 + wsthpr21 +
      DIA_DBV2 + UNDIA_DBV2 + ATR_DBV2 + race_center + work_i02 + sprt_i02 + lisr_i01 +
      Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7 + Factor8 + Factor9 + Factor10 +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
  }
  
  # Run bootstrap 
  boot_res <- bootstrap_cox_monitored(dat, prot, formula_boot, n_boot = 1000,min_events = 15)
  # Store results
  boot_results[[prot]] <- boot_res
}

# Save bootstrap results
save(boot_results, sig_proteins, file = paste0(in.dir, "bootstrap_results_", cancer, ".RData"))
