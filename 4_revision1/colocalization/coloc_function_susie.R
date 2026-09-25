
# ============================================================================
# COLOCALIZATION: PLINK2 .glm pQTL + .tsv.gz GWAS
# ============================================================================

#install.packages("susieR", lib = "/dcs04/nilanjan/data/zwang/tools")
#run this line below in linux
#R_LIBS_USER="/dcs04/nilanjan/data/zwang/tools/susieR" R CMD INSTALL -l "/dcs04/nilanjan/data/zwang/tools" /dcs04/nilanjan/data/zwang/tools/coloc_5.2.3.tar.gz
#install.packages("coloc", lib = "/dcs04/nilanjan/data/zwang/tools")
#library(susieR, lib.loc = "/dcs04/nilanjan/data/zwang/tools")
library(coloc, lib.loc = "/dcs04/nilanjan/data/zwang/tools")
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================
# STEP 1: Load PLINK2 .glm files (pQTL data)
# ============================================================================

load_pqtl_glm <- function(glm_file) {
  
  cat(paste0("Loading pQTL: ", basename(glm_file), "\n"))
  
  pqtl <- fread(glm_file)
  
  # Standardize column names
  if("#CHROM" %in% names(pqtl)) {
    pqtl <- pqtl %>% rename(chr = `#CHROM`)
  }
  
  pqtl_clean <- pqtl %>%
    rename(
      pos = POS,
      rsid = ID,
      ref = REF,
      alt = ALT,
      effect_allele = A1,  # A1 is the tested allele
      beta = BETA,
      se = SE,
      pval = P
    ) %>%
    mutate(
      chr = as.integer(gsub("chr", "", chr)),
      # other_allele is whichever one is NOT the effect allele
      other_allele = ifelse(effect_allele == ref, alt, ref)
    )
  
  # Filter for additive test if present
  if("TEST" %in% names(pqtl_clean)) {
    pqtl_clean <- pqtl_clean %>% filter(TEST == "ADD")
  }
  
  # Keep allele frequency
  if("A1_FREQ" %in% names(pqtl_clean)) {
    pqtl_clean <- pqtl_clean %>% rename(freq = A1_FREQ)
  } else {
    pqtl_clean$freq <- NA
  }
  
  pqtl_clean <- pqtl_clean %>%
    select(rsid, chr, pos, effect_allele, other_allele, beta, se, pval, freq) %>%
    filter(!is.na(beta), !is.na(se), !is.na(pval))
  
  cat(paste0("  Loaded ", nrow(pqtl_clean), " variants\n"))
  
  return(pqtl_clean)
}

# ============================================================================
# Load GWAS data (.tsv.gz)
# ============================================================================

load_gwas_tsv <- function(gwas_file, chr_filter = NULL, 
                          start_pos = NULL, end_pos = NULL) {
  
  
  gwas <- fread(gwas_file)
  
  
  # Auto-detect columns
  col_names <- colnames(gwas)
  
  find_col <- function(patterns) {
    matched <- grep(paste(patterns, collapse = "|"), col_names, ignore.case = TRUE)
    if(length(matched) > 0) return(col_names[matched[1]])
    return(NA)
  }
  
  chr_col <- find_col(c("^chr$$", "^chromosome$$", "^CHR$$", "^#CHROM$$"))
  pos_col <- "base_pair_location"
  rsid_col <- "rsid"
  ea_col <- find_col(c("^effect_allele$$", "^EA$$", "^A1$$", "^ALT$$", "^Allele1$"))
  oa_col <- find_col(c("^other_allele$$", "^OA$$", "^A2$$", "^REF$$", "^Allele2$"))
  beta_col <- find_col(c("^beta$$", "^BETA$$", "^effect$"))
  se_col <- find_col(c("^se$$", "^SE$$", "^standard_error$"))
  pval_col <- find_col(c("^pval$$", "^P$$", "^p_value$"))
  freq_col <- "effect_allele_frequency"
  
  # Check required columns
  if(is.na(chr_col) | is.na(pos_col) | is.na(beta_col) | is.na(se_col) | is.na(pval_col)) {
    stop(paste0("Cannot find required columns.\n",
                "Available: ", paste(col_names, collapse = ", ")))
  }
  
  # Rename columns
  gwas_clean <- gwas
  names(gwas_clean)[names(gwas_clean) == chr_col] <- "chr"
  names(gwas_clean)[names(gwas_clean) == pos_col] <- "pos"
  names(gwas_clean)[names(gwas_clean) == beta_col] <- "beta"
  names(gwas_clean)[names(gwas_clean) == se_col] <- "se"
  names(gwas_clean)[names(gwas_clean) == pval_col] <- "pval"
  
  if(!is.na(rsid_col)) {
    names(gwas_clean)[names(gwas_clean) == rsid_col] <- "rsid"
  } else {
    gwas_clean$rsid <- paste0(gwas_clean$chr, ":", gwas_clean$pos)
    cat("  Note: Creating rsid as chr:pos\n")
  }
  
  if(!is.na(ea_col)) {
    names(gwas_clean)[names(gwas_clean) == ea_col] <- "effect_allele"
  } else {
    gwas_clean$effect_allele <- NA
    warning("Effect allele not found in GWAS")
  }
  
  if(!is.na(oa_col)) {
    names(gwas_clean)[names(gwas_clean) == oa_col] <- "other_allele"
  } else {
    gwas_clean$other_allele <- NA
    cat("  Note: Other allele not found in GWAS\n")
  }
  
  if(!is.na(freq_col)) {
    names(gwas_clean)[names(gwas_clean) == freq_col] <- "freq"
  } else {
    gwas_clean$freq <- NA
  }
  
  gwas_clean <- gwas_clean %>%
    mutate(chr = as.integer(gsub("chr", "", chr))) %>%
    select(rsid, chr, pos, effect_allele, other_allele, beta, se, pval, freq) %>%
    filter(!is.na(beta), !is.na(se), !is.na(pval))
  
  cat(paste0("  Loaded ", nrow(gwas_clean), " variants\n"))
  
  return(gwas_clean)
}

# ============================================================================
# Harmonize data with proper allele checking
# ============================================================================

harmonize_data <- function(pqtl_data, gwas_data, freq_threshold = 0.2) {
  # Merge by rsid first, then by chr:pos for missing rsids
  merged <- merge(pqtl_data, gwas_data, 
                  by.x = c("rsid", "chr", "pos"), 
                  by.y = c("rsid", "chr", "pos"),
                  suffixes = c("_pqtl", "_gwas"))
  
  if(nrow(merged) == 0) {
    warning("No overlapping variants found")
    return(NULL)
  }
  
  # Initialize harmonization columns
  merged$harmonized <- FALSE
  merged$beta_gwas_harm <- merged$beta_gwas
  merged$effect_allele_harm <- merged$effect_allele_pqtl
  
  # Check for palindromic SNPs (A/T or G/C pairs)
  is_palindromic <- (merged$effect_allele_pqtl == "A" & merged$other_allele_pqtl == "T") |
    (merged$effect_allele_pqtl == "T" & merged$other_allele_pqtl == "A") |
    (merged$effect_allele_pqtl == "G" & merged$other_allele_pqtl == "C") |
    (merged$effect_allele_pqtl == "C" & merged$other_allele_pqtl == "G")
  
  # Strand complement function
  complement <- function(allele) {
    chartr("ATGC", "TAGC", allele)
  }
  
  # Check each variant for allele alignment
  for(i in 1:nrow(merged)) {
    ea_pqtl <- merged$effect_allele_pqtl[i]
    oa_pqtl <- merged$other_allele_pqtl[i]
    ea_gwas <- merged$effect_allele_gwas[i]
    oa_gwas <- merged$other_allele_gwas[i]
    freq_pqtl <- merged$freq_pqtl[i]
    freq_gwas <- merged$freq_gwas[i]
    
    # Special handling for palindromic SNPs
    if(is_palindromic[i]) {
      freq_diff_forward <- abs(freq_pqtl - freq_gwas)
      freq_diff_flipped <- abs(freq_pqtl - (1 - freq_gwas))
      
      # Check if forward orientation matches
      if(freq_diff_forward < freq_threshold) {
        merged$harmonized[i] <- TRUE
        # No flipping needed
      }
      # Check if flipped orientation matches
      else if(freq_diff_flipped < freq_threshold) {
        merged$harmonized[i] <- TRUE
        merged$beta_gwas_harm[i] <- (-1) * merged$beta_gwas[i]
        merged$effect_allele_harm[i] <- ea_pqtl
      }
      # If neither orientation matches, exclude this SNP
      else {
        merged$harmonized[i] <- FALSE
      }
    } else {  # Non-palindromic SNPs - standard harmonization 
      # Case 1: Forward match (effect alleles match)
      if(ea_pqtl == ea_gwas & oa_pqtl == oa_gwas) {
        merged$harmonized[i] <- TRUE
        # No flipping needed
      }
      # Case 2: Reverse match (alleles swapped)
      else if(ea_pqtl == oa_gwas & oa_pqtl == ea_gwas) {
        merged$harmonized[i] <- TRUE
        merged$beta_gwas_harm[i] <- (-1) * merged$beta_gwas[i]
        merged$effect_allele_harm[i] <- ea_pqtl
      }
      # Case 3: Strand flip (complement matches)
      else if(complement(ea_pqtl) == ea_gwas & complement(oa_pqtl) == oa_gwas) {
        merged$harmonized[i] <- TRUE
        # No beta flipping, just strand difference
      }
      # Case 4: Strand flip + reverse
      else if(complement(ea_pqtl) == oa_gwas & complement(oa_pqtl) == ea_gwas) {
        merged$harmonized[i] <- TRUE
        merged$beta_gwas_harm[i] <- (-1) * merged$beta_gwas[i]
        merged$effect_allele_harm[i] <- ea_pqtl
      }
    }
  }
  
  harmonized_data <- merged[merged$harmonized, ]
  harmonized_data <- as.data.frame(harmonized_data)
  rownames(harmonized_data) <- NULL
  
  message(sprintf("Harmonized %d of %d variants (%.1f%%), including %d palindromic SNPs", 
                  nrow(harmonized_data), nrow(merged),
                  100 * nrow(harmonized_data) / nrow(merged),
                  sum(is_palindromic & merged$harmonized)))
  
  return(harmonized_data)
  
}


# ============================================================================
# NEW: Calculate LD matrix from genotype data
# ============================================================================

calculate_ld_for_snps <- function(plink_prefix, snp_list, 
                                  plink2_path = "plink2",
                                  temp_dir = NULL) {
  
  # Use specified temp directory or system default
  if(is.null(temp_dir)) {
    temp_dir <- tempdir()
  } else {
    # Create temp directory if it doesn't exist
    if(!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE)
    }
  }
  
  # Create temporary files in specified directory
  temp_snp_file <- file.path(temp_dir, paste0("snplist_", 
                                              format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                              sample(1000:9999, 1), ".txt"))
  
  temp_prefix <- file.path(temp_dir, paste0("ld_calc_", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                            sample(1000:9999, 1)))
  
  # Write SNP list
  write.table(snp_list, temp_snp_file, 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat(sprintf("  Temp files in: %s\n", temp_dir))
  
  # Extract only these SNPs
  cmd_extract <- sprintf("%s --bfile %s --extract %s --make-bed --out %s",
                         plink2_path, plink_prefix, temp_snp_file, temp_prefix)
  system(cmd_extract)
  
  # Calculate LD with --r-unphased (for unphased genotype data)
  cmd_ld <- sprintf("%s --bfile %s --r-unphased square --out %s",
                    plink2_path, temp_prefix, temp_prefix)
  system(cmd_ld)
  
  # FIXED: Check for different possible output file names
  ld_file <- paste0(temp_prefix, ".unphased.vcor1")  # PLINK2 format
  
  if(!file.exists(ld_file)) {
    # Try alternative extension
    ld_file <- paste0(temp_prefix, ".ld")
  }
  
  if(!file.exists(ld_file)) {
    # List what files were actually created
    temp_files <- list.files(temp_dir, pattern = basename(temp_prefix))
    stop(paste0("LD matrix file not created. Files found: ", 
                paste(temp_files, collapse = ", ")))
  }
  
  ld_matrix <- as.matrix(fread(ld_file,header = F))
  
  # Read BIM to get SNP info
  bim <- fread(paste0(temp_prefix, ".bim"))
  colnames(bim) <- c("chr", "rsid", "cm", "pos", "A1", "A2")
  
  snp_info <- data.frame(
    rsid = bim$rsid,
    chr = bim$chr,
    pos = bim$pos,
    ref_allele = bim$A2,
    alt_allele = bim$A1,
    ld_allele = bim$A1  # LD calculated for A1
  )
  rownames(ld_matrix) <- bim$rsid
  colnames(ld_matrix) <- bim$rsid
  
  # Clean up temporary files
  unlink(temp_snp_file)
  unlink(paste0(temp_prefix, ".*"))
  
  cat("  Temp files cleaned up\n")
  
  return(list(ld_matrix = ld_matrix, snp_info = snp_info))
}

# ============================================================================
# NEW: Align summary stats to LD matrix
# ============================================================================

align_to_ld <- function(harm_data, ld_snp_info) {
  # Merge to get LD allele info
  merged <- merge(harm_data, ld_snp_info, by = c("rsid", "chr", "pos"))
  
  # Check if effect_allele matches LD allele
  merged$aligned <- merged$effect_allele_harm == merged$ld_allele
  
  # For misaligned SNPs, flip beta
  merged$beta_pqtl_aligned <- ifelse(merged$aligned, 
                                     merged$beta_pqtl, 
                                     -merged$beta_pqtl)
  
  merged$beta_gwas_aligned <- ifelse(merged$aligned, 
                                     merged$beta_gwas_harm, 
                                     -merged$beta_gwas_harm)
  
  # Update effect allele to match LD
  merged$effect_allele_final <- merged$ld_allele
  
  n_flipped <- sum(!merged$aligned)
  if(n_flipped > 0) {
    message(sprintf("  Flipped beta for %d SNPs to align with LD matrix", n_flipped))
  }
  
  return(merged)
}


# ============================================================================
# Run colocalization
# ============================================================================

run_coloc_analysis <- function(harmonized_data, 
                               protein_name) {
  
  if(nrow(harmonized_data) < 50) {
    warning(paste0("Only ", nrow(harmonized_data), " SNPs for ", protein_name))
  }
  
  # Prepare datasets
  dataset1 <- list(
    beta = harmonized_data$beta_pqtl,
    varbeta = harmonized_data$se_pqtl^2,
    snp = harmonized_data$rsid,
    position = harmonized_data$pos,
    type = "quant",
    MAF = harmonized_data$freq_pqtl,
    sdY = 1
  )

  
  dataset2 <- list(
    beta = harmonized_data$beta_gwas,
    varbeta = harmonized_data$se_gwas^2,
    snp = harmonized_data$rsid,
    position = harmonized_data$pos,
    type = "cc",
    MAF = harmonized_data$freq_gwas
  )

  # Run colocalization
  coloc_result <- coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
  
  # Extract summary
  summary_df <- data.frame(
    protein = protein_name,
    nsnps = coloc_result$summary["nsnps"],
    PP.H0 = coloc_result$summary["PP.H0.abf"],
    PP.H1 = coloc_result$summary["PP.H1.abf"],
    PP.H2 = coloc_result$summary["PP.H2.abf"],
    PP.H3 = coloc_result$summary["PP.H3.abf"],
    PP.H4 = coloc_result$summary["PP.H4.abf"],
    row.names = NULL
  ) %>%
    mutate(
      interpretation = case_when(
        PP.H4 > 0.8 ~ "Strong colocalization",
        PP.H4 > 0.5 ~ "Moderate colocalization",
        PP.H3 > 0.8 ~ "Both associated, different variants",
        TRUE ~ "Weak/unclear"
      )
    )
  
  return(list(
    summary = summary_df,
    detailed = coloc_result
  ))
}


# ============================================================================
# NEW: Run coloc.susie (multiple causal variants)
# ============================================================================

run_coloc_susie_analysis <- function(harm_data_aligned,
                                     N_pqtl = 7213,
                                     N_GWAS = 1865284,
                                     N_case = 3748,
                                     ld_matrix) {
  
  

  # Prepare dataset for pQTL
  dataset_pqtl <- list(
    beta = harm_data_aligned$beta_pqtl_aligned,
    varbeta = harm_data_aligned$se_pqtl^2,
    MAF = harm_data_aligned$freq_pqtl,
    type = "quant",
    snp = harm_data_aligned$rsid,
    LD = ld_matrix,
    position = harm_data_aligned$pos,
    N = N_pqtl,
    sdY = 1
  )
  
  # Prepare dataset for GWAS
  dataset_gwas <- list(
    beta = harm_data_aligned$beta_gwas_aligned,
    varbeta = harm_data_aligned$se_gwas^2,
    type = "cc",
    N = N_GWAS,
    s = N_case/N_GWAS,
    MAF = harm_data_aligned$freq_gwas,
    snp = harm_data_aligned$rsid,
    position = harm_data_aligned$pos,
    LD = ld_matrix
  )
  
  S3=runsusie(dataset_pqtl)
  S4=runsusie(dataset_gwas)
  
  
  # Run coloc.susie
  result <- tryCatch({
    susie.res=coloc.susie(S3,S4)
  }, error = function(e) {
    return(NULL)
  })
  
  return(result)
}

run_coloc_for_proteins <- function(sig_proteins,
                                   pqtl_dir,
                                   gwas_data,
                                   N_pqtl0 = 7213,
                                   N_GWAS0 = 1865284,
                                   N_case0 = 3748,
                                   pqtl_plink_prefix,
                                   output_dir = "coloc_results",
                                   temp_dir = NULL,  # NEW: specify temp directory
                                   method = "both",
                                   plink2_path = "plink2"
 ) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_results <- list()
  failed <- character()
  
  for(i in 1:length(sig_proteins)) {
    
    
    # Find pQTL file
    pqtl_file <- file.path(pqtl_dir, paste0(sig_proteins[i], ".PHENO1.glm.linear"))
    if(!file.exists(pqtl_file)) {
      cat("  ✗ pQTL file not found\n")
      failed <- c(failed, sig_proteins[i])
      next
    }
    
    # Load and process
    pqtl_data <-load_pqtl_glm(pqtl_file)
    
    
    harmonized <- tryCatch({
      harmonize_data(pqtl_data, gwas_data)
    }, error = function(e) {
      cat(paste0("  ✗ Harmonization error: ", e$message, "\n"))
      return(NULL)
    })
    
    coloc_result = list()
  
  
    # Run coloc.abf if requested
    if(method %in% c("abf", "both")) {

      cat("  Running coloc.abf...\n")
      coloc_abf <- tryCatch({
        run_coloc_analysis(harmonized, sig_proteins[i])
      }, error = function(e) {
        cat(paste0("  ✗ Coloc error: ", e$message, "\n"))
        return(NULL)
      })
      
        coloc_result$abf <- coloc_abf

      
    }
    
    # Run coloc.susie if requested
    if(method %in% c("susie", "both")) {
      cat("  Calculating LD matrix...\n")
      
      snp_list <- harmonized$rsid
      
      ld_result <- tryCatch({
        calculate_ld_for_snps(pqtl_plink_prefix, snp_list, 
                              plink2_path, temp_dir)  # Pass temp_dir here
      }, error = function(e) {
        cat(paste0("  ✗ LD calculation error: ", e$message, "\n"))
        return(NULL)
      })
      
      if(!is.null(ld_result)) {
        # Match SNPs between harmonized data and LD
        snps_in_ld <- ld_result$snp_info$rsid
        snps_common <- intersect(snp_list, snps_in_ld)
        
          # Filter and align
          harm_filtered <- harmonized[harmonized$rsid %in% snps_common, ]
          ld_snp_info_filtered <- ld_result$snp_info[ld_result$snp_info$rsid %in% snps_common, ]
          
          # Align to LD
          harm_aligned <- align_to_ld(harm_filtered, ld_snp_info_filtered)
          
          # Reorder
          harm_aligned <- harm_aligned[match(snps_common, harm_aligned$rsid), ]
          idx <- match(snps_common, ld_result$snp_info$rsid)
          ld_matrix <- ld_result$ld_matrix[idx, idx]
          
          cat("  Running coloc.susie...\n")
          coloc_susie <- tryCatch({
            run_coloc_susie_analysis(harm_aligned, N_pqtl = N_pqtl0,
                                     N_GWAS = N_GWAS0,
                                     N_case = N_case0,ld_matrix)
          }, error = function(e) {
            cat(paste0("  ✗ coloc.susie error: ", e$message, "\n"))
            return(NULL)
          })

            coloc_result$susie <- coloc_susie
           

      }
    }
    
   
    if(!is.null(coloc_result)) {
      all_results[[sig_proteins[i]]] <- coloc_result
      saveRDS(coloc_result, file.path(output_dir, paste0(sig_proteins[i], "_coloc.rds")))
      cat(paste0("  ✓ PP.H4 = ", sprintf("%.3f", coloc_result$summary$PP.H4), "\n"))
    } else {
      failed <- c(failed, sig_proteins[i])
    }
  }
  
  cat(paste0("Success: ", length(all_results), " | Failed: ", length(failed), "\n"))

  return(list(results = all_results, failed = failed))

}
