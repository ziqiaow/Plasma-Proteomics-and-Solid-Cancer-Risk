source("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/coloc_function_susie.R")
# 1. Set the path 
folder_path <- "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/EA"

# 2. Get all file names matching the specific pattern
file_names <- list.files(path = folder_path, pattern = "\\.PHENO1\\.glm\\.linear$")

# 3. Strip the suffix to keep only the 'xxx' part
prot <- gsub("\\.PHENO1\\.glm\\.linear$", "", file_names)


#load gwas file once

load_gwas_tsv <- function(gwas_file) {
  
  
  gwas <- fread(gwas_file)
  
  
  # Auto-detect columns
  col_names <- colnames(gwas)
  
  find_col <- function(patterns) {
    matched <- grep(paste(patterns, collapse = "|"), col_names, ignore.case = TRUE)
    if(length(matched) > 0) return(col_names[matched[1]])
    return(NA)
  }
  
  chr_col <- "hm_chrom"
  pos_col <- "hm_pos"
  rsid_col <- "hm_rsid"
  ea_col <- "hm_effect_allele"
  oa_col <- "hm_other_allele"
  beta_col <- "hm_beta"
  se_col <- "standard_error"
  pval_col <- "p_value"
  freq_col <- "hm_effect_allele_frequency"
  # Create a named vector of found columns
  cols_to_keep <- c(
    chr = chr_col,
    pos = pos_col,
    rsid = rsid_col,
    effect_allele = ea_col,
    other_allele = oa_col,
    beta = beta_col,
    se = se_col,
    pval = pval_col,
    freq = freq_col
  )
  
  
  # Subset and rename the dataframe
  gwas_clean <- gwas[, ..cols_to_keep, drop = FALSE]
  colnames(gwas_clean) = names(cols_to_keep)
  
  gwas_clean <- gwas_clean %>%
    mutate(chr = as.integer(gsub("chr", "", chr))) %>%
    select(rsid, chr, pos, effect_allele, other_allele, beta, se, pval, freq) %>%
    filter(!is.na(beta), !is.na(se), !is.na(pval))
  
  cat(paste0("  Loaded ", nrow(gwas_clean), " variants\n"))
  
  return(gwas_clean)
}

gwas_clean = load_gwas_tsv(gwas_file = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/GCST004748.h.tsv.gz")

results <- run_coloc_for_proteins(
  sig_proteins = prot,
  pqtl_dir = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/EA/",
  gwas_data = gwas_clean,
  N_pqtl0 = 7213,
  N_GWAS0 = 85716,
  N_case0 = 29266,
  pqtl_plink_prefix = "/dcs04/legacy-dcs01-arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/plink/ARICEA_TOPMedimputed_maf005.imp30.rsid",
  output_dir = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/coloc_results/EA/lung",
  plink2_path = "/dcs04/nilanjan/data/zwang/tools/20250129/plink2"
 
)
