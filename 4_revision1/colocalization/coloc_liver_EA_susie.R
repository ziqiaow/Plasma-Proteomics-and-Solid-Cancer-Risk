source("/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/coloc_function_susie.R")
# 1. Set the path 
folder_path <- "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/EA"

# 2. Get all file names matching the specific pattern
file_names <- list.files(path = folder_path, pattern = "\\.PHENO1\\.glm\\.linear$")

# 3. Strip the suffix to keep only the 'xxx' part
prot <- gsub("\\.PHENO1\\.glm\\.linear$", "", file_names)


#load gwas file once
gwas_clean = load_gwas_tsv(gwas_file = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/GCST90809296.h.tsv.gz")

results <- run_coloc_for_proteins(
  sig_proteins = prot,
  pqtl_dir = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/data/EA/",
  gwas_data = gwas_clean,
  N_pqtl0 = 7213,
  N_GWAS0 = 1865284,
  N_case0 = 3748,
  pqtl_plink_prefix = "/dcs04/legacy-dcs01-arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/plink/ARICEA_TOPMedimputed_maf005.imp30.rsid",
  output_dir = "/dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/coloc/coloc_results/EA/liver",
  plink2_path = "/dcs04/nilanjan/data/zwang/tools/20250129/plink2"
 
)
