#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-02:30:00

/dcs04/nilanjan/data/zwang/tools/plink \
--bfile /dcs04/legacy-dcs01-arking/ARIC_static/ARIC_Data/GWAS/Black/1_raw/f3v2/B_ARIC_gwas_frz3v2 \
--bmerge /dcs04/legacy-dcs01-arking/ARIC_static/ARIC_Data/GWAS/White/1_raw/f3v2/W_ARIC_gwas_frz3v2 \
--geno 0.05 \
--maf 0.01 \
--mind 0.05 \
--make-bed --out /dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/PCA/ARIC_gwas_merge
