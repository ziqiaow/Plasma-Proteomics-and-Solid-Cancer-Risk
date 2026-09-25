#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-20:30:00

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/HGDP_1000G/GRCh38/mergedplink \
--indep-pairwise 200kb 0.5 \



/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/HGDP_1000G/GRCh38/mergedplink \
--extract /dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/PCA/redo/plink2.prune.in \
--freq counts \
--pca allele-wts vcols=chrom,ref,alt \
--out ref_pcs

#You can then project onto those PCs with
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/PCA/ARIC_gwas_merge_final \
--extract /dcs04/nilanjan/data/zwang/collaboration/ARIC_Cancer/PCA/redo/plink2.prune.in \
--read-freq ref_pcs.acount \
--score ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
               variance-standardize \
--score-col-nums 6-15 \
--out ARIC_projection