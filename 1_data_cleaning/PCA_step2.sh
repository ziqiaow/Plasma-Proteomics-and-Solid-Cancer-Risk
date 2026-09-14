#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0-02:30:00


#module load plink/1.90b

/dcs04/nilanjan/data/zwang/tools/plink \
--merge-list /dcs04/nilanjan/data/zwang/plink/data/mergelist.txt --make-bed --out mergedplink
--indep-pairwise 200kb 0.5 \
