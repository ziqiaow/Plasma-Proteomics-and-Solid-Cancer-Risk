# Analysis code for Identifying Plasma Proteins Associated with Risk of Solid Cancers 
This is the code for the analysis in the manuscript "Identifying Plasma Proteins Associated with Risk of Solid Cancers: A 25-Year Prospective Analysis of 4,712 Circulating Proteins in the ARIC Study."

This repository contains code and documentation for evaluating associations between plasma protein concentrations and the risk of solid cancers in the Atherosclerosis Risk in Communities (ARIC) study.

## Reference
Wang Z, Burk VA, Huang Z, Zahed H, Muller DC, Yarmolinsky J, Lee MA, Joshu CE, Lin ZC, Prizment A, Butler KR, Couper DJ, Smith-Byrne K, Kolijn PM, Vermeulen RCH, Riboli E, Gunter MJ, Coresh J, Chatterjee N, Platz EA. Identifying Plasma Proteins Associated with Risk of Solid Cancers: A 25-Year Prospective Analysis of 4,712 Circulating Proteins in the ARIC Study. medRxiv [Preprint]. 2026 Mar 18:2026.03.16.26348527. doi: 10.64898/2026.03.16.26348527. PMID: 41890985; PMCID: PMC13015667.
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015667/)

## Overview

The primary analyses evaluated prospective associations between plasma protein levels measured at Visit 2 and incident solid cancers using Cox proportional hazards regression models. Analyses were conducted separately for individual cancer types and included adjustment for demographic, behavioral, clinical, genetic ancestry, and technical factors.

Cancer outcomes included:

- Bladder cancer
- Colorectal cancer (Colorectal, colon, rectal)
- Kidney cancer
- Liver cancer
- Lung cancer
- Pancreatic cancer
- Prostate cancer

## Repository structure
.
├── 1_data_cleaning
│   ├── cancer_file.R
│   ├── data_clean_submit.R
│   ├── PCA_step1.sh
│   ├── PCA_step2.sh
│   ├── PCA_step3.sh
│   └── submit_mofa.R
├── 2_main_analyses
│   ├── original_protein_cancer_association
│   │   ├── submit_cancer_minimal_smk.R
│   │   ├── submit_cancer.R
│   │   ├── submit_liver_kidney_lung.R
│   │   └── submit_rectal.R
│   ├── pathway
│   │   └── clusterProfiler_liver.R
│   ├── sensitivity_protein_cancer_association
│   │   ├── submit_cancer_lag5yr.R
│   │   ├── submit_cancer_trunc.R
│   │   ├── submit_kidney_poorfunction.R
│   │   ├── submit_liver_kidney_lung_lag5yr.R
│   │   ├── submit_liver_kidney_lung_trunc.R
│   │   └── submit_liver_poorfunction.R
│   └── time_trajectory
│       ├── protein_correlation_between_visit2.R
│       └── protein_ICC.R
├── 3_main_and_extended_figures
│   ├── Code_Figure2.R
│   ├── Code_Figure3.R
│   ├── Code_Figure4.R
│   ├── Code_Figure5.R
│   └── Code_Supplementary_Figures.R
├── 4_revision1
│   ├── bootstrap
│   │   ├── bootstrap_liver_lung_kidney.R
│   │   └── bootstrap_prostate.R
│   ├── colocalization
│   │   ├── coloc_function_susie.R
│   │   ├── coloc_kidney_EA_susie.R
│   │   ├── coloc_liver_EA_susie.R
│   │   ├── coloc_lung_EA_susie.R
│   │   └── coloc_prostate_EA_susie.R
│   ├── prediction_model
│   │   ├── delta_AUC.R
│   │   └── prediction_UKB.R
│   ├── race_stratified_analysis
│   │   ├── submit_cancer_minimal_smk_black.R
│   │   └── submit_cancer_minimal_smk_white.R
│   ├── raw_protein_cancer_association
│   │   ├── submit_cancer.R
│   │   ├── submit_liver_kidney_lung.R
│   │   ├── submit_rectal.R
│   │   └── submit_summary.R
│   └── Supplementary_Figure_Forest_Plot.R
├── LICENSE
└── README.md
