# Analysis code for Identifying Plasma Proteins Associated with Risk of Solid Cancers 
This is the code for the analysis in the manuscript "Identifying Plasma Proteins Associated with Risk of Solid Cancers: A 25-Year Prospective Analysis of 4,712 Circulating Proteins in the ARIC Study."

This repository contains code and documentation for evaluating associations between plasma protein concentrations and the risk of solid cancers in the Atherosclerosis Risk in Communities (ARIC) study.

## Reference
Wang Z, Burk VA, Huang Z, Zahed H, Muller DC, Yarmolinsky J, Lee MA, Joshu CE, Lin ZC, Prizment A, Butler KR, Couper DJ, Smith-Byrne K, Kolijn PM, Vermeulen RCH, Riboli E, Gunter MJ, Coresh J, Chatterjee N, Platz EA. Identifying Plasma Proteins Associated with Risk of Solid Cancers: A 25-Year Prospective Analysis of 4,712 Circulating Proteins in the ARIC Study. medRxiv [Preprint]. 2026 Mar 18:2026.03.16.26348527. doi: 10.64898/2026.03.16.26348527. PMID: 41890985; PMCID: PMC13015667.
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015667/)


## Overview

The repository includes code for:

- Data cleaning and preparation
- Plasma proteomic data processing
- Protein dimension-reduction and latent-factor analyses
- Protein–cancer association analyses
- Sensitivity analyses
- Protein trajectory and reproducibility analyses
- Pathway enrichment analyses
- Main and supplementary figure generation
- Race-stratified analyses
- Prediction analyses
- Colocalization analyses
- Bootstrap analyses

The cancer outcomes evaluated include:

- Bladder cancer
- Colorectal cancer
- Colon cancer
- Rectal cancer
- Kidney cancer
- Liver cancer
- Lung cancer
- Pancreatic cancer
- Prostate cancer

Individual-level ARIC and proteomic data are not included in this repository because of data-use agreements and participant-confidentiality restrictions.


## Repository Structure

```text
.
├── 1_data_cleaning/
│   ├── cancer_file.R
│   ├── data_clean_submit.R
│   ├── PCA_step1.sh
│   ├── PCA_step2.sh
│   ├── PCA_step3.sh
│   └── submit_mofa.R
│
├── 2_main_analyses/
│   ├── original_protein_cancer_association/
│   │   ├── submit_cancer_minimal_smk.R
│   │   ├── submit_cancer.R
│   │   ├── submit_liver_kidney_lung.R
│   │   └── submit_rectal.R
│   │
│   ├── pathway/
│   │   └── clusterProfiler_liver.R
│   │
│   ├── sensitivity_protein_cancer_association/
│   │   ├── submit_cancer_lag5yr.R
│   │   ├── submit_cancer_trunc.R
│   │   ├── submit_kidney_poorfunction.R
│   │   ├── submit_liver_kidney_lung_lag5yr.R
│   │   ├── submit_liver_kidney_lung_trunc.R
│   │   └── submit_liver_poorfunction.R
│   │
│   └── time_trajectory/
│       ├── protein_correlation_between_visit2.R
│       └── protein_ICC.R
│
├── 3_main_and_extended_figures/
│   ├── Code_Figure2.R
│   ├── Code_Figure3.R
│   ├── Code_Figure4.R
│   ├── Code_Figure5.R
│   └── Code_Supplementary_Figures.R
│
├── 4_revision1/
│   ├── bootstrap/
│   │   ├── bootstrap_liver_lung_kidney.R
│   │   └── bootstrap_prostate.R
│   │
│   ├── colocalization/
│   │   ├── coloc_function_susie.R
│   │   ├── coloc_kidney_EA_susie.R
│   │   ├── coloc_liver_EA_susie.R
│   │   ├── coloc_lung_EA_susie.R
│   │   └── coloc_prostate_EA_susie.R
│   │
│   ├── prediction_model/
│   │   ├── delta_AUC.R
│   │   └── prediction_UKB.R
│   │
│   ├── race_stratified_analysis/
│   │   ├── submit_cancer_minimal_smk_black.R
│   │   └── submit_cancer_minimal_smk_white.R
│   │
│   ├── raw_protein_cancer_association/
│   │   ├── submit_cancer.R
│   │   ├── submit_liver_kidney_lung.R
│   │   ├── submit_rectal.R
│   │   └── submit_summary.R
│   │
│   └── Supplementary_Figure_Forest_Plot.R
│
├── LICENSE
└── README.md
```

## Directory Descriptions

### `1_data_cleaning`

This directory contains scripts for preparing the cancer outcome data, cleaning the analysis datasets, and processing the proteomic data.

- `cancer_file.R`: Prepares cancer-related variables and outcome data.
- `data_clean_submit.R`: Performs data cleaning and prepares analysis datasets.
- `PCA_step1.sh`, `PCA_step2.sh`, and `PCA_step3.sh`: Perform PCA of the genotype data.
- `submit_mofa.R`: Calculates PEER factors for proteomics data.

### `2_main_analyses`

This directory contains the primary protein–cancer association analyses, pathway analyses, sensitivity analyses, and protein trajectory analyses.

#### `original_protein_cancer_association`

These scripts evaluate associations between individual plasma proteins and incident cancers.

#### `pathway`

Performs pathway enrichment analyses for liver cancer.

#### `sensitivity_protein_cancer_association`

These scripts evaluate the robustness of protein–cancer associations under alternative analytic conditions.

- `submit_cancer_lag5yr.R` and `submit_liver_kidney_lung_lag5yr.R`: Performs analyses using a five-year lag period.
- `submit_cancer_trunc.R` and `submit_liver_kidney_lung_trunc.R`: Performs analyses using truncated 10-year follow-up.
- `submit_kidney_poorfunction.R`: Performs kidney cancer analyses after excluding participants with poor kidney function.
- `submit_liver_poorfunction.R`: Performs liver cancer analyses after excluding participants with poor liver function.

#### `time_trajectory`

These scripts evaluate the longitudinal stability and reproducibility of protein measurements.

- `protein_correlation_between_visit2.R`: Calculates correlations between protein measurements at Visits 2,3,5.
- `protein_ICC.R`: Estimates intraclass correlation coefficients for protein measurements.

### `3_main_and_extended_figures`

This directory contains scripts used to generate the main and supplementary figures.

### `4_revision1`

This directory contains analyses added or updated during manuscript revision.

#### `bootstrap`

These scripts perform bootstrap analyses for selected cancer types.

#### `colocalization`

These scripts perform genetic colocalization analyses using SuSiE-based methods.

- `coloc_function_susie.R`: Defines functions used in the colocalization analyses.
- `coloc_kidney_EA_susie.R`
- `coloc_liver_EA_susie.R`
- `coloc_lung_EA_susie.R`
- `coloc_prostate_EA_susie.R`

#### `prediction_model`

These scripts evaluate the predictive performance of protein-based models.

- `delta_AUC.R`: Calculates changes in the area under the receiver operating characteristic curve.
- `prediction_UKB.R`: Performs prediction analyses using UK Biobank data.

#### `race_stratified_analysis`

These scripts perform protein–cancer association analyses separately among Black and White participants.

- `submit_cancer_minimal_smk_black.R`
- `submit_cancer_minimal_smk_white.R`

#### `raw_protein_cancer_association`

These scripts evaluate protein–cancer associations using ANML-normalized protein measurements, without residualization using PC or PEER Factors.


#### Supplementary figure

- `Supplementary_Figure_Forest_Plot.R`: Generates the supplementary forest plot.

## Statistical Analyses

The primary protein–cancer associations were evaluated using Cox proportional hazards regression models. Plasma protein measurements were obtained from blood specimens collected at Visit 2, which was defined as the baseline visit for the proteomic analyses.

The models evaluated associations between individual protein measurements and incident cancer outcomes while adjusting for prespecified demographic, behavioral, clinical, lifestyle, genetic, and technical covariates.

Sensitivity analyses evaluated the robustness of the findings to:

- Five-year lagged follow-up
- Truncated 10-year follow-up
- Exclusion of participants with poor liver/kidney function
- Alternative protein-processing approaches

Race-stratified analyses were conducted separately among Black and White participants. Heterogeneity between the two groups was evaluated using Cochran’s Q test. Because this test may have limited power, particularly when the number of cancer events is small, the absence of statistically significant heterogeneity was interpreted cautiously.

## Protein Trajectory Analyses

Longitudinal protein measurements were evaluated using:

- Correlations between protein measurements at Visits 2, 3, and 5
- Intraclass correlation coefficients to assess measurement reproducibility

## Pathway Analyses

Pathway enrichment analyses were performed using `clusterProfiler` for selected cancer-associated proteins in liver cancer.

## Prediction Analyses

Prediction analyses evaluated the performance of protein-based models and changes in predictive discrimination for liver cancer. These analyses included evaluation of the change in area under the receiver operating characteristic curve and prediction analyses using UK Biobank data.

## Colocalization Analyses

Colocalization analyses assessed whether genetic signals associated with selected proteins and cancer outcomes shared a common causal variant. Single causal variant abf and multiple causal variants SuSiE-based methods were used for the colocalization analyses.


## Reproducibility

To reproduce the analyses:

1. Obtain authorized access to the required ARIC and proteomic datasets.
2. Clone this repository:

   ```bash
   git clone https://github.com/ziqiaow/Plasma-Proteomics-and-Solid-Cancer-Risk.git
   ```

## Contact

For questions about the code or analyses, please contact Ziqiao Wang.
```
