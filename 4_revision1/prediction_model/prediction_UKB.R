#Train a prediction model using the external UK Biobank data, we report a repeated 3-fold cross-validation.
#---------------------------------------------------------------------------------------------------------------------
#for manuscript
library(caret)
library(pROC)

#For liver cancer only
#UKB protein value
prot_ukb=readRDS("./protein/UKB/proteomics_imputed_goodID_EUR.rds")
filter.16 = readRDS("./collaboration/ARIC_Cancer/UKB/cancer_data_cleaned_liverC22.rds")
id = intersect(rownames(prot_ukb),filter.16$eid)
cov_ukb_clean = filter.16[match(id,filter.16$eid),]
table(cov_ukb_clean$liver_incidence_ind)
# 0     1 
# 27161    36 
prot_ukb = prot_ukb[match(id,rownames(prot_ukb)),]
#[1] 27197  2917

anno = readRDS("./protein/UKB/protein_annotation.rds")

#load results for model 2
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_all.RData")
res_model2 = res_final

selected_aric_ids = res_model2$seqid_in_sample[which(res_model2$fdr< 0.05 & res_model2$cancer=="Liver")]
selected_uniprot  <- res_model2$`UniProt ID`[which(res_model2$fdr< 0.05 & res_model2$cancer=="Liver")]

protein_map <- data.frame(
  aric_id = selected_aric_ids,
  uniprot = selected_uniprot,
  ukb_id  = anno$prot_name[match(selected_uniprot, anno$UniProt)],
  stringsAsFactors = FALSE
)

protein_map <- protein_map[
  !is.na(protein_map$ukb_id) &
    protein_map$ukb_id %in% colnames(prot_ukb),
]

cat("Number mapped to UKB/Olink:", nrow(protein_map), "\n")

# extract the test proteins
data_ukb <- prot_ukb[, match(protein_map$ukb_id,colnames(prot_ukb))]
identical(rownames(data_ukb),as.character(cov_ukb_clean$eid))
#[1] TRUE



# Create data
pheno=readRDS("./ukbPheno_032022.rds")
tmp = pheno[match(as.character(cov_ukb_clean$eid),as.character(pheno$IID)),]
cov_ukb_clean = cbind(cov_ukb_clean,tmp)

model_data_test <- data.frame(
  entry=cov_ukb_clean$study.entry.age,exit=cov_ukb_clean$study.exit.age,
  event = cov_ukb_clean$liver_incidence_ind ,
  data_ukb, cov_ukb_clean[,c("sex","BMI","SmokingStatus","AlcoholStatus","age")]
)


model_data_test$event <- factor(
  ifelse(model_data_test$event == 1, "Yes", "No"), 
  levels = c("Yes", "No")
)

model_data_test = model_data_test[complete.cases(model_data_test),]
table(model_data_test$event)
formula_str <- paste("event ~", paste(protein_map$ukb_id, collapse = " + "),"+factor(sex)+BMI+factor(SmokingStatus)+factor(AlcoholStatus)+age")
formula_obj <- as.formula(formula_str)
# Define control settings for 3 folds, repeated 10 times to smooth out variability
cv_control <- trainControl(
  method = "repeatedcv",
  number = 3,                    # Set to 3 folds
  repeats = 10,                  # Increased repeats to stabilize the low event count
  summaryFunction = twoClassSummary, 
  classProbs = TRUE,             
  savePredictions = "final"      
)

# Train the model
set.seed(123)
cv_model <- train(
  formula_obj,
  data = model_data_test, 
  method = "glm", 
  family = "binomial", 
  trControl = cv_control,
  metric = "ROC"                 
)

cv_model$results
# parameter       ROC       Sens      Spec      ROCSD     SensSD       SpecSD
#1      none 0.7522917 0.01439394 0.9999035 0.07763077 0.0327723 9.118612e-05

library(cvAUC)
ci.cvAUC(cv_model$pred$Yes,cv_model$pred$obs,folds = cv_model$pred$Resample)
#$cvAUC
# [1] 0.7522917
# 
# $se
# [1] 0.0142499
# 
# $ci
# [1] 0.7243624 0.7802210

cvAUC(cv_model$pred$Yes,cv_model$pred$obs,folds = cv_model$pred$Resample)
#$cvAUC
#[1] 0.7522917

# ==============================================================================
# 1. EXTRACT INDICES & SET CV CONTROL FROM ORIGINAL MODEL
# ==============================================================================
# Use exact original partitions to ensure identical results
original_indices <- cv_model$control$index

fixed_cv_control <- trainControl(
  method = "repeatedcv",
  index = original_indices,      
  summaryFunction = twoClassSummary, 
  classProbs = TRUE,             
  savePredictions = "final"      
)

# ==============================================================================
# 2. DEFINE THE FORMULAS FOR ALL THREE MODELS
# ==============================================================================
proteins_str <- paste(protein_map$ukb_id, collapse = " + ")
proteins_str_model1 <- paste(protein_map_model1$ukb_id, collapse = " + ")
clinical_str <- "factor(sex) + BMI + factor(SmokingStatus) + factor(AlcoholStatus) + age"

# Model A: Risk Factors Only (The Baseline)
formula_rf <- as.formula(paste("event ~", clinical_str))

# Model B: Proteins in model 2 Only
formula_protein <- as.formula(paste("event ~", proteins_str))

# Model C: Proteins in model 2 + Risk Factors (The Full Model)
formula_protein_rf <- as.formula(paste("event ~", proteins_str, "+", clinical_str))

# ==============================================================================
# 3. TRAIN ALL THREE MODELS ON THE IDENTICAL SPLITS
# ==============================================================================

# Train Model A (Baseline)
cv_model_rf <- train(formula_rf, data = model_data_test, method = "glm", family = "binomial", trControl = fixed_cv_control, metric = "ROC")
ci.cvAUC(cv_model_rf$pred$Yes,cv_model_rf$pred$obs,folds = cv_model_rf$pred$Resample)
#$ci
#[1] 0.6526634 0.7065125
# Train Model B (Proteins in model 2 Only)
cv_model_protein <- train(formula_protein, data = model_data_test, method = "glm", family = "binomial", trControl = fixed_cv_control, metric = "ROC")
ci.cvAUC(cv_model_protein$pred$Yes,cv_model_protein$pred$obs,folds = cv_model_protein$pred$Resample)
#$ci
#[1] 0.6834938 0.7460561


# Train Model C (Full Model)
cv_model_protein_rf <- train(formula_protein_rf, data = model_data_test, method = "glm", family = "binomial", trControl = fixed_cv_control, metric = "ROC")
ci.cvAUC(cv_model_protein_rf$pred$Yes,cv_model_protein_rf$pred$obs,folds = cv_model_protein_rf$pred$Resample)
#$ci
#[1] 0.7243624 0.7802210


cv_model_protein$results$ROC-cv_model_rf$results$ROC #0.03518703
cv_model_protein_rf$results$ROC-cv_model_rf$results$ROC #0.07270376



#--------------------------------------------------------------------------------------------------------
#Report the weights associated with the risk factors for the final risk prediction model of liver cancer
library(survival)
library(caret)
#Create a copy of the dataset to standardize
standardised_data <- model_data_test

for (protein in protein_map$ukb_id) {
  p_mean <- mean(standardised_data[[protein]], na.rm = TRUE)
  p_sd   <- sd(standardised_data[[protein]], na.rm = TRUE)
  standardised_data[[protein]] <- (standardised_data[[protein]] - p_mean) / p_sd
}

standardised_data$sex           <- factor(standardised_data$sex)
standardised_data$SmokingStatus <- factor(standardised_data$SmokingStatus)
standardised_data$AlcoholStatus <- factor(standardised_data$AlcoholStatus)

print(levels(standardised_data$SmokingStatus))
print(levels(standardised_data$AlcoholStatus))

proteins_part <- paste(protein_map$ukb_id, collapse = " + ")
clinical_part <- "sex + BMI + SmokingStatus + AlcoholStatus + age"

formula_str <- paste("event ~", proteins_part, "+", clinical_part)
formula_obj <- as.formula(formula_str)

# Fit the multivariable logistic regression model
final_combined_model <- glm(
  formula_obj, 
  data = standardised_data, 
  family = binomial(link = "logit")
)


model_summary <- summary(final_combined_model)$coefficients
log_odds  <- model_summary[-1, "Estimate"]

reporting_table <- data.frame(
  Risk_factor = rownames(model_summary)[-1],
  Weight_beta  = log_odds
)
reporting_table$protein_uniprot_ID = c(protein_map$uniprot,rep(NA,7))


