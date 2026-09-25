#Train a prediction model using the external UK Biobank data, we report a repeated 3-fold cross-validation.
#For liver cancer only
#load UKB protein values, phenotypes, and liver cancer data
#Note: The steps above are not included in the public GitHub repository due to UKB's strict data access policies

#---------------------------------------------------------------------------------------------------------------------
#for manuscript
library(caret)
library(pROC)


model_data_test <- data.frame(
  event = cov_clean$liver_incidence_ind ,
  protein, cov_clean[,c("sex","BMI","SmokingStatus","AlcoholStatus","age")]
)


model_data_test$event <- factor(
  ifelse(model_data_test$event == 1, "Yes", "No"), 
  levels = c("Yes", "No")
)

model_data_test = model_data_test[complete.cases(model_data_test),]
table(model_data_test$event)
formula_str <- paste("event ~", paste(protein_map$olink_id, collapse = " + "),"+factor(sex)+BMI+factor(SmokingStatus)+factor(AlcoholStatus)+age")
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

library(cvAUC)
ci.cvAUC(cv_model$pred$Yes,cv_model$pred$obs,folds = cv_model$pred$Resample)


cvAUC(cv_model$pred$Yes,cv_model$pred$obs,folds = cv_model$pred$Resample)


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
proteins_str <- paste(protein_map$olink_id, collapse = " + ")
proteins_str_model1 <- paste(protein_map_model1$olink_id, collapse = " + ")
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

# Train Model B (Proteins in model 2 Only)
cv_model_protein <- train(formula_protein, data = model_data_test, method = "glm", family = "binomial", trControl = fixed_cv_control, metric = "ROC")
ci.cvAUC(cv_model_protein$pred$Yes,cv_model_protein$pred$obs,folds = cv_model_protein$pred$Resample)

# Train Model C (Full Model)
cv_model_protein_rf <- train(formula_protein_rf, data = model_data_test, method = "glm", family = "binomial", trControl = fixed_cv_control, metric = "ROC")
ci.cvAUC(cv_model_protein_rf$pred$Yes,cv_model_protein_rf$pred$obs,folds = cv_model_protein_rf$pred$Resample)


cv_model_protein$results$ROC-cv_model_rf$results$ROC #0.03518703
cv_model_protein_rf$results$ROC-cv_model_rf$results$ROC #0.07270376




#--------------------------------------------------------------------------------------------------------
#Report the weights associated with the risk factors for the final risk prediction model of liver cancer
library(survival)
library(caret)
#Create a copy of the dataset to standardize
standardised_data <- model_data_test

for (protein in protein_map$olink_id) {
  p_mean <- mean(standardised_data[[protein]], na.rm = TRUE)
  p_sd   <- sd(standardised_data[[protein]], na.rm = TRUE)
  standardised_data[[protein]] <- (standardised_data[[protein]] - p_mean) / p_sd
}

standardised_data$sex           <- factor(standardised_data$sex)
standardised_data$SmokingStatus <- factor(standardised_data$SmokingStatus)
standardised_data$AlcoholStatus <- factor(standardised_data$AlcoholStatus)

print(levels(standardised_data$SmokingStatus))
print(levels(standardised_data$AlcoholStatus))

proteins_part <- paste(protein_map$olink_id, collapse = " + ")
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

