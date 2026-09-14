library(data.table)
library(pROC)

#functions adapted from the cvAUC package
#cvAUC: Cross-Validated Area Under the ROC Curve Confidence Intervals

############################################################
## Influence-curve contributions for one fold
############################################################

get_auc_ic <- function(pred,
                       outcome,
                       pos,
                       neg,
                       w1,
                       w0) {
  
  outcome <- factor(outcome, levels = c(neg, pos))
  
  n <- length(outcome)
  n_pos <- sum(outcome == pos)
  n_neg <- sum(outcome == neg)
  
  if (n_pos == 0 || n_neg == 0) {
    stop("Each fold must contain both cases and controls.")
  }
  
  auc <- as.numeric(AUC(pred, outcome))
  
  DT <- data.table(
    id = seq_along(pred),
    pred = as.numeric(pred),
    label = outcome
  )
  
  ## Fraction of controls with lower predictions
  DT <- DT[order(pred, -xtfrm(label))]
  DT[, fracNegLabelsWithSmallerPreds :=
       cumsum(label == neg) / n_neg]
  
  ## Fraction of cases with higher predictions
  DT <- DT[order(-pred, label)]
  DT[, fracPosLabelsWithLargerPreds :=
       cumsum(label == pos) / n_pos]
  
  ## Influence contribution
  DT[, icVal :=
       ifelse(
         label == pos,
         w1 * (fracNegLabelsWithSmallerPreds - auc),
         w0 * (fracPosLabelsWithLargerPreds - auc)
       )]
  
  ## Return to original row order
  DT <- DT[order(id)]
  
  list(
    auc = auc,
    ic = DT$icVal
  )
}

############################################################
## Mean fold-specific paired delta-AUC
############################################################
delta_auc_ic_df <- function(data,
                            pred1_col,
                            pred2_col,
                            outcome_col,
                            fold_col,
                            confidence = 0.95) {
  
  data <- as.data.table(data)
  
  required_cols <- c(
    pred1_col,
    pred2_col,
    outcome_col,
    fold_col
  )
  
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  ## Keep complete observations
  data <- data[
    complete.cases(
      data[, ..required_cols]
    )
  ]
  
  ## Convert outcome to a two-level factor
  data[, outcome_tmp := get(outcome_col)]
  
  if (is.factor(data$outcome_tmp)) {
    if (nlevels(data$outcome_tmp) != 2) {
      stop("The outcome must have exactly two levels.")
    }
    
    neg <- levels(data$outcome_tmp)[1]
    pos <- levels(data$outcome_tmp)[2]
    
  } else {
    outcome_values <- sort(unique(data$outcome_tmp))
    
    if (length(outcome_values) != 2) {
      stop("The outcome must have exactly two unique values.")
    }
    
    neg <- as.character(outcome_values[1])
    pos <- as.character(outcome_values[2])
  }
  
  data[, outcome_tmp := factor(
    as.character(outcome_tmp),
    levels = c(neg, pos)
  )]
  
  ## Overall class weights
  n_obs <- nrow(data)
  n_pos <- sum(data$outcome_tmp == pos)
  n_neg <- sum(data$outcome_tmp == neg)
  
  if (n_pos == 0 || n_neg == 0) {
    stop("Both outcome classes must be present.")
  }
  
  w1 <- n_obs / n_pos
  w0 <- n_obs / n_neg
  
  ## Fold labels must uniquely identify each validation fold
  folds <- unique(data[[fold_col]])
  K <- length(folds)
  
  auc1_by_fold <- numeric(K)
  auc2_by_fold <- numeric(K)
  delta_by_fold <- numeric(K)
  ic_delta_by_fold <- vector("list", K)
  
  for (k in seq_along(folds)) {
    
    fold_id <- folds[k]
    
    fold_data <- data[
      data[[fold_col]] == fold_id
    ]
    
    outcome_k <- fold_data$outcome_tmp
    
    if (length(unique(outcome_k)) < 2) {
      stop(
        "Fold '", fold_id,
        "' does not contain both outcome classes."
      )
    }
    
    ic1 <- get_auc_ic(
      pred = fold_data[[pred1_col]],
      outcome = outcome_k,
      pos = pos,
      neg = neg,
      w1 = w1,
      w0 = w0
    )
    
    ic2 <- get_auc_ic(
      pred = fold_data[[pred2_col]],
      outcome = outcome_k,
      pos = pos,
      neg = neg,
      w1 = w1,
      w0 = w0
    )
    
    auc1_by_fold[k] <- ic1$auc
    auc2_by_fold[k] <- ic2$auc
    delta_by_fold[k] <- ic1$auc - ic2$auc
    
    ## Pair the two models within each individual
    ic_delta_by_fold[[k]] <- ic1$ic - ic2$ic
  }
  
  ## Mean fold-specific AUCs and delta-AUC
  auc1 <- mean(auc1_by_fold)
  auc2 <- mean(auc2_by_fold)
  delta_auc <- mean(delta_by_fold)
  
  ## Combine all fold-specific influence contributions
  ic_delta <- unlist(ic_delta_by_fold, use.names = FALSE)
  
  ## Center influence contributions
  ic_delta <- ic_delta - mean(ic_delta)
  
  ## Influence-curve variance and SE
  sigma_delta_sq <- mean(ic_delta^2)
  se_delta <- sqrt(sigma_delta_sq / n_obs)
  
  ## Wald confidence interval
  z <- qnorm(1 - (1 - confidence) / 2)
  
  ci <- c(
    delta_auc - z * se_delta,
    delta_auc + z * se_delta
  )
  
  ## Delta-AUC ranges from -1 to 1
  ci[1] <- max(ci[1], -1)
  ci[2] <- min(ci[2], 1)
  
  list(
    auc_model1 = auc1,
    auc_model2 = auc2,
    auc_model1_by_fold = auc1_by_fold,
    auc_model2_by_fold = auc2_by_fold,
    delta_auc_by_fold = delta_by_fold,
    delta_auc = delta_auc,
    se = se_delta,
    ci = ci,
    confidence = confidence,
    n_obs = n_obs,
    n_folds = K,
    positive_class = pos,
    negative_class = neg,
    influence_delta = ic_delta
  )
}


compare_rf_protein <- data.frame(
  rf = cv_model_rf$pred$Yes,
  rf_prot = cv_model_protein_rf$pred$Yes,
  obs = cv_model_rf$pred$obs,
  fold = cv_model_rf$pred$Resample
)

result12 <- delta_auc_ic_df(data = compare_rf_protein,
                            pred1_col = "rf_prot",
                            pred2_col = "rf",
                            outcome_col = "obs",
                            fold_col = "fold")

result12$ci
#[1] 0.03662457 0.10878296


compare_rf_protein <- data.frame(
  rf = cv_model_rf$pred$Yes,
  rf_prot = cv_model_protein_rf$pred$Yes,
  obs = cv_model_rf$pred$obs,
  fold = cv_model_rf$pred$Resample
)


compare_rf_protein <- data.frame(
  rf = cv_model_rf$pred$Yes,
  rf_prot = cv_model_protein$pred$Yes,
  obs = cv_model_rf$pred$obs,
  fold = cv_model_rf$pred$Resample
)
result13 <- delta_auc_ic_df(data = compare_rf_protein,
                            pred1_col = "rf_prot",
                            pred2_col = "rf",
                            outcome_col = "obs",
                            fold_col = "fold")
result13$ci
#[1] -0.005038331  0.075412383