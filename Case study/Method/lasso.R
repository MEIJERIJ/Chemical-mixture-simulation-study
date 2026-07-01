################################################################################
# 
# Case study: Association between PFAS mixtures and BMI z-score
# 
# Data used from the 3M hotspot study
# 
################################################################################


################################################################################
#                                                                              #
#------------------------------- Lasso regression -----------------------------#
#                                                                              #
################################################################################


### Load packages
#----------------
library(glmnet)


lasso <- function(obs, n_boot = 200, qlow, qhigh) {
  
  
  ### Prepare data
  #---------------
  X <- model.matrix(Y ~ age2 + age3 + gender + income2 + income3 + 
       bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos, data = obs)[, -1]
  Y <- obs$Y
  
  
  ### Penalty factor: 0 for covariates (always in), 1 for exposures (penalized)
  #----------------------------------------------------------------------------
  penalty_factors <- c(rep(0, 5), rep(1, 7))
  exposure_names  <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", "pfos")
  iqr_vec <- qhigh - qlow
  
  
  ### Fit Lasso & Visual Diagnostics
  #---------------------------------
  cv_fit <- cv.glmnet(X, Y, alpha = 1, penalty.factor = penalty_factors, nfolds = 5)
  plot(cv_fit)
  title("Lasso CV (MSE Path)", line = 3)
  
  
  ### Fit final model with lambda.1se (more parsimonious)
  #------------------------------------------------------
  lasso_fit <- glmnet(X, Y, alpha = 1, penalty.factor = penalty_factors, 
                      lambda = cv_fit$lambda.1se)
  
  
  ### Extract point estimates
  #--------------------------
  coef_lasso <- as.vector(coef(lasso_fit))
  coef_full <- coef_lasso[7:13] 
  joint_full <- sum(coef_full * iqr_vec)
  
  
  ### Bootstrap for Confidence Intervals & Selection Frequency (PIP)
  #-----------------------------------------------------------------
  boot_indiv <- matrix(nrow = n_boot, ncol = 7)
  
  for(i in 1:n_boot) {
    idx <- sample(1:nrow(obs), replace = TRUE)
    X_b <- X[idx, ]; Y_b <- Y[idx]
    
    # Fit inside bootstrap
    cv_b <- cv.glmnet(X_b, Y_b, alpha = 1, penalty.factor = penalty_factors, nfolds = 5)
    fit_b <- glmnet(X_b, Y_b, alpha = 1, penalty.factor = penalty_factors, lambda = cv_b$lambda.1se)
    
    b_coefs <- as.vector(coef(fit_b))[7:13]
    boot_indiv[i, ] <- b_coefs
  }
  
  
  ### Summarize Results
  #--------------------
  pip      <- colMeans(boot_indiv != 0)

  
  
  ### Prepare Output
  #-----------------
  res <- data.frame(
    parameter = c("Joint", exposure_names),
    Estimate  = c(joint_full, coef_full),
    lower     = NA,
    upper     = NA,
    Method    = "Lasso regression",
    PIP       = c(NA, pip)
  )
  
  
  par(mfrow = c(1, 1))
  return(res)
}
