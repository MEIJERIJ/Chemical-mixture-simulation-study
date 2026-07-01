################################################################################
# 
# Function to store all the results from Elastic Net regression
#
# 07/01/2026 - Jonas Meijerink
################################################################################


### Load libraries 
#-----------------
library(glmnet)


### Function to store BAS results
#--------------------------------
run_Enet <- function(obs, n_boot = 200, sim_id = 1, qlow, qhigh,  psi, joint_eff){
  ### Input:
  # - obs:        input dataset
  # - n_boot:     bootstrap steps used to construct percentile based confidence intervals
  # - sim_id:     current seed to generate obs
  # - qlow:       vector containing all values for the 25th quantile
  # - qhigh:      vector containing all values for the 75th quantile
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Prepare data
  #---------------
  X <- model.matrix(Y ~ age + gender + Educ2 + Educ3 + bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos, data = obs)[, -1]
  Y <- obs$Y
  penalty_factors <- c(rep(0, 4), rep(1, 7))
  
  
  ### Fit Enet on full data
  #-------------------------
  cv_fit <- cv.glmnet(X, Y, alpha = 0.5, penalty.factor = penalty_factors,
                      type.measure = "mse", nfolds = 5)
  enet_fit <- glmnet(X, Y, alpha = 0.5, penalty.factor = penalty_factors,
                      lambda = cv_fit$lambda.1se)
  
  
  ### Extract coefficients
  #-----------------------
  coef_enet <- as.vector(coef(enet_fit))
  coef_full <- coef_enet[6:12]
  
  
  ### Precompute IQR (Q3 - Q1) per exposure
  #----------------------------------------
  iqr_vec <- qhigh - qlow
  
  
  ### Save joint effect
  #--------------------
  joint_full <- coef_full %*% iqr_vec
  
  
  ### Prepare bootstrap results
  #----------------------------
  chemical_selected <- matrix(FALSE, nrow = n_boot, ncol = 7)
  
  
  ### Run bootstrap sampling
  #-------------------------
  for(i in 1:n_boot){
    
    
    # Bootstrap sample
    boot_indices <- sample(1:nrow(obs), replace = TRUE)
    X_boot <- X[boot_indices, ]
    Y_boot <- Y[boot_indices]
    
    
    # Fit Enet regression using cross-validated lambda
    cv_fit <- cv.glmnet(X_boot, Y_boot, alpha = 0.5, penalty.factor = penalty_factors,
                        type.measure = "mse", nfolds = 5)
    enet_fit <- glmnet(X_boot, Y_boot, alpha = 0.5, penalty.factor = penalty_factors,
                        lambda = cv_fit$lambda.1se)
    
    
    # Extract coefficients
    exposure_estimates <- coef_enet[6:12]
    
    
    # Mark selected if coef != 0
    chemical_selected[i, ] <- exposure_estimates != 0
  }
  
  
  ### Parameter names
  #------------------
  chemicals <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda",  "pfna", "pfos")
  
  
  ### Placeholders for chemical results
  #------------------------------------
  selection_freqs <- numeric(7)
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    # Calculate proportion of times selected (nonzero coef)
    selection_freqs[i] <- mean(chemical_selected[, i])
  }
  
  
  ### Store all results
  #--------------------
  parameters <- c("Joint", chemicals)
  estimates <- c(joint_full, coef_full)
  
  
  ### Store results together
  #-------------------------
  coef_enet <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    parameter = parameters,
    Method = paste0("Elastic Net regression"),
    Size = nrow(obs),
    Power = NA,
    CI_coverage = NA,
    CI_width = NA,
    PIP = c(NA, selection_freqs)
  )
  
  
  return(coef_enet)
}
