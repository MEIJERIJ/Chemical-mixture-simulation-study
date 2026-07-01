################################################################################
#
# Function to store all the results from DART
#
################################################################################


### Load libraries
#-----------------
library(BART)


### Function to run DART
#-----------------------
run_dart <- function(df, q_low , q_high, sim_id = 1, psi, joint_eff) {
  ### Input:
  # - df:         input dataset
  # - q_low:      vector containing all values for the 25th quantile
  # - q_high:     vector containing all values for the 75th quantile
  # - sim_id:     current seed to generate df
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Fit a DART model
  #-------------------
  fit <- gbart(df[-1],df$Y, ntree = 200, nskip = 5000, ndpost = 5000, sparse = TRUE)
  
  
  ### Standardization
  #------------------
  covar_quantiles_low <- df
  covar_quantiles_low[,2:8] <- rep(q_low, each = nrow(df))
  mean_low <- rowMeans(predict(fit, covar_quantiles_low[-1]))
  
  covar_quantiles_high <- df
  covar_quantiles_high[,2:8] <- rep(q_high, each = nrow(df))
  mean_high <- rowMeans(predict(fit, covar_quantiles_high[-1]))
  
  
  ### Difference per simulation
  #----------------------------
  posterior_diff <- mean_high - mean_low
  
  
  ### 95% credible interval
  #------------------------
  ci <- quantile(posterior_diff, probs = c(0.025, 0.975))
  width <- ci[2] - ci[1]
  
  
  ### Calculate coverage & power
  #-----------------------------
  cov <- (ci[1] <= joint_eff * (q_high-q_low) %*% psi) & (ci[2] >=  joint_eff * (q_high-q_low) %*% psi)
  power <- (ci[1] > 0 | ci[2] < 0 )
  
  
  ### Posterior mean difference
  #----------------------------
  joint <- mean(posterior_diff)
  
  
  ### Format output
  #----------------
  coef_dart <- data.frame(
    Estimate = c(joint, rep(NA, 7)),
    Sim = as.factor(sim_id),
    parameter = c("Joint", names(fit$varprob.mean)[1:7]),
    Method = paste0("DART"),
    Size = nrow(df),
    Power = c(power, rep(NA, 7)),
    CI_coverage = c(cov, rep(NA, 7)),
    CI_width = c(width, rep(NA, 7)),
    PIP = c(NA , fit$varcount.mean[1:7])
  )
  
  
  return(coef_dart)
}
