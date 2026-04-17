################################################################################
#
# Function to store all the results from SoftBART
#
################################################################################


### Load libraries
#-----------------
library(SoftBart)


### Function to do standardization with SoftBART
#-----------------------------------------------
run_softbart <- function(df, qlow , qhigh, sim_id = 1, psi, joint_eff) {
  ### Input:
  # - df:         input dataset
  # - qlow:       vector containing all values for the 25th quantile
  # - qhigh:      vector containing all values for the 75th quantile
  # - sim_id:     current seed to generate df
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Fit a SoftBART model
  #-----------------------
  fit <- softbart_regression(Y ~ ., 
                             data = df, 
                             test_data = df,
                             hypers = Hypers(X = df[-1], Y = df$Y ,num_tree = 200),
                             opts = Opts(num_burn = 5000, 
                                         num_save = 5000))
  
  
  ### Extract the PIP
  #------------------
  pip <- posterior_probs(fit) 
  
  
  ### Standardization
  #------------------
  covar_quantiles_low <- df
  covar_quantiles_low[,2:8] <- rep(q_low, each = nrow(df))
  low <- predict.softbart_regression(fit, 
                              newdata = covar_quantiles_low)$mu_mean
  mean_low <- rowMeans(predict.softbart_regression(fit, 
                                         newdata = covar_quantiles_low)$mu)
  covar_quantiles_high <- df
  covar_quantiles_high[,2:8] <- rep(q_high, each = nrow(df))
  high <- predict.softbart_regression(fit, 
                                     newdata = covar_quantiles_high)$mu_mean
  mean_high <- rowMeans(predict.softbart_regression(fit, 
                                                   newdata = covar_quantiles_high)$mu)
  
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
  coef_softbart <- data.frame(
    Estimate = c(joint , rep(NA,7)),
    Sim = as.factor(sim_id),
    parameter = c("Joint", names(pip$post_probs)[1:7]),
    Method = paste0("SoftBART"),
    Size = nrow(df),
    Power = c(power , rep(NA,7)),
    CI_coverage = c(cov , rep(NA,7)),
    CI_width = c(width , rep(NA,7)),
    PIP = c(NA, pip$post_probs[1:7])
  )
  
  
  return(coef_softbart)
}
