################################################################################
#
# Function to store all the results from Bayesian kernel machine regression
#
################################################################################


### Load libraries 
#-----------------
library(bkmr)


### Function to store BKMR results
#---------------------------------
run_bkmr <- function(df, qlow , qhigh, sim_id = 1, psi, joint_eff){
  ### Input:
  # - df:         input dataset
  # - qlow:       vector containing all values for the 25th quantile
  # - qhigh:      vector containing all values for the 75th quantile
  # - sim_id:     current seed to generate df
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Prepare data
  #---------------
  Exp <- df[,2:8]
  y <- df$Y
  Cov <- df[,9:12]
  
  
  ### Fit BKMR
  #-----------
  fit_bkmr <- kmbayes(y = y, Z = Exp, X = Cov, iter = 10000, verbose = TRUE, varsel = TRUE)
  
  
  ### Standardization
  #------------------
  mean_low <- SamplePred(fit_bkmr, Xnew = c(0,0,0,0), Znew = q_low, method="Exact")
  
  mean_high <- SamplePred(fit_bkmr, Xnew = c(0,0,0,0), Znew = q_high, method="Exact")
  
  
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
  
  
  ### Store results together
  #-------------------------
  coef_bkmr <- data.frame(
    Estimate = c(joint , rep(NA,7)),
    Sim = as.factor(sim_id),
    parameter = c("Joint", attributes(fit_bkmr$Z)$dimnames[[2]]),
    Method = paste0("BKMR"),
    Size = nrow(df),
    Power = c(power , rep(NA,7)),
    CI_coverage = c(cov , rep(NA,7)),
    CI_width = c(width , rep(NA,7)),
    PIP = c(NA, ExtractPIPs(fit_bkmr)[,2])
  )
  
  
  return(coef_bkmr)
}
