################################################################################
# Simulation study Master Thesis - Jonas Meijerink
#
# Function to store all the results from OLS regression model
#
################################################################################


run_linear_model <- function(df, qlow , qhigh, sim_id = 1, psi, joint_eff) {
  ### Input:
  # - df:         input dataset
  # - qlow:       vector containing all values for the 25th quantile
  # - qhigh:      vector containing all values for the 75th quantile
  # - sim_id:     current seed to generate df
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Fit linear model
  #-------------------
  linear_model <- lm(Y ~ age + gender + Educ2 + Educ3 + bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos, data = df)
  
  
  ### Extract beta and covariance matrix
  #-------------------------------------
  beta_hat <- coef(linear_model)
  vc <- vcov(linear_model)
  
  
  ### Exposure names
  #-----------------
  exposure_names <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda",  "pfna", "pfos")
  
  
  ### Compute Q25 and Q75 differences for PFAS
  #-------------------------------------------
  step <- qhigh - qlow
  dbar <- numeric(length(beta_hat))
  names(dbar) <- names(beta_hat)
  dbar[exposure_names] <- step
  
  
  ### Calculate variance of weighted sum
  #-------------------------------------
  joint <- sum(dbar * beta_hat)
  se <- sqrt(as.numeric(t(dbar) %*% vc %*% dbar))
  df_resid <- df.residual(linear_model)
  t_crit <- qt(0.975, df_resid)          
  jointci <- c(joint - t_crit * se, joint + t_crit * se)
  
  
  ### Placeholders for chemical results
  #------------------------------------
  chemical_estimates <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    
    chemical_estimates[i] <- coef(linear_model)[exposure_names[i]]
    
    ci <- confint(linear_model)[exposure_names[i],]
    ci_widths[i] <- ci[2] - ci[1]
    empirical_powers[i] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[i] <- (ci[1] <= joint_eff * psi[i]) & (ci[2] >= joint_eff * psi[i])
  }
  
  
  ### Store all results
  #--------------------
  parameters <- c("Joint", exposure_names)
  estimates <- c(joint, chemical_estimates)
  ci_widths_all <- c(jointci[2] - jointci[1], ci_widths)
  ci_power_all <- c((jointci[1] > 0 | jointci[2] < 0), 
                    empirical_powers)
  ci_coverage_all <- c( (jointci[1] <= joint_eff * (qhigh-qlow) %*% psi) & (jointci[2] >=  joint_eff * (qhigh-qlow) %*% psi) , 
                        CI_coverages)
  
  
  ### Prepare output
  #-----------------
  coef_lm <- data.frame(
    Estimate = estimates,
    Sim = as.factor(sim_id),
    Size = nrow(obs),
    parameter = parameters,
    Method = "Multiple pollutant linear model",
    Power = ci_power_all,
    CI_coverage = ci_coverage_all,
    CI_width = ci_widths_all,
    PIP = rep(NA, 8)
  )
  
  
  return(coef_lm)
}
