################################################################################
#
# Function to store all the results from single OLS regression model
#
################################################################################


run_slr <- function(obs, sim_id = 1, psi, joint_eff) {
  ### Input:
  # - obs:        input dataset
  # - sim_id:     current seed to generate obs
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  
  
  ### Placeholders for chemical results
  #------------------------------------
  chemical_estimates <- numeric(7)
  ci_widths <- numeric(7)
  empirical_powers <- numeric(7)
  CI_coverages <- numeric(7)
  exposure_names <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda",  "pfna", "pfos")
  
  
  ### Calculate chemical-specific statistics
  #-----------------------------------------
  for (i in 1:7) {
    
    
    ### Fit single-pollutant linear model
    #------------------------------------
    formula <- as.formula(paste0("Y ~ age + gender + Educ2 + Educ3 + ", exposure_names[i]))
    linear_model <- lm(formula, data = obs)
    
    
    ### Select the estimate
    #----------------------
    chemical_estimates[i] <- coef(linear_model)[exposure_names[i]]
    
    
    ### Select the CI parameters
    #---------------------------
    ci <- confint(linear_model)[exposure_names[i],]
    ci_widths[i] <- ci[2] - ci[1]
    empirical_powers[i] <- (ci[1] > 0 | ci[2] < 0)
    CI_coverages[i] <- (ci[1] <= joint_eff * psi[i]) & (ci[2] >= joint_eff * psi[i])
  }
  
  
  ### Prepare output
  #-----------------
  coef_lm <- data.frame(
    Estimate = chemical_estimates,
    Sim = as.factor(sim_id),
    Size = nrow(obs),
    parameter = exposure_names,
    Method = "Single pollutant linear model",
    Power = empirical_powers,
    CI_coverage = CI_coverages,
    CI_width = ci_widths,
    PIP = rep(NA,7)
  )
  
  return(coef_lm)
}