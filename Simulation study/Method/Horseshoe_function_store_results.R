################################################################################
#
# Store results for Bayesian Horseshoe regression
# 
################################################################################


### Load libraries
#-----------------
library(brms)


### Function to do Bayesian Horseshoe regression
#-----------------------------------------------
run_horseshoe <- function(df, sim_id, psi, joint_eff, q_low, q_high) {
  ### Input:
  # - df:         input dataset
  # - sim_id:     current seed to generate df
  # - psi:        vector containing individual effect sizes
  # - joint_eff:  overall effect of the weighted sum
  # - q_low:      vector containing all values for the 25th quantile
  # - q_high:     vector containing all values for the 75th quantile
  
  
  ### Formula
  #----------
  formula <- bf(
    Y ~ alpha + conf + pfas,      
    alpha ~ 1,                         
    conf ~ 0 + age + gender + Educ2 + Educ3,                             
    pfas ~ 0 + bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos,
    nl = TRUE
  )
  
  
  ### Define prior for all parameters
  #----------------------------------
  prior <- c(
    
    
    ### Normal prior on intercept
    #----------------------------
    set_prior("normal(0, 10)", class = "b", nlpar = "alpha"),
    
    
    ### Normal prior on confounders
    #------------------------------
    set_prior("normal(0, 10)", class = "b", nlpar = "conf"),
    
    
    ### Horseshoe prior on PFAS
    #--------------------------
    set_prior("horseshoe(df = 1)", class = "b", nlpar = "pfas"),
    
    
    ### Prior on the variance
    #------------------------
    set_prior("student_t(3, 0, 2.5)", class = "sigma", lb = 0)
  )

  
  ### Fit the model
  #----------------
  fit <- brm(
    formula = formula,
    data    = df,
    family  = gaussian(),
    prior   = prior,
    seed    = 123,
    chains  = 4, 
    cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4")), 
    iter = 5000,
    warmup = 2000,
    control = list(max_treedepth = 15, adapt_delta = 0.99)
  )
  results <- summary(fit)$fixed
  
  
  ### Get joint effect
  #-------------------
  posterior <- posterior_samples(fit, pars = "^b")[6:12]
  jointpost <- (q_high-q_low) %*% t(posterior_samples(fit, pars = "^b")[6:12])
  joint <- data.frame(
    "Estimate" = mean(jointpost),
    "CI.Upper" = quantile(jointpost, probs = 0.975),
    "CI.Lower" = quantile(jointpost, probs = 0.025)
  )
  
  
  ### Calculate coverage
  #---------------------
  cov_joint <- (joint$CI.Lower <= joint_eff * (q_high-q_low) %*% psi) & (joint$CI.Upper >=  joint_eff * (q_high-q_low) %*% psi)
  cov_pfas <- (results$`l-95% CI`[6:12] <= joint_eff * psi) & (results$`u-95% CI`[6:12] >=  joint_eff * psi)
  
  
  ### Format output
  #----------------
  coef_horseshoe <- data.frame(
    Estimate = c(joint$Estimate ,  results$Estimate[6:12]),
    
    Sim = as.factor(sim_id),
    
    parameter = c("Joint", sub(".*_", "", rownames(results[6:12,]))),
    
    Method = paste0("Horseshoe regression"),
    
    Size = nrow(df),
    
    Power = c(joint$CI.Lower > 0 | joint$CI.Upper < 0  , 
              results$`l-95% CI`[6:12] > 0 | results$`u-95% CI`[6:12] < 0 ),
    
    CI_coverage = c(cov_joint , cov_pfas),
    
    CI_width = c(joint$CI.Upper-joint$CI.Lower,
                 results$`u-95% CI`[6:12]-results$`l-95% CI`[6:12]),
    
    PIP = NA
    
  )
  
  
  return(coef_horseshoe)
}
