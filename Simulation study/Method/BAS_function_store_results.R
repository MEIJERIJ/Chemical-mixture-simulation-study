################################################################################
#
# Function to store all the results from Bayesian additive sampling
#
# 07/01/2026 - Jonas Meijerink
################################################################################


### Load libraries 
#-----------------
library(BAS)


### Function to store BAS results
#--------------------------------
run_bas <- function(obs, sim_id = 1){
  ### Input:
  # - obs:        input dataset
  # - sim_id:     current seed to generate obs
  
  
  ### Fit BAS model
  #----------------
  bas_fit <- bas.lm(
    Y ~ bpfos + pfba + pfda + lbpfhxs + pfna + lbpfoa + pfos, 
    include.always = as.formula(~age + gender + Educ2 + Educ3),
    data  = obs,
    method = "MCMC+BAS", 
    n.models = 2^10, 
    MCMC.iterations = 50000, 
    thin = 10, 
    prior = "JZS",
    modelprior = uniform()
  )
  
  
  ### Extract model-averaged coefficients
  #--------------------------------------
  coef_summary <- coef(bas_fit, estimator = "BMA")
  
  
  ### Store and return as a dataframe
  #----------------------------------
  coef_bayes <- data.frame(
    
    Estimate = NA,
    
    Sim = as.factor(sim_id),
    
    parameter = bas_fit$namesx[2:8],
    
    Method = paste0("Bayesian model averaging"),
    
    Size = nrow(obs),
    
    Power = NA,
    
    CI_coverage = NA,
    
    
    CI_width = NA, 
    
    PIP = c(coef_summary$probne0[-1])
  )
  
  return(coef_bayes)
}
