################################################################################
# 
# Case study: Association between PFAS mixtures and BMI z-score
# 
# Data used from the 3M hotspot study
# 
# 16/03/2026
################################################################################


################################################################################
#                                                                              #
#--------------- Bayesian model averaging and adaptive sampling ---------------#
#
################################################################################


### Load libraries 
#-----------------
library(BAS)


### Function to store BAS results
#--------------------------------
bas <- function(obs){
  
  
  ### Fit BAS model
  #----------------
  bas_fit <- bas.lm(
    Y ~ bpfos + pfba + pfda + lbpfhxs + pfna + lbpfoa + pfos, 
    include.always = as.formula(~age2 + age3 + gender + income2 + income3),
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
  res <- data.frame(
    Estimate = NA,
    parameter = bas_fit$namesx[2:8],
    Method = paste0("Bayesian model averaging"),
    lower     = rep(NA,7),
    upper     = rep(NA,7),
    PIP = c(coef_summary$probne0[-1])
  )
  
  return(res)
}
