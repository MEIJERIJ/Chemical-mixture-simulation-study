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
#--------------------- Bayesian kernel machine regression ---------------------#
#
################################################################################


### Load libraries 
#-----------------
library(bkmr)


### Function to store BKMR results
#--------------------------------
bkmr <- function(df, qlow , qhigh){
 
  
  ### Prepare data
  #---------------
  Exp <- df[,2:8]
  y <- df$Y
  Cov <- df[,9:13]
  
  
  ### Fit BKMR
  #-----------
  fit_bkmr <- kmbayes(y = y, Z = Exp, X = Cov, iter = 10000, verbose = TRUE, varsel = TRUE)
  
  
  ### Standardization
  #------------------
  mean_low <- SamplePred(fit_bkmr, Xnew = c(0,0,0,0,0), Znew = q_low, method="Exact")
  
  mean_high <- SamplePred(fit_bkmr, Xnew = c(0,0,0,0,0), Znew = q_high, method="Exact")
  
  
  ### Difference per simulation
  #----------------------------
  posterior_diff <- mean_high - mean_low
  
  
  ### 95% credible interval
  #------------------------
  ci <- quantile(posterior_diff, probs = c(0.025, 0.975))

  
  ### Posterior mean difference
  #----------------------------
  joint <- mean(posterior_diff)
  
  
  ### Store results together
  #-------------------------
  res <- data.frame(
    Estimate = c(joint , rep(NA,7)),
    parameter = c("Joint", attributes(fit_bkmr$Z)$dimnames[[2]]),
    Method = paste0("BKMR"),
    lower     = c(ci[1], rep(NA,7)),
    upper     = c(ci[2], rep(NA,7)),
    PIP = c(NA, ExtractPIPs(fit_bkmr)[,2])
  )
  
  
  return(res)
}
