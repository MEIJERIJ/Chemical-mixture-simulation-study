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
#------------------- Repeated holdout WQS regression --------------------------#
#                                                                              #
################################################################################


### Load libraries 
#-----------------
library(gWQS)


### Function to store output
#---------------------------
rhwqs <- function(obs){
  
  
  ### Fit a WQS regression model
  #-----------------------------
  resultwqs <- tryCatch({
    suppressWarnings(
      resultwqs <- gwqs(
        Y ~ wqs + age2 + age3 + gender + income2 + income3,
        mix_name = c("bpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
        data = obs,
        q = 4,
        rh = 100,
        b = 100,
        b1_pos = FALSE,
        family = gaussian,
        signal = "one",
      )
    )
    resultwqs
    
  }, error = function(e) {
    message("WQS failed on iteration ", e$message)
    NULL
  })
  
  
  ### Immediately return if the model failed
  #-----------------------------------------
  if (is.null(resultwqs)) {
    return(NULL)
  }
  
  
  ### Joint effect power
  #---------------------
  coefs <- summary(resultwqs)$coefficients
  empirical_power <- as.numeric(coefs[2, "2.5 %"] > 0 | coefs[2, "97.5 %"] < 0)
  
  
  ### Store and return as a dataframe
  #----------------------------------
  res <- data.frame(
    Estimate = c(coefs[2, "Estimate"], rep(NA,7)),
    parameter = c("Joint" , resultwqs$final_weights$mix_name),
    Method = paste0("WQS regression (quartiles)"),
    upper = c(coefs[2, "97.5 %"], rep(NA,7)),
    lower = c(coefs[2, "2.5 %"], rep(NA,7)),
    PIP = c(NA , resultwqs$final_weights$Estimate)
  )
  
  
  return(res)
}
