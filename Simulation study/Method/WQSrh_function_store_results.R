################################################################################
#
# Function to store all the results from WQS repeated holdout regression
#
# 10/01/2025 - Jonas Meijerink
#
# Note that this code only works for gWQS 3.0.5
#
################################################################################


### Load libraries 
#-----------------
library(gWQS)


### Function to store output
#---------------------------
run_WQS_rh <- function(obs, sim_id = 1, q = NULL){
  ### Input:
  # - obs:      input dataset
  # - sim_id:   seed to generate obs
  # - q:        Indicating whether to use quartiles (4), deciles (10) or standardized continuous exposures (NULL) 
  
  
  ### Fit a WQS regression model
  #-----------------------------
  resultwqs <- tryCatch({
    suppressWarnings(
      resultwqs <- gwqs(
        Y ~ wqs + age + gender + Educ2 + Educ3,
        mix_name = c("bpfos" , "pfba" , "pfda" , "lbpfhxs" , "pfna" , "lbpfoa" , "pfos"),
        data = obs,
        q = q,
        rh = 100,
        b = 50,
        b1_pos = FALSE,
        family = gaussian,
        signal = "one",
      )
    )
    resultwqs
    
  }, error = function(e) {
    message("WQS failed on iteration ", sim_id, ": ", e$message)
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
  
  
  ### For printing the method name
  #-------------------------------
  if(is.null(q)){
    q = "standardised continuous"
  }
  else if(q == 4){
    q = "quartiles"
  }
  else if (q == 10){
    q = "deciles"
  }
  else{
    q = q
  }
  
  
  ### Store and return as a dataframe
  #----------------------------------
  coef_WQS <- data.frame(
    Estimate = NA,
    Sim = as.factor(sim_id),
    parameter = c("Joint" , resultwqs$final_weights$mix_name),
    Method = paste0("WQS regression (",q,")"),
    Size = nrow(obs),
    Power = c(empirical_power , rep(NA,7)),
    CI_coverage = NA,
    CI_width = NA,
    PIP = c(NA , resultwqs$final_weights$Estimate)
  )
  
  
  return(coef_WQS)
}
