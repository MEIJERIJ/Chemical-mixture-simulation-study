################################################################################
# Simulation study 
#
# Generate data set  with additive linear dense effect
#
################################################################################


rm(list = ls())


### When using R locally
#-----------------------
setwd("...")


### Load libraries 
#-----------------
library(readxl)


### Load data for a realistic simulation study
#---------------------------------------------
df <- read_excel("~/Dataset_HBM_3M.xlsx")


### Make gender a dummy variable (0/1)
#-------------------------------------
df$GESL <- (df$GESL-1)


### Create data frame
#--------------------
Data <- data.frame("bpfos" = scale(log(df$bpfos_imp)) , "lbpfhxs" = scale(log(df$lbpfhxs_imp)), 
                   "lbpfoa"=scale(log(df$lbpfoa_imp)),
                   "pfba" = scale(log(df$pfba_imp)), "pfda" = scale(log(df$pfda_imp)) , 
                   "pfna" = scale(log(df$pfna_imp)) , "pfos" = scale(log(df$pfos_imp)),
                   "age"=df$OB_Leeftijd, "highest education"=df$ISCED_HH , 
                   "gender"=df$GESL)


### Create dummy variables for education
#---------------------------------------
Data$Educ2 <- ifelse(Data$highest.education==2,1,0)
Data$Educ3 <- ifelse(Data$highest.education==3,1,0)
Data$highest.education <- NULL


### Select continuous variables to simulate
#------------------------------------------
PFAS_clean <- na.omit(Data)


#------------------------------------------------------------------------------#
#
# Function to generate a dataset given a seed based on the 3M hotspot study
#
#------------------------------------------------------------------------------#


### Write a parallel simulation function
#---------------------------------------
Generate_data <-  function(seed = 1234 , sigma = 0.1, psi, joint_eff = -1){
  
  
  ### Define the effect sizes for the confounders
  #----------------------------------------------
  conf_eff <- c(-0.2, 0.3, 0.3 , -0.05)
  
  
  ### Define the outcome
  #---------------------
  signal <- joint_eff * psi %*% t(PFAS_clean[,1:7]) + conf_eff %*% t(PFAS_clean[,8:11])
  
  
  ### Add noise the outcome
  #------------------------
  set.seed(seed)
  eps <- rnorm(300,0, sd = sigma)
  Y <- as.vector(signal) + eps
  
  
  ### Compute actual SNR
  #---------------------
  SNR_actual <- var(as.vector(signal)) / sigma^2
  
  
  ### Return list with data frame and SNR
  #--------------------------------------
  return(list(
    simdata = data.frame(Y = Y, PFAS_clean),
    snr = SNR_actual
  ))
}
