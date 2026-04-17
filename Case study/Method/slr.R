################################################################################
# 
# Case study: Association between PFAS mixtures and BMI z-score
# 
# Data used from the 3M hotspot study
# 
# 16/03/2026
################################################################################


### Load package
#---------------
library(car)


slr <- function(exposure, data) {
  
  
  ### Define and fit the single linear model
  #-----------------------------------------
  form <- as.formula(paste0("Y ~ age2 + age3 + gender + income2 + income3 + ", exposure))
  fit <- lm(form, data = data)
  
  
  ### Model diagnostics
  #--------------------
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  plot(fit, which = 1) 
  plot(fit, which = 2)        
  plot(fit, which = 3)
  plot(data[[exposure]], residuals(fit), 
       xlab = exposure, ylab = "Residuals",
       main = paste("Residuals vs", exposure))
  abline(h = 0, col = "blue", lty = 2)
  lines(lowess(data[[exposure]], residuals(fit)), col = "red") 
  
  
  ### Extract model parameters
  #---------------------------
  res <- data.frame(
    Estimate = coef(fit)[exposure],
    parameter = exposure,
    Method = "Single pollutant linear model",
    lower = confint(fit)[exposure, 1],
    upper = confint(fit)[exposure, 2],
    PIP = NA
  )
  
  
  return(res)
}
