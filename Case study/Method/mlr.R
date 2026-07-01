################################################################################
# 
# Case study: Association between PFAS mixtures and BMI z-score
# 
# Data used from the 3M hotspot study
# 
# 16/03/2026
################################################################################


mlr <- function(df, qlow, qhigh) {
  
  
  ### Fit multiple pollutant linear model
  #--------------------------------------
  exposure_names <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", "pfos")
  formula_str <- paste("Y ~ age2 + age3 + gender + income2 + income3 +", 
                       paste(exposure_names, collapse = " + "))
  linear_model <- lm(as.formula(formula_str), data = df)
  
  
  ### Visual diagnostics
  #---------------------
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  plot(linear_model, which = 1:3, main = "")
  for (exp in exposure_names) {
    par(mfrow = c(1, 1))
    plot(df[[exp]], residuals(linear_model), 
         xlab = exp, ylab = "Residuals", 
         main = paste("Linearity Check:", exp))
    abline(h = 0, col = "blue", lty = 2)
    lines(lowess(df[[exp]], residuals(linear_model)), col = "red")
  }
  par(mfrow = c(1, 1))
  
  
  ### Compute Joint Effect
  #-----------------------
  beta_hat <- coef(linear_model)
  vc <- vcov(linear_model)
  step <- qhigh - qlow
  
  
  # Contrast vector for exposures only
  dbar <- numeric(length(beta_hat))
  names(dbar) <- names(beta_hat)
  dbar[exposure_names] <- step
  
  
  joint_est <- sum(dbar * beta_hat)
  se_joint  <- sqrt(as.numeric(t(dbar) %*% vc %*% dbar))
  t_crit    <- qt(0.975, df.residual(linear_model))
  
  
  joint_lower <- joint_est - t_crit * se_joint
  joint_upper <- joint_est + t_crit * se_joint
  
  
  ### Extract Chemical-Specific results
  #------------------------------------
  chemical_estimates <- beta_hat[exposure_names]
  ci_indiv <- confint(linear_model)[exposure_names, ]
  
  
  ### Prepare output
  #-----------------
  res <- data.frame(
    parameter   = c("Joint", exposure_names),
    Estimate    = c(joint_est, as.numeric(chemical_estimates)),
    lower       = c(joint_lower, ci_indiv[, 1]),
    upper       = c(joint_upper, ci_indiv[, 2]),
    Method      = "Multiple pollutant linear model",
    PIP         = NA
  )
  
  return(res)
}
