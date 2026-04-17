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
#----------------------------- Horseshoe regression ---------------------------#
#                                                                              #
################################################################################


### Load packages
#----------------
library(brms)


horseshoe <- function(df, qlow, qhigh) {
  
  
  ### Formula
  #----------
  formula <- bf(
    Y ~ alpha + conf + pfas,      
    alpha ~ 1,                         
    conf  ~ 0 + age2 + age3 + gender + income2 + income3,                             
    pfas  ~ 0 + bpfos + lbpfhxs + lbpfoa + pfba + pfda + pfna + pfos,
    nl = TRUE
  )
  
  
  ### Priors
  #---------
  prior <- c(
    set_prior("normal(0, 10)", class = "b", nlpar = "alpha"),
    set_prior("normal(0, 10)", class = "b", nlpar = "conf"),
    set_prior("horseshoe(df = 1)", class = "b", nlpar = "pfas"),
    set_prior("student_t(3, 0, 2.5)", class = "sigma", lb = 0)
  )
  
  
  ### Fit the Model
  #----------------
  fit <- brm(
    formula = formula,
    data    = df,
    family  = gaussian(),
    prior   = prior,
    seed    = 123,
    chains  = 2, 
    cores   = 2, 
    iter    = 5000,
    warmup  = 2000,
    control = list(max_treedepth = 15, adapt_delta = 0.995)
  )
  
  
  ### MCMC trace plots
  #-------------------
  plot(fit, pars = "^b_pfas")
  
  
  ### Extract Posteriors
  #---------------------
  draws <- posterior_samples(fit)
  pfas_cols <- grep("^b_pfas_", colnames(draws), value = TRUE)
  pfas_draws <- draws[, pfas_cols]
  
  
  ### Calculate joint effect posterior
  #-----------------------------------
  jointpost <- (q_high-q_low) %*% t(posterior_samples(fit, pars = "^b")[pfas_cols])
  
  
  ### Summary for individual pollutants
  #------------------------------------
  summ <- summary(fit)$fixed
  pfas_summ <- summ[grep("^pfas_", rownames(summ)), ]
  
  
  ### Prepare output
  #-----------------
  exposure_names <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", "pfos")
  
  
  res <- data.frame(
    parameter = c("Joint", exposure_names),
    Estimate  = c(mean(jointpost), pfas_summ$Estimate),
    lower     = c(quantile(jointpost, 0.025), pfas_summ$`l-95% CI`),
    upper     = c(quantile(jointpost, 0.975), pfas_summ$`u-95% CI`),
    Method    = "Horseshoe regression",
    PIP       = NA 
  )
  
  return(res)
}
