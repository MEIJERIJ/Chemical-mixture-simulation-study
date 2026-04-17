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
#--------------------------------- SoftBART -----------------------------------#
#
################################################################################


### Load libraries
#-----------------
library(SoftBart)


softbart <- function(df, q_low, q_high){
  
  
  ### With two dummy columns, SoftBART treats each column as independent
  # Therefore we switch to one categorical variable
  # A sensitivity showed that there are no differences between the two choices.
  df_clean         <- df[, c("Y",
                             "bpfos", "lbpfhxs", "lbpfoa",
                             "pfba",  "pfda",    "pfna",  "pfos",
                             "gender")]
  df_clean$age     <- ifelse(df$age3 == 1, 2,
                             ifelse(df$age2 == 1, 1, 0))
  df_clean$income  <- ifelse(df$income3 == 1, 2,
                             ifelse(df$income2 == 1, 1, 0))
  
  
  ### Fit a SoftBART model
  #-----------------------
  fit <- softbart_regression(Y ~ ., 
                             data = df_clean, 
                             test_data = df_clean,
                             hypers = Hypers(df_clean[-1],
                                             df_clean$Y ,
                                             num_tree = 200),
                             opts = Opts(num_burn = 5000, 
                                         num_save = 5000,
                                         update_s = TRUE)
                             )
  
  
  ### Extract the PIP
  #------------------
  pip <- posterior_probs(fit) 
  
  
  ### Standardization
  #------------------
  covar_quantiles_low <- df_clean
  covar_quantiles_low[,2:8] <- rep(q_low, each = nrow(df))
  low <- predict.softbart_regression(fit, 
                                     newdata = covar_quantiles_low)$mu_mean
  mean_low <- rowMeans(predict.softbart_regression(fit, 
                                                   newdata = covar_quantiles_low)$mu)
  covar_quantiles_high <- df_clean
  covar_quantiles_high[,2:8] <- rep(q_high, each = nrow(df))
  high <- predict.softbart_regression(fit, 
                                      newdata = covar_quantiles_high)$mu_mean
  mean_high <- rowMeans(predict.softbart_regression(fit, 
                                                    newdata = covar_quantiles_high)$mu)
  
  ### Difference per simulation
  #----------------------------
  posterior_diff <- mean_high - mean_low
  
  
  ### 95% credible interval
  #------------------------
  ci <- quantile(posterior_diff, probs = c(0.025, 0.975))
  
  
  ### Posterior mean difference
  #----------------------------
  joint <- mean(posterior_diff)
  
  
  ### Format output
  #----------------
  res <- data.frame(
    Estimate = c(joint , rep(NA,7)),
    parameter = c("Joint", names(pip$post_probs)[1:7]),
    Method = paste0("SoftBART (sparse)"),
    lower     = c(ci[1], rep(NA,7)),
    upper     = c(ci[2], rep(NA,7)),
    PIP = c(NA, pip$post_probs[1:7])
  )
  
  return(res)
}
