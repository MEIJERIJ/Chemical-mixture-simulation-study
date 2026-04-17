################################################################################
# Simulation study
#
# Function to store results form different methods in a single dataframe
#
################################################################################


rm(list = ls())


### Load functions that will be used
#-----------------------------------
source("~/Generate_data.R")

source("~/Method/SoftBART_function_store_results.R")
source("~/Method/SoftBARTsemi_function_store_results.R")
source("~/Method/SoftBART_withoutselection_function_store_results.R")
source("~/Method/BART_function_store_results.R")
source("~/Method/DART_function_store_results.R")
source("~/Method/BKMR_function_store_results.R")
source("~/Method/BKMR_withoutselection_function_store_results.R")
source("~/Method/OLS_function_store_results.R")
source("~/Method/Horseshoe_function_store_results.R")
source("~/Method/regHorseshoe_function_store_results.R")
source("~/Method/BAS_function_store_results.R")
source("~/Method/SLR_function_store_results.R")
source("~/Method/Ridge_function_store_results.R")
source("~/Method/Lasso_function_store_results.R")
source("~/Method/Enet_function_store_results.R")
source("~/Method/WQSrh_function_store_results.R")
source("~/Method/gsoftbart_regression_ortho.R")


### Set working directory to right folder
#----------------------------------------
setwd("~/Dense linear data-generating mechanism")


### Load libraries 
#-----------------
library(data.table)
library(RhpcBLASctl)
library(doParallel)
library(parallel)
library(SoftBart)
library(BART)
library(brms)
library(BAS)
library(glmnet)
library(gWQS)
library(caret)
library(progress)


### Define total number of simulations
#-------------------------------------
m <- 100
n_cores <- 4                 
iterations_per_batch <- n_cores
total_batches <- ceiling(m / iterations_per_batch)


### Define the seed
#------------------
seed <- 85721023 + 0:(m-1)


### Logging and checkpoint
#-------------------------
log_file <- "simulation_log.txt"
checkpoint_file <- "checkpoint_results.rds"


### If file exist check to continue on last batch
#------------------------------------------------
if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  start_batch <- checkpoint$last_batch + 1
  results <- checkpoint$results
  cat("Resuming from batch", start_batch, "\n", file = log_file, append = TRUE)
} else {
  start_batch <- 1
  results <- vector("list", total_batches)
  cat("Starting new simulation\n", file = log_file, append = TRUE)
}


### Create cluster
#-----------------
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)

clusterEvalQ(cluster, {
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  omp_set_num_threads(1)
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  NULL
})
clusterExport(cluster,
              varlist = c("Generate_data", "run_bart", "run_bas","run_dart",
                          "run_Enet", "run_horseshoe", "run_Lasso", "run_linear_model",
                          "run_reghorseshoe", "run_Ridge", "run_slr" , "run_softbart",
                          "run_softbart", "run_softbart_withoutselection", "run_softbart_semi",
                          "run_WQS_rh", "gsoftbart_regression_ortho"),
              envir = environment()
)


### Define your SNRs of interest
#-------------------------------
sigma <- c(0.58 , 0.8 , 1.1 , 2.5)


### Run simulation
#-----------------
timing <- system.time({
  
  for (batch in start_batch:total_batches) {
    
    cat(paste(Sys.time(), "Starting batch", batch, "\n"), file = log_file, append = TRUE)
    
    
    ### Compute seed subset for this batch safely
    #--------------------------------------------
    start_idx <- (batch - 1) * iterations_per_batch + 1
    end_idx <- min(batch * iterations_per_batch, length(seed))
    seed_subset <- seed[start_idx:end_idx]
    iterations_this_batch <- length(seed_subset)
    
    
    ### Store results for all SNRs
    #-----------------------------
    all_results <- list()
    
    
    for (s in sigma) {
      
      cat(paste("  Running SNR =", s, "\n"), file = log_file, append = TRUE)
      
      
      ### Run batch in parallel
      #------------------------
      batch_results <- tryCatch({
        foreach(i = 1:iterations_this_batch, 
                .combine = "list",
                .multicombine = TRUE,
                .packages = c("data.table","SoftBart", "bkmr","BART", 
                              "brms", "BAS", "glmnet", "gWQS", "caret", "progress")) %dopar% {
                  
                  list_results <- vector("list", 16)
                  
                  
                  ### Effect sizes
                  #---------------
                  psi <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4)
                  joint_eff <- -1
                  
                  
                  ### Generate the data with target sigma/SNR
                  #------------------------------------------
                  obs_list <- Generate_data(seed_subset[i], sigma = s, psi = psi, joint_eff = joint_eff)
                  obs <- obs_list$simdata
                  snr_value <- obs_list$snr
                  
                  
                  ### Define the quantiles for standardization
                  #-------------------------------------------
                  q_low  <- sapply(obs[2:8], function(x) quantile(as.numeric(x), 0.25))
                  q_high <- sapply(obs[2:8], function(x) quantile(as.numeric(x), 0.75))
                  
                  
                  
              
                  ### Run the SoftBART function
                  #----------------------------
                  start_time <- Sys.time()
                  cat("Starting SoftBART at ", start_time, file = log_file, append = TRUE)
                  
                  res_softbart <- run_softbart(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_softbart$SNR <- snr_value
                  list_results[[1]] <- res_softbart
                  
                  end_time <- Sys.time()
                  cat("Finished SoftBART at ", end_time,
                          " (", round(difftime(end_time, start_time, units = "secs"), 3), " sec)",
                      file = log_file, append = TRUE)
                  
                  
                  ### Run the BART function
                  #------------------------
                  res_bart <- run_bart(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_bart$SNR <- snr_value
                  list_results[[2]] <- res_bart
                  
                  
                  ### Run the linear model for joint effect
                  #----------------------------------------
                  res_lm <- run_linear_model(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_lm$SNR <- snr_value
                  list_results[[3]] <- res_lm
                  
                  
                  ### Run the DART function
                  #------------------------
                  res_dart <- run_dart(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_dart$SNR <- snr_value
                  list_results[[4]] <- res_dart
                  
                  
                  ### Run the SoftBART function
                  #----------------------------
                  res_softbart_sl <- run_softbart_withoutselection(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_softbart_sl$SNR <- snr_value
                  list_results[[5]] <- res_softbart_sl
                  
                  
                  ### Run the Horseshoe function
                  #-----------------------------
                  start_time <- Sys.time()
                  cat("Starting Horseshoe at ", start_time, file = log_file, append = TRUE)
                  
                  list_results[[6]] <- tryCatch({
                    tmp <- run_horseshoe(
                      df = obs,
                      sim_id = seed_subset[i],
                      psi = psi,
                      joint_eff = joint_eff,
                      q_low, q_high
                    )
                    tmp$SNR <- snr_value
                    tmp
                  }, error = function(e) {
                    cat("Horseshoe failed: ", conditionMessage(e), file = log_file, append = TRUE)
                    NULL
                  })
                  
                  end_time <- Sys.time()
                  cat("Finished Horseshoe at ", end_time,
                          " (", round(difftime(end_time, start_time, units = "secs"), 3), " sec)",
                      file = log_file, append = TRUE)
                  
                  
                  
                  ### Run the Regularized Horseshoe function
                  #-----------------------------------------
                  list_results[[7]] <- tryCatch({
                    tmp <- run_reghorseshoe(
                      df = obs,
                      sim_id = seed_subset[i],
                      psi = psi,
                      joint_eff = joint_eff,
                      q_low, q_high
                    )
                    tmp$SNR <- snr_value
                    tmp
                  }, error = function(e) {
                    cat("RegHorseshoe failed: ", conditionMessage(e),file = log_file, append = TRUE)
                    NULL
                  })
                  
                  
                  ### Run the BAS function
                  #-----------------------
                  res_bas <- run_bas(obs, sim_id = seed_subset[i])
                  res_bas$SNR <- snr_value
                  list_results[[8]] <- res_bas
                  
                  
                  ### Run the SLR function
                  #-----------------------
                  res_slr <- run_slr(obs, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                  res_slr$SNR <- snr_value
                  list_results[[9]] <- res_slr
                  
                  
                  ### Run the Ridge function
                  #-------------------------
                  res_ridge <- run_Ridge(obs , n_boot = 500, sim_id = seed_subset[i] , 
                                         qlow = q_low, qhigh = q_high, psi = psi, joint_eff = joint_eff)
                  res_ridge$SNR <- snr_value
                  list_results[[10]] <- res_ridge
                  
                  
                  ### Run the Lasso function
                  #-------------------------
                  res_lasso <- run_Lasso(obs , n_boot = 500, sim_id = seed_subset[i] , 
                                         qlow = q_low, qhigh = q_high, psi = psi, joint_eff = joint_eff)
                  res_lasso$SNR <- snr_value
                  list_results[[11]] <- res_lasso
                  
                  
                  ### Run the Elastic Net function
                  #-------------------------------
                  res_enet <- run_Enet(obs , n_boot = 500, sim_id = seed_subset[i] , 
                                         qlow = q_low, qhigh = q_high, psi = psi, joint_eff = joint_eff)
                  res_enet$SNR <- snr_value
                  list_results[[12]] <- res_enet
                  
                  
                  ### Run the WQS with standardized continuous exposures
                  #-----------------------------------------------------
                  start_time <- Sys.time()
                  cat("Starting WQS standardised continuous at ", start_time, file = log_file, append = TRUE)
                  
                  list_results[[13]] <- tryCatch({
                    tmp <- run_WQS_rh(obs = obs, sim_id = seed_subset[i], q = NULL)
                    tmp$SNR <- snr_value
                    tmp
                  }, error = function(e) {
                    cat("WQS_cont failed: ", conditionMessage(e),file = log_file, append = TRUE)
                    NULL
                  })
                  
                  end_time <- Sys.time()
                  cat("Finished WQS standardised continuous at ", end_time,
                          " (", round(difftime(end_time, start_time, units = "secs"), 3), " sec)",
                      file = log_file, append = TRUE)
                  
                  
                  ### Run the WQS with quartiles
                  #-----------------------------
                  list_results[[14]] <- tryCatch({
                    tmp <- run_WQS_rh(obs = obs, sim_id = seed_subset[i], q = 4)
                    tmp$SNR <- snr_value
                    tmp
                  }, error = function(e) {
                    cat("WQS_quartile failed: ", conditionMessage(e), file = log_file, append = TRUE)
                    NULL
                  })
                  
              
                  ### Run WQS with deciles
                  #-----------------------
                  list_results[[15]] <- tryCatch({
                      tmp <- run_WQS_rh(obs = obs, sim_id = seed_subset[i], q = 10)
                      tmp$SNR <- snr_value
                      tmp
                    }, error = function(e) {
                      cat("WQS_decile failed: ", conditionMessage(e), file = log_file, append = TRUE)
                      NULL
                    })
                  
                  
                  ### Combine all results into one data frame
                  #------------------------------------------
                  list_results <- list_results[!sapply(list_results, is.null)]
                  data.table::rbindlist(list_results, fill = TRUE)
                  
                }
      }, error = function(e) {
        cat("Error in batch", batch, "SNR", s, ":", conditionMessage(e), "\n", file = log_file, append = TRUE)
        NULL
      })
      
      
      ### Add Batch and SNR columns
      #----------------------------
      if (!is.null(batch_results)) {
        batch_results$Batch <- batch
        
        
        ### Name for list element
        #------------------------
        list_name <- paste0("Batch_", batch, "_SNR_", s)
        
        
        ### Store in top-level results list
        #----------------------------------
        results[[list_name]] <- batch_results
        
        
        ### Save batch results separately if desired
        #-------------------------------------------
        cat("Batch", batch, "SNR", s, "results saved\n", file = log_file, append = TRUE)
      }
      
    }
    
    ### Save checkpoint
    #------------------
    saveRDS(list(last_batch = batch, results = results), file = checkpoint_file)
    cat("Checkpoint saved at batch", batch, "\n", file = log_file, append = TRUE)
  }
})


### Stop cluster
#---------------
stopCluster(cluster)
closeAllConnections()
