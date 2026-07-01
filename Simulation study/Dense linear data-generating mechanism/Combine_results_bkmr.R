################################################################################
# Simulation study
#
# Function to store results form different methods in a single dataframe
#
################################################################################


rm(list = ls())


### Load functions that will be used
#-----------------------------------
source("~/Generate_data_local.R")
source("~/Method/BKMR_function_store_results.R")


setwd("~/Dense linear data-generating mechanism")


### Load libraries 
#-----------------
library(data.table)
library(RhpcBLASctl)
library(doParallel)
library(parallel)
library(bkmr)


### Define total number of simulations
#-------------------------------------
m <- 100
n_cores <- 6                 
iterations_per_batch <- n_cores
total_batches <- ceiling(m / iterations_per_batch)


### Define the seed
#------------------
seed <- 85721023 + 0:(m-1)


### Logging and checkpoint
#-------------------------
log_file <- "simulation_log_bkmr.txt"
checkpoint_file <- "checkpoint_results_bkmr.rds"


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
cluster <- makeCluster(n_cores ,  type = "PSOCK")
registerDoParallel(cluster)
clusterExport(cluster,
              varlist = c("Generate_data", "run_bkmr"),
              envir = environment()
)


### Define your SNRs of interest
#-------------------------------
sigma <- c(0.58 , 0.8 , 1.1 , 2.5)


### Run simulation
#-----------------
timing <- system.time({
  
  for (batch in start_batch:total_batches) {
    
    cat("Running batch", batch, "\n")
    
    
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
      
      cat("Running batch for SNR =", s, "\n")
      
      
      ### Run batch in parallel
      #------------------------
      batch_results <- tryCatch({
        foreach(i = 1:iterations_this_batch, 
                .combine = "list",
                .multicombine = TRUE,
                .packages = c("data.table","bkmr")) %dopar% {
                                
                                list_results <- vector("list", 1)
                                
                                
                                ### Effect sizes
                                #---------------
                                psi <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4)
                                joint_eff <- -1
                                
                                
                                ### Generate the data with target sigma/SNR
                                #------------------------------------------
                                obs_list <- Generate_data(seed_subset[i], sigma = s, psi = psi, joint_eff = joint_eff)
                                obs <- obs_list$simdata
                                snr_value <- obs_list$snr
                                
                                
                                ### Define the quantiles for standardisation
                                #-------------------------------------------
                                q_low  <- sapply(obs[2:8], function(x) quantile(as.numeric(x), 0.25))
                                q_high <- sapply(obs[2:8], function(x) quantile(as.numeric(x), 0.75))
                                
                                
                                start_time <- Sys.time()
                                cat("Starting BKMR at ", start_time, file = log_file, append = TRUE)
                                
                                
                                ### Run the BKMR function
                                #------------------------
                                res_bkmr <- run_bkmr(obs, q_low, q_high, sim_id = seed_subset[i], psi = psi, joint_eff = joint_eff)
                                res_bkmr$SNR <- snr_value
                                list_results[[1]] <- res_bkmr
                                
                                
                                end_time <- Sys.time()
                                cat("Finished BKMR at ", end_time,
                                    " (", round(difftime(end_time, start_time, units = "secs"), 3), " sec)",
                                    file = log_file, append = TRUE)
                                
                                
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

