# rm(list = ls())

########################################
# Load calibration related functions
source("Calibration/2_Workflow/util.R")

########################################
#Read-in save results
calibration_results_1 <- readRDS(paste0("Calibration/3_Results/calib_Combine_1_7114a36_50-40_SFA.rds"))
calibration_results_2 <- readRDS(paste0("Calibration/3_Results/calib_Combine_2_7114a36_50-40_SFA.rds"))
calibration_results_3 <- readRDS(paste0("Calibration/3_Results/calib_Combine_3_7114a36_50-40_SFA.rds"))

combine_top_calibration <- function(..., top_prop = 0.50) {
  
  calibration_results <- list(...)
  
  if (top_prop <= 0 || top_prop > 1) {
    stop("top_prop must be > 0 and <= 1.")
  }
  
  selected <- lapply(calibration_results, function(x) {
    
    burned_lls <- x$burned_lls
    burned_samples <- x$burned_samples
    
    if (length(burned_lls) != length(burned_samples)) {
      stop("burned_lls and burned_samples must have the same length.")
    }
    
    # Number of observations to keep
    n_keep <- ceiling(length(burned_lls) * top_prop)
    
    # Identify which observations have the highest likelihoods
    # (this determines WHICH observations are selected)
    top_idx <- order(burned_lls, decreasing = TRUE)[seq_len(n_keep)]
    
    # Sort the selected indices back into their ORIGINAL order
    top_idx <- sort(top_idx)
    
    list(
      samples = burned_samples[top_idx],
      lls = burned_lls[top_idx]
    )
  })
  
  # Combine while preserving the original order within each calibration result
  combined_samples <- unlist(
    lapply(selected, `[[`, "samples"),
    recursive = FALSE
  )
  
  combined_lls <- unlist(
    lapply(selected, `[[`, "lls")
  )
  
  list(
    combined_samples = combined_samples,
    combined_lls = combined_lls
  )
}

combined_results <- combine_top_calibration(
  calibration_results_1,
  calibration_results_2,
  calibration_results_3,
  top_prop = 0.50
)

posterior_samples_combined <- resample_using_importance_weights(combined_results$combined_samples,
                                                         combined_results$combined_lls,
                                                         ndraws = 100)
target_result<- run_model_with_sample(posterior_samples_combined, calpar.info, vparameters, 
                                      vlist, calib.target, valid.target, 
                                      mode = "target", num_targets = num_targets)
calibration_results <- list(posterior_samples_combined = posterior_samples_combined,
                            target_model = target_result$model,
                            target_obs = target_result$obs,
                            outa       = target_result$outa,
                            all_lls    = combined_results$combined_lls)
saveRDS(calibration_results, paste0("Calibration/3_Results/calib_test_",CITY,".rds"))
