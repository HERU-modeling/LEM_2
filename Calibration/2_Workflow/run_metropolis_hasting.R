library(foreach)
library(doParallel)
library(tictoc)



#' Metropolis-Hasting accepting criterion
#'
#' @param init_ll likelihood of the initial sample
#' @param new_ll likelihood of the newly proposed sample
#' @param init_prior_weight likelihood of obtaining the initial prior sample
#' @param new_prior_weight likelihood of obtaining the new prior sample
#'
#' @return True if new proposal is accepted, False if not
#' @export
#'
#' @examples
mh_accept <- function(init_prior_weight, init_ll, new_prior_weight, new_ll){
	alpha <- (new_ll + rowSums(new_prior_weight)) - (init_ll + rowSums(init_prior_weight))
	alpha <- exp(alpha)
	u <- runif(1, 0, 1)
	return(u < alpha)
}


#' Metropolis-Hasting algorithm.
#'
#' @param calpar list of free parameters
#' @param calpar.info summary of the list of free parameters
#' @param vlist unstructured info of all parameters
#' @param vparameters structured parameters
#' @param calib.target calibration targets
#' @param valid.target validation targets
#' @param num_mh_iter number of Metropolis Hasting iterations
#' @return a list with accepted_samples, accepted log-likelihood, posterior samples,
#' predicted target, and observation
#' @export
#'
#' @examples
metropolis_hasting <- function(calpar, calpar.info, vlist, vparameters, 
															 calib.target, valid.target, num_mh_iter = 200,
															 num_resamples = 100, num_targets = 144,
															 get_target = TRUE){
	# Saving results
	npar <- nrow(calpar)
	accepted_samples <- list()
	accepted_ll <- c()

	# Initial run, serve as comparison for the first acceptance/rejection
	init_sample_output <- sample_dist_wrapper(1, calpar, prior_weight = TRUE)
	init_prior_samples <- init_sample_output$randsp
	init_prior_weight <- init_sample_output$sample_weights

	init_ll <- run_model_with_sample(init_prior_samples, calpar.info, vparameters, vlist,
																	 calib.target, valid.target, mode = "ll", num_targets = num_targets,
																	 verbal = FALSE)
	
	pb <- progress_bar$new(format = "  Running Metropolis Hasting [:bar] :percent, eta: :eta",
													 total = num_mh_iter)

	# MH iteration
	for(i in 1:num_mh_iter){
		pb$tick()
		
		# Number of param selected to be updated range from half of total calib param to total calib param -1
	  # if npar = 10, no_selected_param range from 5 to 9
		no_selected_param <- if (npar != 1) {sample(floor(npar/2):(npar - 1), 1)} else 1 #Reminder npar considered each strata of a calib param as different param, so 1 param stratify by 3 races only contribute 3 to npar
		
		idx_selected_param <- sample(1:npar, no_selected_param) 

		# Update new sample (Tri: "Jump" dist) - note that this is ok if the priors of the free parameters
		# are independent with each other. If later on we establish connections between
		# parameters, then it's important to not only update the selected parameter,
		# but also parameters that depend on that parameter. A way to decide the impacted
		# parameter is to topologically sort the parameters
		
		new_prior_samples <- init_prior_samples
		new_prior_weight <- init_prior_weight
		
		new_sample_output <- sample_dist_wrapper(1, calpar[idx_selected_param, ], prior_weight = TRUE)
		
		new_prior_samples[idx_selected_param] <- new_sample_output$randsp
		new_prior_weight[idx_selected_param] <- new_sample_output$sample_weights
		
		# Run with new sample
		new_ll <- run_model_with_sample(new_prior_samples, calpar.info, vparameters, vlist,
																		calib.target, valid.target, mode = "ll", num_targets = num_targets,
																		verbal = FALSE)
		
		accept <- mh_accept(init_prior_weight, init_ll, new_prior_weight, new_ll)

		if (is.na(accept)) {
			cat("NA accept at iteration:", i, "\n")
			browser()
		}
		if(accept){
			accepted_samples[[i]] <- new_prior_samples
			accepted_ll <- c(accepted_ll, new_ll)
			
			init_prior_samples <- new_prior_samples
			
			init_prior_weight <- new_prior_weight
			init_ll <- new_ll
			
		}else{
		  accepted_samples[[i]] <- init_prior_samples
			accepted_ll <- c(accepted_ll, init_ll)
		}
	}
	

	if(get_target){
		# Resample
		posterior_samples <- resample_using_importance_weights(accepted_samples, accepted_ll, 
																													 ndraws = num_resamples)
		
		# Create output table: marginalized posterior
		output_df <- create_posterior_table(posterior_samples, nrow(calpar), calpar.info)
		
		# If we want to predict targets with the posterior distribution
	  target_result<- run_model_with_sample(posterior_samples, calpar.info, vparameters, 
	  																			vlist, calib.target, valid.target, 
	  																			mode = "target", num_targets = num_targets) # Default verbal = TRUE
	  
	}else{
		posterior_samples <- NULL
		output_df <- NULL
	  target_result <- NULL
	}

	list(accepted_samples = accepted_samples,
			 accepted_ll = accepted_ll,
			 output_df = output_df,
			 posterior_samples = posterior_samples,
			 target_model = target_result$model,
			 target_obs = target_result$obs
			 ,outa       = target_result$outa
				) #Output for diagnostic plot, optional for calibration
}




#' Multichain metropolis hasting
#'
#' @param num_chains number of metropolis-hasting chains to run
#' @param num_cores number of cores in the CPU for parallel processing
#' @inheritParams metropolis_hasting
#' @return posterior samples, model prediction, observation, all accepted samples,
#' and the likelihood for all the accepted samples
#' @export
#'
#' @examples
metropolis_hasting_multichain_combine <- function(calpar, calpar.info, vlist, vparameters, 
                                          calib.target, valid.target, num_mh_iter = 200,
                                          num_resamples = 100, num_targets = 144, 
                                          num_chains = 1, num_cores = 2){
  all_lls <- list()
  registerDoParallel(num_cores)
  
  # Create a cluster
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)  # Register the cluster for parallel computing
  clusterSetRNGStream(cl, iseed = 1111)
  
  # Export entire global environment
  clusterExport(cl, varlist = ls(globalenv()), envir = globalenv())
  
  mh_result_size <- 7 # Add outa to to the above list
  mh_ll_result_indices <- (1:num_chains)*mh_result_size - 5
  mh_sample_result_indices <- (1:num_chains)*mh_result_size - 6

  # Run MH in parallel
  tic()
    multichain_result <- foreach(i=1:num_chains, .combine = c, .packages = c("mc2d","openxlsx","deSolve","progress","posterior")) %dopar% {
    # multichain_result <- foreach(i=1:num_chains, .combine = c, .packages = c("mc2d","openxlsx","deSolve","progress","posterior")) %do% {
    metropolis_hasting(calpar, calpar.info, vlist, vparameters,
                                      calib.target, valid.target, num_mh_iter,
                                      num_resamples, num_targets, get_target = FALSE)
  }
  toc()
  # Stop the cluster to free resources
  stopCluster(cl)
  # Process MH output
  all_lls     <- multichain_result[mh_ll_result_indices] # This would be where I can adjust to create a burn-in period
  all_samples <- multichain_result[mh_sample_result_indices]
  
  # Burn-in
  burned <- create_burnin(all_samples, all_lls, burn_frac = 0.5)
  
  # Diagnostics
  diagnostics <- run_diagnostics(
  	burned_samples = burned$all_samples,
  	burned_lls     = burned$all_lls
  )
  
  # Process MH output
  # all_lls     <- unname(unlist(multichain_result[mh_ll_result_indices])) 
  # all_samples <- unname(unlist(multichain_result[mh_sample_result_indices]))
  burned_lls <- unname(unlist(burned$all_lls)) 
  burned_samples <- unname(unlist(burned$all_samples, recursive = FALSE))
  
  # posterior_samples <- resample_using_importance_weights(all_samples, all_lls, 
  #                                                        ndraws = num_resamples)
  posterior_samples <- resample_using_importance_weights(burned_samples, burned_lls,
  	ndraws = num_resamples)
  
  target_result<- run_model_with_sample(posterior_samples, calpar.info, vparameters, 
  																			vlist, calib.target, valid.target, 
  																			mode = "target", num_targets = num_targets) # Default verbal = TRUE

  
  return(list(posterior_samples = posterior_samples,
              target_model = target_result$model,
              target_obs = target_result$obs,
  						outa       = target_result$outa, #Output for diagnostic plot, optional for calibration
  						all_lls  = all_lls,
              all_samples = all_samples,
              burned_lls = burned_lls,
              burned_samples = burned_samples
							,diagnostics = diagnostics
  ))
}

