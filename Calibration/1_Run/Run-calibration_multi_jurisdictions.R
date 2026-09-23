#Step 0: Clear R environment
rm(list = ls())
#Step 1: Importing R functions & set seeds
library(openxlsx)
library(mc2d)
library(posterior)
library(tictoc)
# set.seed(1111)


Jurisdictions <- c("MAR") #"MAR","SFA", "LOR", "SAC", "SDI", "DUV","TAM","BMP", "ORF", "ATL", "COO","MIN", "EBR", "ORL", "BAL", "SUF","WAY", "CLA",  "NYC", "NJE", "MEC","CUY", "FRA", "HAM","PHI","SHE", "BEX","DFW","HAR", "TRV", "KIN", "DCM"


for (CITY in Jurisdictions) {
	tic(paste("Model run for", CITY))

################### Start Copy step 3,4,5 from calibration-demo-multichain-mcmc #####

	#Step 3: Importing LEM related packages
	# Importing baseline environment
	source("ModelCoreModules/Parameter_input/LoadBaselineWorkspace.R")
	
	#  Load vparameters
	# source("ModelCoreModules/Parameter_input/_targets.R") #Update Target as needed - Part of par_pop.R 
	source("ModelCoreModules/Parameter_input/Module-Data.input.R")
	
	# Load calibration related functions
	source("Calibration/2_Workflow/util.R")
	source("Calibration/2_Workflow/run_metropolis_hasting.R")
	source("Calibration/2_Workflow/objective_fn.R")
	
	# Load additional data not covered by current code
	source("Calibration/2_Workflow/util_intervention.R")
	source("Calibration/4_Visualization/plots.R")
	
	########################################
	#Step 4: Calibation parameter pre-processing
	#This step select the free parameters to calibration. 
	calib_param <- c("prop.HIV.aware","noG","nsG","s" ,"ass.eO" ,"ass.eS","psi" ,"phi1_" ,"alpha1_","mix.prob.YO", "mu_adj","d", "uio", "uis","v2","v3","T1_O","O_T")#
	param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
	
	# update newly loaded parameters
	calpar <- param_info$calpar
	calpar.info <- param_info$calpar.info
	
	# Update vparameters with service updates (on psi) & intervention
	vparameters <- add_intervention_to_param(vparameters, "No interventions")

	########################################
	#Step 5:
	num_targets = sum(sapply(calib.target, length)) #Reduce hard-code, original 72
	calib_name = "get_target_true"
	num_mh_iter = 100
	num_resamples = 100
	# seed = 3333
	git_commit <- system("git rev-parse --short HEAD", intern = TRUE)
	# git_commit <- "04cff6d"
	
	seed_list          <- c(1111,2222,3333,4444,5555)
	# Store results from each seed
	calibration_results_all <- vector("list", length(seed_list))
	
	for (i in seq_along(seed_list)) {

	  seed <- seed_list[i]
	#Step 5.1: Run multiple chain Metropolis Hasting for each seed
	calibration_results_multiple <- metropolis_hasting_multichain_combine(calpar, calpar.info, vlist, vparameters,
	                                                                      calib.target, valid.target, num_mh_iter = num_mh_iter,
	                                                                      num_resamples = num_resamples, num_targets = num_targets, get_target = TRUE,
	                                                                      num_chains = 4, num_cores = 4, seed = seed)

	#Add meta-data
	calibration_results_multiple$calib_param <- calib_param
	saveRDS(calibration_results_multiple, paste0("Calibration/3_Results/calib_",calib_name,"_",seed,"_",git_commit,"_",num_mh_iter,"-",num_resamples,"_",CITY,".rds"))
	# # Store result for later combination
	# calibration_results_all[[i]] <- calibration_results_multiple
	# Remove large result from RAM before next iteration
	rm(calibration_results_multiple)

	# Force R to release unused memory
	gc()
	}
	
	calibration_results_1 <- readRDS(paste0(
	  "Calibration/3_Results/calib_", calib_name, "_1111_", git_commit, "_",
	  num_mh_iter, "-", num_resamples, "_", CITY, ".rds"
	))
	calibration_results_2 <- readRDS(paste0(
	  "Calibration/3_Results/calib_", calib_name, "_2222_", git_commit, "_",
	  num_mh_iter, "-", num_resamples, "_", CITY, ".rds"
	))
	calibration_results_3 <- readRDS(paste0(
	  "Calibration/3_Results/calib_", calib_name, "_3333_", git_commit, "_",
	  num_mh_iter, "-", num_resamples, "_", CITY, ".rds"
	))
	calibration_results_4 <- readRDS(paste0(
	  "Calibration/3_Results/calib_", calib_name, "_4444_", git_commit, "_",
	  num_mh_iter, "-", num_resamples, "_", CITY, ".rds"
	))
	calibration_results_5 <- readRDS(paste0(
	  "Calibration/3_Results/calib_", calib_name, "_5555_", git_commit, "_",
	  num_mh_iter, "-", num_resamples, "_", CITY, ".rds"
	))
	
	calibration_results_all <- list(
	  calibration_results_1,
	  calibration_results_2,
	  calibration_results_3,
	  calibration_results_4,
	  calibration_results_5
	)
	# Step 5.2: Combine top results across seeds
	# num_targets = sum(sapply(calib.target, length)) #Reduce hard-code, original 72
	calib_name = "get_target_true_combined"
	num_resamples = 100
	set.seed(1111)
	
	combined_results <- do.call( combine_top_calibration,  c(calibration_results_all, list(top_prop = 0.50)))

	# Resample
	posterior_samples_combined <- resample_using_importance_weights(combined_results$combined_samples,
	                                                                combined_results$combined_lls,
	                                                                ndraws = num_resamples) #num_resamples
	target_result <- run_model_with_sample(posterior_samples_combined, calpar.info, vparameters, 
	                                      vlist, calib.target, valid.target, 
	                                      mode = "target", num_targets = num_targets)
	
	calibration_results <- list(posterior_samples_combined = posterior_samples_combined,
	                            target_model = target_result$model,
	                            target_obs = target_result$obs,
	                            outa       = target_result$outa,
	                            all_lls    = combined_results$combined_lls)
	
	saveRDS(calibration_results, paste0("Calibration/3_Results/calib_",calib_name,"_combined_",git_commit,"_",num_mh_iter,"-",num_resamples,"_",CITY,".rds"))
###################END   Copy step 3,4,5 from calibration-demo-multichain-mcmc#####

cat("Completed model run for", CITY, "\n")
# --- FORCED CLEANUP TO AVOID SOCKET LEAKS ---
# Close *all* nonessential connections except stdin/out/err (0,1,2)
open_con <- showConnections(all = TRUE)
non_essential <- as.integer(rownames(open_con)[as.integer(rownames(open_con)) > 2])
for (con in non_essential) {
	try(close(con), silent = TRUE)
}

# Reclaim memory after closing connections
gc()
}
