#Step 0: Clear R environment
rm(list = ls())
#Step 1: Importing R functions & set seeds
library(openxlsx)
library(mc2d)
library(posterior)
library(tictoc)
set.seed(1111)


Jurisdictions <- c("SFA") #"MAR","SFA", "LOR", "SAC", "SDI", "DUV","TAM","BMP", "ORF", "ATL", "COO","MIN", "EBR", "ORL", "BAL", "SUF","WAY", "CLA",  "NYC", "NJE", "MEC","CUY", "FRA", "HAM","PHI","SHE", "BEX","DFW","HAR", "TRV", "KIN", "DCM"
#

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
	calib_name = "Parallel3"
	num_mh_iter = 50
	num_resamples = 40
	git_commit <- system("git rev-parse --short HEAD", intern = TRUE)
	
	#Step 5.1: Run single chain Metropolis Hasting
	# calibration_results_single <- metropolis_hasting(calpar, calpar.info, vlist, vparameters,
	# 																							 calib.target, valid.target, num_mh_iter = 8,
	# 																						 num_resamples = 8, num_targets = num_targets)
	
	#saveRDS(calibration_results_single, paste0("Model_development/7.Red_strata/Calibration/Results/calibration_single-chain_results-", CITY, ".rds"))

	
	#Step 5.2: Run multiple chain Metropolis Hasting
	calibration_results_multiple <- metropolis_hasting_multichain_combine(calpar, calpar.info, vlist, vparameters,
																																				calib.target, valid.target, num_mh_iter = num_mh_iter,
																																				num_resamples = num_resamples, num_targets = num_targets, num_chains = 4,
																																				num_cores = 4)
	
	#Add meta-data
	calibration_results_multiple$calib_param <- calib_param

	
	saveRDS(calibration_results_multiple, paste0("Calibration/3_Results/calib_",calib_name,"_",git_commit,"_", num_mh_iter,"-", num_resamples,"_", CITY, ".rds"))

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
