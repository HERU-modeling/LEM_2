#Step 0: Clear R environment
rm(list = ls())
#Step 1: Importing R functions & set seeds
library(openxlsx)
library(mc2d)
library(posterior)


######### 1. Set up ############
#Step 2: Choosing city
CITY <- "MAR" # MIA1, NYC1 for original, KIN, BAL, PHI


#Step 3: Importing LEM related packages
# Importing baseline environment
source("ModelCoreModules/Parameter_input/LoadBaselineWorkspace.R")

#  Load vparameters
source("ModelCoreModules/Parameter_input/Par_targets.R") #Update Target as needed - Part of par_pop.R 
source("ModelCoreModules/Parameter_input/Module-Data.input.R")

# Load calibration related functions
source("Calibration/2_Workflow/util.R")
source("Calibration/2_Workflow/run_metropolis_hasting.R")
source("Calibration/2_Workflow/objective_fn.R")

# Load additional data not covered by current code
source("Calibration/2_Workflow/util_intervention.R")
source("Calibration/4_Visualization/plots.R")


#Step 4: Calibation parameter pre-processing
#This step select the free parameters to calibration. 
calib_param <- c("prop.HIV.aware","noG","nsG","s" ,"ass.eO" ,"ass.eS","psi" ,"phi1_" ,"alpha1_","mix.prob.YO", "mu_adj","d", "uio", "uis","v2","v3","T1_O","O_T")#
param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 

# update newly loaded parameters
calpar <- param_info$calpar
calpar.info <- param_info$calpar.info

# Update vparameters with service updates (on psi) & intervention
vparameters <- add_intervention_to_param(vparameters, "No interventions")

#Step 5:
num_targets = sum(sapply(calib.target, length)) #Reduce hard-code, original 72
calib_name = "Seed"
num_mh_iter = 10
num_resamples = 5
git_commit <- system("git rev-parse --short HEAD", intern = TRUE)

#########Read-in save results####
calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_Combine_1_7114a36_50-40_SFA.rds"))
#calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_iter100_resam40_chains1_cores1_2_7114a36_100-40_SFA.rds"))
calib_param <- calibration_results_multiple$calib_param
#Re-Run toget meta data
param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
calpar <- param_info$calpar
calpar.info <- param_info$calpar.info


posterior_samples <- resample_using_importance_weights(calibration_results_multiple$burned_samples,
																											 calibration_results_multiple$burned_lls,
																											 ndraws = 1000)
prior_samples <- sample_dist_wrapper(1000, calpar)

# Examine the PE, Upper, Lower & 
colnames_all_samples <- paste0(calpar$par, "_", calpar$gender, "_", calpar$ethnicity, "_",
															 calpar$risk, "_", calpar$sex.risk.age.group)
colnames(posterior_samples) <- c(colnames_all_samples)
colnames(prior_samples)     <- c(colnames_all_samples)

posterior_stats <- descriptive_stats(posterior_samples)
prior_stats     <- descriptive_stats(prior_samples)
#posterior_stats[posterior_stats$name == "noG_m_b_mwid_Young", ]
#prior_stats[prior_stats$name == "psi_f_b_het_na", ]
#print(posterior_stats[grepl("^prop.HIV.aware", posterior_stats$name), ], n = Inf)
#print(prior_stats[grepl("^noG", posterior_stats$name), ], n = Inf)
# posterior_stats_long <- expand_posterior_stats(posterior_stats)

# Generate plot
# plot_prior_posterior_comparison(prior_samples, posterior_samples, calpar, adjust = 2)
# Comment for ggsave()
# What you should be aiming for. TODO: this is not very useful, good to show, 
# if prior/posterior similar, then dont put this in the paper. what should you expect.

#Some literature argue that the similarity between prior and posterior is a sign for 
#non-identifiability, because updating prior with data does not seem to alter
#the distribution of it. This could be true, but it can also be because the posterior
#should indeed be very similar to the prior, because of good data quality.

#Therefore I suggest to include this comparison plot in the paper only if the 
#prior and posterior differs, as a way to avoid argument and confusion.

#Prior and posterior are considered similar if they overlap more than 30% (suggested
#by the person who created the Latent Class Identifiability Display)


#6.2 Plot MCMC traceplots
#plot_multichain_trace(calibration_results_single$accepted_ll)
plot_multichain_trace(calibration_results_multiple$all_lls)


#Stationary: the traces are stablized and not going upward or downward.
#Mixing: all four chains are mixed together, no one chain dominate all the time

#6.3 Check collinearity result
# Check structure (just to confirm again)
str(posterior_samples)
n_cols <- ncol(posterior_samples)
chunk_size <- 21

# All
png(paste0("Model_development/7.Red_strata/Calibration/posterior_distribution_marginal-", CITY, ".png"),
		width = 30, height = 30, units = "in", res = 300)  # make it larger for 14 parameters

psych::pairs.panels(posterior_samples)

dev.off()

# Loop through posterior_samples in chunks Or...
for (i in seq(1, n_cols, by = chunk_size)) {
	# Define the end of the chunk, making sure we don't exceed 121
	end <- min(i + chunk_size - 1, n_cols)
	
	# Construct the filename
	filename <- paste0("Model_development/7.Red_strata/Calibration/posterior_distribution_marginal-",
										 CITY, "_cols_", i, "_to_", end, ".png")
	
	# Save plot
	png(filename, width = 12, height = 12, units = "in", res = 300)
	psych::pairs.panels(posterior_samples[, i:end])
	dev.off()
}
# ... Pick specific parameters to check collinearlirty 
colnames_all_samples
par_of_interest <- c("prop.HIV.aware", "mor_S")
# Create a regex pattern like "^prop.HIV.aware|^mor_S"
pattern <- paste0("^", par_of_interest, collapse = "|")
par_ind <- colnames_all_samples[grep(pattern, colnames_all_samples)]
# Collapse into a single string for the filename
par_str <- paste(par_of_interest, collapse = "_")
png(paste0("Model_development/7.Red_strata/Calibration/posterior_distribution_marginal-", CITY,"-for-", par_str, ".png"),
		width = 12, height = 12, units = "in", res = 300)  
psych::pairs.panels(posterior_samples[,par_ind])
dev.off()

#
library(psych)
# Calculate correlation matrix only (without plotting)
#cor_mat <- cor(posterior_samples[, par_ind], use = "pairwise.complete.obs")
cor_mat <- cor(posterior_samples, use = "pairwise.complete.obs")
# Get upper triangle of correlation matrix (excluding diagonal)
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
# Convert to long format to filter
library(reshape2)
cor_df <- melt(cor_mat, na.rm = TRUE)
# Filter for strong correlations (e.g., > 0.4)
high_corr <- subset(cor_df, value > 0.4)
# View result
print(high_corr)



#######Obs Vs Target####
#NEW Plot: based on plot_target_series
  target_name <- c("diag.W","diag.B","diag.H","ndiag.M_HET","ndiag.F_HET","ndiag.Young","ndiag.Old","ndiag.PWID","ndiag.W","ndiag.B","ndiag.H",
  								 "p.diag_late.W","p.diag_late.B","p.diag_late.H","p.viral_sup.W","p.viral_sup.B","p.viral_sup.H")
  


  lem2_targets_output <- plot_target_2(results = calibration_results_single, 
  																	 target_name= target_name, y_from_zero = TRUE)
  lem2_targets_output <- plot_target_2_median(results = calibration_results_multiple, 
  																	 target_name= target_name, y_from_zero = TRUE)
  
  lem3_targets_output <- plot_target_distribution(
  	results = calibration_results_multiple,
  	target_name = target_name,
  	y_from_zero = TRUE
  )
#Table Results
  target_obs_table <- matrix(calibration_results_multiple$target_obs,
  													 nrow = 17, ncol = 6, byrow = TRUE)
  rownames(target_obs_table) <- target_name

  # Median 50% quantile
  target_model_table_median <- matrix(apply(calibration_results_multiple$target_model, 2, quantile, probs = 0.5, na.rm = TRUE)
  														 , nrow = 17, ncol = 6, byrow = TRUE) 
  # Lower 2.5% quantile
  target_model_table_lower <- matrix(apply(calibration_results_multiple$target_model, 2, quantile, probs = 0.025, na.rm = TRUE)
  														 , nrow = 17, ncol = 6, byrow = TRUE) 
  # Upper 97.5% quantile
  target_model_table_upper <- matrix(apply(calibration_results_multiple$target_model, 2, quantile, probs = 0.975, na.rm = TRUE)
  														 , nrow = 17, ncol = 6, byrow = TRUE) 
  # Mean
  target_model_table_mean <- matrix(colMeans(calibration_results_multiple$target_model),
  														 nrow = 17, ncol = 6, byrow = TRUE) 
  rownames(target_model_table_median) <- rownames(target_model_table_lower) <- rownames(target_model_table_upper) <- rownames(target_model_table_mean) <- target_name
  
source("Model_development/Optional_util_diagnostic_code/Utility_functions.R")
  p <- plot_pct_diff_heatmap(
  	target_obs_table,
  	target_model_table_mean,
  	output_dir = file.path(getwd(), "Model_development", "Optional_util_diagnostic_code"),
  	cap_value = 100
  )
 
  ###Loop all 5 Cities
  calib_NYC1 <- readRDS(paste0("Model_development/7.Red_strata/Calibration/Results/calib_comb-Base_v2v3_f23ecaf4-251126_100-100_NYC1.rds"))
  calib_MIA1 <- readRDS(paste0("Model_development/7.Red_strata/Calibration/Results/calib_comb-Base_v2v3_f23ecaf4-251126_100-100_MIA1.rds"))
  calib_BAL  <- readRDS(paste0("Model_development/7.Red_strata/Calibration/Results/calib_comb-Base_v2v3_f23ecaf4-251126_100-100_BAL.rds"))
  calib_PHI  <- readRDS(paste0("Model_development/7.Red_strata/Calibration/Results/calib_comb-Base_v2v3_f23ecaf4-251126_100-100_PHI.rds"))
  calib_KIN  <- readRDS(paste0("Model_development/7.Red_strata/Calibration/Results/calib_comb-Base_v2v3_f23ecaf4-251126_100-100_KIN.rds"))
  
  avg_target_obs <- rowMeans(
  	cbind(calib_NYC1$target_obs, calib_MIA1$target_obs, calib_BAL$target_obs, calib_PHI$target_obs, calib_KIN$target_obs),
  	na.rm = TRUE)
  target_obs_table_average <- matrix(avg_target_obs,
  													 nrow = 17, ncol = 6, byrow = TRUE)
  rownames(target_obs_table_average) <- target_name
  avg_target_model <- rowMeans(
  	cbind(colMeans(calib_NYC1$target_model), colMeans(calib_MIA1$target_model), colMeans(calib_BAL$target_model), colMeans(calib_PHI$target_model), colMeans(calib_KIN$target_model)),
  	na.rm = TRUE)  
  target_model_table_average <- matrix(avg_target_model,
  														 nrow = 17, ncol = 6, byrow = TRUE) 
  
  rownames(target_model_table_average) <- target_name
  p <- plot_pct_diff_heatmap(
  	target_obs_table_average,
  	target_model_table_average,
  	output_dir = file.path(getwd(), "Model_development", "Optional_util_diagnostic_code"),
  	cap_value = 100
  )
  #Best Measurement 1:
  pct_diff_average <- (target_model_table_average - target_obs_table_average) / target_obs_table_average * 100
  mean(abs(pct_diff_average))
  colMeans(abs(pct_diff_average))
  rowMeans(abs(pct_diff_average))
  #Best Measurement 2:
  c(mean(calib_NYC1$all_lls),mean(calib_MIA1$all_lls),mean(calib_BAL$all_lls),mean(calib_PHI$all_lls),mean(calib_KIN$all_lls))
  mean(c(mean(calib_NYC1$all_lls),mean(calib_MIA1$all_lls),mean(calib_BAL$all_lls),mean(calib_PHI$all_lls),mean(calib_KIN$all_lls)))

###########Extract accepted samples list#####
  # Extract accepted samples list
  # samples <- calibration_results_multiple[["all_samples"]][["accepted_samples"]]
  samples <- calibration_results_multiple$burned_samples
  # Extract value at param position (ex, psi at 51) from each of the 2000 samples
  values <- sapply(samples, function(x) x[1, 49])
  
  # Create dataframe for plotting
  df_plot <- data.frame(
  	sample = seq_along(values),
  	value  = values
  )
  
  # Plot
  library(ggplot2)
  
  ggplot(df_plot, aes(x = sample, y = value)) +
  	geom_line() +
  	labs(
  		x = "Sample Index",
  		y = "Value",
  		title = "Trace Plot of Parameter"
  	) +
  	theme_bw()
####### Output Posteriors for multiple results ####
source("Model_development/Optional_util_diagnostic_code/Utility_functions.R") 
  param <- "All" #"v2v3","uio","uis","O_T","d" - Note to match save results names
  prefix    <- paste0("Model_development/7.Red_strata/Calibration/Results/calib_",param,"_10000-5000-jun23_")
  save_path <- paste0("Model_development/7.Red_strata/Calibration/Results/",param,"_posterior_values_output.xlsx")
 Jurisdictions <- all.cities #all.cities #c("BAL", "PHI", "KIN", "NYC1", "MIA1") #Default is all jurisdictions
  posterior_results <- output_posterior_val(Jurisdictions, prefix, save_path)  # If Error, re-run read-in save results with correct param
  # Change its to more readable form
  posterior_results_long <- expand_posterior_stats(posterior_results)
  write_xlsx(posterior_results_long, save_path)
  
####### Examine Posterior Distribution ####
  library(dplyr)
  library(ggplot2)
  
  filter_combinations <- posterior_results_long %>%
  	distinct(par, gender, risk, sex.risk.age.group)
  
  df_plot <- posterior_results_long %>%
  	filter(
  		par == "O_T",
  		gender == "na",
  		risk == "na",
  		sex.risk.age.group == "na"
  	)
  
  # build dynamic title from the data (safe even if you change filters later)
  title_text <- df_plot %>%
  	summarise(
  		txt = paste0(
  			"par: ", first(par),
  			" | gender: ", first(gender),
  			" | risk: ", first(risk),
  			" | age group: ", first(sex.risk.age.group)
  		)
  	) %>%
  	pull(txt)
  
  ggplot(df_plot, aes(x = mean)) +
  	geom_histogram(
  		binwidth = 0.001, #0.001
  		fill = "steelblue",
  		color = "black",     # makes bars clearly separated
  		linewidth = 0.3
  	) +
  	facet_grid(ethnicity ~ par) +
  	labs(
  		title = title_text,
  		x = "Mean",
  		y = "Count"
  	) +
  	theme_minimal() +
  	theme(
  		panel.spacing = unit(1, "lines"),
  		strip.text = element_text(face = "bold"),
  		plot.title = element_text(hjust = 0.5)
  	)

####### New posterior sample from saved results ####  
 # Get calibration_results_multiple
  
  num_resamples = 6888
  
  burned_lls    <- calibration_results_multiple$burned_lls
  burned_samples<- calibration_results_multiple$burned_samples
  
  #browser()
  # posterior_samples <- resample_using_importance_weights(all_samples, all_lls, 
  #                                                        ndraws = num_resamples)
  posterior_samples <- resample_using_importance_weights(burned_samples, burned_lls,
  																											 ndraws = num_resamples)
  
  #browser()
  target_result<- run_model_with_sample(posterior_samples, calpar.info, vparameters, 
  																			vlist, calib.target, valid.target, 
  																			mode = "target", num_targets = num_targets) # Default verbal = TRUE
  

####### Exaime join posterior distribution from posterior samples ####  
  ##: Run 1.Set-up
# Here is the idea, Look at 2 parameters (have to be fully stratified) posterior values 
# Plot them to see how they cluster
  
  # Extract Seed
  calibration_results_1 <- readRDS(paste0("Calibration/3_Results/calib_Seed_1111_590c3bb_1000-500_KIN.rds"))
  calibration_results_2 <- readRDS(paste0("Calibration/3_Results/calib_Seed_2222_590c3bb_1000-500_KIN.rds"))
  calibration_results_3 <- readRDS(paste0("Calibration/3_Results/calib_Seed_3333_590c3bb_1000-500_KIN.rds"))
  calibration_results_4 <- readRDS(paste0("Calibration/3_Results/calib_Seed_4444_590c3bb_1000-500_KIN.rds"))
  
  #Re-Run to get meta data
  calib_param <- calibration_results_4$calib_param
  param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
  calpar <- param_info$calpar
  calpar.info <- param_info$calpar.info
  
  # Re-sample to get the posterior sampkes
  posterior_samples <- resample_using_importance_weights(calibration_results_1$burned_samples,
                                                         calibration_results_1$burned_lls,
                                                         ndraws = 1000)
  colnames_all_samples <- paste0(calpar$par, "_", calpar$gender, "_", calpar$ethnicity, "_",
                                 calpar$risk, "_", calpar$sex.risk.age.group)
  colnames(posterior_samples) <- c(colnames_all_samples)

  
  plot_posterior_pair <- function(posterior_samples, param_x, param_y) {
    
    x <- posterior_samples[, param_x]
    y <- posterior_samples[, param_y]
    
    plot(
      x, y,
      xlab = param_x,
      ylab = param_y,
      pch = 16,
      cex = 0.7
    )
  }
  plot_posterior_pair(
    posterior_samples,
    param_x = "nsG_na_w_msm_old",
    param_y = "psi_na_w_na_na"
  )
  