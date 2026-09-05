#Step 0: Clear R environment
rm(list = ls())
#Step 1: Importing R functions & set seeds
library(openxlsx)
library(mc2d)
library(posterior)
# set.seed(1111)

#Setting seed for reproducibility. Note that setting seeds doesn't guarantee the 
#same output on different computer therefore you might be seeing different results
#than I do. Setting seeds only guarantees the same output on the same computer.

########################################
#Step 2: Choosing city
CITY <- "MAR" # MIA1, NYC1 for original, KIN, BAL, PHI


########################################
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
calib_name = "Seed"
num_mh_iter = 10
num_resamples = 5
# seed = 3333
git_commit <- system("git rev-parse --short HEAD", intern = TRUE)

########################################
#Read-in save results
calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_Combine_1_7114a36_50-40_SFA.rds"))
#calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_iter100_resam40_chains1_cores1_2_7114a36_100-40_SFA.rds"))
calib_param <- calibration_results_multiple$calib_param
#Re-Run toget meta data
param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
calpar <- param_info$calpar
calpar.info <- param_info$calpar.info

#External test
#calibration_results_multiple <- readRDS(paste0("Model_development/external_rds/Azzam-0911_multiple-chain_results-", CITY, ".rds"))

#Step 6: Diagnoses on calibration result

#Diagnoses in the format of plot are all of the object `ggplot2`. Meaning that you
#could assign a variable name to the plot and use `+` to customize plot, just
#like with function `ggplot()`


#6.1 Check parameter prior & marginalized posterior

# Resample samples with their likelihood
# posterior_samples <- resample_using_importance_weights(calibration_results_multiple$all_samples, 
# 																											 calibration_results_multiple$all_lls, 
# 																											 ndraws = 1000)
posterior_samples <- resample_using_importance_weights(calibration_results_multiple$burned_samples,
																											 calibration_results_multiple$burned_lls,
																											 ndraws = 1000)
# posterior_samples <- resample_using_importance_weights(calibration_results_single$accepted_samples, 
# 																											 calibration_results_single$accepted_ll, 
# 																											 ndraws = 1000)

# Obtain prior samples - this ideally comes from the proposal distribution,
# it's ok for now since priors are independent
prior_samples <- sample_dist_wrapper(1000, calpar)
# gamma (55.14, 0.0127)

#The fractional ESS should be within a certain range, for example, [30%, 80%]. A
#value too low means that only a few of the prior samples gave reasonable likelihood
#between model output and observation, a value too high means all prior samples
#gave similar likelihood, which may reflect the lack of variability in parameters.

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


#6.4 Create LEM 1.0 output & plot comparison
# Sourcing to reset calpar to the Morris' calpar & run to generate model output
#source("ModelCoreModules/Parameter_input/CascadeCEA-Model-1-Module-Data.input.R")
#source("ModelCoreModules/Calibration/load_weights.R")
#calpar.info <- create_calpar_info(calpar)
#load(paste0("Data_Files/Result/LEM1.0/LEM_1.0_output/", "outputs_", lem_1.0_output_file, ".RData"))
#lem1.0_result <- run_model_with_sample(t(nmcl), calpar.info, vparameters, vlist,
#																			 calib.target, valid.target, mode = "target", 
#																			 num_targets = 72)
#lem1.0_result$target_model <- lem1.0_result$model
#lem1.0_result$target_obs <- lem1.0_result$obs

#This code reset the free parameters to the free parameters selected based on Morris
#Method, and run sampling to get the prediction for new diagnoses



# Generate comparison
#plot_lem_comparison(lem1.0_result, calibration_results, gp18.gn, names18)

#Customizability of ggplot
#g <- plot_lem_comparison(lem1.0_result, calibration_results, gp18.gn, names18)
#g + theme_gray()


  
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
  

####### Stacked Plot ####
  outa <- apply(calibration_results_multiple$outa[ , , , ], c(2,3,4), mean) #Average the value of all the simulation across 3 dimension: Time x state x strata
  outa <- out_sens
  outa <- simul.in
  
  #
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  
  # All 30 strata
#  n.strata    <- apply(outa[ , , ], c(1,3), sum) # across 2 dimension: Time x strata
  n.strata    <- apply(outa[ , c(1:16), ], c(1,3), sum) # across 2 dimension: Time x strata; forget that outa also include other, so look at pop have to sub set c(1:16)
	  #
  plot_stacked_area <- function(time, 
  															data, 
  															name, 
  															CITY = CITY, 
  															xmin = 0, xmax = 348,
  															y_label = "Total population: 15-64 (in millions)",
  															start_year = 2016,
  															color_values = NULL) {
  	
  	# Prepare dataframe
  	df <- as.data.frame(cbind(time, data))
  	colnames(df) <- c("time", name)
  	
  	# Optional zoom in time
  	df <- df[xmin:xmax, ]
  	
  	# Pivot to long format
  	df_long <- df %>%
  		pivot_longer(
  			cols = -time,
  			names_to = "variable",
  			values_to = "value"
  		)
# browser()	
		var_order <- df_long %>%
			distinct(variable) %>%
			pull(variable)
  	# Compute cumulative stacked heights
  	df_long <- df_long %>%
  		group_by(time) %>%
  		# arrange(variable) %>%  # arrange alphabetically
  		arrange(match(variable, var_order), .by_group = TRUE) %>%
  		mutate(
  			cum_value = cumsum(value),
  			mid_value = cum_value - value/2,
  			total_value = sum(value),
  			pct = round((value / total_value) * 100, 1)
  		) %>%
  		ungroup()
# browser()  	
  	# Compute year corresponding to each time
  	df_long$year <- start_year + (df_long$time) / 12
  	
  	# Subset the labels "value" at min and max time point
 # 	label_data <- df_long %>% filter(time %in% c(xmin, xmax-1)) #miss-match between time & index
  	# Labels only at first and last year
  	label_data <- df_long %>%
  		filter(time == min(time) | time == max(time)) %>%
  		filter(pct >= 1) %>%   # hide tiny groups
  		mutate(
  			hjust = ifelse(time == min(time), 0, 1),
  			year  = ifelse(time == min(time), year + 1.6, year - 1.6)
  		) #%>%
  		# group_by(time) %>%
  		# mutate(
  		# 	label_y = mid_value + (row_number() - mean(row_number())) * 0.05 * max(mid_value)
  		# ) %>%
  		# ungroup()
# browser()  	
  	# Reverse factor levels for correct stacking
  	df_long$variable <- factor(df_long$variable, levels = rev(unique(df_long$variable)))
  	label_data$variable <- factor(label_data$variable, levels = rev(unique(label_data$variable)))
 # browser() 	
  	# Find ymax automatically
  	ymax_auto <- max(df_long$total_value, na.rm = TRUE)
  	# Create plot
  	g <- ggplot(df_long, aes(x = year, y = value, fill = variable)) +
  		geom_area(color = "black", alpha = 0.8) +
  		geom_text(
  			data = label_data,
  			aes(x = year, y = mid_value, label = paste0(round(value)," (",pct,"%)")),
  			# aes(x = year, y = label_y, label = paste0(round(value)," (",pct,"%)")),
  			color = "black",
  			size = 4,
  			position = "identity"
  		) +
  		coord_cartesian(ylim = c(0, ymax_auto)) +  # <-- Only ymin = 0 and ymax from data
  		theme_minimal(base_size = 14) +
  		theme(legend.position = "bottom") +
  		labs(
  			# title = paste0("Model population (", CITY, ") by group"),
  			x = "Year",
  			y = y_label,
  			fill = NULL #"Group"
  		) +
  		# scale_y_continuous(labels = scales::comma_format(scale = 1e-6, suffix = " M")) +
  		scale_y_continuous(labels = scales::comma_format(scale = 1e-6)) +
  		scale_x_continuous(
  			breaks = seq(2015, 2045, by = 5),
  			labels = as.character(seq(2015, 2045, by = 5))
  		) +
  		guides(fill = guide_legend(reverse = FALSE))
  	
  	if (!is.null(color_values)) {
  		g <- g + scale_fill_manual(values = color_values)
  	}
  	# else: do nothing → ggplot default (same as before)
  	
  	return(g)
  }
  
  ## Compare Pop Proportion across strata
  # 30 strata
  g1 <- plot_stacked_area(time = vt[-1], data=n.strata, name = vparameters$names.gp, CITY=CITY) #time = vt[-1] for str[outa] =349, =vt for 348
		
  # 18 strata
  group18 = function(group30){
  	g18 = matrix(0, length(vt[c(-1)]), 18) #length(vt[c(-1)]) or length(vt)
  	for (i in 1:18){
  		ind = which(rname %in% vparameters$names18[i])
  		g18[ ,i] = rowSums(group30[ ,ind, drop = FALSE]) #must include drop= FALSE
  		colnames(g18) <- vparameters$names18 
  	}
  	return(g18)
  }
  n.strata.18    <- group18(n.strata)
  g2 <- plot_stacked_area(time = vt[-1], data=n.strata.18, name = vparameters$names18, CITY) #time = vt[-1]
  # 18 strata - Exclude HET
  n.strata.18.exHET   <-n.strata.18[,1:12]
  g2.ex <- plot_stacked_area(time = vt[-1], data=n.strata.18.exHET , name = vparameters$names18[1:12], CITY) #time = vt[-1]
  
  # Calibration-level strata
  group.calibration.exp = function(group30){
  	g.cal = matrix(0,length(vt[c(-1)]),7)           #length(vt[c(-1)]) or length(vt)
  	g.cal[ ,1] = rowSums(group30[ ,gp.gn$L.het.f])  #F/HET Low
  	g.cal[ ,2] = rowSums(group30[ ,gp.gn$H.het.f])  #F/HET high
  	g.cal[ ,3] = rowSums(group30[ ,gp.gn$L.het.m])  #M/HET Low
  	g.cal[ ,4] = rowSums(group30[ ,gp.gn$H.het.m])  #M/HET High
  	g.cal[ ,5] = rowSums(group30[ ,gp.gn$young])   #15-34/MSMMWID
  	g.cal[ ,6] = rowSums(group30[ ,gp.gn$old])  #35-34/MSMWID
  	g.cal[ ,7] = rowSums(group30[ ,gp.gn$pwid])   #PWID
  	colnames(g.cal) <- c("Female:Het/Low","Female:Het/High","Male:Het/Low","Male:Het/High","Male:MSM/MWID 15-34","Male:MSM/MWID 35-64","PWID")
  	return(g.cal)
  }
  group.calibration = function(group30){
 	g.cal = matrix(0,length(vt[c(-1)]),5)           #length(vt[c(-1)]) or length(vt)
 	g.cal[ ,1] = rowSums(group30[ ,gp.gn$het.f])  #F/HET 
 	g.cal[ ,2] = rowSums(group30[ ,gp.gn$het.m])  #M/HET
 	g.cal[ ,3] = rowSums(group30[ ,gp.gn$young])   #15-34/MSM,not include MWID
 	g.cal[ ,4] = rowSums(group30[ ,gp.gn$old])  #35-34/MSM,not include MWID
 	g.cal[ ,5] = rowSums(group30[ ,gp.gn$pwid])   #PWID/MWID
 	colnames(g.cal) <- c("Female:Het","Male:Het","Male:MSM/MWID 15-34","Male:MSM/MWID 35-64","PWID")
  	# g.cal = matrix(0,length(vt),2)           #length(vt[c(-1)]) or length(vt)
  	# g.cal[ ,1] = rowSums(group30[ ,gp.gn$young])   #15-34/MSM,not include MWID
  	# g.cal[ ,2] = rowSums(group30[ ,gp.gn$old])  #35-34/MSM,not include MWID
  	# colnames(g.cal) <- c("Young","Old")
  	return(g.cal)
  }
  n.strata.cal    <- group.calibration(n.strata)
  g3 <- plot_stacked_area(time = vt[-1], data=n.strata.cal , name= colnames(n.strata.cal), CITY) #time = vt[-1] or vt
  
  # Race-level strata
  group.race = function(group30){
  	g.race = matrix(0,length(vt[c(-1)]),3)   #length(vt[c(-1)]) or length(vt)
  	g.race[ ,1] = rowSums(group30[ ,gp.gn$white])  #
  	g.race[ ,2] = rowSums(group30[ ,gp.gn$hisp])   #
  	g.race[ ,3] = rowSums(group30[ ,gp.gn$black])  #
  	colnames(g.race) <- c("White/Other","Hispanic","Black") # In this order for Black to be on top
  	return(g.race)
  }
  n.strata.race    <- group.race(n.strata)
  g4 <- plot_stacked_area(time = vt[-1], data=n.strata.race, name = colnames(n.strata.race), CITY) #time = vt[-1] or vt
  
  
  ## Compare Pop Proportion across states
  # All 24 states
  n.state    <- apply(outa[ , , ], c(1,2), sum) # across 2 dimension: Time x states
  n.state.16    <- n.state[,1:16]
# n.state.16    <- n.state[,1:19] #LEm 1.0
	#
  df.state   <- as.data.frame(cbind(vt[-1],n.state.16)) #time = vt[-1] or vt
#  df.state   <- as.data.frame(cbind(vt,n.state.16))
  colnames(df.state)    <- c("time",vparameters$state.name[1:16])
  
  # All states except S
  g5 <- plot_stacked_area(
  		time = vt,  #time = vt[-1] or vt
# 		data = df.state[, c("S1","S2","Ia","I1", "I2","I3","Da","D1","D2","D3","T1","T2","T3","O1","O2","O3")],  # << pick 3 columns (S1, S2, Ia)
#  		names.gp = c("S1","S2","Ia","I1", "I2","I3","Iap","Ip","Da","D1","D2","D3","T1","T2","T3","O1","O2","O3"),
  		data = df.state[, c("Ia","I1", "I2","I3")],  
  		name = c("Ia","I1", "I2","I3"),
  		CITY, start_year=2017)
  
  # Race-level strata, except S
  group.state = function(group30){
  	g.state = matrix(0,length(vt[c(-1)]),4) #time = vt[-1] or vt
#  	g.state = matrix(0,length(vt),4)
  	g.state[ ,1] = rowSums(df.state[ ,c("Ia","I1","I2","I3")])  
  	g.state[ ,2] = rowSums(df.state[ ,c("Da","D1","D2","D3")])  
  	g.state[ ,3] = rowSums(df.state[ ,c("T1","T2","T3")]) 
  	g.state[ ,4] = rowSums(df.state[ ,c("O1","O2","O3")])
  	colnames(g.state) <- c("Infected, undiagnosed","Diagnosed, not on ART","On ART","Off ART")
  	return(g.state)
  }
  n.state.4    <- group.state(df.state)
#  g6 <- plot_stacked_area(time = vt, data = n.state.4, name = colnames(n.state.4), CITY)
 g6 <- plot_stacked_area(time = vt, data = n.state.4, name = colnames(n.state.4), CITY, start_year = 2017, xmax =144)  #time = vt[-1] or vt, remove xmax to see full time period

 # Race-level strata, among I, D, T
 df.state[,c("Ia","I1","I2","I3")]
 g6.I <- plot_stacked_area(time = vt[-1], data = df.state[,c("Ia","I1","I2","I3")], name = colnames(df.state[,c("Ia","I1","I2","I3")]), CITY)
# g6.I <- plot_stacked_area(time = vt[-1], data = df.state[,c("Ia","I1","I2","I3")], name = colnames(df.state[,c("Ia","I1","I2","I3")]), CITY, start_year = 2017, xmax =144)
 g6.D <- plot_stacked_area(time = vt[-1], data = df.state[,c("Da","D1","D2","D3")], name = colnames(df.state[,c("Da","D1","D2","D3")]), CITY)
 g6.T <- plot_stacked_area(time = vt[-1], data = df.state[,c("T1","T2","T3")], name = colnames(df.state[,c("T1","T2","T3")]), CITY)
 
 ## Compare Pop Proportion across states (For I only for new diagnoses)
	# Het only 
	# str(outa[ ,c("Ia","I1","I2","I3"), gp.gn$het])
  n.strata   <- apply(outa[ ,c("Ia","I1","I2","I3"), gp.gn$pwid], c(1,3), sum) # across 2 dimension: Time x strata
  plot_stacked_area(time = vt, data = n.strata, name = colnames(n.strata), CITY, start_year = 2017, xmax = 72)
  n.state    <- apply(outa[ ,c("Ia","I1","I2","I3"), gp.gn$pwid], c(1,2), sum) # across 2 dimension: Time x states
  plot_stacked_area(time = vt, data = n.state, name = colnames(n.state), CITY, start_year = 2017, xmax = 72)
  # str(outa[ ,c("Ia","I1","I2","I3"), gp.gn$all.msm])
  n.strata   <- outa[ ,c("death"), ] # across 2 dimension: Time x strata
  plot_stacked_area(time = vt, data = n.strata, name = colnames(n.strata), CITY, start_year = 2017)
  
#######Line Plot####  		

  # 
  plot_line <- function(data, name, TITLE = na, CITY = CITY, start_year = 2017, y_min = NULL) {
  	
  	# Convert to dataframe and set column names
  	df <- as.data.frame(data)
  	colnames(df) <- name
  	
  	# Add Year column
  	df$Year <- seq(from = start_year, by = 1, length.out = nrow(df))
  	
  	# Pivot to long format
  	df_long <- df %>%
  		pivot_longer(
  			cols = name,
  			names_to = "Legend",
  			values_to = "Value"
  		)
  	
  	# Subset the labels "value" at min and max time point
  	label_data <- df_long %>%	filter(Year == min(Year) | Year == max(Year))
  	
  	# Generate plot
  	g <- ggplot(df_long, aes(x = Year, y = Value, color = Legend)) +
  		geom_line(size = 1.2) +
  		geom_text(
  			data = label_data,
  			aes(label = round(Value, 1)),  # label = value rounded to 1 decimal
  			vjust = -0.5,
  			size = 3.5,
  			show.legend = FALSE) +
  		theme_minimal(base_size = 14) +
  		labs(
  			title = paste(TITLE ," ", CITY),
  			x = "Year",
  			y = TITLE,
  			color = "Legend"
  		) +
  		scale_color_brewer(palette = "Set1") +
  		theme(
  			legend.position = "right",
  			plot.title = element_text(hjust = 0.5)
  		)
  	
  	# Add y-axis minimum if specified
  	if (!is.null(y_min)) {
  		g <- g + expand_limits(y = y_min)
  	}
  	
  	print(g)
  }
  
  ### Target : New infection
  # Get yearly number of infections 
  n.inf        <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$white ]), 2, diff)
  n.inf.w      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$white ]), c(1), sum) #inc_bo,inc_bs,inc_g
  n.inf.w      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$white ]),diff(n.inf.w))
  n.inf.b      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$black ]), c(1), sum)
  n.inf.b      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$black ]),diff(n.inf.b))
  n.inf.h      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$hisp ]), c(1), sum)
  n.inf.h      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$hisp ]),diff(n.inf.h))
  n.inf        <- cbind(n.inf.w, n.inf.b, n.inf.h)
  plot_line(data=n.inf, name = c("White","Black","Hisp"), TITLE = "Number of infections", CITY = CITY )
# plot_line(data=n.inf, name = c("White","Black","Hisp"), TITLE = "Number of infections", CITY = CITY )
  #LEM 1.0
  n.inf.msm      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$msm ]), c(1), sum) #inc_bo,inc_bs,inc_g
  n.inf.msm      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$msm ]),diff(n.inf.msm))
  n.inf.mwid      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$mwid]), c(1), sum) #inc_bo,inc_bs,inc_g
  n.inf.mwid      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$mwid]),diff(n.inf.mwid))
  n.inf.pwid      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$pwid ]), c(1), sum)
  n.inf.pwid      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$pwid]),diff(n.inf.pwid))
  n.inf.het      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$het ]), c(1), sum)
  n.inf.het      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$het ]),diff(n.inf.het))
  n.inf.risk        <- cbind(n.inf.msm, n.inf.mwid, n.inf.pwid, n.inf.het)
  plot_line(data=n.inf.risk , name = c("msm","mwid","pwid","HET"), TITLE = "Number of new infections", CITY = CITY) 
  #LEM 2.0
  n.inf.Young      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$young ]), c(1), sum) #inc_bo,inc_bs,inc_g
  n.inf.Young      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$young ]),diff(n.inf.Young))
  n.inf.Old      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$old]), c(1), sum) #inc_bo,inc_bs,inc_g
  n.inf.Old      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$old]),diff(n.inf.Old))
  n.inf.PWID      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$pwid ]), c(1), sum)
  n.inf.PWID      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$pwid ]),diff(n.inf.PWID))
  n.inf.het      <- apply((outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"),gp.gn$het ]), c(1), sum)
  n.inf.het      <- c(sum(outa[12, c("inc_bo","inc_bs","inc_g"),gp.gn$het ]),diff(n.inf.het))
  n.inf.risk        <- cbind(n.inf.Young, n.inf.Old, n.inf.PWID, n.inf.het)
  plot_line(data=n.inf.risk , name = c("Young","Old","PWID","HET"), TITLE = "Number of new infections", CITY = CITY)
  
  # Get yearly Rate of infections per 100,000
  #LEM 1.0
  r.inf.w <- n.inf.w/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$white], c(1), sum) *100000
  r.inf.b <- n.inf.b/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$black], c(1), sum) *100000
  r.inf.h <- n.inf.h/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$hisp], c(1), sum) *100000
  #LEM 2.0
  r.inf.w <- n.inf.w/apply(outa[end_yr_ind, c("S1","S2"), gp.gn$white], c(1), sum) *100000
  r.inf.b <- n.inf.b/apply(outa[end_yr_ind, c("S1","S2"), gp.gn$black], c(1), sum) *100000
  r.inf.h <- n.inf.h/apply(outa[end_yr_ind, c("S1","S2"), gp.gn$hisp], c(1), sum) *100000
  #
  r.inf   <- cbind(r.inf.w, r.inf.b, r.inf.h)
  plot_line(data=r.inf, name = c("White","Black","Hisp"), TITLE = "Rate of infections by 100,000", CITY = CITY)
# plot_line(data=r.inf, name = c("White","Black","Hisp"), TITLE = "Rate of infections by 100,000", CITY = CITY, start_year =2012)
  #LEM 1.0
  r.inf.msm  <- n.inf.w/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$msm], c(1), sum) *100000
  r.inf.mwid <- n.inf.b/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$mwid], c(1), sum) *100000
  r.inf.pwid <- n.inf.h/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$pwid], c(1), sum) *100000
  r.inf.het  <- n.inf.h/apply(outa[end_yr_ind, c("S1","S2","Sp"), gp.gn$het], c(1), sum) *100000
  r.inf   <- cbind(r.inf.msm, r.inf.mwid, r.inf.pwid, r.inf.het)
  plot_line(data=r.inf, name = c("msm","mwid","pwid","HET"), TITLE = "Rate of infections by 100,000", CITY = CITY)
  #LEM 2.0
  r.inf.Young   <- n.inf.Young  /apply(outa[end_yr_ind, c("S1","S2"), gp.gn$young], c(1), sum) *100000
  r.inf.Old   <- n.inf.Old  /apply(outa[end_yr_ind, c("S1","S2"), gp.gn$old], c(1), sum) *100000
  r.inf.PWID <- n.inf.PWID/apply(outa[end_yr_ind, c("S1","S2"), gp.gn$pwid], c(1), sum) *100000
  r.inf.het    <- n.inf.het  /apply(outa[end_yr_ind, c("S1","S2"), gp.gn$het], c(1), sum) *100000
  r.inf.risk   <- cbind(r.inf.Young, r.inf.Old, r.inf.PWID, r.inf.het)
  plot_line(data=r.inf.risk , name = c("Young","Old","PWID","HET"), TITLE = "Rate of infections by 100,000", CITY = CITY)
  
  # Opposite sex infection
  inc_bo     =  rbind(outa[12, "inc_bo", ],  apply(outa[end_yr_ind, "inc_bo", ], 2, diff)) 
  inc_bo.Young       = rowSums(inc_bo[ , gp.gn$young])
  inc_bo.Old       = rowSums(inc_bo[ , gp.gn$old])
  inc_bo.PWID     = rowSums(inc_bo[ , gp.gn$pwid])
  inc_bo.het        = rowSums(inc_bo[ , gp.gn$het])
  inc_bo.risk    = cbind(inc_bo.Young,inc_bo.Old ,inc_bo.PWID, inc_bo.het)
  plot_line(data=inc_bo.risk, name = c("Young","Old","PWID","HET"), TITLE = "inc_bo.risk", CITY = CITY)
  # Same sex infection
  inc_bs     =  rbind(outa[12, "inc_bs", ],  apply(outa[end_yr_ind, "inc_bs", ], 2, diff)) 
  inc_bs.Young       = rowSums(inc_bs[ , gp.gn$young])
  inc_bs.Old       = rowSums(inc_bs[ , gp.gn$old])
  inc_bs.PWID     = rowSums(inc_bs[ , gp.gn$pwid])
  inc_bs.het        = rowSums(inc_bs[ , gp.gn$het])
  inc_bs.risk    = cbind(inc_bs.Young,inc_bs.Old ,inc_bs.PWID, inc_bs.het)
  plot_line(data=inc_bs.risk, name = c("Young","Old","PWID","HET"), TITLE = "inc_bs.risk", CITY = CITY) 
  # Shared needle infection
  inc_g     =  rbind(outa[12, "inc_g", ],  apply(outa[end_yr_ind, "inc_g", ], 2, diff)) 
  inc_g.Young       = rowSums(inc_g[ , gp.gn$young])
  inc_g.Old       = rowSums(inc_g[ , gp.gn$old])
  inc_g.PWID     = rowSums(inc_g[ , gp.gn$pwid])
  inc_g.het        = rowSums(inc_g[ , gp.gn$het])
  inc_g.risk    = cbind(inc_g.Young,inc_g.Old ,inc_g.PWID, inc_g.het)
  plot_line(data=inc_g.risk, name = c("Young","Old","PWID","HET"), TITLE = "inc_g.risk", CITY = CITY) 
  
  # Target 1: Total Diagnosed PLHIV
#  end_yr_ind   = c(c(12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,300,312,324,336)+1)
#  end_yr_ind_plus   = c(1,end_yr_ind)
  diag         = apply(outa[end_yr_ind, 7:16, ], c(1,3), sum)
#  diag         = apply(outa[, 7:16, ], c(1,3), sum)
  diag.all     = rowSums(diag)
  # 3 strata
  diag.w       = rowSums(diag[ , gp.gn$white])
  diag.b       = rowSums(diag[ , gp.gn$black])
  diag.h       = rowSums(diag[ , gp.gn$hisp])
  diag.race    = cbind(diag.w,diag.b,diag.h)
  plot_line(data= diag.race, name = c("White","Black","Hisp"), TITLE = "Number of Total diagnoses", CITY = CITY)
  # Target
  vparameters$pop.diag # First year 2016
  source("Model_development/7.Red_strata/Calibration/7.load_targets.R")
  calib.target$diag.obs # Remaining years 2016
 
  new_row <- data.frame(
  	W = sum(vparameters$pop.diag[gp18.gn$white]),
  	B = sum(vparameters$pop.diag[gp18.gn$black]),
  	H = sum(vparameters$pop.diag[gp18.gn$hisp])
  )
  diag.obs <- rbind(new_row, calib.target$diag.obs)
  rownames(diag.obs) <- 2016:(2016 + nrow(calib.target$diag.obs))
  plot_line(data= cbind(diag.race[c(1:7),1],diag.obs[,1]), name = c("White-Model","White-Target"), TITLE = "Number of Total diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(diag.race[c(1:7),2],diag.obs[,2]), name = c("Black-Model","Black-Target"), TITLE = "Number of Total diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(diag.race[c(1:7),3],diag.obs[,3]), name = c("Hisp-Model","Hisp-Target"), TITLE = "Number of Total diagnoses", CITY = CITY, y_min=0)
  #
  diag.Young   = rowSums(diag[ , gp.gn$young])
  diag.Old     = rowSums(diag[ , gp.gn$old])
  diag.PWID    = rowSums(diag[ , gp.gn$pwid])
  diag.het.m    = rowSums(diag[ , gp.gn$het.m])
  diag.het.f    = rowSums(diag[ , gp.gn$het.f])
  diag.risk   = cbind(diag.Young, diag.Old, diag.PWID, diag.het.m, diag.het.f)
  plot_line(data= diag.risk, name = c("Young","Old","PWID","HET.m","HET.f"), TITLE = "Number of Total diagnoses", CITY = CITY)
  
  ### Target 2: New diagnosed among infected
# ndiag        = rbind(outa[12, "diag", ],  apply(outa[end_yr_ind, "diag", ], 2, diff))
  ndiag        = rbind(outa[12, "diag", ],  apply(outa[end_yr_ind_calibration, "diag", ], 2, diff))
  ndiag.w      = rowSums(ndiag[ , gp.gn$white])
  ndiag.b      = rowSums(ndiag[ , gp.gn$black])
  ndiag.h      = rowSums(ndiag[ , gp.gn$hisp])
  ndiag.race    = cbind(ndiag.w,ndiag.b,ndiag.h)
  plot_line(data=ndiag.race, name = c("White","Black","Hisp"), TITLE = "Number of new diagnoses", CITY = CITY)
  # Target
  ndiag.obs     = calib.target$ndiag.obs 
#  ndiag.obs <- cbind(
#  	ndiag.obs,
#  	Young = ndiag.obs[, "X15.34.MSM"] + ndiag.obs[, "X15.34.MWID"],
#  	Old   = ndiag.obs[, "X35.64.MSM"] + ndiag.obs[, "X35.64.MWID"]
#  ) 	# add new cols
#  ndiag.obs <- ndiag.obs[, c("M.HET", "F.HET", "Young", "Old", "PWID", "W", "B", "H")]	# drop originals MSM & MWID
  
  plot_line(data= cbind(ndiag.race[c(1:6),"ndiag.w"],ndiag.obs[,"W"]), name = c("White-Model","White-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.race[c(1:6),"ndiag.b"],ndiag.obs[,"B"]), name = c("Black-Model","Black-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.race[c(1:6),"ndiag.h"],ndiag.obs[,"H"]), name = c("Hisp-Model","Hisp-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  #
  ndiag.Young   = rowSums(ndiag[ , gp.gn$young])
  ndiag.Old     = rowSums(ndiag[ , gp.gn$old])
  ndiag.PWID    = rowSums(ndiag[ , gp.gn$pwid])
  ndiag.het.m    = rowSums(ndiag[ , gp.gn$het.m])
  ndiag.het.f    = rowSums(ndiag[ , gp.gn$het.f])
  ndiag.het    = rowSums(ndiag[ , gp.gn$het])
  ndiag.risk   = cbind(ndiag.Young,ndiag.Old,ndiag.PWID,ndiag.het.m,ndiag.het.f,ndiag.het)
  plot_line(data=ndiag.risk, name = c("Young","Old","PWID","HET.m","HET.f","HET"), TITLE = "Number of new diagnoses", CITY = CITY)
  # Target
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.Young"]  ,ndiag.obs[,3]), name = c("Young-Model","Young-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.Old"]  ,ndiag.obs[,4]), name = c("Old-Model","Old-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.PWID"],ndiag.obs[,5]), name = c("PWID-Model","PWID-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.het.m"] ,ndiag.obs[,1]), name = c("Het.m-Model","Het.m-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.het.f"] ,ndiag.obs[,2]), name = c("Het.f-Model","Het.f-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  plot_line(data= cbind(ndiag.risk[c(1:6),"ndiag.het"] ,rowSums(ndiag.obs[,c(1:2)])), name = c("Het-Model","Het-Target"), TITLE = "Number of New diagnoses", CITY = CITY, y_min=0)
  # Get yearly Rate of new diagnoses per 100,000
  r.ndiag.w    = ndiag.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.ndiag.b    = ndiag.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.ndiag.h    = ndiag.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.ndiag   <- cbind(r.ndiag.w, r.ndiag.b, r.ndiag.h)
  plot_line(data=r.ndiag, name = c("White","Black","Hisp"), TITLE = "Rate of new diagnoses by 100,000", CITY = CITY)
  r.ndiag.Young    = ndiag.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.ndiag.Old    = ndiag.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.ndiag.PWID  = ndiag.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.ndiag.het     = ndiag.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.ndiag.risk   <- cbind(r.ndiag.Young, r.ndiag.Old, r.ndiag.PWID,r.ndiag.het)
  plot_line(data=r.ndiag.risk, name = c("Young","Old","PWID","HET"), TITLE = "Rate of new diagnoses by 100,000", CITY = CITY)
  
  ### Target : Normal initiation
  D_T        = rbind(outa[12, "D_T", ],  apply(outa[end_yr_ind, "D_T", ], 2, diff))
  D_T.w      = rowSums(D_T[ , gp.gn$white])
  D_T.b      = rowSums(D_T[ , gp.gn$black])
  D_T.h      = rowSums(D_T[ , gp.gn$hisp])
  D_T.race    = cbind(D_T.w,D_T.b,D_T.h)
  plot_line(data = D_T.race, name = c("White","Black","Hisp"), TITLE = "Number of D_T", CITY = CITY)
# plot_line(data = D_T.race, name = c("White","Black","Hisp"), TITLE = "Number of D_T", CITY = CITY, start_year =2012) #LEM 1.0  
  D_T.Young   = rowSums(D_T[ , gp.gn$young])
  D_T.Old   = rowSums(D_T[ , gp.gn$old])
  D_T.PWID = rowSums(D_T[ , gp.gn$pwid])
  D_T.het    = rowSums(D_T[ , gp.gn$het])
  D_T.risk   = cbind(D_T.Young,D_T.Old,D_T.PWID,D_T.het)
  plot_line(data=D_T.risk, name = c("Young","Old","PWID","HET"), TITLE = "Number of D_T", CITY = CITY)
  # Get yearly Rate of new diagnoses per 100,000
  r.D_T.w    = D_T.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.D_T.b    = D_T.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.D_T.h    = D_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.D_T     <- cbind(r.D_T.w, r.D_T.b, r.D_T.h)
  plot_line(data=r.D_T, name = c("White","Black","Hisp"), TITLE = "Rate of D_T by 100,000", CITY = CITY)
  r.D_T.Young    = D_T.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.D_T.Old    = D_T.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.D_T.PWID  = D_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.D_T.het     = D_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.D_T.risk   <- cbind(r.D_T.Young, r.D_T.Old, r.D_T.PWID,r.D_T.het)
  plot_line(data=r.D_T.risk, name = c("Young","Old","PWID","HET"), TITLE = "Rate of D_T by 100,000", CITY = CITY)
  
  ### Target : Immediate initiation
  I_T        = rbind(outa[12, "I_T", ],  apply(outa[end_yr_ind, "I_T", ], 2, diff))
  I_T.w      = rowSums(I_T[ , gp.gn$white])
  I_T.b      = rowSums(I_T[ , gp.gn$black])
  I_T.h      = rowSums(I_T[ , gp.gn$hisp])
  I_T.race    = cbind(I_T.w,I_T.b,I_T.h)
  plot_line(data = I_T.race, name = c("White","Black","Hisp"), TITLE = "Number of Immediate initation", CITY = CITY)
  I_T.Young   = rowSums(I_T[ , gp.gn$young])
  I_T.Old   = rowSums(I_T[ , gp.gn$old])
  I_T.PWID = rowSums(I_T[ , gp.gn$pwid])
  I_T.het    = rowSums(I_T[ , gp.gn$het])
  I_T.risk   = cbind(I_T.Young,I_T.Old,I_T.PWID,I_T.het)
  plot_line(data=I_T.risk, name = c("Young","Old","PWID","HET"), TITLE = "Number of Immediate initation", CITY = CITY)
  # Get yearly Rate of new diagnoses per 100,000
  r.I_T.w    = I_T.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.I_T.b    = I_T.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.I_T.h    = I_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.I_T     <- cbind(r.I_T.w, r.I_T.b, r.I_T.h)
  plot_line(data=r.I_T, name = c("White","Black","Hisp"), TITLE = "Rate of Immediate initation by 100,000", CITY = CITY)
  r.I_T.Young    = I_T.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.I_T.Old    = I_T.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.I_T.PWID  = I_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.I_T.het     = I_T.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.I_T.risk   <- cbind(r.I_T.Young, r.I_T.Old, r.I_T.PWID,r.I_T.het)
  plot_line(data=r.I_T.risk, name = c("Young","Old","PWID","HET"), TITLE = "Rate of Immediate initation by 100,000", CITY = CITY)
  
  ### Target : I to D
  I_D        = rbind(outa[12, "I_D", ],  apply(outa[end_yr_ind, "I_D", ], 2, diff))
  I_D.w      = rowSums(I_D[ , gp.gn$white])
  I_D.b      = rowSums(I_D[ , gp.gn$black])
  I_D.h      = rowSums(I_D[ , gp.gn$hisp])
  I_D.race    = cbind(I_D.w,I_D.b,I_D.h)
  plot_line(data = I_D.race, name = c("White","Black","Hisp"), TITLE = "Number of I to D", CITY = CITY)
# plot_line(data = I_D.race, name = c("White","Black","Hisp"), TITLE = "Number of I to D", CITY = CITY, start_year = 2012) #LEM 1.0
  I_D.Young   = rowSums(I_D[ , gp.gn$young]) #gp.gn$Y.msm
  I_D.Old     = rowSums(I_D[ , gp.gn$old])   #gp.gn$O.msm
  I_D.PWID   = rowSums(I_D[ , gp.gn$pwid]) #gp.gn$all.pwid
  I_D.het    = rowSums(I_D[ , gp.gn$het])
  I_D.risk   = cbind(I_D.Young,I_D.Old,I_D.PWID,I_D.het)
  plot_line(data=I_D.risk, name = c("Young","Old","PWID","HET"), TITLE = "Number of I to D", CITY = CITY)
  # Get yearly Rate of new diagnoses per 100,000
  r.I_D.w    = I_D.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.I_D.b    = I_D.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.I_D.h    = I_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.I_D     <- cbind(r.I_D.w, r.I_D.b, r.I_D.h)
  plot_line(data=r.I_D, name = c("White","Black","Hisp"), TITLE = "Rate of I to D by 100,000", CITY = CITY)
  r.I_D.Young   = I_D.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.I_D.Old     = I_D.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.I_D.PWID    = I_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.I_D.het     = I_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.I_D.risk   <- cbind(r.I_D.young, r.I_D.Old, r.I_D.PWID,r.I_D.het)
  plot_line(data=r.I_D.risk, name = c("young","Old","PWID","HET"), TITLE = "Rate of I to D by 100,000", CITY = CITY)
  
  ### Target 3: Death among infected
  death        = rbind(outa[12, "death", ],  apply(outa[end_yr_ind, "death", ], 2, diff))
  death.w      = rowSums(death[ , gp.gn$white])
  death.b      = rowSums(death[ , gp.gn$black])
  death.h      = rowSums(death[ , gp.gn$hisp])
  death.race   = cbind(death.w,death.b,death.h)
  plot_line(data=death.race, name = c("White","Black","Hisp"), TITLE = "Number of death among PLHIV", CITY = CITY)
  #
  death.young  = rowSums(death[ , gp.gn$young])
  death.Old    = rowSums(death[ , gp.gn$old])
  death.PWID   = rowSums(death[ , gp.gn$pwid])
  death.het    = rowSums(death[ , gp.gn$het])
  death.risk   = cbind(death.Young,death.Old,death.PWID,death.het)
  plot_line(data=death.risk, name = c("Young","Old","PWID","HET"), TITLE = "Number of death among PLHIV", CITY = CITY)
  # Get yearly Rate of new diagnoses per 100,000
  r.death.w    = death.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.death.b    = death.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.death.h    = death.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.death   <- cbind(r.death.w, r.death.b, r.death.h)
  plot_line(data=r.death, name = c("White","Black","Hisp"), TITLE = "Rate of death among PLHIV by 100,000", CITY = CITY)
  r.death.Young   = death.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.death.Old     = death.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.death.PWID  = death.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.death.het     = death.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.death.risk   <- cbind(r.death.Young, r.death.Old, r.death.PWID,r.death.het)
  plot_line(data=r.death.risk, name = c("Young","Old","PWID","HET"), TITLE = "Rate of death among PLHIV by 100,000", CITY = CITY)
  
  ### Target 3: Death among D only
  death_D        = rbind(outa[12, "death_D", ],  apply(outa[end_yr_ind, "death_D", ], 2, diff))
  death_D.w      = rowSums(death_D[ , gp.gn$white])
  death_D.b      = rowSums(death_D[ , gp.gn$black])
  death_D.h      = rowSums(death_D[ , gp.gn$hisp])
  death_D.race   = cbind(death_D.w,death_D.b,death_D.h)
  plot_line(data=death_D.race, name = c("White","Black","Hisp"), TITLE = "Number of death_D among PLHIV", CITY = CITY)
# plot_line(data=death_D.race, name = c("White","Black","Hisp"), TITLE = "Number of death_D among PLHIV", CITY = CITY, start_year = 2012)
  death_D.Young  = rowSums(death_D[ , gp.gn$young])
  death_D.Old    = rowSums(death_D[ , gp.gn$old])
  death_D.PWID   = rowSums(death_D[ , gp.gn$pwid])
  death_D.het    = rowSums(death_D[ , gp.gn$het])
  death_D.risk   = cbind(death_D.Young,death_D.Old,death_D.PWID,death_D.het)
  plot_line(data=death_D.risk, name = c("Young","Old","PWID","HET"), TITLE = "Number of death_D among PLHIV", CITY = CITY)
  # Get yearly Rate of new diagnoses per 100,000
  r.death_D.w    = death_D.w/apply(outa[end_yr_ind, 1:16, gp.gn$white], c(1), sum) *100000
  r.death_D.b    = death_D.b/apply(outa[end_yr_ind, 1:16, gp.gn$black], c(1), sum) *100000
  r.death_D.h    = death_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$hisp ], c(1), sum) *100000
  r.death_D   <- cbind(r.death_D.w, r.death_D.b, r.death_D.h)
  plot_line(data=r.death_D, name = c("White","Black","Hisp"), TITLE = "Rate of death_D among PLHIV by 100,000", CITY = CITY)
  r.death_D.Young   = death_D.w/apply(outa[end_yr_ind, 1:16, gp.gn$young], c(1), sum) *100000
  r.death_D.Old     = death_D.b/apply(outa[end_yr_ind, 1:16, gp.gn$old], c(1), sum) *100000
  r.death_D.PWID    = death_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$pwid], c(1), sum) *100000
  r.death_D.het     = death_D.h/apply(outa[end_yr_ind, 1:16, gp.gn$het], c(1), sum) *100000
  r.death_D.risk   <- cbind(r.death_D.Young, r.death_D.Old, r.death_D.PWID,r.death_D.het)
  plot_line(data=r.death_D.risk, name = c("Young","Old","PWID","HET"), TITLE = "Rate of death_D among PLHIV by 100,000", CITY = CITY)
  
  # Target 4: (Proportion) Yearly late diagnoses: diag_late <- psi[6] * I3 / diagm: psi[3] * Ia + psi[4] * I1 + psi[5] * I2 + psi[6] * I3 (Target)
  ldiag  <- rbind(outa[12, "diag_late", ],  apply(outa[end_yr_ind_calibration, "diag_late", ], 2, diff)) # Exactly same as New diag (only that diag late is for I3 only), then divided by ndiag
  p.diag_late.all = ldiag/ndiag
  
  p.diag_late.w   = rowSums(ldiag[, gp.gn$white])/rowSums(ndiag[, gp.gn$white])
  p.diag_late.b   = rowSums(ldiag[, gp.gn$black])/rowSums(ndiag[, gp.gn$black])
  p.diag_late.h   = rowSums(ldiag[, gp.gn$hisp ])/rowSums(ndiag[, gp.gn$hisp ])
  p.diag_late.cal.model   = cbind(p.diag_late.w,p.diag_late.b,p.diag_late.h)
  
  W_diag_late <- plot_line(cbind(p.diag_late.cal.model[1:6,1], calib.target$late_diag.obs[,"W"]),
  												 c("Sens","Target"),
  												 TITLE="late diagnoses W", CITY=CITY, y_min=0) + theme(legend.position = "none")        
  B_diag_late <- plot_line(cbind(p.diag_late.cal.model[1:6,2], calib.target$late_diag.obs[,"B"]),
  												 c("Sens","Target"),
  												 TITLE="late diagnoses B", CITY=CITY, y_min=0) + theme(legend.position = "none")
  H_diag_late <- plot_line(cbind(p.diag_late.cal.model[1:6,3], calib.target$late_diag.obs[,"H"]),
  												 c("Sens","Target"),
  												 TITLE="late diagnoses H", CITY=CITY, y_min=0) + theme(legend.position = "none")
  
  # Target 5: (Proportion) Yearly Viral suppression 2017 to 2021
  p.viral_sup.all = apply(outa[end_yr_ind_calibration, c(11,12), ], c(1,3), sum) / apply(outa[end_yr_ind_calibration, c(7:16), ], c(1,3), sum) #T1, T2
  p.viral_sup.w   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$white], c(1,3), sum))/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$white], c(1,3), sum))
  p.viral_sup.b   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$black], c(1,3), sum))/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$black], c(1,3), sum))
  p.viral_sup.h   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$hisp ], c(1,3), sum))/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$hisp ], c(1,3), sum))
  p.viral_sup.cal.model   = cbind(p.viral_sup.w,p.viral_sup.b,p.viral_sup.h)
  ### Check First 12 months
  p.viral_sup.w   = rowSums(apply(outa[c(1:12), c(11,12), gp.gn$white], c(1,3), sum))/rowSums(apply(outa[c(1:12), c(7:16), gp.gn$white], c(1,3), sum))
  p.viral_sup.b   = rowSums(apply(outa[c(1:12), c(11,12), gp.gn$black], c(1,3), sum))/rowSums(apply(outa[c(1:12), c(7:16), gp.gn$black], c(1,3), sum))
  p.viral_sup.h   = rowSums(apply(outa[c(1:12), c(11,12), gp.gn$hisp ], c(1,3), sum))/rowSums(apply(outa[c(1:12), c(7:16), gp.gn$hisp ], c(1,3), sum))
  