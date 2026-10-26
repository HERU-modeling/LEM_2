#Step 0: Clear R environment
rm(list = ls())
#Step 1: Importing R functions & set seeds
library(openxlsx)
library(mc2d)
library(posterior)


######### 1. Set up ############
#Step 2: Choosing city
CITY <- "SHE" # MIA1, NYC1 for original, KIN, BAL, PHI


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
calibration_results_multiple1 <- readRDS(paste0("Calibration/3_Results/calib_fix_scen2_five-seed_job3917938_combined_156a386_2000-2000_SHE.rds"))
calibration_results_multiple2 <- readRDS(paste0("Calibration/3_Results/calib_fix_scen2_reduced-h-nsg_job4283360_combined_156a386_2000-2000_SHE.rds"))
get_calpar_info <- readRDS(paste0("Calibration/3_Results/calib_All_10000-5000-jun23_SHE.rds"))
calib_param <- get_calpar_info$calib_param
#Re-Run toget meta data
param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
calpar <- param_info$calpar
calpar.info <- param_info$calpar.info


# posterior_samples <- resample_using_importance_weights(calibration_results_multiple$burned_samples,
# 																											 calibration_results_multiple$burned_lls,
# 																											 ndraws = 1000)
posterior_samples_1 <- calibration_results_multiple1$posterior_samples_combined
posterior_samples_2 <- calibration_results_multiple2$posterior_samples_combined
prior_samples <- sample_dist_wrapper(2000, calpar)

# Examine the PE, Upper, Lower & 
colnames_all_samples <- paste0(calpar$par, "_", calpar$gender, "_", calpar$ethnicity, "_",
															 calpar$risk, "_", calpar$sex.risk.age.group)
colnames(posterior_samples_1) <- c(colnames_all_samples)
colnames(posterior_samples_2) <- c(colnames_all_samples)
colnames(prior_samples)     <- c(colnames_all_samples)

posterior_stats_1 <- descriptive_stats(posterior_samples_1)
posterior_stats_2 <- descriptive_stats(posterior_samples_2)
prior_stats     <- descriptive_stats(prior_samples)

posterior_stats_1[c(15:20,34:36),]
posterior_stats_2[c(15:20,34:36),]
prior_stats[c(15:20,34:36),]
#posterior_stats[posterior_stats$name == "noG_m_b_mwid_Young", ]
#prior_stats[prior_stats$name == "psi_f_b_het_na", ]
#print(posterior_stats[grepl("^prop.HIV.aware", posterior_stats$name), ], n = Inf)
#print(prior_stats[grepl("^noG", posterior_stats$name), ], n = Inf)
# posterior_stats_long <- expand_posterior_stats(posterior_stats)

# Generate plot
# plot_prior_posterior_comparison(prior_samples, posterior_samples, calpar, adjust = 2)
plot_prior_posterior_comparison_2(
  prior_samples,
  posterior_samples_1,
  posterior_samples_2 = posterior_samples_2,
  calpar,
  par_selected = c("nsG","psi"), Jurisdiction_names = CITY,
  adjust =10
)
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
  
  source("Calibration/6_Others/Utility_functions.R")
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
  source("Calibration/6_Others/Utility_functions.R")
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
  
  ####### Examine: pop.change.30  ####  
  ###############1. Figure - 1 strata Panel #########
  
  # Extract Seed
calibration_results <- readRDS(paste0("Calibration/3_Results/calib_prop.adj_TRUE_1234_714fdc9_100-100_SHE.rds"))
  # str(calibration_results$pop.change.30): 1000 samples * 348 times * 30 strata
  pop.change.30 <-  aperm(calibration_results$pop.change.30, c(1,3,2))

plot_strata_ci <- function(arr, strata_name) {
  
  strata_idx <- which(dimnames(arr)$group == strata_name)
  
  tmp <- arr[, , strata_idx]
  
  df <- data.frame(
    month = seq_len(ncol(tmp)),
    median = apply(tmp, 2, median),
    lower  = apply(tmp, 2, quantile, 0.025),
    upper  = apply(tmp, 2, quantile, 0.975)
  )
#browser()  
  ggplot(df, aes(month, median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3) +
    geom_line(size = 0.8) +
    labs(
      title = strata_name,
      x = "Month",
      y = "Value"
    ) +
    theme_bw()
}
plot_strata_ci(
  calibration_results$pop.change.30/12,
  "MSM/young/black"
) 

###############1. Figure - 1 strata Panel #########


###############2. Figure - Multiple strata Panel #########
## Plot monthly population entrants (pop.change.30 / 12) by stratum:
## posterior median + 95% interval, one figure per city, 30 panels (3 columns x 10 rows).
## Run from the project root.

dir <- "prop.adj_FALSE1_1234_1775acb_100-100"
Jurisdictions  <- c("SHE")

for (j in seq_along(Jurisdictions)) {
  calibration_results <- readRDS(
    paste0(
      "Calibration/3_Results/calib_",
      dir, "_", Jurisdictions[j], ".rds"
    )
  )
  pc <-  calibration_results$pop.change.30 / 12
  # f <- list.files("Calibration/3_Results", full.names = TRUE,
  #                 pattern = paste0(run_tag, "_combined_.*_", city, "\\.rds$"))
  # 
  # pc <- readRDS(f)$pop.change.30 / 12           # [sample x month x group], monthly entrants
  
  q <- apply(pc, c(2, 3), quantile, probs = c(0.025, 0.5, 0.975))   # [quantile x month x group]
  groups <- dimnames(pc)$group
  df <- data.frame(
    month = rep(seq_len(dim(pc)[2]), times = length(groups)),
    group = factor(rep(groups, each = dim(pc)[2]), levels = groups),
    lo    = as.vector(q[1, , ]),
    med   = as.vector(q[2, , ]),
    hi    = as.vector(q[3, , ])
  )
 
  g <- ggplot(df, aes(month)) +
    geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2a78d6", alpha = 0.25) +
    geom_line(aes(y = med), colour = "#2a78d6", linewidth = 0.5) +
    facet_wrap(~ group, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 348, 60)) +
    labs(title    = paste0(Jurisdictions[j], ": monthly population entrants by stratum"),
         subtitle = paste0("Posterior median (line) and 95% interval (band), ", dim(pc)[1],
                           " draws. Negative = net exits."),
         x = "Month (1 = Jan 2017)", y = "People entering per month") +
    theme_minimal(base_size = 10) +
    theme(plot.background  = element_rect(fill = "#fcfcfb", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.25),
          strip.text       = element_text(colour = "#0b0b0b", face = "bold", hjust = 0),
          axis.text        = element_text(colour = "#898781", size = 7),
          axis.title       = element_text(colour = "#52514e"),
          plot.title       = element_text(colour = "#0b0b0b", face = "bold"),
          plot.subtitle    = element_text(colour = "#52514e"))
  
  out <- paste0("Calibration/5_Plots/pop_change_entrants_", Jurisdictions[j], ".png")
  ggsave(out, g, width = 9, height = 16, dpi = 200)
  cat("Saved", out, "\n")
  rm(pc); gc()
}

###############3. Figure - Stacked Number/Proportion #########
# pop.change sum across White:
pop.change.W <- apply(calibration_results$pop.change.30[,,gp.gn$white], c(1, 2), sum) # Same across sample & time
# Pop.change sum across MSM
pop.change.msm   <- apply(calibration_results$pop.change.30[,,gp.gn$msm], c(1, 2), sum)
pop.change.mwid  <- apply(calibration_results$pop.change.30[,,gp.gn$mwid], c(1, 2), sum)
pop.change.pwid.m<- apply(calibration_results$pop.change.30[,,gp.gn$pwid.m], c(1, 2), sum)
pop.change.het.m <- apply(calibration_results$pop.change.30[,,gp.gn$het.m], c(1, 2), sum)
pop.change.pwid.f<- apply(calibration_results$pop.change.30[,,gp.gn$pwid.f], c(1, 2), sum)
pop.change.het.f <- apply(calibration_results$pop.change.30[,,gp.gn$het.f], c(1, 2), sum)
# Check: pop.change.msm+pop.change.mwid+pop.change.pwid.m+pop.change.het.m +pop.change.pwid.f+pop.change.het.f
nsamp <- nrow(pop.change.msm)
nmonth <- ncol(pop.change.msm)
pop.change <- array(
  NA_real_,
  dim = c(nsamp, nmonth, 6),
  dimnames = list(
    sample = NULL,
    month  = colnames(pop.change.msm),
    group  = c("msm", "mwid", "pwid.m", "het.m", "pwid.f", "het.f")
  )
)

pop.change[, , "msm"]    <- pop.change.msm
pop.change[, , "mwid"]   <- pop.change.mwid
pop.change[, , "pwid.m"] <- pop.change.pwid.m
pop.change[, , "het.m"]  <- pop.change.het.m
pop.change[, , "pwid.f"] <- pop.change.pwid.f
pop.change[, , "het.f"]  <- pop.change.het.f

plot_strata_ci(pop.change,"het.f") 


library(ggplot2)

# median across samples
med  <- apply(pop.change, c(2, 3), median)
med <- apply(pop.change, c(2, 3), mean)
df <- data.frame(
  month = rep(seq_len(dim(med)[1]), times = dim(med)[2]),
  group = factor(
    rep(colnames(med), each = dim(med)[1]),
    levels = colnames(med)
  ),
  median = c(med)
)
browser()
ggplot(df,
       aes(x = month,
           y = median,
           fill = group)) +
  geom_area(position = "stack") +
  labs(
    x = "Month",
    y = "Population",
    fill = "Group"
  ) +
  theme_bw()

###############4. Calculate the MEAN proportion for pop.change.18 #########
# Run Prepare Environment from calibration 1 run first with file contains pop.change
# Create empty object before loop
all_pop_props <- data.frame()

for (j in seq_along(Jurisdictions)) {
  
  message("GEt Pop.change.18.prop ", Jurisdictions[j])
  
  ## ---- Load Jurisdiction-specific results ----
  calibration_results <- readRDS(
    paste0(
      "Calibration/3_Results/calib_",
      dir, "_", Jurisdictions[j], ".rds"
    )
  )
  
# prop_1:MSM/white
#      rowSums(calibration_results$pop.change.30[1,,intersect(gp.gn$white,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[1,,])  # values          across t of 1st sample set
# mean(rowSums(calibration_results$pop.change.30[1,,intersect(gp.gn$white,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[1,,])) # mean (of value) across t of 1st sample set; 1st value of  the next row
#      rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,])  # mean across 1000 sample sets
prop_1 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_2:MSM/black
prop_2 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_3:MSM/hisp 
prop_3 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,]))
##   mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,]))
##  +mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,]))
##  +mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,])) 
##  =mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$msm])/rowSums(calibration_results$pop.change.30[ ,,]))
  
# prop_4:MSM/PWID/white
prop_4 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$mwid)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_5:MSM/PWID/black
prop_5 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$mwid)])/rowSums(calibration_results$pop.change.30[ ,,]))  
# prop_6:MSM/PWID/hisp
prop_6 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$mwid)])/rowSums(calibration_results$pop.change.30[ ,,]))
##  =mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$mwid])/rowSums(calibration_results$pop.change.30[ ,,]))  
  
# prop_7:PWID/male/white
prop_7 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$pwid.m)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_8:PWID/male/black
prop_8 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$pwid.m)])/rowSums(calibration_results$pop.change.30[ ,,]))  
# prop_9:PWID/male/hisp
prop_9 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$pwid.m)])/rowSums(calibration_results$pop.change.30[ ,,]))  
  
# prop_10:Het/male/white
prop_10 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$het.m)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_11:Het/male/black
prop_11 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$het.m)])/rowSums(calibration_results$pop.change.30[ ,,]))  
# prop_12:Het/male/hisp
prop_12 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$het.m)])/rowSums(calibration_results$pop.change.30[ ,,]))  
## =mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$het.m])/rowSums(calibration_results$pop.change.30[ ,,])) 
  
# prop_13:PWID/female/white
prop_13 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$pwid.f)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_14:PWID/female/black
prop_14 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$pwid.f)])/rowSums(calibration_results$pop.change.30[ ,,]))  
# prop_15:PWID/female/hisp
prop_15 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$pwid.f)])/rowSums(calibration_results$pop.change.30[ ,,]))    
  
# prop_16:het/female/white
prop_16 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$white,gp.gn$het.f)])/rowSums(calibration_results$pop.change.30[ ,,]))
# prop_17:het/female/black
prop_17 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$het.f)])/rowSums(calibration_results$pop.change.30[ ,,]))  
# prop_18:het/female/hisp
prop_18 <- mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$hisp ,gp.gn$het.f)])/rowSums(calibration_results$pop.change.30[ ,,])) 

# Check:  
# all proportion add up to 1: 1 = mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$msm])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$mwid])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$pwid.m])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$pwid.f])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$het.m])/rowSums(calibration_results$pop.change.30[ ,,])) + (rowSums(calibration_results$pop.change.30[ ,,gp.gn$het.f])/rowSums(calibration_results$pop.change.30[ ,,])) 
# Add up proportion for male: mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$male])/rowSums(calibration_results$pop.change.30[ ,,])) = mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$msm])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$mwid])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$pwid.m])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$het.m])/rowSums(calibration_results$pop.change.30[ ,,]))
# Add up proportion for black:mean(rowSums(calibration_results$pop.change.30[ ,,gp.gn$black])/rowSums(calibration_results$pop.change.30[ ,,]))= mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$msm)])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$mwid)])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$pwid.m)])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$het.m)])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$pwid.f)])/rowSums(calibration_results$pop.change.30[ ,,])) + mean(rowSums(calibration_results$pop.change.30[ ,,intersect(gp.gn$black,gp.gn$het.f)])/rowSums(calibration_results$pop.change.30[ ,,]))
pop.change.18.prop <- c(prop_1,prop_2,prop_3,prop_4,prop_5,prop_6,prop_7,prop_8,prop_9,prop_10,prop_11,prop_12,prop_13,prop_14,prop_15,prop_16,prop_17,prop_18)


# add jurisdiction name as first column
all_pop_props <- rbind(
  all_pop_props,
  cbind(
    Jurisdiction = Jurisdictions[j],
    as.data.frame(t(pop.change.18.prop))
  )
)
}

# Set column names
colnames(all_pop_props) <- c(
  "Jurisdiction",
  paste0("prop_", 1:18)
)

# Write to Excel
library(openxlsx)

write.xlsx(
  all_pop_props,
  file = "Calibration/3_Results/pop_change_18_prop.xlsx",
  rowNames = FALSE
)
