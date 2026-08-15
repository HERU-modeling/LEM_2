rm(list = ls())
library(openxlsx)
library(mc2d)
library(posterior)
library(ggplot2)
library(patchwork)
library(cowplot)
set.seed(1111)

################ [Run this] Prepare Environment ###############
#Step 2: Choosing city
CITY <- "MAR" # MIA1, NYC1 for original, KIN, BAL, PHI -Pick any random cities, just need to get this running


#Step 3: Importing LEM related packages
# Importing baseline environment
# Importing baseline environment
source("ModelCoreModules/Parameter_input/LoadBaselineWorkspace.R")

#  Load vparameters
source("Calibration/2_Workflow/load_targets.R") #Update Target as needed - Part of par_pop.R - Need to move up here since diag18.obs is needed for initial pop
source("ModelCoreModules/Parameter_input/Module-Data.input.R")

# Load calibration related functions
source("Calibration/2_Workflow/util.R")
source("Calibration/2_Workflow/run_metropolis_hasting.R")
source("Calibration/2_Workflow/objective_fn.R")

# Load additional data not covered by current code
source("Calibration/2_Workflow/util_intervention.R")
source("Calibration/4_Visualization/plots.R")


#Step 6:
#Read-in save results
calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_check_af921fc_10-10_MAR.rds"))
calib_param <- calibration_results_multiple$calib_param
#Re-Run toget meta data
param_info <- add_param(names = calib_param, par_info = par_info, WB = WB, CITY = CITY) 
calpar <- param_info$calpar
calpar.info <- param_info$calpar.info

#SOme stuffs
target_name_extended <- c("diag.W","diag.B","diag.H","ndiag.M_HET","ndiag.F_HET","ndiag.Young","ndiag.Old","ndiag.PWID","ndiag.W","ndiag.B","ndiag.H",
													"p.diag_late.W","p.diag_late.B","p.diag_late.H","p.viral_sup.W","p.viral_sup.B","p.viral_sup.H","ndiag.MSM","ndiag.total")
years = c(2017, 2018, 2019, 2020, 2021, 2022)
y_from_zero = TRUE

# Add labels to 3 city code
Jurisdiction_labels <- c(
	MAR = "Phoenix, AZ",      # MAR - Phoenix, AZ
	SFA = "San Francisco, CA",# SFA - San Francisco, CA
	LOR = "Los Angeles, CA",  # LOR - Los Angeles, CA
	SAC = "Sacramento, CA",   # SAC - Sacramento, CA
	SDI = "San Diego, CA",    # SDI - San Diego, CA
	DUV = "Jacksonville, FL", # DUV - Jacksonville, FL
	TAM = "Tampa Bay, FL",    # TAM - Tampa Bay, FL
	BMP = "Miami, FL",        # BMP - Miami, FL
	ORF = "Orlando, FL",      # ORF - Orlando, FL
	ATL = "Atlanta, GA",      # ATL - Atlanta, GA
	COO = "Chicago, IL",      # COO - Chicago, IL
	MIN = "Indianapolis, IN", # MIN - Indianapolis, IN
	EBR = "Baton Rouge, LA",  
	ORL = "New Orleans, LA",
	BAL = "Baltimore, MD",
	SUF = "Boston, MA",
	WAY = "Detroit, MI",
	CLA = "Las Vegas, NV",
	NYC = "New York City, NY",
	NJE = "Newark–Jersey City, NJ",
	MEC = "Charlotte, NC",
	CUY = "Cleveland, OH",
	FRA = "Columbus, OH",
	HAM = "Cincinnati, OH",
	PHI = "Philadelphia, PA",
	SHE = "Memphis, TN",
	BEX = "San Antonio, TX",
	DFW = "Dallas–Fort Worth, TX",
	HAR = "Houston, TX",
	TRV = "Austin, TX",
	KIN = "Seattle, WA",
	DCM = "Washington, DC"
)
target_labels <- c(
	diag.W = "Prevalence | White",
	diag.B = "Prevalence | Black",
	diag.H = "Prevalence | Hispanic",
	ndiag.M_HET = "New diagnoses | Male, HET",
	ndiag.F_HET = "New diagnoses | Female, HET",
	ndiag.Young = "New diagnoses | MSMWID, 15-34",
	ndiag.Old   = "New diagnoses | MSMWID, 35-64",
	ndiag.PWID  = "New diagnoses | PWID",
	ndiag.W = "New diagnoses | White",
	ndiag.B = "New diagnoses | Black",
	ndiag.H = "New diagnoses | Hispanic",
	p.diag_late.W = "Pct late diagnoses | White",
	p.diag_late.B = "Pct late diagnoses | Black",
	p.diag_late.H = "Pct late diagnoses | Hispanic",
	p.viral_sup.W = "Pct pVL suppressed | White",
	p.viral_sup.B = "Pct pVL suppressed | Black",
	p.viral_sup.H = "Pct pVL suppressed | Hispanic",
	ndiag.MSM = "New diagnoses | MSMWID",
	ndiag.total = "New diagnoses | Total"
	)
#All Functions
plot_mean_sd <- function(input_mean, input_sd, ci_mult = 1.96,
												 time0 = 2016,	strata_labels = NULL, y_label = "Y_name", legend_title = "Legend_Name") {
	
	#  Coerce vectors to 1-column matrices 
	if (is.vector(input_mean)) input_mean <- matrix(input_mean, ncol = 1)
	if (is.vector(input_sd))   input_sd   <- matrix(input_sd,   ncol = 1)
	
	stopifnot(all(dim(input_mean) == dim(input_sd)))
	
	n_time   <- nrow(input_mean)
	n_strata <- ncol(input_mean)
	
	#  Optional strata labels 
	if (!is.null(strata_labels)) {
		stopifnot(length(strata_labels) == n_strata)
		strata_levels <- as.character(strata_labels)
	} else {
		strata_levels <- as.character(seq_len(n_strata))
	}
	
	#  Map index -> calendar year 
	years <- time0 + seq_len(n_time)
	
	df_plot <- expand.grid(
		year   = years,
		strata = factor(seq_len(n_strata), levels = seq_len(n_strata),
										labels = strata_levels)
	)
	#browser()	
	df_plot$mean <- as.vector(input_mean)
	df_plot$sd   <- as.vector(input_sd)
	
	df_plot$upper <- df_plot$mean + ci_mult * df_plot$sd
	df_plot$lower <- df_plot$mean - ci_mult * df_plot$sd
	
	ggplot(df_plot) +
		geom_line(aes(x = year, y = mean, group = strata, colour = factor(strata))) +
		geom_ribbon(aes(x = year, ymin = lower,	ymax = upper,	fill = factor(strata),group = strata),alpha = 0.15) +
		labs(title = Jurisdictions[j],	x = "Year",	y = y_label, colour = legend_title,	fill = legend_title) +
		theme_bw()+
		scale_colour_hue(aesthetics = c("colour", "fill"))
}

plot_median_quantile_old <- function(input_array,
																 time0 = 2016,
																 strata_labels = NULL,
																 y_label = "Y_name",
																 legend_title = "Legend_Name") {
	
	# Coerce matrix to 3D array (samples × time × strata) 
	if (is.matrix(input_array)) {
		input_array <- array(input_array,
												 dim = c(nrow(input_array),
												 				ncol(input_array),
												 				1))
	}
	
	stopifnot(length(dim(input_array)) == 3)
	
	n_time   <- dim(input_array)[2]
	n_strata <- dim(input_array)[3]
	
	# Hard-coded quantiles 
	median_vals <- apply(input_array, c(2, 3), quantile, probs = 0.50)
	lower_vals  <- apply(input_array, c(2, 3), quantile, probs = 0.025)
	upper_vals  <- apply(input_array, c(2, 3), quantile, probs = 0.975)
	#browser()	
	# Optional strata labels 
	if (!is.null(strata_labels)) {
		stopifnot(length(strata_labels) == n_strata)
		strata_levels <- as.character(strata_labels)
	} else {
		strata_levels <- as.character(seq_len(n_strata))
	}
	
	# Map index -> calendar year 
	years <- time0 + seq_len(n_time)
	
	df_plot <- expand.grid(
		year   = years,
		strata = factor(seq_len(n_strata),
										levels = seq_len(n_strata),
										labels = strata_levels)
	)
	
	df_plot$median <- as.vector(median_vals)
	df_plot$lower  <- as.vector(lower_vals)
	df_plot$upper  <- as.vector(upper_vals)
	
	#  Plot 
	ggplot(df_plot) +
		geom_line(
			aes(x = year, y = median, group = strata, colour = strata)
		) +
		geom_ribbon(
			aes(x = year, ymin = lower, ymax = upper,
					fill = strata, group = strata),
			alpha = 0.15
		) +
		labs(title = Jurisdictions[j],
				 x = "Year",
				 y = y_label,
				 colour = legend_title,
				 fill   = legend_title
		) +
		theme_bw() +
		scale_colour_hue(aesthetics = c("colour", "fill"))
}

plot_median_quantile <- function(input_array_model,
																 input_array_target = NULL,
																 time0 = 2016,
																 skip_time = 0,
																 dashed_year = NULL,             # Vertical Dashed line
																 y_label = "Y_name",
																 y_limit = NULL,
																 show_legend = TRUE,
																 strata_labels = NULL,           # Labels of group within Legend
																 legend_title = "Legend_Name",   # Specify the Legend label
																 show_linetype_legend = FALSE,   # show_linetype_legend = TRUE; shows solid = Model Projection (or custom), dotted = Target Data (or custom)
																 color_values = NULL,            # Specify the color of the legend (as oppsoed to linetype_legend below) 
																 source_label = NULL   ) {       # Specify the label for linetype_legend, take vector of 2; default/Null "Model Projection" "Target Data"; show when show_linetype_legend = TRUE
 
	

	# Helper: force to 3D array

	to_3d <- function(x) {
		if (is.matrix(x)) {
			array(x, dim = c(nrow(x), ncol(x), 1))
		} else {
			x
		}
	}

	input_array_model  <- to_3d(input_array_model)
	if (!is.null(input_array_target)) {
		input_array_target <- to_3d(input_array_target)
	}
	
	stopifnot(length(dim(input_array_model)) == 3)
	
	# Skip initial months if requested
	if (skip_time > 0) {
		input_array_model <- input_array_model[, -(seq_len(skip_time)), , drop = FALSE]
		
		if (!is.null(input_array_target)) {
			input_array_target <- input_array_target[, -(seq_len(skip_time)), , drop = FALSE]
		}
	}
# browser()	
	n_time   <- dim(input_array_model)[2]
	n_strata <- dim(input_array_model)[3]
	
	# QUANTILES — MODEL
	median_vals_model <- apply(input_array_model, c(2,3), quantile, probs = 0.50)
	# median_vals_model <- apply(input_array_model, c(2,3), mean)
	lower_vals_model  <- apply(input_array_model, c(2,3), quantile, probs = 0.025)
	upper_vals_model  <- apply(input_array_model, c(2,3), quantile, probs = 0.975)
	
	# STRATA LABELS
	if (!is.null(strata_labels)) {
		stopifnot(length(strata_labels) == n_strata)
		strata_levels <- as.character(strata_labels)
	} else {
		strata_levels <- as.character(seq_len(n_strata))
	}
	
	# YEAR INDEX
	years <- time0 + skip_time/12 + seq_len(n_time)
	

	# DATAFRAME — MODEL
	df_plot_model <- expand.grid(
		year   = years,
		strata = factor(seq_len(n_strata),
										levels = seq_len(n_strata),
										labels = strata_levels)
	)
	
	df_plot_model$median <- as.vector(median_vals_model)
	df_plot_model$lower  <- as.vector(lower_vals_model)
	df_plot_model$upper  <- as.vector(upper_vals_model)
	
	# DATAFRAME — TARGET
	# Allow target and model to have different
	# numbers of time points
	if (!is.null(input_array_target)) {
		
		target_n_time   <- dim(input_array_target)[2]
		target_n_strata <- dim(input_array_target)[3]
		
		# only require same number of strata
		stopifnot(target_n_strata == n_strata)
		
		target_vals <- apply(
			input_array_target,
			c(2,3),
			median
		)
		
		# target years can be shorter than model years
		target_years <- time0 + skip_time/12 + seq_len(target_n_time)
		
		df_plot_target <- expand.grid(
			year   = target_years,
			strata = factor(
				seq_len(n_strata),
				levels = seq_len(n_strata),
				labels = strata_levels
			)
		)
		
		df_plot_target$value <- as.vector(target_vals)
	}
	

	# SOURCE LABELS
	if (is.null(source_label)) {
		source_model  <- "Model Projection"
		source_target <- "Target Data"
	} else {
		stopifnot(length(source_label) == 2)
		source_model  <- source_label[1]
		source_target <- source_label[2]
	}
	
	# PLOT
	
	# Tag model source
	df_plot_model$source <- source_model
	
	p <- ggplot() +
		
		# Model ribbon
		geom_ribbon(
			data = df_plot_model,
			aes(x = year, ymin = lower, ymax = upper,
					fill = strata, group = strata),
			alpha = 0.15
		) +
		
		# Model median
		geom_line(
			data = df_plot_model,
			aes(x = year, y = median, group = strata,
					colour = strata, linetype = source),
			linewidth = 1
		)
	
	# Add target if present
	if (!is.null(input_array_target)) {
		
		df_plot_target$source <- source_target 
		
		p <- p +
			geom_line(
				data = df_plot_target,
				aes(x = year, y = value, group = strata,
						colour = strata, linetype = source),
				linewidth = 1
			)
	}
	# Add optional vertical dashed line
	if (!is.null(dashed_year)) {
		p <- p +
			geom_vline(
				xintercept = dashed_year,
				linetype = "dashed",
				colour = "black"
			)
	}
	# Add optional ability to specify color
	if (!is.null(color_values)) {
		p <- p +		labs(
			title  = Jurisdictions[j],
			x      = "Year",
			y      = y_label,
			colour = if (show_legend) legend_title else NULL,   # hide legend title when legend is off
			fill   = if (show_legend) legend_title else NULL,   # hide legend title when legend is off
			linetype = if (show_linetype_legend) "" else NULL
		) +
			scale_linetype_manual(
				values = setNames(
					c("solid", "dotted"),
					c(source_model, source_target)
				),
				guide = if (show_linetype_legend) "legend" else "none"   # "legend" or"none" are values for guide
			) +
			scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
			{ if(!is.null(y_limit)){coord_cartesian(ylim = c(0, y_limit))} }+
			theme_bw(base_size = 11) +
			theme(
				text = element_text(family = "Arial", colour = "black"),
				plot.title = element_text(face = "bold")
			) +
			# scale_colour_manual(values = color_values, guide = if (show_legend) "legend" else "none")  +
			# scale_fill_manual(values = color_values,	guide = if (show_legend) "legend" else "none")  
			scale_colour_manual(values = color_values, guide = guide_legend()) +
			scale_fill_manual(values = color_values, guide = guide_legend()) +
			theme(
				legend.position = if (show_legend) "bottom" else "none"
			)
	} else {
		p <- p +		labs(
			title  = Jurisdictions[j],
			x      = "Year",
			y      = y_label,
			colour = if (show_legend) legend_title else NULL,   # hide legend title when legend is off
			fill   = if (show_legend) legend_title else NULL,   # hide legend title when legend is off
			linetype = if (show_linetype_legend) "" else NULL
		) +
			scale_linetype_manual(
				values = setNames(
					c("solid", "dotted"),
					c(source_model, source_target)
				),
				guide = if (show_linetype_legend) "legend" else "none"   # "legend" or"none" are values for guide
			) +
			scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
			{ if(!is.null(y_limit)){coord_cartesian(ylim = c(0, y_limit))} }+
			theme_bw(base_size = 11) +
			theme(
				text = element_text(family = "Arial", colour = "black"),
				plot.title = element_text(face = "bold")
			) +
			# scale_colour_hue(aesthetics = c("colour", "fill"),guide = if (show_legend) "legend" else "none" )
			scale_colour_hue(
				aesthetics = c("colour", "fill"),
				guide = "none"
			) +
			theme(
				legend.position = if (show_linetype_legend) "bottom" else "none"
			)
	}
}


################ Adapt plot_target_2 ###############

#' Figure 1: Plot ndiag.MSM & ndiag.total w CI: Plot targets - Adapt from plot_target_2 & plot_one_time_serie
#'
#' @param results results from inference
#' @param target_name calibration targets differ between LEM 1.0 and 2.0
#' @param gp18.gn index of subgroups within the 18 targets
#' @param years years of calibration targets
#' @param plot_dir directory to save plots
#'
#' @return
#' @export
#'
#' @examples
plot_target_2_median <- function(results,target_name,
		years = c(2017, 2018, 2019, 2020, 2021, 2022),
		plot_dir = "Calibration/5_Plots/timeplots",
		comparator = NULL,
		y_from_zero = FALSE,
		y_axis_zoom = 0.25) {
	plot_names <- target_name
	timespan   <- length(years)
	n_targets  <- length(target_name)
	
	## Helper: column indices for each target 
	get_target_indices <- function(i) {
		rep(i * timespan, timespan) + (-timespan + 1):0
	}
	plots <- list()
	for (i in seq_len(n_targets)) {
		
		ith_target_cols <- get_target_indices(i)
# browser()		
		##  Model samples: (n_samples × n_years) 
		ith_target_mat <- results$target_model[, ith_target_cols, drop = FALSE]
		ith_target_obs <- results$target_obs[ith_target_cols]
		
		##  Median & quantiles (across samples) 
		# ith_median <- apply(ith_target_mat, 2, quantile, probs = 0.50)
		ith_median <- apply(ith_target_mat, 2, mean) #Check
		ith_lower  <- apply(ith_target_mat, 2, quantile, probs = 0.025)
		ith_upper  <- apply(ith_target_mat, 2, quantile, probs = 0.975)
	
		df_to_plot <- data.frame(
			time   = years,
			obs    = ith_target_obs,
			median = ith_median,
			lower  = ith_lower,
			upper  = ith_upper
		)
		
		plot_name <- gsub("/", "_", plot_names[i])
		save_name <- file.path(
			plot_dir, "all_groups", paste0(plot_name, ".jpg")
		)
browser()		# get df_to_plot
		g <- ggplot(df_to_plot) +
			geom_line(
				aes(x = time, y = obs, colour = "Observation"),
				linewidth = 0.8
			) +
			geom_line(
				aes(x = time, y = median, colour = "Prediction"),
				linewidth = 0.8
			) +
			geom_ribbon(
				aes(x = time, ymin = lower, ymax = upper),
				alpha = 0.15,
				fill  = "blue"
			) +
			labs(
				title = plot_name,
				x     = "Year",
				y     = "Persons",
				colour = NULL
			) +
			theme_bw(base_size = 11) +
			theme(
				text = element_text(family = "Arial", colour = "black"),
				plot.title = element_text(
					face = "bold",
					size = 40)
			)
		
		if (y_from_zero) {
			g <- g +
				scale_y_continuous(
					limits = c(0, NA)
				)
		} else {
			g <- g +
				scale_y_continuous(
					expand = expansion(mult = y_axis_zoom) # Adjust mult to data range as padding on both ends of the axis and preserves ggplot's automatic scaling
				)
		}
		ggsave(save_name, g, width = 6, height = 4, dpi = 300)
		plots[[plot_name]] <- g
	}
	
	return(plots)
}


#However, want to have the ability to get specific transformation from plot_target_2

calibration_results_multiple$target_model
model_pe <- apply(calibration_results_multiple$target_model, 2, mean) # Sasme as colMeans(calibration_results_multiple$target_model) 
model_sd <- apply(calibration_results_multiple$target_model, 2, sd)   # Do this cause there is no sd() equivalent for colMeans
obs_pe   <- calibration_results_multiple$target_obs

target_name <- c("Prevalence W","Prevalence B","Prevalence H","New diag M_HET","New diag.F_HET","New diag MSM/MWID 15-34","New diag MSM/MWID 35-64","New diag PWID","New diag W","New diag B","New diag H",
								 "Prop late diag W","Prop late diag B","Prop late diag H","Prop Viral Suppression W","Prop Viral Suppression B","Prop Viral Suppression H")
years = c(2017, 2018, 2019, 2020, 2021, 2022)

#Table Results
obs_pe_table <- matrix(obs_pe,	 nrow = 17, ncol = 6, byrow = TRUE)
rownames(obs_pe_table) <- target_name
colnames(obs_pe_table) <- years
model_pe_table <- matrix(model_pe, nrow = 17, ncol = 6, byrow = TRUE) 
rownames(model_pe_table) <- target_name
colnames(model_pe_table) <- years
model_sd_table <- matrix(model_sd, nrow = 17, ncol = 6, byrow = TRUE) 
rownames(model_sd_table) <- target_name
colnames(model_sd_table) <- years

targets <- c("ndiag.Young", "ndiag.Old")

df_to_plot <- do.call(rbind, lapply(targets, function(tgt) {
	
	pred_mean <- as.numeric(model_pe_table[tgt, ])
	obs       <- as.numeric(obs_pe_table[tgt, ])
	sd_vals   <- as.numeric(model_sd_table[tgt, ])
	
	data.frame(
		target    = tgt,
		time      = years,
		obs       = obs,
		pred_mean = pred_mean,
		upper     = pred_mean + 2 * sd_vals,   # 95% CI
		lower     = pred_mean - 2 * sd_vals
	)
}))
#Single Target
g <- ggplot(df_to_plot[df_to_plot$target == "ndiag.Young", ])+
	geom_line(aes(x = time, y = obs, colour = "Observation"))+
	geom_line(aes(x = time, y = pred_mean, colour = "Prediction"))+
	geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
							alpha = 0.1, fill = "blue") +   # adjust color to match predictions
	labs(title = "ndiag.Young", x = "Year", y = "Person")+
	theme_bw()
# Add y-axis starting at 0 only if requested
if (y_from_zero) {
	g <- g + scale_y_continuous(limits = c(0, NA))
}
print(g)
#2 or more targets
g <- ggplot(df_to_plot, aes(x = time, group = target)) +
	geom_line(aes(y = obs, colour = paste(target, "Observation"))) +
	geom_line(aes(y = pred_mean, colour = paste(target, "Prediction"))) +
	geom_ribbon(aes(ymin = lower, ymax = upper, fill = target), alpha = 0.1) +
	labs(title = "Predictions vs Observations", x = "Year", y = "Person") +
	theme_bw()

if (y_from_zero) {
	g <- g + scale_y_continuous(limits = c(0, NA))
}

print(g)

################ Separate plot_target_2 ###############
#### Function 1: Extract the plotting data ( so that I can subset only the target I need)
build_target_plot_data <- function(
		results,
		target_name,
		target_label,
		years = c(2017, 2018, 2019, 2020, 2021, 2022)
) {
	
	timespan  <- length(years)
	n_targets <- length(target_name)
	
	get_target_indices <- function(i) {
		rep(i * timespan, timespan) + (-timespan + 1):0
	}
	
	df_to_plot <- vector("list", n_targets)
	
	for (i in seq_len(n_targets)) {
		
		ith_target_cols <- get_target_indices(i)
		
		ith_target_mat <- results$target_model[
			,
			ith_target_cols,
			drop = FALSE
		]
		
		ith_target_obs <- results$target_obs[
			ith_target_cols
		]
		
		ith_median <- apply(ith_target_mat, 2, mean)
		ith_lower  <- apply(ith_target_mat, 2, quantile, probs = 0.025)
		ith_upper  <- apply(ith_target_mat, 2, quantile, probs = 0.975)
		
		df_to_plot[[i]] <- data.frame(
			target_name  = target_name[i],
			target_label = unname(target_label[i]),
			plot_name    = gsub("/", "_", target_name[i]),
			n_targets    = n_targets,
			time         = years,
			obs          = ith_target_obs,
			median       = ith_median,
			lower        = ith_lower,
			upper        = ith_upper,
			stringsAsFactors = FALSE
		)
	}
	
	dplyr::bind_rows(df_to_plot)
}

#### Function 2: Plot from prepared data
build_target_plots <- function(
		df_to_plot,
		target_name = NULL,
		target_order = NULL,
		y_from_zero = FALSE,
		y_axis_zoom = 0.25
) {
	
	## 1. Filter by target_name
	if (!is.null(target_name)) {
		df_to_plot <- df_to_plot[
			df_to_plot$target_name %in% target_name,
		]
	}
	
	## 2. Split by target_name (stable key)
	g_list <- split(df_to_plot, df_to_plot$target_name)
	
	g_list <- lapply(g_list, function(df) {
		
		ggplot(df) +
			
			geom_line(
				aes(time, obs, colour = "Observation"),
				linewidth = 0.8
			) +
			
			geom_line(
				aes(time, median, colour = "Prediction"),
				linewidth = 0.8
			) +
			
			geom_ribbon(
				aes(time, ymin = lower, ymax = upper),
				alpha = 0.15,
				fill = "blue"
			) +
			
			labs(
				title = unique(df$target_label),   # 👈 IMPORTANT: label comes from df
				x = "Year",
				y = "Persons",
				colour = NULL
			) +
			
			theme_bw(base_size = 11) +
			
			theme(
				plot.title = element_text(
					size = 8,
					lineheight = 0.9,
					margin = margin(b = 2)
				)
			) +
			
			if (y_from_zero) {
				scale_y_continuous(limits = c(0, NA))
			} else {
				scale_y_continuous(
					expand = expansion(mult = y_axis_zoom)
				)
			}
	})
	
	## 3. Ordering (optional)
	if (!is.null(target_order)) {
		g_list <- g_list[target_order]
	}
	
	g_list
}
################ [Archive] Target - Multiple City ###############
#Run: plot_target_2_median 

#Plan to get code for 1 Jurisdictions then loop across the 5 cities
library(patchwork)
library(ggplot2)

Jurisdictions <- c("MAR", "SFA", "LOR", "SAC")
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 

# all.cities
dir <- "All_10000-5000-smdm"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

target_name <- target_name_extended[-c(18,19)]
target_label  <- target_labels[target_name_extended[-c(18,19)]]

target_name_ndiag <- target_name_extended[c(4:11)]
target_label_ndiag  <- target_labels[target_name_extended[c(4:11)]]
target_name_other <- target_name_extended[c(1:3,12:17)]
target_label_other  <- target_labels[target_name_extended[c(1:3,12:17)]]

## Store plots as a list of rows 
row_plots <- list()

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	## Load results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
browser()
	## Generate median plots (returns named list) 
	g_list <- plot_target_2_median(
		results     = calibration_results_multiple,
		target_name = target_label,
		y_from_zero = FALSE
	) 
	
	#  FIX Label & Title
	g_list <- lapply(g_list, function(p) {
		p +
			theme(
				plot.title = element_text(
					size = 8,
					lineheight = 0.9,
					margin = margin(b = 2)
				)
			) +
			labs(title = stringr::str_wrap(p$labels$title, width = 18)) #How Much character in the Title
	})
	if (length(g_list) == 0) {
		stop("No plots returned for ", Jurisdictions[j])
	}
	## Enforce target order explicitly 
  g_list <- g_list[target_label]

	## Add Jurisdiction label (row strip)
  title_plot <- ggplot() +
  	theme_void() +
  	labs(title = Jurisdiction_names[j]) +
  	theme(
  		plot.title = element_text(
  			hjust = 0,
  			face = "bold",
  			size = 14
  		),
  		plot.margin = margin(0, 0, 0, 5)
  	)
  
  row_plots[[Jurisdictions[j]]] <-
  	title_plot /
  	wrap_plots(g_list, ncol = length(target_name)) +
  	plot_layout(heights = c(0.08, 1))
 # browser()	
	# Save individual Jurisdictions
	
	# Create temporary combined plot for saving
	save_plot <- wrap_plots(g_list, ncol = 6) +
		plot_layout(guides = "collect") & 
		theme(legend.position = "bottom")   # or "right"
	
	plot_dir <- "Model_development/7.Red_strata/Calibration/Results/Plots/Individual"
	
	if (!dir.exists(plot_dir)) {
		dir.create(plot_dir, recursive = TRUE)
	}
	
	ggsave(
		filename = file.path(
			plot_dir,
			paste0("CalibTargets_", Jurisdictions[j], ".jpg")
		),
		plot  = save_plot,
		width = 20,
		height = 10,
		dpi = 300
	)
	## Remove large temporary objects
	rm(calibration_results_multiple,
		 g_list,
		 save_plot,
		 title_plot)
	
	gc()
}

## Stack all Jurisdictions (rows) 
combined_plot <-
	wrap_plots(row_plots, ncol = 1) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")

## Save
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/CalibTargets_all_Jurisdictions.jpg",
	plot     = combined_plot,
	width    = 2.5 * length(target_name),
	height   = 2.5 * length(Jurisdictions),
	dpi      = 150
)

################ [Run this] Target - Multiple City ###############
#Run: 

library(patchwork)
library(ggplot2)

Jurisdictions <- c("MAR")
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC")
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM")

# all.cities
dir <- "check_af921fc_10-10"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

target_name <- target_name_extended[c(1:17)]
target_label  <- target_labels[target_name_extended[c(1:17)]]
target_name_ndiag <- target_name_extended[c(4:11)]
target_label_ndiag  <- target_labels[target_name_extended[c(4:11)]]
target_name_other <- target_name_extended[c(1:3,12:17)]
target_label_other  <- target_labels[target_name_extended[c(1:3,12:17)]]

## Store plots as a list of rows 
row_plots <- list()

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	calibration_results_multiple <- readRDS(
		paste0(
			"Calibration/3_Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	
	df_to_plot <- build_target_plot_data(
		results      = calibration_results_multiple,
		target_name  = target_name_extended[c(1:17)],
		target_label  = target_labels[target_name_extended[c(1:17)]],
		years        = years
	)
# browser()	
	g_list_diag <- build_target_plots(
		df_to_plot,
		target_name  = target_name_extended[c(1:3)],
		target_order = target_name_extended[c(1:3)],
		y_from_zero  = FALSE,
		y_axis_zoom  = 2
	)
	g_list_ndiag <- build_target_plots(
		df_to_plot,
		target_name  = target_name_extended[c(4:11)],
		target_order = target_name_extended[c(4:11)],
		y_from_zero  = TRUE
	)
	g_list_ldiag <- build_target_plots(
		df_to_plot,
		target_name  = target_name_extended[c(12:14)],
		target_order = target_name_extended[c(12:14)],
		y_from_zero  = FALSE,
		y_axis_zoom  = 0.8
	)
	g_list_vs <- build_target_plots(
		df_to_plot,
		target_name  = target_name_extended[c(15:17)],
		target_order = target_name_extended[c(15:17)],
		y_from_zero  = FALSE,
		y_axis_zoom  = 2
	)
	g_list <- c(
		g_list_diag,
		g_list_ndiag,
		g_list_ldiag,
		g_list_vs
	)
 # browser()	
	title_plot <- ggplot() +
		theme_void() +
		labs(title = Jurisdiction_names[j]) +
		theme(
			plot.title = element_text(
				hjust = 0,
				face = "bold",
				size = 14
			),
			plot.margin = margin(
				0, 0, 0, 5
			)
		)
	
	row_plots[[Jurisdictions[j]]] <-
		title_plot /
		wrap_plots(
			g_list,
			ncol = length(target_name)
		) +
		plot_layout(
			heights = c(0.08, 1)
		)
	
	# Create temporary combined plot for saving
	save_plot <- wrap_plots(g_list, ncol = 3) +
		plot_layout(guides = "collect") & 
		theme(legend.position = "bottom")   # or "right"
	
	plot_dir <- "Calibration/5_Plots/Individual"
	
	if (!dir.exists(plot_dir)) {
		dir.create(plot_dir, recursive = TRUE)
	}
	save_plot <-
		wrap_plots(
			g_list,
			ncol = 3
		) +
		plot_layout(
			guides = "collect"
		) &
		theme(
			legend.position = "bottom"
		)
	
	ggsave(
		filename = file.path(
			plot_dir,
			paste0(
				"CalibTargets_",
				Jurisdictions[j],
				".jpg"
			)
		),
		plot = save_plot,
		width = 10,
		height = 20,
		dpi = 300
	)
	
	rm(
		calibration_results_multiple,
		df_to_plot,
		g_list,
		save_plot,
		title_plot
	)
	
	gc()
}

## Stack all Jurisdictions (rows) 
combined_plot <-
	wrap_plots(row_plots, ncol = 1) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")

## Save
ggsave(
	filename = "Calibration/5_Plots/CalibTargets_all_Jurisdictions.jpg",
	plot     = combined_plot,
	width    = 2.5 * length(target_name),
	height   = 2.5 * length(Jurisdictions),
	dpi      = 150
)

################ [Run this] Target - Multiple City - %dev Obs Vs Target Table ###############
# Setup

Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

# target_name <- target_name_extended[-c(18,19)]
target_label  <- target_labels[target_name_extended]

# target_name_extended <- unique(c("Prevalence W","Prevalence B","Prevalence H",
# 												"New diag M_HET","New diag.F_HET","New diag MSM_MWID 15-34","New diag Old MSM_MWID 35-64","New diag PWID",
# 												"New diag W","New diag B","New diag H",
# 												"Prop late diag W","Prop late diag B","Prop late diag H",
# 												"Prop Viral Suppression W","Prop Viral Suppression B","Prop Viral Suppression H"))

# Create empty tables (numeric only)
model_table <- matrix(NA_real_, 
											nrow = length(Jurisdictions), 
											ncol = length(target_name_extended),
											dimnames = list(Jurisdictions, target_name_extended))

obs_table <- matrix(NA_real_, 
										nrow = length(Jurisdictions), 
										ncol = length(target_name_extended),
										dimnames = list(Jurisdictions, target_name_extended))


# Loop over Jurisdictions 
for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	## Load results 
	res <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)

	## Derived targets (samples × years) - Average across 6 yrs, comment out as needed
	target_model_diag_W <- res$target_model[, 1:6]
	target_model_diag_B <- res$target_model[, 7:12]
	target_model_diag_H <- res$target_model[, 13:18]
	target_model_ndiag_M.Het <- res$target_model[, 19:24]
	target_model_ndiag_F.Het <- res$target_model[, 25:30]
	target_model_ndiag_Young <- res$target_model[, 31:36]
	target_model_ndiag_Old   <- res$target_model[, 37:42]
	target_model_ndiag_PWID  <- res$target_model[, 43:48]
	target_model_ndiag_W     <- res$target_model[, 49:54]
	target_model_ndiag_B     <- res$target_model[, 55:60]
	target_model_ndiag_H     <- res$target_model[, 61:66]
	target_model_ldiag_W     <- res$target_model[, 67:72]
	target_model_ldiag_B     <- res$target_model[, 73:78]
	target_model_ldiag_H     <- res$target_model[, 79:84]
	target_model_vs_W     <- res$target_model[, 85:90]
	target_model_vs_B     <- res$target_model[, 91:96]
	target_model_vs_H     <- res$target_model[, 97:102]
	target_model_ndiag_MSM <- res$target_model[, 31:36] + res$target_model[, 37:42]
	target_model_ndiag_total <- res$target_model[, 49:54] +		res$target_model[, 55:60] +		res$target_model[, 61:66]

	target_obs_diag_W <- res$target_obs[1:6]
	target_obs_diag_B <- res$target_obs[7:12]
	target_obs_diag_H <- res$target_obs[13:18]
	target_obs_ndiag_M.Het <- res$target_obs[19:24]
	target_obs_ndiag_F.Het <- res$target_obs[25:30]
	target_obs_ndiag_Young <- res$target_obs[31:36]
	target_obs_ndiag_Old   <- res$target_obs[37:42]
	target_obs_ndiag_PWID  <- res$target_obs[43:48]
	target_obs_ndiag_W     <- res$target_obs[49:54]
	target_obs_ndiag_B     <- res$target_obs[55:60]
	target_obs_ndiag_H     <- res$target_obs[61:66]
	target_obs_ldiag_W     <- res$target_obs[67:72]
	target_obs_ldiag_B     <- res$target_obs[73:78]
	target_obs_ldiag_H     <- res$target_obs[79:84]
	target_obs_vs_W     <- res$target_obs[85:90]
	target_obs_vs_B     <- res$target_obs[91:96]
	target_obs_vs_H     <- res$target_obs[97:102]
	target_obs_ndiag_MSM <- res$target_obs[31:36] + res$target_obs[37:42]
	target_obs_ndiag_total <- res$target_obs[49:54] +		res$target_obs[55:60] +		res$target_obs[61:66]

	## Derived targets (samples × years) - Average across 2017,2018
	# target_model_diag_W <- res$target_model[, 1:2]
	# target_model_diag_B <- res$target_model[, 7:8]
	# target_model_diag_H <- res$target_model[, 13:14]
	# target_model_ndiag_M.Het <- res$target_model[, 19:20]
	# target_model_ndiag_F.Het <- res$target_model[, 25:26]
	# target_model_ndiag_Young <- res$target_model[, 31:32]
	# target_model_ndiag_Old   <- res$target_model[, 37:38]
	# target_model_ndiag_PWID  <- res$target_model[, 43:44]
	# target_model_ndiag_W     <- res$target_model[, 49:50]
	# target_model_ndiag_B     <- res$target_model[, 55:56]
	# target_model_ndiag_H     <- res$target_model[, 61:62]
	# target_model_ldiag_W     <- res$target_model[, 67:68]
	# target_model_ldiag_B     <- res$target_model[, 73:74]
	# target_model_ldiag_H     <- res$target_model[, 79:80]
	# target_model_vs_W     <- res$target_model[, 85:86]
	# target_model_vs_B     <- res$target_model[, 91:92]
	# target_model_vs_H     <- res$target_model[, 97:98]
	# target_model_ndiag_MSM <- res$target_model[, 31:32] + res$target_model[, 37:38]
	# target_model_ndiag_total <- res$target_model[, 49:50] +		res$target_model[, 55:56] +		res$target_model[, 61:62]
	# 
	# target_obs_diag_W <- res$target_obs[1:2]
	# target_obs_diag_B <- res$target_obs[7:8]
	# target_obs_diag_H <- res$target_obs[13:14]
	# target_obs_ndiag_M.Het <- res$target_obs[19:20]
	# target_obs_ndiag_F.Het <- res$target_obs[25:26]
	# target_obs_ndiag_Young <- res$target_obs[31:32]
	# target_obs_ndiag_Old   <- res$target_obs[37:48]
	# target_obs_ndiag_PWID  <- res$target_obs[43:44]
	# target_obs_ndiag_W     <- res$target_obs[49:50]
	# target_obs_ndiag_B     <- res$target_obs[55:56]
	# target_obs_ndiag_H     <- res$target_obs[61:62]
	# target_obs_ldiag_W     <- res$target_obs[67:68]
	# target_obs_ldiag_B     <- res$target_obs[73:74]
	# target_obs_ldiag_H     <- res$target_obs[79:80]
	# target_obs_vs_W     <- res$target_obs[85:86]
	# target_obs_vs_B     <- res$target_obs[91:96]
	# target_obs_vs_H     <- res$target_obs[97:98]
	# target_obs_ndiag_MSM <- res$target_obs[31:32] + res$target_obs[37:38]
	# target_obs_ndiag_total <- res$target_obs[49:50] +		res$target_obs[55:56] +		res$target_obs[61:62]
	
	## Derived targets (samples × years) - Average across 2021,2022
	# target_model_diag_W <- res$target_model[, 5:6]
	# target_model_diag_B <- res$target_model[, 11:12]
	# target_model_diag_H <- res$target_model[, 17:18]
	# target_model_ndiag_M.Het <- res$target_model[, 23:24]
	# target_model_ndiag_F.Het <- res$target_model[, 29:30]
	# target_model_ndiag_Young <- res$target_model[, 35:36]
	# target_model_ndiag_Old   <- res$target_model[, 41:42]
	# target_model_ndiag_PWID  <- res$target_model[, 47:48]
	# target_model_ndiag_W     <- res$target_model[, 53:54]
	# target_model_ndiag_B     <- res$target_model[, 59:60]
	# target_model_ndiag_H     <- res$target_model[, 65:66]
	# target_model_ldiag_W     <- res$target_model[, 71:72]
	# target_model_ldiag_B     <- res$target_model[, 77:78]
	# target_model_ldiag_H     <- res$target_model[, 83:84]
	# target_model_vs_W     <- res$target_model[, 89:90]
	# target_model_vs_B     <- res$target_model[, 95:96]
	# target_model_vs_H     <- res$target_model[, 101:102]
	# target_model_ndiag_MSM <- res$target_model[, 35:36] + res$target_model[, 41:42]
	# target_model_ndiag_total <- res$target_model[, 53:54] +		res$target_model[, 59:60] +		res$target_model[, 65:66]
	# 
	# target_obs_diag_W <- res$target_obs[5:6]
	# target_obs_diag_B <- res$target_obs[11:12]
	# target_obs_diag_H <- res$target_obs[17:18]
	# target_obs_ndiag_M.Het <- res$target_obs[23:24]
	# target_obs_ndiag_F.Het <- res$target_obs[29:30]
	# target_obs_ndiag_Young <- res$target_obs[35:36]
	# target_obs_ndiag_Old   <- res$target_obs[41:42]
	# target_obs_ndiag_PWID  <- res$target_obs[47:48]
	# target_obs_ndiag_W     <- res$target_obs[53:54]
	# target_obs_ndiag_B     <- res$target_obs[59:60]
	# target_obs_ndiag_H     <- res$target_obs[65:66]
	# target_obs_ldiag_W     <- res$target_obs[71:72]
	# target_obs_ldiag_B     <- res$target_obs[77:78]
	# target_obs_ldiag_H     <- res$target_obs[83:84]
	# target_obs_vs_W     <- res$target_obs[89:90]
	# target_obs_vs_B     <- res$target_obs[95:96]
	# target_obs_vs_H     <- res$target_obs[101:102]
	# target_obs_ndiag_MSM <- res$target_obs[35:36] + res$target_obs[41:42]
	# target_obs_ndiag_total <- res$target_obs[53:54] +		res$target_obs[59:60] +		res$target_obs[65:66]
	
	## Compute mean across years median
	model_table[j, "ndiag.MSM"]    <- mean(apply(target_model_ndiag_MSM, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.total"]  <- mean(apply(target_model_ndiag_total, 2, quantile, probs = 0.5))
	model_table[j, "diag.W"] <- mean(apply(target_model_diag_W, 2, quantile, probs = 0.5))
	model_table[j, "diag.B"] <- mean(apply(target_model_diag_B, 2, quantile, probs = 0.5))
	model_table[j, "diag.H"] <- mean(apply(target_model_diag_H, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.M_HET"] <- mean(apply(target_model_ndiag_M.Het, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.F_HET"] <- mean(apply(target_model_ndiag_F.Het, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.Young"] <- mean(apply(target_model_ndiag_Young, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.Old"]   <- mean(apply(target_model_ndiag_Old, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.PWID"]  <- mean(apply(target_model_ndiag_PWID, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.W"] <- mean(apply(target_model_ndiag_W, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.B"] <- mean(apply(target_model_ndiag_B, 2, quantile, probs = 0.5))
	model_table[j, "ndiag.H"] <- mean(apply(target_model_ndiag_H, 2, quantile, probs = 0.5))
	model_table[j, "p.diag_late.W"] <- mean(apply(target_model_ldiag_W, 2, quantile, probs = 0.5))
	model_table[j, "p.diag_late.B"] <- mean(apply(target_model_ldiag_B, 2, quantile, probs = 0.5))
	model_table[j, "p.diag_late.H"] <- mean(apply(target_model_ldiag_H, 2, quantile, probs = 0.5))
	model_table[j, "p.viral_sup.W"] <- mean(apply(target_model_vs_W, 2, quantile, probs = 0.5))
	model_table[j, "p.viral_sup.B"] <- mean(apply(target_model_vs_B, 2, quantile, probs = 0.5))
	model_table[j, "p.viral_sup.H"] <- mean(apply(target_model_vs_H, 2, quantile, probs = 0.5))
	
	
	obs_table[j, "ndiag.MSM"]    <- mean(target_obs_ndiag_MSM)
	obs_table[j, "ndiag.total"]  <- mean(target_obs_ndiag_total)
	obs_table[j, "diag.W"] <- mean(target_obs_diag_W)
	obs_table[j, "diag.B"] <- mean(target_obs_diag_B)
	obs_table[j, "diag.H"] <- mean(target_obs_diag_H)
	obs_table[j, "ndiag.M_HET"] <- mean(target_obs_ndiag_M.Het)
	obs_table[j, "ndiag.F_HET"] <- mean(target_obs_ndiag_F.Het)
	obs_table[j, "ndiag.Young"] <- mean(target_obs_ndiag_Young)
	obs_table[j, "ndiag.Old"]   <- mean(target_obs_ndiag_Old)
	obs_table[j, "ndiag.PWID"]  <- mean(target_obs_ndiag_PWID)
	obs_table[j, "ndiag.W"] <- mean(target_obs_ndiag_W)
	obs_table[j, "ndiag.B"] <- mean(target_obs_ndiag_B)
	obs_table[j, "ndiag.H"] <- mean(target_obs_ndiag_H)
	obs_table[j, "p.diag_late.W"] <- mean(target_obs_ldiag_W)
	obs_table[j, "p.diag_late.B"] <- mean(target_obs_ldiag_B)
	obs_table[j, "p.diag_late.H"] <- mean(target_obs_ldiag_H)
	obs_table[j, "p.viral_sup.W"] <- mean(target_obs_vs_W)
	obs_table[j, "p.viral_sup.B"] <- mean(target_obs_vs_B)
	obs_table[j, "p.viral_sup.H"] <- mean(target_obs_vs_H)
	}
	# Output
	model_table
	obs_table
source("Model_development/Optional_util_diagnostic_code/Utility_functions.R")
p <- plot_pct_diff_heatmap(
	obs_table,
	model_table,
	output_dir = file.path(getwd(), "Model_development", "Optional_util_diagnostic_code"),
	cap_value = 100,
	jurisdiction_labels = Jurisdiction_labels
) 

## Compute Median and quantile and see how many obs are within the range of 95% quantile
# Create empty tables (numeric only)
yrs <- c(2017:2022)
model_table_lower <- matrix(NA_real_, 
											nrow = length(target_name_extended), 
											ncol = length(yrs), 
											dimnames = list(target_name_extended,yrs))
model_table_upper <- matrix(NA_real_, 
											nrow = length(target_name_extended), 
											ncol = length(yrs), 
											dimnames = list(target_name_extended,yrs))

obs_table <- matrix(NA_real_, 
										nrow = length(target_name_extended), 
										ncol = length(yrs), 
										dimnames = list(target_name_extended,yrs))

summary_table <- matrix(NA_real_, 
												nrow = length(Jurisdictions), 
												ncol = length(target_name_extended),
												dimnames = list(Jurisdictions, target_name_extended))
# Loop over Jurisdictions
for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	## Load results 
	res <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	
	## Derived targets (samples × years) 
	target_model_diag_W <- res$target_model[, 1:6]
	target_model_diag_B <- res$target_model[, 7:12]
	target_model_diag_H <- res$target_model[, 13:18]
	target_model_ndiag_M.Het <- res$target_model[, 19:24]
	target_model_ndiag_F.Het <- res$target_model[, 25:30]
	target_model_ndiag_Young <- res$target_model[, 31:36]
	target_model_ndiag_Old   <- res$target_model[, 37:42]
	target_model_ndiag_PWID  <- res$target_model[, 43:48]
	target_model_ndiag_W     <- res$target_model[, 49:54]
	target_model_ndiag_B     <- res$target_model[, 55:60]
	target_model_ndiag_H     <- res$target_model[, 61:66]
	target_model_ldiag_W     <- res$target_model[, 67:72]
	target_model_ldiag_B     <- res$target_model[, 73:78]
	target_model_ldiag_H     <- res$target_model[, 79:84]
	target_model_vs_W     <- res$target_model[, 85:90]
	target_model_vs_B     <- res$target_model[, 91:96]
	target_model_vs_H     <- res$target_model[, 97:102]
	target_model_ndiag_MSM <- res$target_model[, 31:36] + res$target_model[, 37:42]
	target_model_ndiag_total <- res$target_model[, 49:54] +		res$target_model[, 55:60] +		res$target_model[, 61:66]
	
	target_obs_diag_W <- res$target_obs[1:6]
	target_obs_diag_B <- res$target_obs[7:12]
	target_obs_diag_H <- res$target_obs[13:18]
	target_obs_ndiag_M.Het <- res$target_obs[19:24]
	target_obs_ndiag_F.Het <- res$target_obs[25:30]
	target_obs_ndiag_Young <- res$target_obs[31:36]
	target_obs_ndiag_Old   <- res$target_obs[37:42]
	target_obs_ndiag_PWID  <- res$target_obs[43:48]
	target_obs_ndiag_W     <- res$target_obs[49:54]
	target_obs_ndiag_B     <- res$target_obs[55:60]
	target_obs_ndiag_H     <- res$target_obs[61:66]
	target_obs_ldiag_W     <- res$target_obs[67:72]
	target_obs_ldiag_B     <- res$target_obs[73:78]
	target_obs_ldiag_H     <- res$target_obs[79:84]
	target_obs_vs_W     <- res$target_obs[85:90]
	target_obs_vs_B     <- res$target_obs[91:96]
	target_obs_vs_H     <- res$target_obs[97:102]
	target_obs_ndiag_MSM <- res$target_obs[31:36] + res$target_obs[37:42]
	target_obs_ndiag_total <- res$target_obs[49:54] +		res$target_obs[55:60] +		res$target_obs[61:66]

	## Compute median and quantile for that area (keep year)
	model_table_lower["ndiag.MSM", c(1:6)]    <- apply(target_model_ndiag_MSM, 2, quantile, probs = 0.025) 
	model_table_lower["ndiag.total", c(1:6)]  <- apply(target_model_ndiag_total, 2, quantile, probs = 0.025)  
	model_table_lower["diag.W", c(1:6)] <- apply(target_model_diag_W, 2, quantile, probs = 0.025)  
	model_table_lower["diag.B", c(1:6)] <- apply(target_model_diag_B, 2, quantile, probs = 0.025) 
	model_table_lower["diag.H", c(1:6)] <- apply(target_model_diag_H, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.M_HET", c(1:6)] <- apply(target_model_ndiag_M.Het, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.F_HET", c(1:6)] <- apply(target_model_ndiag_F.Het, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.Young", c(1:6)] <- apply(target_model_ndiag_Young, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.Old", c(1:6)]   <- apply(target_model_ndiag_Old, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.PWID", c(1:6)]  <- apply(target_model_ndiag_PWID, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.W", c(1:6)] <- apply(target_model_ndiag_W, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.B", c(1:6)] <- apply(target_model_ndiag_B, 2, quantile, probs = 0.025)  
	model_table_lower["ndiag.H", c(1:6)] <- apply(target_model_ndiag_H, 2, quantile, probs = 0.025)  
	model_table_lower["p.diag_late.W", c(1:6)] <- apply(target_model_ldiag_W, 2, quantile, probs = 0.025)  
	model_table_lower["p.diag_late.B", c(1:6)] <- apply(target_model_ldiag_B, 2, quantile, probs = 0.025)  
	model_table_lower["p.diag_late.H", c(1:6)] <- apply(target_model_ldiag_H, 2, quantile, probs = 0.025)  
	model_table_lower["p.viral_sup.W", c(1:6)] <- apply(target_model_vs_W, 2, quantile, probs = 0.025)  
	model_table_lower["p.viral_sup.B", c(1:6)] <- apply(target_model_vs_B, 2, quantile, probs = 0.025)  
	model_table_lower["p.viral_sup.H", c(1:6)] <- apply(target_model_vs_H, 2, quantile, probs = 0.025)  
	
	model_table_upper["ndiag.MSM", c(1:6)]    <- apply(target_model_ndiag_MSM, 2, quantile, probs = 0.975) 
	model_table_upper["ndiag.total", c(1:6)]  <- apply(target_model_ndiag_total, 2, quantile, probs = 0.975)  
	model_table_upper["diag.W", c(1:6)] <- apply(target_model_diag_W, 2, quantile, probs = 0.975)  
	model_table_upper["diag.B", c(1:6)] <- apply(target_model_diag_B, 2, quantile, probs = 0.975) 
	model_table_upper["diag.H", c(1:6)] <- apply(target_model_diag_H, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.M_HET", c(1:6)] <- apply(target_model_ndiag_M.Het, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.F_HET", c(1:6)] <- apply(target_model_ndiag_F.Het, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.Young", c(1:6)] <- apply(target_model_ndiag_Young, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.Old", c(1:6)]   <- apply(target_model_ndiag_Old, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.PWID", c(1:6)]  <- apply(target_model_ndiag_PWID, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.W", c(1:6)] <- apply(target_model_ndiag_W, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.B", c(1:6)] <- apply(target_model_ndiag_B, 2, quantile, probs = 0.975)  
	model_table_upper["ndiag.H", c(1:6)] <- apply(target_model_ndiag_H, 2, quantile, probs = 0.975)  
	model_table_upper["p.diag_late.W", c(1:6)] <- apply(target_model_ldiag_W, 2, quantile, probs = 0.975)  
	model_table_upper["p.diag_late.B", c(1:6)] <- apply(target_model_ldiag_B, 2, quantile, probs = 0.975)  
	model_table_upper["p.diag_late.H", c(1:6)] <- apply(target_model_ldiag_H, 2, quantile, probs = 0.975)  
	model_table_upper["p.viral_sup.W", c(1:6)] <- apply(target_model_vs_W, 2, quantile, probs = 0.975)  
	model_table_upper["p.viral_sup.B", c(1:6)] <- apply(target_model_vs_B, 2, quantile, probs = 0.975)  
	model_table_upper["p.viral_sup.H", c(1:6)] <- apply(target_model_vs_H, 2, quantile, probs = 0.975) 
	
	obs_table["ndiag.MSM", c(1:6)]    <- target_obs_ndiag_MSM
	obs_table["ndiag.total", c(1:6)]  <- target_obs_ndiag_total
	obs_table["diag.W", c(1:6)] <- target_obs_diag_W
	obs_table["diag.B", c(1:6)] <- target_obs_diag_B
	obs_table["diag.H", c(1:6)] <- target_obs_diag_H
	obs_table["ndiag.M_HET", c(1:6)] <- target_obs_ndiag_M.Het
	obs_table["ndiag.F_HET", c(1:6)] <- target_obs_ndiag_F.Het
	obs_table["ndiag.Young", c(1:6)] <- target_obs_ndiag_Young
	obs_table["ndiag.Old", c(1:6)]   <- target_obs_ndiag_Old
	obs_table["ndiag.PWID", c(1:6)]  <- target_obs_ndiag_PWID
	obs_table["ndiag.W", c(1:6)] <- target_obs_ndiag_W
	obs_table["ndiag.B", c(1:6)] <- target_obs_ndiag_B
	obs_table["ndiag.H", c(1:6)] <- target_obs_ndiag_H
	obs_table["p.diag_late.W", c(1:6)] <- target_obs_ldiag_W
	obs_table["p.diag_late.B", c(1:6)] <- target_obs_ldiag_B
	obs_table["p.diag_late.H", c(1:6)] <- target_obs_ldiag_H
	obs_table["p.viral_sup.W", c(1:6)] <- target_obs_vs_W
	obs_table["p.viral_sup.B", c(1:6)] <- target_obs_vs_B
	obs_table["p.viral_sup.H", c(1:6)] <- target_obs_vs_H
	
	within_bounds <- obs_table >= model_table_lower &	obs_table <= model_table_upper
  counts <- rowSums(within_bounds, na.rm = TRUE)
  summary_table[j, names(counts)] <- counts
}
#Visualize it
df <- as.data.frame(as.table(summary_table)) # Need to convert matrix to table first before use as.dataframe to convert to long format
colnames(df) <- c("scenario", "indicator", "count")

p <- ggplot(df, aes(indicator, scenario, fill = count)) +
	geom_tile() +
	geom_text(aes(label = count)) +
	scale_fill_gradient(
		low = "red",
		high = "white",
		limits = c(0, 6)
	) +
	theme_bw(base_size = 11) +
	theme(
		text = element_text(family = "Arial", colour = "black"),
		plot.title = element_text(face = "bold")
	)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Optional relabeling
if (!is.null(Jurisdiction_labels)) {
	p <- p + scale_y_discrete(labels = Jurisdiction_labels)
} 
  p <- p + scale_x_discrete(labels = target_label)
p


################ [Run this] Figure : Projected new ndiagnosis &  Rate by races (w,b,h) - Multiple City ###############

Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

plots_ndiag <- vector("list", length(Jurisdictions))
plots_ndiagr  <- vector("list", length(Jurisdictions))

table <- matrix(NA, nrow = length(Jurisdictions), ncol = 10)
colnames(table) <- c("Median_2025_num", "Median_2045_num","Median_2025_r","Median_2045_r","Median_2025_r_W","Median_2045_r_W","Median_2025_r_B","Median_2045_r_B","Median_2025_r_H","Median_2045_r_H")
rownames(table) <- Jurisdiction_names
years_idx <- c(
	"2025 Estimates" = 9,
	#	"2025 Estimates" = 9,
	"2045 Estimates" = 29
)

compute_ndiag_total <- function(outa, ndiag_inds, end_yr_ind, start_yr_ind = 12) {
	
	# ndiag differences
	x_diff <- apply(
		outa[, end_yr_ind, c(20) , ndiag_inds, drop = FALSE], 
		c(1, 3, 4),
		diff
	)    # Outa: sample*time*state*strata; apply diff over time(2) dimension , keep other (1,3,4) same -> x_diff here (before aperm): time*sample*state*strata
	#browser()	
	x_diff <- aperm(x_diff, c(2, 1, 3, 4)) # x_diff (after aperm): sample*time*state*strata (again)
	
	n_sample  <- dim(x_diff)[1]
	n_state   <- dim(x_diff)[3]
	n_strata  <- dim(x_diff)[4]
	
	# First year
	x0 <- outa[, start_yr_ind, c(20) , ndiag_inds, drop = FALSE]     # x0: sample*time*state*strata, drop = False to maintain dim
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))
	
	# Full time series
	x_full <- abind(x0, x_diff, along = 2)                  # x_full: sample*time*state*strata, bind along "time" dimension (2)
	
	# Sum over strata + state
	apply(x_full, c(1, 2), sum)	                            # sum over everything while keeping sample (1) & time (2) dimensions
	#browser()
}

for (j in seq_along(Jurisdictions)) {
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	outa <- calibration_results_multiple$outa
	
	# Incidence index mapping
	ndiag_inds <- c(
		White = gp.gn$white,
		Black = gp.gn$black,
		Hisp  = gp.gn$hisp
	)
	
	# Compute totals
	ndiag_totals <- lapply(
		ndiag_inds,
		compute_ndiag_total,
		outa        = outa,
		end_yr_ind  = end_yr_ind,
		start_yr_ind = 12
	)
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
	ndiag_totals <- list(
		black = Reduce(`+`, ndiag_totals[grep("^Black", names(ndiag_totals))]),
		hisp  = Reduce(`+`, ndiag_totals[grep("^Hisp",  names(ndiag_totals))]),
		white = Reduce(`+`, ndiag_totals[grep("^White", names(ndiag_totals))])
	) #ndiag_sum <- Reduce(`+`, ndiag_totals) -> will give the same value for calculate from total above
	
	ndiag_array <- abind::abind(ndiag_totals, along = 3)
	
	g_ndiag <- plot_median_quantile(ndiag_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new HIV diagnoses", legend_title = "Race/Ethnicity",
																	color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle(Jurisdiction_names[j])

	plots_ndiag[[j]] <- g_ndiag
		
	#Population sums 
	pop_sum_white   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1,2),	sum) # Per Xiao, denominator for ndiag is also All pop
	pop_sum_black   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1,2),	sum)
	pop_sum_hisp    <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1,2),	sum)
	
	
	##IRs 
	ndiagR_sum_white        <- ndiag_totals$white        / pop_sum_white      * 100000
	ndiagR_sum_black        <- ndiag_totals$black        / pop_sum_black      * 100000
	ndiagR_sum_hisp         <- ndiag_totals$hisp         / pop_sum_hisp       * 100000
	
	ndiagR_array <- abind::abind(ndiagR_sum_black,ndiagR_sum_hisp,ndiagR_sum_white, along = 3)
	g_ndiagr <- plot_median_quantile(ndiagR_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV diagnoses per 100,000", legend_title = "Race/Ethnicity",
																	 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle(Jurisdiction_names[j])
	plots_ndiagr[[j]] <- g_ndiagr
# browser()
	## Diagnostic Table
	table[j,c(1:2)] <- apply( ndiag_array[,years_idx,],2, median)
	table[j,c(3:4)] <- apply(ndiagR_array[,years_idx,],2, median)
	table[j,c(5:6)]  <- apply(ndiagR_sum_white[,years_idx],2, median)
	table[j,c(7:8)]  <- apply(ndiagR_sum_black[,years_idx],2, median)
	table[j,c(9:10)] <- apply( ndiagR_sum_hisp[,years_idx],2, median)
}
## Diagnostic table
file_name <- "Model_development/7.Red_strata/Calibration/Results/Table_ndiag.rds"
# saveRDS(table, file_name)
table <- readRDS(file_name)
table <- as.data.frame(table)
## Change in rate
table$change_r <- table$Median_2045_r - table$Median_2025_r
table[order(table$change_r, decreasing = TRUE),]
table$pct_change_r <- table$Median_2045_r/table$Median_2025_r
table[order(table$pct_change_r, decreasing = TRUE),]

for (j in seq_along(plots_ndiag)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiag_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ndiag[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_ndiagr)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiagR_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ndiagr[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}

# Now Remove y-axis from individual plots again to apply a common plots later
g <- ggplotGrob(plots_ndiagr[[1]] +
									theme(legend.text = element_text(size = 16),
												legend.title = element_text(size = 18))) # Option to change legend size
guide_positions <- which(sapply(g$grobs, function(x) x$name == "guide-box"))
legend <- g$grobs[[guide_positions[1]]]

plots_ndiagr <- lapply(plots_ndiagr, function(p) p + labs(y = NULL))
plots_ndiagr <- lapply(plots_ndiagr, function(p) p + theme(legend.position = "none"))
combined_plot_ndiagr_clean <- wrap_plots(plots_ndiagr, ncol = 8)

#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot_ndiag <- wrap_plots(plots_ndiag, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
# combined_plot_ndiagr  <- wrap_plots(plots_ndiagr, ncol = 8) +
# 	plot_layout(guides = "collect") &
# 	theme(legend.position = "bottom")
combined_plot_ndiagr <- ggdraw() +
	draw_plot(combined_plot_ndiagr_clean , x = 0.005, y = 0.005, width = 0.95, height = 1) +
	draw_label("Rate of new HIV diagnoses per 100,000",
						          x = 0.005  , y = 0.5, angle = 90  , vjust = 0.5) +
	draw_plot(legend,   x = 0.5, y = 0.005, width = 0.10, height = 0.005)

#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/ndiag_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_ndiag,
	width = 24,
	height = 12,
	dpi = 300
)
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/ndiagR_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_ndiagr,
	width = 24,
	height = 12,
	dpi = 300
)

################ [Run this] Figure : ndiagnosis rates (target) & incidence Rate (projected) All pop - Multiple City ###############
Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]
end_yr_ind_plus6 <- c(end_yr_ind_calibration,84,96,108,120,132,144,156,168,180,192)

plots_number <- vector("list", length(Jurisdictions))
plots_rate  <- vector("list", length(Jurisdictions))

table <- matrix(NA, nrow = length(Jurisdictions), ncol = 5)
colnames(table) <- c("Target_2016_rate","Target_2022_rate","Median_2016_rate","Median_2022_rate", "Median_2032_rate")
rownames(table) <- Jurisdiction_names
years_idx <- c(
	"2022 Estimates" = 6,
	#	"2025 Estimates" = 9,
	"2032 Estimates" = 16
)

compute_strata_total <- function(outa, strata_inds, end_yr_ind, start_yr_ind = 12) {
	
	if (missing(strata_inds) || is.null(strata_inds)) {
		strata_inds <- seq_len(dim(outa)[4])
	}
	
	# Differences (ndiag or inc)
	x_diff <- apply(
		outa[, end_yr_ind, c(17:19) , strata_inds, drop = FALSE], # incidence all source 17:19
		c(1, 3, 4),
		diff
	)    # Outa: sample*time*state*strata; apply diff over time(2) dimension , keep other (1,3,4) same -> x_diff here (before aperm): time*sample*state*strata
	#browser()	
	x_diff <- aperm(x_diff, c(2, 1, 3, 4)) # x_diff (after aperm): sample*time*state*strata (again)
	
	n_sample  <- dim(x_diff)[1]
	n_state   <- dim(x_diff)[3]
	n_strata  <- dim(x_diff)[4]
	
	# First year  (ndiag or inc)
	x0 <- outa[, start_yr_ind, c(17:19) , strata_inds, drop = FALSE]     # x0: sample*time*state*strata, drop = False to maintain dim
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))
	
	# Full time series
	x_full <- abind(x0, x_diff, along = 2)                  # x_full: sample*time*state*strata, bind along "time" dimension (2)
	
	# Sum over strata + state
	apply(x_full, c(1, 2), sum)	                            # sum over everything while keeping sample (1) & time (2) dimensions
	#browser()
}

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	outa <- calibration_results_multiple$outa
	
	
	# strata_inds <- c(
	# 	# White = gp.gn$white,
	# 	# Black = gp.gn$black,
	# 	all.msm  = gp.gn$all.msm
	# )
	strata_inds <- list() # All pop
	
	# Compute totals
	totals <- lapply(
		if (length(strata_inds)) strata_inds else list(NULL),
		compute_strata_total,
		outa = outa,
		# end_yr_ind = end_yr_ind_plus6, #end_yr_ind_calibration If you only want to do calib period, they are exactly the same
		end_yr_ind = end_yr_ind,
		start_yr_ind = 12
	)
#browser()	#median_totals <- lapply(totals, function(x) apply(x, 2, median)) # Use this to check with the other version later
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
	# ndiag_totals <- list(
	# 	white = Reduce(`+`, ndiag_totals[grep("^White", names(ndiag_totals))]),
	# 	black = Reduce(`+`, ndiag_totals[grep("^Black", names(ndiag_totals))]),
	# 	hisp  = Reduce(`+`, ndiag_totals[grep("^Hisp",  names(ndiag_totals))])
	# ) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
	# 
	input_array_model <- abind::abind(totals, along = 3)
	
	### Target
	# diag.obs  = data.matrix(subset(data.frame(read.xlsx("Model_development/7.Red_strata/Calibration/251022_calibration_targets.xlsx", sheet="TotalDiagnosedPLHIV",  colNames =T)), city == Jurisdictions[j] & year != 2016)[ , c("W","B","H")])
	ndiag.obs = data.matrix(subset(data.frame(read.xlsx("Model_development/7.Red_strata/Calibration/251022_calibration_targets.xlsx", sheet="NewDiagnoses", colNames =T)), city == Jurisdictions[j] & year != 2016)[ , c("M.HET","F.HET","X15.34.MSM","X35.64.MSM","X15.34.MWID","X35.64.MWID","PWID","W","B","H")])
	# obs.inc.all = subset(data.frame(read.xlsx("Model_development/7.Red_strata/Calibration/251022_calibration_targets.xlsx", sheet="obs.inc.all", colNames =T)), city == Jurisdictions[j])[ , c("value","low","high")]
	# obs.inc.msm = subset(data.frame(read.xlsx("Model_development/7.Red_strata/Calibration/251022_calibration_targets.xlsx", sheet="obs.inc.msm", colNames =T)), city == Jurisdictions[j])[ , c("value","low","high")]

	target_vec <- rowSums(ndiag.obs[, 8:10]) # sum of W,B,H give all
	input_array_target <- array(
		matrix(target_vec, nrow = nrow(totals[[1]]), ncol = length(target_vec), byrow = TRUE),
		dim = c(nrow(totals[[1]]), length(target_vec), 1) # convert input into 3 dimension arrayc
	)
	#browser()	
	###Back to plot
	g_number <- plot_median_quantile(input_array_model  =input_array_model ,
																	 input_array_target =input_array_target,
																	 time0 = 2016, 
																	 dashed_year = 2022,
																	 strata_labels = c("Entire Population"), 
																	 y_label = " HIV cases", 
																	 legend_title = "Legend") + ggtitle(Jurisdiction_names[j])
	
	plots_number[[j]] <- g_number
	
	#Population sums 
	# pop_sum_all   <- apply(	calibration_results_multiple$outa[ , end_yr_ind_plus6, 1:16, ],	c(1,2),	sum) # Per Xiao, denominator for ndiag is also All pop
	pop_sum_all   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, ],	c(1,2),	sum)
	
# browser()	
	##IRs 
	rate_sum_all_model       <- totals[[1]]       / pop_sum_all      * 100000
	# ndiagR_sum_black        <- ndiag_totals$black        / pop_sum_black      * 100000
	# ndiagR_sum_hisp         <- ndiag_totals$hisp         / pop_sum_hisp       * 100000
	rate_sum_all_target <- sweep(	pop_sum_all[,c(1:length(target_vec))],2,	target_vec,
																			FUN = function(pop, obs) obs / pop * 100000	)
	
	rate_array_model  <- abind::abind(rate_sum_all_model , along = 3)
	rate_array_target <- abind::abind(rate_sum_all_target, along = 3)
	
	g_rate <- plot_median_quantile(input_array_model=rate_array_model, 
																 input_array_target =rate_array_target, 
																 dashed_year = 2022,
																 strata_labels = c("Entire Population"), 
																 y_label = "HIV cases rates, per 100,000", 
																 y_limit = 60,
																 show_legend = FALSE,
																 show_linetype_legend = TRUE,
																 source_label = c("New incidence HIV infections (model)","New diagnoses (obs/target)")
																 ) + ggtitle(Jurisdiction_names[j])
						
	plots_rate[[j]] <- g_rate
	
## Diagnostic Table
table[j,c(1:2)]  <- apply(rate_array_target[,c(1,6)   ,],2, median)
table[j,c(3:5)]  <- apply( rate_array_model[,c(1,6,16),],2, median)
}
file_name <- "Model_development/7.Red_strata/Calibration/Results/Table_ndiag_inc_all.rds"
saveRDS(table, file_name)
table <- readRDS(file_name)
table <- as.data.frame(table)
## Diagnostic
increase_both        <- table[	table$Target_2016_rate < table$Target_2022_rate &		table$Median_2022_rate < table$Median_2032_rate,]
target_down_model_up <- table[	table$Target_2016_rate > table$Target_2022_rate &		table$Median_2022_rate < table$Median_2032_rate,]

for (j in seq_along(plots_number)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiag-inc_all_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_number[[j]]+ ggtitle(""),
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_rate)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiagR-incR_all_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_rate[[j]]+ ggtitle(""),
		width = 5,
		height = 4,
		dpi = 300
	)
}
# Now Remove y-axis from individual plots again to apply a common plots later
g <- ggplotGrob(plots_rate[[1]] +
									theme(legend.text = element_text(size = 16),
												legend.title = element_text(size = 18))) # Option to change legend size
guide_positions <- which(sapply(g$grobs, function(x) x$name == "guide-box"))
legend <- g$grobs[[guide_positions[1]]]

plots_rate <- lapply(plots_rate, function(p) p + labs(y = NULL))
plots_rate <- lapply(plots_rate, function(p) p + theme(legend.position = "none"))
combined_plot_rate_clean <- wrap_plots(plots_rate, ncol = 8)

#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot_number <- wrap_plots(plots_number, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
# combined_plot_rate  <- wrap_plots(plots_rate, ncol = 8) +
# 	plot_layout(guides = "collect") &
# 	theme(legend.position = "bottom")
combined_plot_rate <- ggdraw() +
	draw_plot(combined_plot_rate_clean , x = 0.01, y = 0.02, width = 0.98, height = 0.98) +
	draw_label("Rate of new HIV cases per 100,000",
						 x = 0.005  , y = 0.5, angle = 90  , vjust = 0.5) +
	draw_plot(legend,   x = 0.5, y = 0.01, width = 0.10, height = 0.005)

#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/ndiag-inc_all_Jurisdictions.jpg",
	plot = combined_plot_number ,
	width = 24,
	height = 12,
	dpi = 300
)
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/ndiagR-incR_all_Jurisdictions.jpg",
	plot = combined_plot_rate,
	width = 24,
	height = 12,
	dpi = 300
)

################ [        ] Figure 2: Incidence Rate - Multiple City ###############
Jurisdictions <- c("COO") #all.cities #c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL")
dir <- "Phase2_c98714e_50-50"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

plots <- vector("list", length(Jurisdictions))
table <- matrix(NA, nrow = length(Jurisdictions), ncol = 2)
colnames(table) <- c("Median_2017", "Median_2045")
#colnames(table) <- c("Median_2020", "Median_2045")
rownames(table) <- Jurisdictions

for (j in seq_along(Jurisdictions)) {
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	
	#Incidence differences (RECOMPUTED EACH LOOP)
	x_diff <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 17:19, ], 
		c(1, 3, 4),
		diff
	)    # Outa: sample*time*state*strata; apply diff over time(2) dimension , keep other (1,3,4) same -> x_diff here (before aperm): time*sample*state*strata

	x_diff <- aperm(x_diff, c(2, 1, 3, 4)) # x_diff (after aperm): sample*time*state*strata (again)
		n_sample <- dim(x_diff)[1]
		n_state <- dim(x_diff)[3]
		n_strata<- dim(x_diff)[4]
	x0 <- calibration_results_multiple$outa[, 12, 17:19, ]                  # x0: sample*time*state*strata 
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))                #Convert 1st year 3-dim array into 4 dim, second dim = 1 due to being 1 yr
	
	x_full <- abind(x0, x_diff, along = 2) # x_full: sample*time*state*strata, bind along "time" dimension (2)
	
	#Incidence sums 
	#	inc_sum_w <- apply(x_full[,,, gp.gn$white], c(1, 2), sum)
	#	inc_sum_b <- apply(x_full[,,, gp.gn$black], c(1, 2), sum)
	#	inc_sum_h <- apply(x_full[,,, gp.gn$hisp ], c(1, 2), sum)
	inc_sum_total <- apply(x_full, c(1, 2), sum) #2017 (end_yr_ind=12) to 2045 (end_yr_ind=348), check: length(2017:2045)=length(end_yr_ind)
	inc_sum_total_2020onward <- inc_sum_total[,-(1:3)] #sample * time =length(2020:2045) 
	
browser()	#median_inc <- apply(inc_sum_total, 2, median)
	#Population sums 
	#	pop_sum_w <- apply(
	#		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$white],
	#		c(1, 2), sum
	#	)
	#	pop_sum_b <- apply(
	#		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$black],
	#		c(1, 2), sum
	#	)
	#	pop_sum_h <- apply(
	#		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$hisp ],
	#		c(1, 2), sum
	#	)
	pop_sum_total <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 1:16, ],
		c(1, 2),
		sum
	)
	pop_sum_total_2020onward <- pop_sum_total[,-(1:3)] #sample * time =length(2020:2045)
	
	##IRs 
	#	IR_mean_total <- inc_mean_total / pop_mean_total
	#	IR_sd_total <- IR_mean_total * sqrt(
	#		(inc_sd_total / inc_mean_total)^2 +
	#		(pop_sd_total / pop_mean_total)^2
	#	)
	IR_sum_total            <- inc_sum_total            / pop_sum_total            * 100000
	IR_sum_total_2020onward <- inc_sum_total_2020onward / pop_sum_total_2020onward * 100000
		
	#	plot_mean_sd(
	#		IR_mean_total,
	#		IR_sd_total,
	#		ci_mult = 1.96,
	#		strata_labels = c("IR Total Pop"),
	#		y_label = "IR",
	#		legend_title = "Legend"
	#	)
	
	g <- plot_median_quantile(
		IR_sum_total,	time0 = 2017, #all years
#		IR_sum_total_2020onward,	time0 = 2020,
		strata_labels = c("Incidence Rate (IR) Total Population"),
		y_label = "IR per 100,000",
		legend_title = "Legend"
	) + ggtitle(Jurisdiction_names[j])
	
	if (y_from_zero) {
		g <- g + scale_y_continuous(limits = c(0, NA))
	}
	
	plots[[j]] <- g
	
	first_year <- IR_sum_total[, 1]
#	first_year <- IR_sum_total_2020onward[, 1]
	last_year  <- IR_sum_total[, ncol(IR_sum_total)]
#	last_year  <- IR_sum_total_2020onward[, ncol(IR_sum_total_2020onward)]
	
	table[j, ] <- c(
		quantile(first_year, 0.5, na.rm = TRUE),
		quantile(last_year,  0.5, na.rm = TRUE)
	)
}
print(table)
Reduction_pct <- (table[, "Median_2017"] - table[, "Median_2045"]) /	table[, "Median_2017"] * 100
Reduction_pct <- round(Reduction_pct, 1)
table <- cbind(table, Reduction_pct)
print(table)


#3 - Save individual plots
dir.create("Plots/Individual", showWarnings = FALSE)

for (j in seq_along(plots)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/IR_total_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}

#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot <- wrap_plots(plots, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")

#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/IR_total_all_Jurisdictions.jpg",
	plot = combined_plot,
	width = 24,
	height = 12,
	dpi = 300
)

################ [Run this] Figure : Incidence &  Rate by races (w,b,h) - Multiple City ###############
Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

plots_inc <- vector("list", length(Jurisdictions))
plots_ir  <- vector("list", length(Jurisdictions))

table <- matrix(NA, nrow = length(Jurisdictions), ncol = 15)
colnames(table) <- c("Median_2017_num","Median_2025_num", "Median_2045_num",
										 "Median_2017_r","Median_2025_r","Median_2045_r",
										 "Median_2017_r_W","Median_2025_r_W","Median_2045_r_W",
										 "Median_2017_r_B","Median_2025_r_B","Median_2045_r_B",
										 "Median_2017_r_H","Median_2025_r_H","Median_2045_r_H")
rownames(table) <- Jurisdiction_names
years_idx <- c(
	"2017 Estimates" = 1,
	"2025 Estimates" = 9,
	"2045 Estimates" = 29
)

compute_inc_total <- function(outa, inc_inds, end_yr_ind, start_yr_ind = 12) {
	
	# Incidence differences
	x_diff <- apply(
		outa[, end_yr_ind, c(17:19) , inc_inds, drop = FALSE], 
		c(1, 3, 4),
		diff
	)    # Outa: sample*time*state*strata; apply diff over time(2) dimension , keep other (1,3,4) same -> x_diff here (before aperm): time*sample*state*strata
	
	x_diff <- aperm(x_diff, c(2, 1, 3, 4)) # x_diff (after aperm): sample*time*state*strata (again)
	
	n_sample  <- dim(x_diff)[1]
	n_state   <- dim(x_diff)[3]
	n_strata  <- dim(x_diff)[4]
	
	# First year
	x0 <- outa[, start_yr_ind, c(17:19) , inc_inds, drop = FALSE]     # x0: sample*time*state*strata, drop = False to maintain dim
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))
	
	# Full time series
	x_full <- abind(x0, x_diff, along = 2)                  # x_full: sample*time*state*strata, bind along "time" dimension (2)
	
	# Sum over strata + state
	apply(x_full, c(1, 2), sum)	                            # sum over everything while keeping sample (1) & time (2) dimensions
# browser()
	}

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	outa <- calibration_results_multiple$outa
	
	# Incidence index mapping
	inc_inds <- c(
		White = gp.gn$white,
		Black = gp.gn$black,
		Hisp  = gp.gn$hisp
	)
#browser()	
	# Compute totals
	inc_totals <- lapply(
		inc_inds,
		compute_inc_total,
		outa        = outa,
		end_yr_ind  = end_yr_ind,
		start_yr_ind = 12
	)
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
	inc_totals <- list(
		black = Reduce(`+`, inc_totals[grep("^Black", names(inc_totals))]),
		hisp  = Reduce(`+`, inc_totals[grep("^Hisp",  names(inc_totals))]),
		white = Reduce(`+`, inc_totals[grep("^White", names(inc_totals))])
	) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
	inc_all = inc_totals$white +	inc_totals$black +	inc_totals$hisp
	
	inc_array <- abind::abind(inc_totals, along = 3)
	
	g_inc <- plot_median_quantile(inc_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new infections", legend_title = "Race/Ethnicity",
																color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle(Jurisdiction_names[j])
	plots_inc[[j]] <- g_inc
	
	#Population sums 
	pop_sum_white   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1,2),	sum) # Per Xiao, denominator for inc is All pop
	pop_sum_black   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1,2),	sum)
	pop_sum_hisp    <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1,2),	sum)
	pop_sum_all     <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16,  ],	c(1,2),	sum)

	
	##IRs 
	IR_sum_white        <- inc_totals$white        / pop_sum_white      * 100000
	IR_sum_black        <- inc_totals$black        / pop_sum_black      * 100000
	IR_sum_hisp         <- inc_totals$hisp         / pop_sum_hisp       * 100000
	IR_sum_all          <- inc_all                 / pop_sum_all        * 100000

	IR_array <- abind::abind(IR_sum_black,IR_sum_hisp,IR_sum_white, along = 3)

	g_ir <- plot_median_quantile(IR_array,skip_time = 0,y_limit = 100 ,dashed_year = 2025 ,strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV incidence per 100,000", legend_title = "Race/Ethnicity",
															 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle(Jurisdiction_names[j])
	
	plots_ir[[j]] <- g_ir

	## Diagnostic Table
	table[j,c(1:3)] <- apply( inc_all[,years_idx],2, median)
	table[j,c(4:6)] <- apply( IR_sum_all[,years_idx],2, median)
	table[j,c(7:9)]  <- apply(IR_sum_white[,years_idx],2, median)
	table[j,c(10:12)]  <- apply(IR_sum_black[,years_idx],2, median)
	table[j,c(13:15)] <- apply( IR_sum_hisp[,years_idx],2, median)
}
## Diagnostic table
file_name <- "Model_development/7.Red_strata/Calibration/Results/Table_inc_race.rds"
saveRDS(table, file_name)
table <- readRDS(file_name)
table <- as.data.frame(table)
## Change in rate
table$change_r <- table$Median_2045_r - table$Median_2025_r
table[order(table$change_r, decreasing = TRUE),]
table$pct_change_r <- table$Median_2045_r/table$Median_2025_r
table[order(table$pct_change_r, decreasing = TRUE),]
## Average increase per Races
#mean(table$Median_2045_r_H - table$Median_2025_r_H)
# table[table$Median_2025_r_H < table$Median_2025_r_B & table$Median_2045_r_H > table$Median_2045_r_B, ]
## Increase in Hisp 17 -> 45 causing decrease in Black or White 17 -> 25 
## followed by subsequent increase 25 -> 45
# table[
# 	table$Median_2045_r_H > table$Median_2017_r_H &
# 		(
# 			(	table$Median_2025_r_B < table$Median_2017_r_B &
# 				table$Median_2045_r_B > table$Median_2025_r_B) |
# 			(	table$Median_2025_r_W < table$Median_2017_r_W &
# 				table$Median_2045_r_W > table$Median_2025_r_W)
# 		),
# ]
# OR both Black and White populations must show this reversal pattern <- report this currently
table[
	table$Median_2045_r_H > table$Median_2017_r_H &
		table$Median_2025_r_B < table$Median_2017_r_B &
		table$Median_2045_r_B > table$Median_2025_r_B &
		table$Median_2025_r_W < table$Median_2017_r_W &
		table$Median_2045_r_W > table$Median_2025_r_W,
]
for (j in seq_along(plots_inc)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/nInc_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_inc[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_ir)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/IR_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ir[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
# Now Remove y-axis from individual plots again to apply a common plots later
g <- ggplotGrob(plots_ir[[1]] +
		theme(legend.text = element_text(size = 16),
			    legend.title = element_text(size = 18))) # Option to change legend size
guide_positions <- which(sapply(g$grobs, function(x) x$name == "guide-box"))
legend <- g$grobs[[guide_positions[1]]]

plots_ir <- lapply(plots_ir, function(p) p + labs(y = NULL))
plots_ir <- lapply(plots_ir, function(p) p + theme(legend.position = "none"))
combined_plot_ir_clean <- wrap_plots(plots_ir, ncol = 8)

#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot_inc <- wrap_plots(plots_inc, ncol = 4) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
# combined_plot_ir  <- wrap_plots(plots_ir, ncol = 4) +
# 	plot_layout(guides = "collect") &
# 	theme(legend.position = "bottom")
combined_plot_ir <- ggdraw() +
	draw_plot(combined_plot_ir_clean , x = 0.01, y = 0.02, width = 0.98, height = 0.98) +
	draw_label("Rate of new HIV incidence per 100,000",
						 x = 0.005  , y = 0.5, angle = 90  , vjust = 0.5) +
	draw_plot(legend,   x = 0.5, y = 0.01, width = 0.10, height = 0.005)

#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/nInc_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_inc,
	width = 24,
	height = 12,
	dpi = 300
)
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/IR_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_ir,
	width = 24,
	height = 12,
	dpi = 300
)

################ [        ] Figure : Incidence &  Rate by risk group (young, old, het.m, het.f, pwid) - Multiple City ###############
table <- matrix(NA, nrow = length(Jurisdictions), ncol = 10)
colnames(table) <- c("Median_2025_r_old","Median_2045_r_old","Median_2025_r_young","Median_2045_r_young","Median_2025_r_het.m","Median_2045_r_het.m","Median_2025_r_het.f","Median_2045_r_het.f","Median_2025_r_het.high","Median_2045_r_het.high")
rownames(table) <- Jurisdiction_names
years_idx <- c(
	"2025 Estimates" = 9,
	#	"2025 Estimates" = 9,
	"2045 Estimates" = 29
)

# Run everything up to "compute_inc_total" from "Figure : Incidence &  Rate by races (w,b,h) - Multiple City" again
for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	outa <- calibration_results_multiple$outa
	
	# Incidence index mapping
	inc_inds <- c(
		old   = gp.gn$old,
		young = gp.gn$young,
		pwid  = gp.gn$pwid,
		het.m = gp.gn$het.m,
		het.f = gp.gn$het.f,
		het.high = gp.gn$high
	)
	#browser()	
	# Compute totals
	inc_totals <- lapply(
		inc_inds,
		compute_inc_total,
		outa        = outa,
		end_yr_ind  = end_yr_ind,
		start_yr_ind = 12
	)
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, ! H
	inc_totals <- list(
		old   = Reduce(`+`, inc_totals[grep("^old"  , names(inc_totals))]),
		young = Reduce(`+`, inc_totals[grep("^young", names(inc_totals))]),
		pwid  = Reduce(`+`, inc_totals[grep("^pwid" , names(inc_totals))]),
		het.m = Reduce(`+`, inc_totals[grep("^het.m", names(inc_totals))]),
		het.f = Reduce(`+`, inc_totals[grep("^het.f" , names(inc_totals))]),
		het.high = Reduce(`+`, inc_totals[grep("^het.high" , names(inc_totals))])
	) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
	
	inc_array <- abind::abind(inc_totals, along = 3)
	
	# g_inc <- plot_median_quantile(inc_array, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new infections", legend_title = "Race/Ethnicity",
	# 															color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle(Jurisdiction_names[j])
	g_inc <- plot_median_quantile(inc_array, strata_labels = c("old","young","pwid","het.m","het.f","het.high"), y_label = "Incidence", legend_title = "Legend",
																	color_values = c("old" = "#F8766D",	"young" = "#619CFF","pwid" = "#00BA38",	"het.m" = "yellow","het.f"="purple", "het.high"="black")) + ggtitle(Jurisdiction_names[j])
	plots_inc[[j]] <- g_inc
	
	#Population sums 
	pop_sum_old     <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$old  ],	c(1,2),	sum) # Per Xiao, denominator for inc is All pop
	pop_sum_young   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$young],	c(1,2),	sum)
	pop_sum_pwid    <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$pwid ],	c(1,2),	sum)
	pop_sum_het.m   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$het.m],	c(1,2),	sum)
	pop_sum_het.f   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$het.f ],	c(1,2),	sum)
	pop_sum_het.high<- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$high ],	c(1,2),	sum)

	##IRs 
	IR_sum_old          <- inc_totals$old          / pop_sum_old        * 100000
	IR_sum_young        <- inc_totals$young        / pop_sum_young      * 100000
	IR_sum_pwid         <- inc_totals$pwid         / pop_sum_pwid       * 100000
	IR_sum_het.m        <- inc_totals$het.m        / pop_sum_het.m      * 100000
	IR_sum_het.f        <- inc_totals$het.f        / pop_sum_het.f      * 100000
	IR_sum_het.high     <- inc_totals$het.high     / pop_sum_het.high   * 100000
	
	IR_array <- abind::abind(IR_sum_old, IR_sum_young, IR_sum_pwid, IR_sum_het.m, IR_sum_het.f,IR_sum_het.high, along = 3)
	g_ir <- plot_median_quantile(IR_array ,strata_labels = c("old","young","pwid","het.m","het.f","het.high"), y_label = "IR per 100,000", legend_title = "Legend",
															 color_values = c("old" = "#F8766D",	"young" = "#619CFF","pwid" = "#00BA38",	"het.m" = "yellow","het.f"="purple", "het.high"="black")) + ggtitle(Jurisdiction_names[j])
	plots_ir[[j]] <- g_ir
	# browser()
	## Diagnostic Table
	table[j,c(1:2)] <- apply(IR_sum_old[,years_idx],2, median)
	table[j,c(3:4)] <- apply(IR_sum_young[,years_idx],2, median)
	table[j,c(5:6)] <- apply(IR_sum_het.m[,years_idx],2, median)
	table[j,c(7:8)] <- apply(IR_sum_het.f[,years_idx],2, median)
	table[j,c(9:10)]<- apply(IR_sum_het.high[,years_idx],2, median)
}
## Diagnostic table
file_name <- "Model_development/7.Red_strata/Calibration/Results/Table_inc_risk.rds"
saveRDS(table, file_name)
table <- readRDS(file_name)
table <- as.data.frame(table)
## Change in rate
# table$change_r <- table$Median_2045_r - table$Median_2025_r
# table[order(table$change_r, decreasing = TRUE),]
# table$pct_change_r <- table$Median_2045_r/table$Median_2025_r
# table[order(table$pct_change_r, decreasing = TRUE),]

for (j in seq_along(plots_inc)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/nInc_separate_risk_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_inc[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_ir)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/IR_separate_risk_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ir[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot_inc <- wrap_plots(plots_inc, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
combined_plot_ir  <- wrap_plots(plots_ir, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/nInc_total_separate_risk_Jurisdictions.jpg",
	plot = combined_plot_inc,
	width = 24,
	height = 12,
	dpi = 300
)
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/IR_total_separate_risk_Jurisdictions.jpg",
	plot = combined_plot_ir,
	width = 24,
	height = 12,
	dpi = 300
)

################ [Run this] Figure 3: IRR - Multiple City ###############
# STEP 0: Setup
library(ggpubr)

Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

years_idx <- c(
	"2025 Estimates" = 9,
#	"2030 Estimates" = 14,
	"2045 Estimates" = 29
)


df_list <- list()   # only summaries live outside
#### Step 1: Loop across Jurisdiction
for (j in seq_along(Jurisdictions)) {
	
	## ---- Load Jurisdiction-specific results ----
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	## ---- Incidence differences (RECOMPUTED EACH LOOP) ----
	x_diff <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 17:19, ],
		c(1, 3, 4),
		diff
	)      # 28 time x samples x 3 states x 30 strata
	x_diff <- aperm(x_diff, c(2, 1, 3, 4))
	n_sample <- dim(x_diff)[1]
	n_state  <- dim(x_diff)[3]
	n_strata <- dim(x_diff)[4]
	
	x0 <- calibration_results_multiple$outa[, 12, 17:19, ]
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))                #Convert 1st year 3-dim array into 4 dim, second dim = 1 due to being 1 yr
	
	x_full <- abind(x0, x_diff, along = 2)
	
	## ---- Incidence sums ----
	inc_sum_w <- apply(x_full[,,, gp.gn$white], c(1, 2), sum)
	inc_sum_b <- apply(x_full[,,, gp.gn$black], c(1, 2), sum)
	inc_sum_h <- apply(x_full[,,, gp.gn$hisp ], c(1, 2), sum)
	
	## ---- Population sums ----
	pop_sum_w <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$white], 
		c(1, 2), sum
	)
	pop_sum_b <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$black],
		c(1, 2), sum
	)
	pop_sum_h <- apply(
		calibration_results_multiple$outa[, end_yr_ind, 1:16, gp.gn$hisp ],
		c(1, 2), sum
	)
	
	## ---- IRs and IRRs ----
	IR_w <- inc_sum_w / pop_sum_w * 100000
	IR_b <- inc_sum_b / pop_sum_b * 100000
	IR_h <- inc_sum_h / pop_sum_h * 100000
	
	IRR_bw <- IR_b / IR_w
	IRR_hw <- IR_h / IR_w
# browser()	
	## ---- Extract summary stats ONLY ----
#	for (yr in c(4, 14, 29)) {
	for (yr in c(9,     29)) {		
#		scenario <- c("2020 Estimates","2030 Estimates","2045 Estimates")[match(yr, c(4,14,29))]
		scenario <- c("2025 Estimates",                 "2045 Estimates")[match(yr, c(9,   29))]
		
		df_list[[length(df_list) + 1]] <- data.frame(
			city = Jurisdiction_names[j],
			scenario = scenario,
			comparison = "Black vs White",
			irr   = quantile(IRR_bw[, yr], 0.50),
			lower = quantile(IRR_bw[, yr], 0.025),
			upper = quantile(IRR_bw[, yr], 0.975)
		)
	
		df_list[[length(df_list) + 1]] <- data.frame(
			city = Jurisdiction_names[j],
			scenario = scenario,
			comparison = "Hispanic vs White",
			irr   = quantile(IRR_hw[, yr], 0.50),
			lower = quantile(IRR_hw[, yr], 0.025),
			upper = quantile(IRR_hw[, yr], 0.975)
		)
	}
}

df <- dplyr::bind_rows(df_list)

#### Step 2: Add spacer rows + row type
df <- df %>%
	arrange(city, scenario, comparison) %>%
	mutate(row_type = "data")

df_spacer <- df %>%
	distinct(city) %>%
	mutate(
		scenario = "CITY_HEADER",
		comparison = NA,
		irr = NA,
		lower = NA,
		upper = NA,
		row_type = "spacer"
	)

df_plot <- bind_rows(df_spacer, df) %>%
	mutate(
		scenario = factor(
			scenario,
			levels = c("CITY_HEADER",
								 "2025 Estimates",
#								 "2030 Estimates",
								 "2045 Estimates")
		)
	) %>%
	arrange(city, scenario)

#### Step 3: Factor ordering 
df_plot <- df_plot %>%
	mutate(
		city_scenario = factor(
			paste(city, scenario, sep = " — "),
			levels = rev(unique(paste(city, scenario, sep = " — ")))
		)
	)

#### City separator lines
spacer_y <- df_plot %>%
	filter(row_type == "spacer") %>%
	mutate(y = as.numeric(city_scenario)) %>%
	pull(y)

#### Step 4: Forest plots function (fixed)
forest_panel <- function(data, title) {
	ggplot(data, aes(x = irr, y = city_scenario, color = scenario)) +
		## City separator lines
		geom_hline(
			yintercept = spacer_y,
			color = "grey70",
			linewidth = 0.4
		) +
		# Only draw points/error bars for data rows
		geom_point(data = subset(data, row_type == "data"), size = 2, na.rm = TRUE) +
		## 🔹 IRR value labels
		geom_text(
			data = subset(data, row_type == "data"),
			aes(label = paste0("IRR=", sprintf("%.2f", irr))),
			hjust = -0.3, vjust = 0.3,
			size = 4,
			show.legend = FALSE
		) +
		##
		geom_errorbarh(
			data = subset(data, row_type == "data"),
			aes(xmin = lower, xmax = upper),
			height = 0.2,
			linetype = "dotted",
			alpha = 0.85
		) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
		scale_color_manual(values = c("#F8766D","#619CFF","#2E7D32")) +
		coord_cartesian(xlim = c(0, 15)) +
		# scale_y_discrete(drop = FALSE) +
		scale_y_discrete(
			drop = FALSE,
			labels = function(x) {
				sapply(x, function(lbl) {
					if (grepl("CITY_HEADER$", lbl)) {
						as.expression(
							bquote(bold(.(sub(" — CITY_HEADER", "", lbl))))
						)
					} else if (grepl("2025 Estimates$", lbl)) {
						"2025 IRR"
					} else if (grepl("2045 Estimates$", lbl)) {
						"2045 IRR"
					} else {""}
				})
			})+
		labs(
			title = title,
			x = NULL,
			y = NULL,
			color = NULL
		) +
		theme_bw(base_size = 18) +
		theme(
			# axis.text.y  = element_blank(),
			# axis.ticks.y = element_blank(),
			axis.text.y  = element_text(size = 12),
			axis.ticks.y = element_line(),
			panel.grid.major.y = element_blank(),
			panel.grid.minor = element_blank(),
			plot.title = element_text(face = "bold"),
			legend.position = "none" 
		)
}

# #### Step 5: Left panel (labels)
# # Left labels panel aligned with forest plots
# p0 <- ggplot(df_plot, aes(y = city_scenario)) +
# 	geom_text(
# 		aes(
# 			x = 0,
# 			label = case_when(
# 				scenario == "CITY_HEADER"    ~ city,
# 				scenario == "2025 Estimates" ~ "– 2025 IRR",
# #				scenario == "2030 Estimates" ~ "– 2030 IRR",
# 				scenario == "2045 Estimates" ~ "– 2045 IRR",
# 				TRUE ~ ""
# 			),
# 			fontface = ifelse(scenario == "CITY_HEADER", "bold", "plain")
# 		),
# 		hjust = 0
# 	) +
# 	# Make the plot "invisible" except text
# 	scale_x_continuous(limits = c(0, 1)) +
# 	scale_y_discrete(labels = Jurisdiction_labels, drop = FALSE) +        # keep all factor levels for alignment
# 	labs(
# 		title = "City / Scenario",
# 		x = NULL,
# 		y = NULL
# 	) +
# 	theme_minimal(base_size = 10) +
# 	theme(
# 		axis.text.y = element_blank(),    # y labels handled by geom_text
# 		axis.ticks.y = element_blank(),
# 		panel.grid.major.y = element_blank(),
# 		panel.grid.minor = element_blank(),
# 		plot.title = element_text(face = "bold"),
# 		plot.background = element_blank(),
# 		panel.background = element_blank(),
# 		axis.title.y = element_blank(),
# 		axis.text.x = element_text(size = 10)
# 	)
#### Step 5: Left panel (labels)

# Dummy dataframe for shared legend
legend_df <- data.frame(
	scenario = c("2025 Estimates", "2045 Estimates"),
	x = c(0, 0),
	y = c(1, 1)
)

p0 <- ggplot(df_plot, aes(y = city_scenario)) +
	
	# Dummy points ONLY for legend
	geom_point(
		data = legend_df,
		aes(x = x, y = y, color = scenario),
		alpha = 0
	) +
	
	geom_text(
		aes(
			x = 0,
			label = case_when(
				scenario == "CITY_HEADER"    ~ city,
				scenario == "2025 Estimates" ~ "– 2025 IRR",
				#				scenario == "2030 Estimates" ~ "– 2030 IRR",
				scenario == "2045 Estimates" ~ "– 2045 IRR",
				TRUE ~ ""
			),
			fontface = ifelse(scenario == "CITY_HEADER", "bold", "plain")
		),
		hjust = 0
	) +
	
	scale_color_manual(
		values = c(
			"2025 Estimates" = "#F8766D",
			"2045 Estimates" = "#619CFF"
		),
		name = NULL
	) +
	
	# Make the plot "invisible" except text
	scale_x_continuous(limits = c(0, 1)) +
	scale_y_discrete(labels = Jurisdiction_labels, drop = FALSE) +
	
	labs(
		title = "City / Scenario",
		x = NULL,
		y = NULL
	) +
	
	theme_minimal(base_size = 20) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold"),
		plot.background = element_blank(),
		panel.background = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(),
		
		# KEEP legend here
		legend.position = "right"
	)

#### Step 6: Forest plots
p1 <- forest_panel(
	df_plot %>% filter(comparison == "Black vs White"),
	"Black versus White"
)

p2 <- forest_panel(
	df_plot %>% filter(comparison == "Hispanic vs White"),
	"Hispanic versus White"
)

#### Step 7: Combine panels
final_plot <- ggarrange(
# p0, p1, p2,
	p1, p2,
# ncol = 3,
	ncol = 2,
# widths = c(0.45, 1, 1),
  widths = c( 1, 1),
	common.legend = TRUE,
	legend = "bottom",#"right",
	align = "v"
)

final_plot

#### Step 8: Save
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/Figure_3_IRR.jpg",
	plot = final_plot,
	width = 14,
	height = 14,
	units = "in",
	dpi = 600
)
#### Step 9: Diagnostic code
file_name <- "Model_development/7.Red_strata/Calibration/Results/df_IRR.rds"
saveRDS(df, file = file_name)
df <- readRDS("Model_development/7.Red_strata/Calibration/Results/df_IRR.rds")

wide <- reshape(
	df[, c("city", "scenario", "comparison", "irr")],
	idvar = c("city", "comparison"),
	timevar = "scenario",
	direction = "wide"
)
## How many is increased from 2020-> 2045
wide$Increased <-	wide$`irr.2045 Estimates` > wide$`irr.2025 Estimates`
tapply(
	wide$Increased,
	wide$comparison,
	sum,
	na.rm = TRUE
)
## Lowest and highest in 2045
wide %>%
	group_by(comparison) %>%
	summarise(
		lowest_city  = city[which.min(`irr.2045 Estimates`)],
		lowest_irr   = min(`irr.2045 Estimates`, na.rm = TRUE),
		highest_city = city[which.max(`irr.2045 Estimates`)],
		highest_irr  = max(`irr.2045 Estimates`, na.rm = TRUE),
		.groups = "drop"
	)
# 1.Cities where Black vs White increased, ranked largest to smallest increase
bw_increase <- wide %>%
	filter(comparison == "Black vs White",
				 Increased) %>%
	mutate(change = `irr.2045 Estimates` - `irr.2025 Estimates`) %>%
	arrange(desc(change))

bw_increase

bw_decrease <- wide %>%
	filter(comparison == "Black vs White",
				 !Increased) %>%
	mutate(change = (`irr.2025 Estimates`-`irr.2045 Estimates`)/`irr.2025 Estimates`) %>%
	arrange(change)

# 2. Cities where Hispanic vs White increased, ranked largest to smallest increase
hw_increase <- wide %>%
	filter(comparison == "Hispanic vs White",
				 Increased) %>%
	mutate(change = `irr.2045 Estimates` - `irr.2025 Estimates`) %>%
	arrange(desc(change))

hw_increase

hw_decrease <- wide %>%
	filter(comparison == "Hispanic vs White",
				 !Increased) %>%
	mutate(change = (`irr.2025 Estimates`-`irr.2045 Estimates`)/`irr.2025 Estimates`) %>%
	arrange(change)
# 3. Hispanic vs White > Black vs White in 2025, but Black vs White > Hispanic vs White in 2045
tmp <- wide %>%
	select(city, comparison,
				 irr2025 = `irr.2025 Estimates`,
				 irr2045 = `irr.2045 Estimates`) %>%
	pivot_wider(
		names_from = comparison,
		values_from = c(irr2025, irr2045)
	)

#Hispanic > Black in 2025, but Black > Hispanic in 2045
tmp %>%
	filter(
		`irr2025_Hispanic vs White` > `irr2025_Black vs White`,
		`irr2045_Black vs White` > `irr2045_Hispanic vs White`
	) %>%
	select(city)
#Black > Hispanic in 2025, but Hispanic > Black in 2045
tmp %>%
	filter(
		`irr2025_Black vs White` > `irr2025_Hispanic vs White`,
		`irr2045_Hispanic vs White` > `irr2045_Black vs White`
	) %>%
	select(city)

# 4. IRR for Black & Hisp decrease 2025 to 2045
wide %>%
	filter(
		comparison %in% c("Black vs White", "Hispanic vs White")
	) %>%
	group_by(city) %>%
	summarise(
		both_decrease = all(`irr.2045 Estimates` < `irr.2025 Estimates`),
		.groups = "drop"
	) %>%
	filter(both_decrease) %>%
	select(city)

################ [Run this] Figure : Demographic & cascade of care change ###############

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
group18 = function(group30){
	g18 = matrix(0, length(vt[c(-1)]), 18) #length(vt[c(-1)]) or length(vt)
	for (i in 1:18){
		ind = which(rname %in% vparameters$names18[i])
		g18[ ,i] = rowSums(group30[ ,ind, drop = FALSE]) #must include drop= FALSE
		colnames(g18) <- vparameters$names18 
	}
	return(g18)
}
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
group.race = function(group30){
	g.race = matrix(0,length(vt[c(-1)]),3)   #length(vt[c(-1)]) or length(vt)
	g.race[ ,1] = rowSums(group30[ ,gp.gn$white])  #
	g.race[ ,2] = rowSums(group30[ ,gp.gn$hisp])   #
	g.race[ ,3] = rowSums(group30[ ,gp.gn$black])  #
	colnames(g.race) <- c("White/Other","Hispanic","Black") # In this order for Black to be on top
	return(g.race)
}
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
#Plan to get code for 1 Jurisdictions then loop across the 5 cities
library(patchwork)
library(ggplot2)

#Run the following functions:
	#plot_stacked_area
  #group18  #use vt[c(-1)]
	#group.calibration
	#group.race
	#group.state
Jurisdictions <- c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

## Store plots as a list of rows 
row_plots <- list()

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	## Load results 
	res <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)

	# Average the sample result
	outa <- apply(res$outa[ , , , ], c(2,3,4), mean)
	
	# All 30 strata
	# All 30 strata
	#  n.strata    <- apply(outa[ , , ], c(1,3), sum) # across 2 dimension: Time x strata
	n.strata    <- apply(outa[ , c(1:16), ], c(1,3), sum) # across 2 dimension: Time x strata
	
	# g1 <- plot_stacked_area(time = vt[-1], data=n.strata, name = vparameters$names.gp, CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j])
	
	# 18 strata
	# n.strata.18    <- group18(n.strata)
	# g2 <- plot_stacked_area(time = vt[-1], data=n.strata.18, name = vparameters$names18, CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j]) #time = vt[-1]
	# 18 strata - Exclude HET
	# n.strata.18.exHET   <-n.strata.18[,1:12]
	# g2.ex <- plot_stacked_area(time = vt[-1], data=n.strata.18.exHET , name = vparameters$names18[1:12], CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j]) #time = vt[-1]
	# Risk group
	n.strata.cal    <- group.calibration(n.strata)
	g3 <- plot_stacked_area(time = vt[-1], data=n.strata.cal , name= colnames(n.strata.cal), CITY=Jurisdictions[j])# + ggtitle(Jurisdiction_names[j])
	## Risk group - Exclude Sus, and HET include Low/High
  # n.strata.ex.S    <- apply(outa[ , c(14:16), ], c(1,3), sum) # S: c(1:2), I: c(3:6), D: c(7:10), T: c(11:13), O: c(14:16)
  # n.strata.cal.ex.S <- group.calibration.exp(n.strata.ex.S)
  # g3.ex.S <- plot_stacked_area(time = vt[-1], data=n.strata.cal.ex.S , name= colnames(n.strata.cal.ex.S), CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j])
  # g3.ex.S
	
	# Race-level strata
	n.strata.race    <- group.race(n.strata)
	g4 <- plot_stacked_area(time = vt[-1], data=n.strata.race, name = colnames(n.strata.race), CITY=Jurisdictions[j],  
													color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	))# + ggtitle(Jurisdiction_names[j]) #time = vt[-1] or vt
# browser()		
	## Compare Pop Proportion across states
	# All 24 states
	n.state    <- apply(outa[ , , ], c(1,2), sum) # across 2 dimension: Time x states
	n.state.16    <- n.state[,1:16]
	df.state   <- as.data.frame(cbind(vt[-1],n.state.16)) #time = vt[-1] or vt
	colnames(df.state)    <- c("time",vparameters$state.name[1:16])
	n.state.4    <- group.state(df.state)
	g6 <- plot_stacked_area(time = vt[-1], 
													data = n.state.4, 
													name = colnames(n.state.4), 
													CITY=Jurisdictions[j],
													y_label = "Model-estimated population of PLHIV: 15-64 (in millions)")# + ggtitle(Jurisdiction_names[j])
	
	# Race-level strata, among I, D, T
	# df.state[,c("Ia","I1","I2","I3")]
	# g6.I <- plot_stacked_area(time = vt[-1], data = df.state[,c("Ia","I1","I2","I3")], name = colnames(df.state[,c("Ia","I1","I2","I3")]), CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j])
	# g6.D <- plot_stacked_area(time = vt[-1], data = df.state[,c("Da","D1","D2","D3")], name = colnames(df.state[,c("Da","D1","D2","D3")]), CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j])
	# g6.T <- plot_stacked_area(time = vt[-1], data = df.state[,c("T1","T2","T3")], name = colnames(df.state[,c("T1","T2","T3")]), CITY=Jurisdictions[j]) + ggtitle(Jurisdiction_names[j])
	

	plot_dir <- "Model_development/7.Red_strata/Calibration/Results/Plots/Individual"
	
	if (!dir.exists(plot_dir)) {
		dir.create(plot_dir, recursive = TRUE)
	}
	# Save individual 18-level plot (g2)
	ggsave(
		filename = file.path(
			plot_dir,
			paste0("pop_risk_", Jurisdictions[j], ".jpg")
		),
		plot  = g3,        
		width = 12,
		height = 9,
		dpi = 300
	)
	# Save individual race-level plot (g4)
	ggsave(
		filename = file.path(
			plot_dir,
			paste0("pop_race_", Jurisdictions[j], ".jpg")
		),
		plot  = g4,        
		width = 12,
		height = 9,
		dpi = 300
	)
	# Save individual state-level plot (g6)
	ggsave(
		filename = file.path(
			plot_dir,
			paste0("pop_state_", Jurisdictions[j], ".jpg")
		),
		plot  = g6,        
		width = 12,
		height = 9,
		dpi = 300
	)
	
## Combine plots for this Jurisdiction 
g_list <- list(g3, g4, g6) #g1, g2, g2.ex, g3, g4, g6, g6.I, g6.D, g6.T

	row_plots[[j]] <- wrap_plots(
		list(g3, g4, g6), #g1, g2, g2.ex, g3, g4, g6, g6.I, g6.D, g6.T
		ncol = 3 #9
	) +
		plot_annotation(
			subtitle = CITY,
			theme = theme(
				plot.subtitle = element_text(face = "bold", size = 2)
			)
		)
}

## Stack all Jurisdictions 
combined_plot <-
	wrap_plots(row_plots, ncol = 1) +   #ncol = 1, for when include Multiple plot ypes
	plot_layout(guides = "collect", widths = 1) &   
	theme(
		legend.position = "bottom",
		legend.box = "horizontal"
	)


## Save 
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/CalibTargets_all_DemographicCascadeofCare.jpg",
	plot     = combined_plot,
	width    = 49,
	height   = 49, #10
	dpi      = 150
)

#########[Run this] Combined Pop figures ####
library(ggplot2)
library(ggpubr)
library(cowplot)
Jurisdictions <- all.cities

plot_dir <- "Model_development/7.Red_strata/Calibration/Results/Plots/Individual/"

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	pA_file <- paste0(
		plot_dir,
		"pop_race_",
		Jurisdictions[j],
		".jpg"
	)
	
	pB_file <- paste0(
		plot_dir,
		"pop_risk_",
		Jurisdictions[j],
		".jpg"
	)
	
	## image panels
	img_A <- ggdraw() +
		draw_image(pA_file)
	
	img_B <- ggdraw() +
		draw_image(pB_file)
	
	## separate title strips
	title_A <- ggdraw() +
		draw_label(
			"Panel A: Population projection by race/ethnicity",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	title_B <- ggdraw() +
		draw_label(
			"Panel B: Population projection by HIV risk group",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	## combine title + image for each panel
	panel_A <- ggarrange(
		title_A,
		img_A,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	panel_B <- ggarrange(
		title_B,
		img_B,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	## final figure: 1 column, 2 rows
	final_plot <- ggarrange(
		panel_A,
		panel_B,
		ncol = 1,
		nrow = 2
	)
	
	ggsave(
		filename = paste0(
			plot_dir,
			"pop_race_risk_",
			Jurisdictions[j],
			".jpg"
		),
		plot = final_plot,
		width = 7,
		height = 10,
		dpi = 300
	)
}

#### 4 panels figures of
library(ggpubr)
library(cowplot)

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	pA_file <- paste0(
		plot_dir,
		"ndiag_separate_race_",
		Jurisdictions[j],
		".jpg"
	)
	
	pB_file <- paste0(
		plot_dir,
		"ndiagR_separate_race_",
		Jurisdictions[j],
		".jpg"
	)
	
	pC_file <- paste0(
		plot_dir,
		"nInc_separate_race_",
		Jurisdictions[j],
		".jpg"
	)
	
	pD_file <- paste0(
		plot_dir,
		"IR_separate_race_",
		Jurisdictions[j],
		".jpg"
	)
	
	## Images
	img_A <- ggdraw() + draw_image(pA_file)
	img_B <- ggdraw() + draw_image(pB_file)
	img_C <- ggdraw() + draw_image(pC_file)
	img_D <- ggdraw() + draw_image(pD_file)
	
	## Panel titles
	title_A <- ggdraw() +
		draw_label(
			"Panel A: New HIV diagnoses: Total number",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	title_B <- ggdraw() +
		draw_label(
			"Panel B: New HIV diagnoses: Rate per 100,000 population",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	title_C <- ggdraw() +
		draw_label(
			"Panel C: HIV incidence: Total number",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	title_D <- ggdraw() +
		draw_label(
			"Panel D: HIV incidence: Rate per 100,000 population",
			x = 0,
			hjust = 0,
			fontface = "bold"
		)
	
	## Combine title + image
	panel_A <- ggarrange(
		title_A,
		img_A,
		ncol = 1,
		heights = c(0.05, 1)
	)
	
	panel_B <- ggarrange(
		title_B,
		img_B,
		ncol = 1,
		heights = c(0.05, 1)
	)
	
	panel_C <- ggarrange(
		title_C,
		img_C,
		ncol = 1,
		heights = c(0.05, 1)
	)
	
	panel_D <- ggarrange(
		title_D,
		img_D,
		ncol = 1,
		heights = c(0.05, 1)
	)
	
	## Final 4-panel figure
	final_plot <- ggarrange(
		panel_A,
		panel_B,
		panel_C,
		panel_D,
		ncol = 2,
		nrow = 2
	)
	
	ggsave(
		filename = paste0(
			plot_dir,
			"ndiag_ndiagR_inc_IR_combined_",
			Jurisdictions[j],
			".jpg"
		),
		plot = final_plot,
		width = 8,
		height = 8,
		dpi = 300
	)
}

#########[Run this] Combined Inc ndiag figures ####
library(ggplot2)
library(ggpubr)
library(cowplot)


Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

plots_ndiag <- vector("list", length(Jurisdictions))
plots_ndiagr  <- vector("list", length(Jurisdictions))
plots_inc <- vector("list", length(Jurisdictions))
plots_ir  <- vector("list", length(Jurisdictions))

# Run
## compute_inc_total
## compute_ndiag_total

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	outa <- calibration_results_multiple$outa
	
	# Incidence index mapping
	ndiag_inds <- c(
		White = gp.gn$white,
		Black = gp.gn$black,
		Hisp  = gp.gn$hisp
	)
	
	# Compute totals
	ndiag_totals <- lapply(
		ndiag_inds,
		compute_ndiag_total,
		outa        = outa,
		end_yr_ind  = end_yr_ind,
		start_yr_ind = 12
	)
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
	ndiag_totals <- list(
		black = Reduce(`+`, ndiag_totals[grep("^Black", names(ndiag_totals))]),
		hisp  = Reduce(`+`, ndiag_totals[grep("^Hisp",  names(ndiag_totals))]),
		white = Reduce(`+`, ndiag_totals[grep("^White", names(ndiag_totals))])
	) #ndiag_sum <- Reduce(`+`, ndiag_totals) -> will give the same value for calculate from total above
	
	ndiag_array <- abind::abind(ndiag_totals, along = 3)
	
	g_ndiag <- plot_median_quantile(ndiag_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new HIV diagnoses", legend_title = "Race/Ethnicity",
																	color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	))+ ggtitle("")
	
	plots_ndiag[[j]] <- g_ndiag
	
	#Population sums 
	pop_sum_white   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1,2),	sum) # Per Xiao, denominator for ndiag is also All pop
	pop_sum_black   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1,2),	sum)
	pop_sum_hisp    <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1,2),	sum)
	
	
	##IRs 
	ndiagR_sum_white        <- ndiag_totals$white        / pop_sum_white      * 100000
	ndiagR_sum_black        <- ndiag_totals$black        / pop_sum_black      * 100000
	ndiagR_sum_hisp         <- ndiag_totals$hisp         / pop_sum_hisp       * 100000
	
	ndiagR_array <- abind::abind(ndiagR_sum_black,ndiagR_sum_hisp,ndiagR_sum_white, along = 3)
	g_ndiagr <- plot_median_quantile(ndiagR_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV diagnoses per 100,000", legend_title = "Race/Ethnicity",
																	 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle("")
	plots_ndiagr[[j]] <- g_ndiagr
	# browser()

	# Incidence index mapping
	inc_inds <- c(
		White = gp.gn$white,
		Black = gp.gn$black,
		Hisp  = gp.gn$hisp
	)
	#browser()	
	# Compute totals
	inc_totals <- lapply(
		inc_inds,
		compute_inc_total,
		outa        = outa,
		end_yr_ind  = end_yr_ind,
		start_yr_ind = 12
	)
	# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
	inc_totals <- list(
		black = Reduce(`+`, inc_totals[grep("^Black", names(inc_totals))]),
		hisp  = Reduce(`+`, inc_totals[grep("^Hisp",  names(inc_totals))]),
		white = Reduce(`+`, inc_totals[grep("^White", names(inc_totals))])
	) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
	inc_all = inc_totals$white +	inc_totals$black +	inc_totals$hisp
	
	inc_array <- abind::abind(inc_totals, along = 3)
	
	g_inc <- plot_median_quantile(inc_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new infections", legend_title = "Race/Ethnicity",
																color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle("")# + ggtitle(Jurisdiction_names[j])
	plots_inc[[j]] <- g_inc
	
	#Population sums 
	pop_sum_white   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1,2),	sum) # Per Xiao, denominator for inc is All pop
	pop_sum_black   <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1,2),	sum)
	pop_sum_hisp    <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1,2),	sum)
	pop_sum_all     <- apply(	calibration_results_multiple$outa[ , end_yr_ind, 1:16,  ],	c(1,2),	sum)
	
	
	##IRs 
	IR_sum_white        <- inc_totals$white        / pop_sum_white      * 100000
	IR_sum_black        <- inc_totals$black        / pop_sum_black      * 100000
	IR_sum_hisp         <- inc_totals$hisp         / pop_sum_hisp       * 100000
	IR_sum_all          <- inc_all                 / pop_sum_all        * 100000
	
	IR_array <- abind::abind(IR_sum_black,IR_sum_hisp,IR_sum_white, along = 3)
	
	g_ir <- plot_median_quantile(IR_array,skip_time = 0,y_limit = 100 ,dashed_year = 2025 ,strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV incidence per 100,000", legend_title = "Race/Ethnicity",
															 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) + ggtitle("")# + ggtitle(Jurisdiction_names[j])
	
	plots_ir[[j]] <- g_ir
	
	########### NEW
	
	# Extract one shared legend
	g_leg <- ggplotGrob(g_ir)
	guide_positions <- which(
		sapply(g_leg$grobs, function(x) x$name == "guide-box")
	)
	legend <- g_leg$grobs[[guide_positions[1]]]
	
	# Remove legends from all four plots
	g_ndiag_clean  <- g_ndiag  + theme(legend.position = "none")
	g_ndiagr_clean <- g_ndiagr + theme(legend.position = "none")
	g_inc_clean    <- g_inc    + theme(legend.position = "none")
	g_ir_clean     <- g_ir     + theme(legend.position = "none")
	
	# Add panel titles
	panel_A <- ggarrange(
		ggdraw() +
			draw_label(
				"Panel A: New HIV diagnoses: \nTotal number",
				x = 0, hjust = 0, fontface = "bold",
				size = 12
			),
		g_ndiag_clean,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	panel_B <- ggarrange(
		ggdraw() +
			draw_label(
				"Panel B: New HIV diagnoses: \nRate per 100,000 population",
				x = 0, hjust = 0, fontface = "bold",
				size = 12
			),
		g_ndiagr_clean,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	panel_C <- ggarrange(
		ggdraw() +
			draw_label(
				"Panel C: HIV incidence: \nTotal number",
				x = 0, hjust = 0, fontface = "bold",
				size = 12
			),
		g_inc_clean,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	panel_D <- ggarrange(
		ggdraw() +
			draw_label(
				"Panel D: HIV incidence: \nRate per 100,000 population",
				x = 0, hjust = 0, fontface = "bold",
				size = 12
			),
		g_ir_clean,
		ncol = 1,
		heights = c(0.08, 1)
	)
	
	# 2 x 2 layout
	panel_grid <- ggarrange(
		panel_A, panel_B,
		panel_C, panel_D,
		ncol = 2,
		nrow = 2
	)
	
	# Add one shared legend at bottom
	combined_plot_ndiag_ndiagR_inc_IR <- ggdraw() +
		draw_plot(
			panel_grid,
			x = 0,
			y = 0.08,
			width = 1,
			height = 0.92
		) +
		draw_grob(
			legend,
			x = 0.25,
			y = 0.00,
			width = 0.50,
			height = 0.08
		)
	
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiag_ndiagR_inc_IR_combined_",
			Jurisdictions[j],
			".jpg"
		),
		plot = combined_plot_ndiag_ndiagR_inc_IR,
		width = 12,
		height = 10,
		dpi = 300
	)
}

for (j in seq_along(plots_ndiag)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiag_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ndiag[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_ndiagr)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/ndiagR_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ndiagr[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_inc)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/nInc_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_inc[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}
for (j in seq_along(plots_ir)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/IR_separate_race_",
			Jurisdictions[j],
			".jpg"
		),
		plot = plots_ir[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}

################ [Run this] Figure : MCMC Traceplots ###############
source("Model_development/7.Red_strata/Calibration/7.plots.R", local = .GlobalEnv)

Jurisdictions <- all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
dir <- "All_10000-5000-jun23"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

## Store plots as a list of rows 
row_plots <- list()

for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	res <- readRDS(
		paste0(
			"Model_development/7.Red_strata/Calibration/Results/calib_",
			dir, "_", Jurisdictions[j], ".rds"
		)
	)
	
	trace_plot <- plot_multichain_trace(res$all_lls) + 
		ggtitle(Jurisdiction_names[j])
	
	row_plots[[j]] <- trace_plot
}
for (j in seq_along(row_plots)) {
	ggsave(
		filename = paste0(
			"Model_development/7.Red_strata/Calibration/Results/Plots/Individual/traceplot_",
			Jurisdictions[j],
			".jpg"
		),
		plot = row_plots[[j]],
		width = 5,
		height = 4,
		dpi = 300
	)
}

#4 - Combine into 4 × 8 grid with ONE shared legend
combined_plot_trace <- wrap_plots(row_plots, ncol = 8) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")

#5 - Save combined plot
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/traceplot_all.jpg",
	plot = combined_plot_trace,
	width = 24,
	height = 12,
	dpi = 300
)

################ Plots all target the deterministic ###############

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
# Load from check_prior_on_target
	# extract_outputs
	# plot_line
source("Model_development/Optional_util_diagnostic_code/Utility_functions.R") 	
source("Model_development/7.Red_strata/7.Red_strata-LoadBaselineWorkspace.R")
source("Model_development/7.Red_strata/Calibration/7.load_targets.R")
source("Model_development/7.Red_strata/calibration/7.plots.R")


# Target
vparameters$pop.diag # First year 2016
calib.target$diag.obs # Remaining years 2016
new_row <- data.frame(
	W = sum(vparameters$pop.diag[gp18.gn$white]),
	B = sum(vparameters$pop.diag[gp18.gn$black]),
	H = sum(vparameters$pop.diag[gp18.gn$hisp])
)
diag.obs <- rbind(new_row, calib.target$diag.obs)


Jurisdictions <- c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# all.cities
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

## Store plots as a list of rows
vparameters_orig <- vparameters
row_plots <- list()

interv <- c(17) # No intervention
cc <- interv[1]


# Need to be part of Jurisdictions loop
for (j in seq_along(Jurisdictions)) {
	
	message("Processing ", Jurisdictions[j])
	
	vparameters <- vparameters_orig   # reset each iteration
	CITY <- all.cities[[j]]
	
	source("Model_development/7.Red_strata/7.Red_strata-Module-Data.input.R")
	
	out_base <- intervention.model.combination(
		input.parameters = vparameters,
		current.int = interventions[cc]
	)
	base <- extract_outputs(out_base)
	
	# --- Sensitivity change ---
	vlist$mix.prob.YO$pe <- 0.54
	
	source("Model_development/7.Red_strata/Calibration/7.Update_vlist_to_vparameters.R")
	vparameters <- par_fix(vparameters, vlist)
	
	result      <- par_pop(vparameters, vlist)
	vparameters <- result$vparameters
	x           <- result$x
	
	out_sens <- intervention.model.combination(
		input.parameters = vparameters,
		current.int = interventions[cc]
	)
	sens <- extract_outputs(out_sens)
	
	g_all <- make_calibration_plots(base, sens, CITY)
	
	if (is.null(g_all)) stop("Plot failed for ", Jurisdictions[j])
	
	row_plots[[Jurisdictions[j]]] <- g_all
}

combined_plot <-
	wrap_plots(row_plots, ncol = 1) +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")
## Save
ggsave(
	filename = "Model_development/7.Red_strata/Calibration/Results/Plots/Check_priors.jpg",
	plot     = combined_plot,
	width    = 2.2 * length(target_name),
	height   = 1.8 * length(Jurisdictions),
	dpi      = 150
)

################ Population growth ###############
Jurisdictions <-  all.cities
# c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") #BEX, "DFW"

# Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

## Store plots as a list of rows 
myFile.name <- c("Model_development/7.Red_strata/7.Red_strata-Evidence-Inputs-Master_01292026")
WB <- openxlsx::loadWorkbook(paste0(myFile.name, ".xlsx"))

# for (j in seq_along(Jurisdictions)) {
# 		 vlist[[i]] <- subset(read.xlsx(WB, sheet = "1.1.1-PopTotal"), city == Jurisdictions)[, -1]
# }

pop_total <- read.xlsx(WB, sheet = "1.1.1-PopTotal")

library(dplyr)
# Cities with overall population increase
city_growth <- pop_total %>%
	mutate(year = as.numeric(year)) %>%
	filter(year %in% c(2016, 2045)) %>%
	group_by(city, year) %>%
	summarise(pop = sum(pe), .groups = "drop") %>%
	tidyr::pivot_wider(
		names_from = year,
		values_from = pop,
		names_prefix = "y"
	) %>%
	mutate(
		change = y2045 - y2016,
		pct_change = 100 * change / y2016
	)

city_growth<-city_growth %>%
	filter(change > 0) %>%
	arrange(desc(pct_change)) # arrange(desc(change)) #print(city_growth, n = Inf)

#Cities with population increase by ethnicity (White, Black, Hispanic)
eth_growth <- pop_total %>%
	mutate(year = as.numeric(year)) %>%
	filter(
		year %in% c(2016, 2045),
		ethnicity %in% c("w", "b", "h")
	) %>%
	group_by(city, ethnicity, year) %>%
	summarise(pop = sum(pe), .groups = "drop") %>%
	tidyr::pivot_wider(
		names_from = year,
		values_from = pop,
		names_prefix = "y"
	) %>%
	mutate(
		change = y2045 - y2016,
		pct_change = 100 * change / y2016
	)
#3. Which cities increased for each ethnicity?
eth_growth %>%
	filter(change > 0) %>%
	arrange(ethnicity, desc(change))
#4. Summary table (increase = TRUE/FALSE)
eth_growth1 <- eth_growth %>%
	transmute(
		city,
		ethnicity,
		increase = change > 0,
		pct_change = round(pct_change, 1)
	) %>%
	tidyr::pivot_wider(
		names_from = ethnicity,
		values_from = c(increase, pct_change)
	) #print(eth_growth1 %>% arrange(desc(pct_change_w)), n = Inf)

################ VALIDATING LEM PREMISE ###############
rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
CITY <- Jurisdictions <- "MAR" # MIA1, NYC1 for original, KIN, BAL, PHI -Pick any random cities, just need to get this running
current.int = c("No interventions")
j = 1

source("Model_development/7.Red_strata/7.Red_strata-LoadBaselineWorkspace.R")
source("Model_development/7.Red_strata/Calibration/7.load_targets.R") #Update Target as needed - Part of par_pop.R - Need to move up here since diag18.obs is needed for initial pop
source("Model_development/7.Red_strata/7.Red_strata-Module-Data.input.R")
intervention.model.combination <- function(input.parameters,
																					 current.int) {
	# Set TRUE/FALSE intervention parameters
	
	tf.ls <- true.false.interventions.combination(current.int)
	
	
	new.inputs <- c(input.parameters, tf.ls,
									current.int = list(current.int),
									int.start = 0, int.end = 0
	)
	
	
	out_euler <- deSolve::euler(y = x, times = vt, func = ode_model, parms = new.inputs)[, -1] # add parameter names to map
	#  outa <- array(out_euler[-1 , ], dim = c(n, 
	outa <- array(out_euler[   , ], dim = c(n+1, #If want t=0, to compared number of diagnosed PLHIV matched target
																					length(input.parameters$state.name),
																					length(input.parameters$names.gp)
	)) # out_euler[ , ] & n+1, if want to keep initial value 
	
	#  dimnames(outa)[[1]] <- 1:n # 1:(n+1)
	dimnames(outa)[[1]] <- 1:(n+1)
	dimnames(outa)[[2]] <- input.parameters$state.name
	dimnames(outa)[[3]] <- input.parameters$names.gp
	
	return(outa)
}

## 3. evaluate the model and generating the result matrix
outa <- intervention.model.combination(input.parameters = vparameters,
																						 current.int = interventions[16])
n.strata    <- apply(outa[ , c(1:16), ], c(1,3), sum)
n.strata    <- n.strata[,1:6]
# Population Plot
plot_stacked_area(time = vt, data=n.strata , name= colnames(n.strata), CITY="MAR")
# Incidence Projection plot
## Add 1 dim to match existing code:
outa <- array(
	outa,
	dim = c(1, dim(outa)),
	dimnames = c(list("1"), dimnames(outa))
)
### Run compute_inc_total,plot_target_2_median
# Incidence index mapping
inc_inds <- c(
	White = gp.gn$white,
	Black = gp.gn$black,
	Hisp  = gp.gn$hisp
)
#browser()	
# Compute totals
inc_totals <- lapply(
	inc_inds,
	compute_inc_total,
	outa        = outa,
	end_yr_ind  = end_yr_ind,
	start_yr_ind = 12
)
# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
inc_totals <- list(
	black = Reduce(`+`, inc_totals[grep("^Black", names(inc_totals))]),
	hisp  = Reduce(`+`, inc_totals[grep("^Hisp",  names(inc_totals))]),
	white = Reduce(`+`, inc_totals[grep("^White", names(inc_totals))])
) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
inc_all = inc_totals$white +	inc_totals$black +	inc_totals$hisp

inc_array <- abind::abind(inc_totals, along = 3)

g_inc <- plot_median_quantile(inc_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new infections", legend_title = "Race/Ethnicity",color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) 
g_inc

#Population sums 
pop_sum_white   <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1),	sum) # Per Xiao, denominator for inc is All pop
pop_sum_black   <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1),	sum)
pop_sum_hisp    <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1),	sum)
pop_sum_all     <- apply(	outa[ , end_yr_ind, 1:16,  ],	c(1),	sum)


##IRs 
IR_sum_white        <- inc_totals$white        / pop_sum_white      * 100000
IR_sum_black        <- inc_totals$black        / pop_sum_black      * 100000
IR_sum_hisp         <- inc_totals$hisp         / pop_sum_hisp       * 100000
IR_sum_all          <- inc_all                 / pop_sum_all        * 100000

IR_array <- abind::abind(IR_sum_black,IR_sum_hisp,IR_sum_white, along = 3)

g_ir <- plot_median_quantile(IR_array,skip_time = 0,y_limit = 100 ,dashed_year = 2025 ,strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV incidence per 100,000", legend_title = "Race/Ethnicity",
														 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	))

g_ir
