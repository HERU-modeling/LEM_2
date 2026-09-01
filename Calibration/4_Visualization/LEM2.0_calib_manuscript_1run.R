rm(list = ls())
library(openxlsx)
library(mc2d)
library(posterior)
library(ggplot2)
library(patchwork)
library(cowplot)
set.seed(1111)

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

################ CHANGE THIS AS NEEDED ###############
Jurisdictions <- c("MAR", "SFA") #, "LOR", "SAC", "SDI"
Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]
dir <- "check_156a386_10000-5000-aug26"

################ [Run this] Prepare Environment ###############
#Step 2: Choosing city
CITY <- "NJE" # -Pick any random cities, just need to get this running
#Read-in save results
calibration_results_multiple <- readRDS(paste0("Calibration/3_Results/calib_check_2c91081_10-10_MAR.rds")) # Any files in the results

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


#Step 6:

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


################ NEW plot_target_2 ###############
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
				title = unique(df$target_label),   #  IMPORTANT: label comes from df
				x = "Year",
				y = "Persons",
				colour = NULL
			) +
			
			theme_bw(base_size = 20) + #11 -  IMPORTANT: Affect "Persons" & numbers on y-axis, "Year" on x-axis & numbers
			
			theme(
				plot.title = element_text(
					size = 20,              #8 -  IMPORTANT: Affect Targets on Title, i.e., New  diagnoses, Prevalence,... & Observation/Prediction line in legend
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
################ [Run this] Target - Multiple City ###############
#Run: 

library(patchwork)
library(ggplot2)




target_name <- target_name_extended[c(1:17)]
target_label  <- target_labels[target_name_extended[c(1:17)]]
target_name_ndiag <- target_name_extended[c(4:11)]
target_label_ndiag  <- target_labels[target_name_extended[c(4:11)]]
target_name_other <- target_name_extended[c(1:3,12:17)]
target_label_other  <- target_labels[target_name_extended[c(1:3,12:17)]]

## Store plots as a list of rows 
row_plots <- list()

for (j in seq_along(Jurisdictions)) {
	
	message("Processing Target plots ", Jurisdictions[j])
	
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
#browser()	
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

	title_plot <- ggplot() +
		theme_void() +
		labs(title = Jurisdiction_names[j]) +
		theme(
			plot.title = element_text(
				hjust = 0,
				face = "bold",
				size = 14                         #14 -  IMPORTANT: Affect Jurisdictions/ City Names (i.e., Phoenix, AZ); affect CalibTargets_all_jurisdictions, not individual/report plots
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
#browser()	 #save_plot
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
		width = 20,           #Width 10 & Height 20 - Impacting Individual city report & Prevent title overlap (Increase width)
		height = 25,
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



################ [Run this] Figure : Incidence &  Rate by races (w,b,h) - Multiple City ###############
# Jurisdictions <- c("MAR", "SFA", "LOR", "SAC", "SDI", "TAM", "BMP")
# # c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL") 
# # c("MAR", "SFA", "LOR", "SDI", "DUV", "TAM", "ATL") # "SAC", "BMP", "ORF", 
# # c("COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC") 
# # c("COO", "EBR", "ORL", "BAL", "NYC") # "MIN","SUF", "WAY", "CLA", "NJE", "MEC"
# # c("CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM") 
# # all.cities
# dir <- "check_156a386_10000-5000-aug26"#"All_10000-5000-jun23"
# Jurisdiction_names <- Jurisdiction_labels[Jurisdictions]

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
	
	message("Processing Incidence plots", Jurisdictions[j])
	
	#Load Jurisdiction-specific results 
	calibration_results_multiple <- readRDS(
		paste0(
			"Calibration/3_Results/calib_",
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
file_name <- "Calibration/3_Results/Table_inc_race.rds"
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
			"Calibration/5_Plots/Individual/nInc_separate_race_",
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
			"Calibration/5_Plots/Individual/IR_separate_race_",
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
combined_plot_ir_clean <- wrap_plots(plots_ir, ncol = 4)#8

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
	filename = "Calibration/5_Plots/nInc_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_inc,
	width = 24,
	height = 12,
	dpi = 300
)
ggsave(
	filename = "Calibration/5_Plots/IR_total_separate_race_Jurisdictions.jpg",
	plot = combined_plot_ir,
	width = 24,
	height = 12,
	dpi = 300
)

