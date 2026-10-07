# Utility_functions

##Function to save a matrix to a CSV file with an option to choose the file location####
#install.packages("writexl")
# library(writexl)
# save_matrix_to_excel <- function(matrix_data) {
# 	
# 	# Construct file path
# 	save_path <- file.path("Calibration/5_Plots/Optional_util_diagnostic_code", "saved_matrix.xlsx")
# 	# Convert matrix to a data frame since writexl works with data frames
# 	matrix_df <- as.data.frame(matrix_data)
# 	# Export matrix to the chosen CSV file
# 	write_xlsx(matrix_df, path = save_path)
# 	
# 	# Return a message indicating the file has been saved
# 	cat("Matrix has been saved to", save_path, "\n")
# }


##Function to Takes two tables, Computes the % difference, Plots a blue–white–red heatmap,Excludes or caps values where |%diff| > 100 from the color scale (to prevent extreme values from skewing colors)####
plot_pct_diff_heatmap <- function(target_obs_table,
																	target_model_table,
																	output_dir = getwd(),
																	filename = "heatmap_model_vs_obs_pctdiff.png",
																	cap_value = 100,   jurisdiction_labels = NULL) {
	library(tidyverse)
	
	# ---- 1. Compute % difference ----
	pct_diff <- (target_model_table - target_obs_table) / target_obs_table * 100
	df <- as.data.frame(pct_diff)
	df$parameter <- rownames(pct_diff)
	
	# ---- 2. Cap extreme values for coloring ----
	df_long <- df %>%
		pivot_longer(
			cols = -parameter,
			names_to = "target",
			values_to = "pct_diff"
		) %>%
		mutate(
			pct_diff_capped = ifelse(pct_diff > cap_value, cap_value,
															 ifelse(pct_diff < -cap_value, -cap_value, pct_diff))
		)
	
	# ---- 3. Preserve parameter order ----
	df_long$parameter <- factor(df_long$parameter, levels = rev(rownames(pct_diff)))
	
	# ---- 4. Create heatmap ----
	p <- ggplot(df_long, aes(x = target, y = parameter, fill = pct_diff_capped)) +
		geom_tile(color = "lightgray", size = 0.2) +
		geom_text(aes(label = sprintf("%.1f%%", pct_diff)), size = 1.5) +
		scale_fill_gradient2(
			low = "#d7191c",
			mid = "white",
			high = "#2c7bb6",
			midpoint = 0,
			limits = c(-cap_value, cap_value),
			# name = sprintf("%% diff (Model - Obs)\nCapped at ±%d%%", cap_value) 
			name = "% Difference:\nModel estimate – observation" 
		) +
		labs(x = NULL, y = NULL) +
		theme_bw(base_size = 8) +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			plot.background  = element_rect(fill = "white", color = NA),
			panel.background = element_rect(fill = "white", color = NA),
			axis.text.y      = element_text(size = 6),
			axis.text.x      = element_text(size = 6, angle = 45, hjust = 1),
			legend.position  = "right"
		)
	if (!is.null(jurisdiction_labels)) {
		p <- p + scale_y_discrete(labels = jurisdiction_labels)
		p <- p + scale_x_discrete(labels = target_label)
	}
	
	# ---- 5. Save output ----
	ggsave(
		filename = file.path(output_dir, filename),
		plot     = p,
		width    = 6,
		height   = 6, #nrow(df) * 0.15
		dpi      = 300
	)
# browser()	
	# ---- 4b. Summary diagnostics ----
	
	# By parameter
	param_summary <- df_long %>%
	  group_by(parameter) %>%
	  summarise(
	    mean_abs_pct_diff = mean(abs(pct_diff), na.rm = TRUE),
	    min_abs_pct_diff  = min(abs(pct_diff), na.rm = TRUE),
	    max_abs_pct_diff  = max(abs(pct_diff), na.rm = TRUE),
	    .groups = "drop"
	  ) %>%
	  arrange(match(
	    as.character(parameter),
	    unique(as.character(df_long$parameter))
	  ))
	
	param_max <- param_summary %>% slice_max(mean_abs_pct_diff, n = 1)
	param_min <- param_summary %>% slice_min(mean_abs_pct_diff, n = 1)
	
	# By target
	target_summary <- df_long %>%
	  group_by(target) %>%
	  summarise(
	    mean_abs_pct_diff = mean(abs(pct_diff), na.rm = TRUE),
	    min_abs_pct_diff  = min(abs(pct_diff), na.rm = TRUE),
	    max_abs_pct_diff  = max(abs(pct_diff), na.rm = TRUE),
	    .groups = "drop"
	  ) %>%
	  arrange(match(
	    as.character(target),
	    unique(as.character(df_long$target))
	  ))
	
	target_max <- target_summary %>% slice_max(mean_abs_pct_diff, n = 1)
	target_min <- target_summary %>% slice_min(mean_abs_pct_diff, n = 1)
	
	# ---- Messages ----
	message("Parameter with HIGHEST avg |pct_diff|: ",
					param_max$parameter,
					" | mean=", round(param_max$mean_abs_pct_diff, 2),
					"% | min=", round(param_max$min_abs_pct_diff, 2),
					"% | max=", round(param_max$max_abs_pct_diff, 2), "%")
	
	message("Parameter with LOWEST avg |pct_diff|: ",
					param_min$parameter,
					" | mean=", round(param_min$mean_abs_pct_diff, 2),
					"% | min=", round(param_min$min_abs_pct_diff, 2),
					"% | max=", round(param_min$max_abs_pct_diff, 2), "%")
	
	message("Target with HIGHEST avg |pct_diff|: ",
					target_max$target,
					" | mean=", round(target_max$mean_abs_pct_diff, 2),
					"% | min=", round(target_max$min_abs_pct_diff, 2),
					"% | max=", round(target_max$max_abs_pct_diff, 2), "%")
	
	message("Target with LOWEST avg |pct_diff|: ",
					target_min$target,
					" | mean=", round(target_min$mean_abs_pct_diff, 2),
					"% | min=", round(target_min$min_abs_pct_diff, 2),
					"% | max=", round(target_min$max_abs_pct_diff, 2), "%")
	message(
	  "Average Parameter Dev: ",
	  round(mean(param_summary$mean_abs_pct_diff, na.rm = TRUE), 4)
	)
	message("✅ Heatmap saved to: ", file.path(output_dir, filename))
return(list(
  plot = p,
  param_summary = param_summary,
  target_summary = target_summary
))
}

##Function to output Posterior values across cities####
output_posterior_val <- function(Jurisdictions, file_prefix, save_path) {
	
	library(readr)
	library(dplyr)
	library(writexl)
	set.seed(1111)
	
	output_mat <- NULL
	
	for (j in Jurisdictions) {
		
		# ---- Build calibration filename ----
		calib_file <- paste0(file_prefix, j, ".rds")
		
		calibration_results_multiple <- readRDS(calib_file)
		
		# ---- Resample posterior ----
		# posterior_samples <- resample_using_importance_weights(calibration_results_multiple$all_samples, 
		# 																											 calibration_results_multiple$all_lls, 
		# 																											 ndraws = 1000)
		posterior_samples <- resample_using_importance_weights(calibration_results_multiple$burned_samples,
																													 calibration_results_multiple$burned_lls,
																													 ndraws = 1000)
# browser()		
		# ---- Parameter names ----
		colnames_all_samples <- paste0(
			calpar$par, "_",
			calpar$gender, "_",
			calpar$ethnicity, "_",
			calpar$risk, "_",
			calpar$sex.risk.age.group
		)
		colnames(posterior_samples) <- colnames_all_samples
		
		# ---- Descriptive statistics ----
		posterior_stats <- descriptive_stats(posterior_samples)
# browser()			
		# ---- Add jurisdiction, reorder ----
		posterior_stats$Jurisdiction <- j
		posterior_stats <- posterior_stats %>%
			select(Jurisdiction, everything())
		
		output_mat <- bind_rows(output_mat, posterior_stats)
	}
# browser()		
	# ---- Write to Excel in user-specified location ----
	write_xlsx(output_mat, save_path)
	
	return(output_mat)
}

##Function to expanded name of posterior stats####
expand_posterior_stats <- function(posterior_stats) {
	
	parts <- strsplit(as.character(posterior_stats$name), "_")
	
	expanded <- do.call(
		rbind,
		lapply(parts, function(x) {
			n <- length(x)
			tibble::tibble(
				par                = paste(x[1:(n - 4)], collapse = "_"),
				gender             = x[n - 3],
				ethnicity          = x[n - 2],
				risk               = x[n - 1],
				sex.risk.age.group = x[n]
			)
		})
	)
	
	tibble::as_tibble(
		cbind(
			posterior_stats[, c("Jurisdiction")],
			expanded,
			posterior_stats[, c("mean", "sd","median", "lower", "upper")]
		)
	)
}
