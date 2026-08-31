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
compute_inc_total <- function(outa, inc_inds, end_yr_ind, start_yr_ind = 12) {
	
	# Incidence differences
	x_diff <- apply(
		outa[, end_yr_ind, c(18) , inc_inds, drop = FALSE], # for PREMISE col 17 (of out/outa) is equal 18 in LEM, i.e., inc_bs
		c(1, 3, 4),
		diff
	)    # Outa: sample*time*state*strata; apply diff over time(2) dimension , keep other (1,3,4) same -> x_diff here (before aperm): time*sample*state*strata
	
	x_diff <- aperm(x_diff, c(2, 1, 3, 4)) # x_diff (after aperm): sample*time*state*strata (again)
	
	n_sample  <- dim(x_diff)[1]
	n_state   <- dim(x_diff)[3]
	n_strata  <- dim(x_diff)[4]
	
	# First year
	x0 <- outa[, start_yr_ind, c(18) , inc_inds, drop = FALSE]     # x0: sample*time*state*strata, drop = False to maintain dim
	x0 <- array(x0, dim = c(n_sample, 1, n_state, n_strata))
	
	# Full time series
	x_full <- abind(x0, x_diff, along = 2)                  # x_full: sample*time*state*strata, bind along "time" dimension (2)
	
	# Sum over strata + state
	apply(x_full, c(1, 2), sum)	                            # sum over everything while keeping sample (1) & time (2) dimensions
	# browser()
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
