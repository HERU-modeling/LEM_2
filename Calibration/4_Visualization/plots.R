#' Plot predicted targets
#'
#' @param target_model targets generated from the model
#' @param target_obs target observed in the real world
#' @param names18 name of targets
#' @param n_target_to_plot number of targets to plots
#'
#' @return
#' @export
#'
#' @examples
plot_targets <- function(result, names18, n_target_to_plot = 144,
												 plot_dir){
	
	names18_death <- unlist(lapply(names18, function(x) paste0(x, "_death")))
	plot_names <- c(names18, names18_death)
	plot_names <- repeat_strings(strings_list = plot_names,
															 repeats_list = rep(4, n_target_to_plot))

	for(i in 1:n_target_to_plot){ # each of the 11 targets across 4 years
		bound <- c(min(result$target_model[, i]), max(result$target_model[, i]), result$target_obs[i])
		name <- plot_names[i]
		name <- gsub("/","_", name) # get rid of `/` because os treat it as folder
		
		g <- ggplot(data.frame(x = result$target_model[, i]))+
			geom_histogram(aes(x = x), bins = 100)+
			geom_vline(aes(xintercept = result$target_obs[i]), color = "red")+
			coord_cartesian(xlim = c(min(bound), max(bound)))+
			theme_bw()
		
		plot_name <- file.path(plot_dir, paste0("diag_", name, i, ".png"))
		ggsave(plot_name)
	}
}


#' Generate the Latent Class Identifiability Display
#'
#' @param calpar list of free parameters
#' @param calpar.info info of the list of free parameters
#' @param posterior_samples posterior samples of the free parameters
#' @param plot_name name of the plot to save to
#'
#' @return
#' @export
#'
#' @examples
plot_lcid <- function(calpar, calpar.info, posterior_samples,
											plot_name = "Data_Files/Result/Calibration/lcid.png"){
	
	prior_samples <- sample_dist_wrapper(500, calpar)
	
	colnames_df <- append_los_with_index(repeat_strings(calpar.info$names, calpar.info$plength))
	
	df_post_samp <- data.frame(posterior_samples)
	colnames(df_post_samp) <- colnames_df
	df_post_samp$where <- "post"
	
	df_prior_samp <- data.frame(prior_samples)
	colnames(df_prior_samp) <- colnames_df
	df_prior_samp$where <- "prior"
	
	df <- rbind(df_post_samp, df_prior_samp)
	
	df <- df %>%
		pivot_longer(!where, names_to = "parameter", values_to = "value")
	
	lcid <- ggplot(df, aes(x = value))+
		facet_wrap(vars(parameter), scales = "free")+
		geom_density(aes(x = value, y = after_stat(density), fill = where), alpha = 0.7)
	
	ggsave(plot_name, width = 20, height = 20)
}





#' Plot targets
#'
#' @param results results from inference
#' @param names18 names of the 18 targets
#' @param gp18.gn index of subgroups within the 18 targets
#' @param which_group if `na`, plot all subgroups; if it's an element in `gp18.gn`,
#' combine corresponding subgroups
#' @param years years of calibration targets
#' @param plot_dir directory to save plots
#'
#' @return
#' @export
#'
#' @examples
plot_target_series <- function(results, names18, gp18.gn, which_group,
															 years = c(2012, 2013, 2014, 2015),
															 plot_dir = "Model_development/7.Red_strata/Calibration/timeplots",
															 comparator = NULL){
	
	type <- match.arg(which_group, choices = c("na", names(gp18.gn)))
	
	# Death counts not included
	plot_names <- c(names18)
	timespan <- length(years)
	
	# Death counts not included
	n_targets <- ifelse(type == "na", length(names18), 1)
	get_target_indices <- function(i) rep(i*timespan, timespan) + (-timespan+1):0
	
	for(i in 1:n_targets){
		
		if (type != "na"){ # For a particular subgroup
			get_all_target_indices <- function(is) t(sapply(is, get_target_indices))
			
			target_indices_matrix <- matrix(get_all_target_indices(unlist(gp18.gn[which_group])),
																			nrow = length(unlist(gp18.gn[which_group])),
																			ncol = 4)
			
			combine_targets <- function(i){
				m <- matrix(i[target_indices_matrix], nrow = nrow(target_indices_matrix),
										ncol = ncol(target_indices_matrix))
				return(colSums(m))
			}
			
			ith_target_mat <- t(apply(results$target_model, MARGIN = 1, FUN = combine_targets))
			ith_target_obs <- apply(t(matrix(results$target_obs)), MARGIN = 1, FUN = combine_targets)
		
			
			save_name <- file.path(plot_dir, "subgroups", paste0(type, ".png"))
			plot_name <- type
		
		} else { # Else plot all subgroups
			ith_target_cols <- get_target_indices(i)
			ith_target_mat <- results$target_model[, ith_target_cols]
			ith_target_obs <- results$target_obs[ith_target_cols]
			save_name <- file.path(plot_dir, "all_groups", gsub("/","_", paste0(plot_names[i], ".png")))
			plot_name <- gsub("/","_", plot_names[i])
		}
		
		ith_mean <- apply(ith_target_mat, MARGIN = 2, FUN = mean)
		ith_sd <- apply(ith_target_mat, MARGIN = 2, FUN = sd)
		
		df_to_plot <- data.frame(time = years,
														 obs = ith_target_obs,
														 pred_mean = ith_mean,
														 upper = ith_mean + 2*ith_sd,
														 lower = ith_mean - 2*ith_sd)
		
		plot_one_time_serie(df_to_plot, plot_name, save_name, years)
	}
	return(df_to_plot)
}

#' Plot targets - LEM 2.0 Only
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
plot_target_2  <- function(results, target_name ,
												 #names18, gp18.gn, which_group,
												 years = c(2017, 2018, 2019, 2020, 2021, 2022),
												 plot_dir = "Model_development/7.Red_strata/Calibration/timeplots",
												 comparator = NULL, y_from_zero = FALSE){
	
	#  	type <- match.arg(which_group, choices = c("na", names(gp18.gn)))
	
	# Death counts not included
	#  	plot_names <- c(names18)
	plot_names <- c(target_name) 
	timespan <- length(years)
	#browser() 	
	# Death counts not included
	#  	n_targets <- ifelse(type == "na", length(names18), 1)
	n_targets <- length(target_name)
	get_target_indices <- function(i) rep(i*timespan, timespan) + (-timespan+1):0
	
	for(i in 1:n_targets){
		
		# 		if (type != "na"){ # For a particular subgroup
		# 			get_all_target_indices <- function(is) t(sapply(is, get_target_indices))
		
		# 			target_indices_matrix <- matrix(get_all_target_indices(unlist(gp18.gn[which_group])),
		# 																			nrow = length(unlist(gp18.gn[which_group])),
		# 																			ncol = 4)
		
		# 			combine_targets <- function(i){
		#  				m <- matrix(i[target_indices_matrix], nrow = nrow(target_indices_matrix),
		#  										ncol = ncol(target_indices_matrix))
		#  				return(colSums(m))
		#  			}
		
		#  			ith_target_mat <- t(apply(results$target_model, MARGIN = 1, FUN = combine_targets))
		#  			ith_target_obs <- apply(t(matrix(results$target_obs)), MARGIN = 1, FUN = combine_targets)
		
		
		#  			save_name <- file.path(plot_dir, "subgroups", paste0(type, ".png"))
		#  			plot_name <- type
		
		#  		} else { # Else plot all subgroups
		ith_target_cols <- get_target_indices(i)
		ith_target_mat <- results$target_model[, ith_target_cols]
		ith_target_obs <- results$target_obs[ith_target_cols]
		save_name <- file.path(plot_dir, "all_groups", gsub("/","_", paste0(plot_names[i], ".png")))
		plot_name <- gsub("/","_", plot_names[i]) #plot_names <- c(names18)
		# 		}
		
		ith_mean <- apply(ith_target_mat, MARGIN = 2, FUN = mean)
		ith_sd <- apply(ith_target_mat, MARGIN = 2, FUN = sd)
		
		df_to_plot <- data.frame(time = years,
														 obs = ith_target_obs,
														 pred_mean = ith_mean,
														 upper = ith_mean + 2*ith_sd,
														 lower = ith_mean - 2*ith_sd)
		
		plot_one_time_serie(df_to_plot, plot_name, save_name, years, y_from_zero)
	}
	return(df_to_plot)
}


#' Plot time serie plot of one dataframe
#'
#' @param df_to_plot dataframe to plot. Column consists of years, observed targets,
#' predicted sd and mean of targets, 95% confidence interval of the predicted targets
#' @param plot_name title of the plot
#' @param save_name where the plot should be saved
#' @param years years of calibration targets
#'
#' @return
#' @export
#'
#' @examples
plot_one_time_serie <- function(df_to_plot, plot_name, save_name, years,
																y_from_zero = FALSE){
	g <- ggplot(df_to_plot)+
		geom_line(aes(x = time, y = obs, colour = "Observation"))+
		geom_line(aes(x = time, y = pred_mean, colour = "Prediction"))+
		geom_errorbar(aes(x = time, ymin = lower, ymax = upper), width=0.2)+
		scale_color_manual(name = "Target Type", 
											 values = c("Observation" = "darkblue", "Prediction" = "red"))+
		labs(title = plot_name, x = "Year", y = "Person")+
		theme_bw()
	# Add y-axis starting at 0 only if requested
	if (y_from_zero) {
		g <- g + scale_y_continuous(limits = c(0, NA))
	}
	print(g)
	ggsave(save_name)
}



#' Plot multi chain trace for multiprocessing
#'
#' @param multichain_results output of the multichain metropolis hasting function
#'
#' @return multi-chain trace plot
#' @export
#'
#' @examples
plot_multichain_trace <- function(multichain_results){
  traces <- data.frame(multichain_results)
  traces$iterations <- 1:nrow(traces)
  
  p <- traces %>%
    pivot_longer(!iterations, names_to = "chain", values_to = "value") %>%
    ggplot()+
    geom_line(aes(x = iterations, y = value, group = chain, color = chain))+
    theme_bw()
  
  p
}



#' Generate the plot comparing parameters' prior and posterior distributions
#'
#' @param prior_samples prior samples, output from Likelihood Weighting & Metropolis-Hasting
#' algorithm
#' @param posterior_samples posterior samples, output from Likelihood Weighting & Metropolis-Hasting
#' algorithm
#' @param calpar data frame of free parameters
#' @param adjust adjust the smoothness of the posterior plots
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_prior_posterior_comparison <- function(prior_samples, posterior_samples, calpar,
                                            adjust = 1){
#browser()  
  posterior_samples <- data.frame(posterior_samples)
  posterior_samples$type <- "posterior"
  prior_samples <- data.frame(prior_samples)
  prior_samples$type <- "prior"
  
  all_samples <- rbind(prior_samples, posterior_samples)
  
  colnames_all_samples <- paste0(calpar$par, "_", calpar$gender, "_", calpar$ethnicity, "_",
                                 calpar$risk, "_", calpar$sex.risk.age.group)
  colnames(all_samples) <- c(colnames_all_samples, "type")
  
  g1 <- all_samples %>%
    pivot_longer(!type, names_to = "Parameters", values_to = "Value")%>%
    ggplot()+
    facet_wrap(vars(Parameters), scales = "free")+
    geom_density(aes(x = Value, color = type), adjust = adjust)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_blank())
  
  g2 <- all_samples %>%
  	pivot_longer(!type, names_to = "Parameters", values_to = "Value") %>%
  	ggplot(aes(x = Value, color = type)) +
  	facet_wrap(vars(Parameters), scales = "free") +
  	geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  	theme_bw() +
  	theme(
  		strip.text = element_text(size = 10),
  		axis.text.x = element_blank()
  	)
  return(list(g1 = g1, g2 = g2))
}


#' Plot multichain traceplot for Gibbs within Metropolis-Hasting
#'
#' @param multichain_results output of `metropolis_hasting_multichain_combine`
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_multichain_trace <- function(multichain_results){
  traces <- data.frame(multichain_results)
  traces$iterations <- 1:nrow(traces)
  
  p <- traces %>%
    pivot_longer(!iterations, names_to = "chain", values_to = "value") %>%
    ggplot()+
    geom_line(aes(x = iterations, y = value, group = chain, color = chain))+
    theme_bw()
  
  p
}



#' Plot new diagnoses comparison between LEM 1.0 and LEM 2.0
#'
#' @param lem1_cal_result result of LEM 1.0
#' @param lem2_cal_result result of LEM 2.0
#' @param gp18.gn index of the 18 strata
#' @param names18 name of the 18 strata
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_lem_comparison <- function(lem1_cal_result, lem2_cal_result, 
                                gp18.gn, names18){
  gp18.gn$all <- 1:18
  
  lem2_targets_output <- plot_target_series(results = lem2_cal_result, 
                                            names18 = names18, gp18.gn = gp18.gn, 
                                            which_group = "all")
  
  lem1_targets_output <- plot_target_series(results = lem1_cal_result, 
                                            names18 = names18, gp18.gn = gp18.gn, 
                                            which_group = "all")
  
  lem1_targets_output <- lem1_targets_output %>%
    rename("time_gv1" = "time")%>%
    rename("obs_gv1" = "obs")%>%
    rename("Cali1.0_pred_mean" = "pred_mean")%>%
    rename("Cali1.0_upper" = "upper") %>%
    rename("Cali1.0_lower" = "lower")
  
  g <- ggplot(cbind(lem1_targets_output, lem2_targets_output))+
    geom_point(aes(x = time, y = obs)) + 
    geom_line(aes(x = time, y = pred_mean), color = "#A31223") + 
    geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.1, fill = "#A31223")+
    geom_line(aes(x = time, y = Cali1.0_pred_mean), color = "#87AEDC") + 
    geom_ribbon(aes(x = time, ymin = Cali1.0_lower, ymax = Cali1.0_upper), alpha = 0.1, fill = "#87AEDC")+
    theme_bw()+
    labs(title = "Calibration 1.0 (blue) vs Calibration 2.0 (red)", x = "Year", y = "Person")
  
  return(g)
}

plot_one_time_serie_distribution <- function(
		df_to_plot,
		plot_name,
		save_name,
		years,
		y_from_zero = FALSE
){
	
	# Distribution 1:
	# Violin + boxplot
	
	g1 <- ggplot(df_to_plot) +
		
		# violin distribution
		geom_violin(
			aes(x = factor(time), y = value),
			fill = "lightblue",
			alpha = 0.5,
			trim = FALSE
		) +
		
		# median point
		stat_summary(
			aes(x = factor(time), y = value),
			fun = median,
			geom = "point",
			colour = "blue",
			size = 3
		) +
		
		# 95% credible interval
		stat_summary(
			aes(x = factor(time), y = value),
			fun.data = function(x){
				
				data.frame(
					y    = median(x),
					ymin = quantile(x, 0.025),
					ymax = quantile(x, 0.975)
				)
				
			},
			geom = "errorbar",
			width = 0.5,
			colour = "blue"
		) +
		
		# observed value
		geom_point(
			aes(x = factor(time), y = obs),
			colour = "red",
			size = 3
		) +
		
		labs(
			title = paste0(plot_name, " - Distribution 1"),
			x = "Year",
			y = "Value"
		) +
		
		theme_bw()
	
	# Add y-axis starting at 0 only if requested
	if (y_from_zero) {
		g1 <- g1 + scale_y_continuous(limits = c(0, NA))
	}
	
	print(g1)
	
	# save violin/boxplot
	ggsave(
		filename = save_name,
		plot = g1,
		width = 8,
		height = 5
	)
	
	
	# Distribution 2:
	# Histogram by year
	
	g2 <- ggplot(df_to_plot, aes(x = value)) +
		
		geom_histogram(
			bins = 40,
			fill = "lightblue",
			colour = "black"
		) +
		
		facet_wrap(
			~time,
			scales = "free"
		) +
		
		labs(
			title = paste0(plot_name, " - Distribution 2"),
			x = "Value",
			y = "Count"
		) +
		
		theme_bw()
	
	print(g2)
	
	# histogram save name
	save_name_hist <- gsub(
		"distribution1",
		"distribution2",
		save_name
	)
	
	# save histogram
	ggsave(
		filename = save_name_hist,
		plot = g2,
		width = 10,
		height = 6
	)
}



plot_target_distribution <- function(
		results,
		target_name,
		years = c(2017, 2018, 2019, 2020, 2021, 2022),
		plot_dir = "Model_development/7.Red_strata/Calibration/timeplots",
		y_from_zero = FALSE
){
	
	plot_names <- c(target_name)
	
	timespan  <- length(years)
	n_targets <- length(target_name)
	
	get_target_indices <- function(i) {
		rep(i * timespan, timespan) + (-timespan + 1):0
	}
	
	for(i in 1:n_targets){
		
		ith_target_cols <- get_target_indices(i)
		
		# 5000 x 6 matrix
		ith_target_mat <- results$target_model[, ith_target_cols]
		
		# observed values
		ith_target_obs <- results$target_obs[ith_target_cols]
		
		# save name for distribution1
		save_name <- file.path(
			plot_dir,
			"all_groups",
			gsub(
				"/",
				"_",
				paste0(plot_names[i], "_distribution1.png")
			)
		)
		
		plot_name <- gsub("/", "_", plot_names[i])
		
		# convert matrix to long format
		df_to_plot <- data.frame(
			time  = rep(years, each = nrow(ith_target_mat)),
			value = as.vector(ith_target_mat)
		)
		
		# observed values
		df_obs <- data.frame(
			time = years,
			obs  = ith_target_obs
		)
		
		# merge observed values
		df_to_plot <- merge(df_to_plot, df_obs, by = "time")
		
		# create plots
		plot_one_time_serie_distribution(
			df_to_plot  = df_to_plot,
			plot_name   = plot_name,
			save_name   = save_name,
			years       = years,
			y_from_zero = y_from_zero
		)
	}
	
	return(df_to_plot)
}
