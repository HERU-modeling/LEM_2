library(tictoc)
library(ggplot2)
library(tidyverse)
library(progress)
library(mc2d)
library(posterior)

library(openxlsx)
library(abind)
library(deSolve)

library(tictoc)


#' Create information table for free parameters
#'
#' @param calpar free parameters
#'
#' @return `calpar.info`, summary of the free parameters
#' @export
#'
#' @examples
create_calpar_info <- function(calpar){
	calpar.info = as.list(2)      #Contains information for free.par
	calpar.info$names   = unique(calpar$par)
	calpar.info$plength = as.numeric(table(factor(calpar$par, levels=calpar.info$names)))
	return(calpar.info)
}

#' Wrapper for the function to sample from distribution
#'
#' @param nsample number of samples to take for each parameter
#' @param calpar dataframe of all free parameters of a given city
#'
#' @return prior samples which is a nsample by npar matrix
#' @export
#'
#' @examples
sample_dist_wrapper <- function(nsample, calpar, prior_weight = FALSE){
  npar <- nrow(calpar)
  randsp <- matrix(0, nrow = nsample, ncol = npar)
  sample_weights <- matrix(0, nsample, npar)
#browser()
  for (i in 1:npar){
# browser(expr= {i >= 24}) #print(calpar, n = Inf)
#browser()  	
  	sample_output <- sample_dist(calpar$dist[i], calpar$no.par[i], 
                                calpar[c('par1', 'par2', 'par3', 'lower', 'upper')][i, ],
                                nsample, prior_weight = prior_weight)
  	
  	if (any(is.nan(sample_output$weight) | is.infinite(sample_output$weight) | sample_output$weight == 0)) {   
  		cat("sample_output$weight is NaN, Inf, or numerically underflowed to 0 — entering debug mode\n")  
  		browser()
  	}  
  	
  	# adding this try catch
  	tryCatch({randsp[, i] <- sample_output$sample},
  	         error = function(e){
  	           print("error occured!!!")
  	           print(calpar[i,])
  	           print("$$$")
  	           print(sample_output)
  	         })

  	if(prior_weight) sample_weights[, i] <- sample_output$weight
  }

  if(prior_weight){
    return(list(randsp = randsp, sample_weights = sample_weights))
  } else {
  	return(randsp)
  }
}

#' Sample from distribution
#'
#' @param param_dist distribution of random variable to sample
#' @param param_npar number of parameters that parametrize the prior distribution
#' of the free parameter
#' @param param_pars the row corresponding to the free parameter to sample from
#' @param nsample number of samples to sample from the distribution
#'
#' @return prior samples of the target free parameter
#' @export
#'
#' @examples
sample_dist <- function(param_dist, param_npar, param_pars, nsample, prior_weight = FALSE){
  output <- c()
  suppressWarnings({param_pars <- lapply(param_pars, as.numeric)}) # par3 is often NA, therefore this gives a warning
  if (param_npar != 3){
    if (param_dist == "beta"){
#browser()       	
      output = rbeta(nsample, param_pars$par1, param_pars$par2)
      if(prior_weight) weight <- dbeta(output, param_pars$par1, param_pars$par2, log = TRUE)
    }
    else if (param_dist == "pert") {
#browser()    	
      output = rbetagen(nsample, param_pars$par1, param_pars$par2, 
      									param_pars$lower, param_pars$upper)
      ### Clamp samples to stay strictly within bounds (avoid boundary values that cause infinite densities
      eps <- 1e-12
      output <- pmin(pmax(output, param_pars$lower + eps),
      							 param_pars$upper - eps) 
      ###
      if(prior_weight) weight <- dbetagen(output, 
      																		param_pars$par1, 
      																		param_pars$par2, 
      																		param_pars$lower, param_pars$upper, log = TRUE)
    }
    else if (param_dist == "1/pert") {
#browser()    	
    	sample <- rbetagen(nsample, param_pars$par1, param_pars$par2, 
    										 1/param_pars$upper, 1/param_pars$lower)
      output = 1/sample
      if(prior_weight) weight <- dbetagen(sample, param_pars$par1, param_pars$par2, 
      																		1/param_pars$upper, 1/param_pars$lower, log = TRUE)
    }
    else if (param_dist == "gamma") {
#browser()    	
      output = rgamma(nsample, param_pars$par1, 1/param_pars$par2)
      if(prior_weight) weight <- dgamma(output, param_pars$par1, 1/param_pars$par2, log = TRUE)
    }
    else if (param_dist == "uni") {
#browser()    	
      output = runif(nsample, param_pars$lower, param_pars$upper)
      if(prior_weight) weight <- dunif(output, param_pars$lower, param_pars$upper, log = TRUE)
#browser()      
    }
    else if (param_dist == "1/uni") {
# browser()    	
    	sample <- runif(nsample, 1/param_pars$upper, 1/param_pars$lower)
      output = 1/(sample)
      if(prior_weight) weight <- dunif(sample, 1/param_pars$upper, 1/param_pars$lower, log = TRUE)
    }
    else if (param_dist == "ln") {
#browser()
      output = rlnorm(nsample, param_pars$par1, param_pars$par2)
      if(prior_weight) weight <- dlnorm(output, param_pars$par1, param_pars$par2, log = TRUE)
    }
    else if (param_dist == "1-lognormal") {
#browser()    	
    	sample <- rlnorm(nsample, as.numeric(param_pars$par1), as.numeric(param_pars$par2))
      output = 1- sample
      if(prior_weight) weight <- dlnorm(sample, as.numeric(param_pars$par1), as.numeric(param_pars$par2), log = TRUE)
    }
  	else if (param_dist == "1 - lognormal") {
#browser()  		
  		sample <- rlnorm(nsample, as.numeric(param_pars$par1), as.numeric(param_pars$par2))
  		output = 1- sample
  		if(prior_weight) weight <- dlnorm(sample, as.numeric(param_pars$par1), as.numeric(param_pars$par2), log = TRUE)
  	}
    else if (param_dist == "1/ln") {
#browser()    	
    	sample <- rlnorm(nsample, as.numeric(param_pars$par1), as.numeric(param_pars$par2))
    	output = 1- sample
    	if(prior_weight) weight <- dlnorm(sample, as.numeric(param_pars$par1), as.numeric(param_pars$par2), log = TRUE)
    }
#    else if (param_dist == "poisson") {
    else if (param_dist == "pois") {
#browser()
    	sample <- rpois(nsample, param_pars$par1)
      output = sample / param_pars$par2
      if(prior_weight) weight <- dpois(sample, param_pars$par1, log = TRUE)
    }
  }
#  else if (param_dist == "dirichlet_sp"){
  else if (param_dist == "dir"){
#browser()  	
    dirsm        = rdirichlet(nsample, as.numeric(c(param_pars$par1, param_pars$par2, param_pars$par3)))
    output = dirsm[,1]*0.3 + dirsm[ ,2]*0 + dirsm[ ,3]*1 # what's with the weight?
#    if(prior_weight) weight <- ddirichlet(dirsm, as.numeric(c(param_pars$par1, param_pars$par2, param_pars$par3)), log = TRUE)
    if(prior_weight) weight <- log(ddirichlet(dirsm, as.numeric(c(param_pars$par1, param_pars$par2, param_pars$par3))))
  }
  
  if(prior_weight){
  	return(list(sample = output, weight = weight))
  }else{
  	return(list(sample = output))
  }
 
}


#' Resample using importance weights
#'
#' @param samples samples drawn from prior of free parameters
#' @param weights corresponding weights for each prior sample
#' @param ndraws number of draws to resample
#'
#' @return samples that estimate the probability density function (PDF) of 
#' posterior distribution of free parameters
#' @export
#'
#' @examples
resample_using_importance_weights <- function(samples, weights, ndraws = 300, burn_frac = 0.5){
  if(typeof(samples) == "list"){
    nsamples <- length(samples)
  }else if(typeof(samples) == "double"){
    nsamples <- nrow(samples)
  }
# #browser()	
  weights <- exp(weights)
  weights <- weights/sum(weights)

  # calculate effective sample size: Effective Sample Size for Importance Sampling based on discrepancy measures L. Martino
  nweights <- length(weights)
  ESS_IS <- 1/sum(weights^2)
  print(paste0("ESS_IS: ", ESS_IS))
  print(paste0("Fraction ESS_IS: ", ESS_IS/nweights))
 
  # Resampling Based on Weights
  indices <- posterior::resample_draws(as_draws(as.matrix(1:nsamples)), weights = weights,
                            method = "simple", ndraws = ndraws)
#browser()   
  if(typeof(samples) == "list"){
    resamples <- samples[indices]
    npar <- length(resamples[[1]])
    resamples <- t(matrix(unlist(resamples), nrow = npar, ncol =ndraws))
  }else if(typeof(samples) == "double"){
    resamples <- samples[indices, ]
  }
  return(resamples)
}


#' Repeat a list of strings a given number of times for each string in the list
#'
#' @param strings_list list of strings
#' @param repeats_list number of times each string should repeat
#'
#' @return
#' @export
#'
#' @examples repeat_strings(c("r1, r3, r5"), c(1, 3, 5))
repeat_strings <- function(strings_list, repeats_list) {
  repeated_strings <- character()
  
  for (i in seq_along(strings_list)) {
    repeated_strings <- c(repeated_strings, rep(strings_list[i], repeats_list[i]))
  }
  
  return(repeated_strings)
}

append_los_with_index <- function(los){
  sapply(1:length(los), FUN = function(x) paste0(x, "_ ",  los[x]))
}


#' Create table for the marginalized posterior distributions
#'
#' @param new_samples samples that estimate the posterior distribution
#'
#' @return
#' @export
#'
#' @examples
create_posterior_table <- function(new_samples, npar, names){
  output_matrix <- matrix(nrow = npar, ncol = 5) # mean, sd, 0.05, median, 0.95
  for(i in 1:npar){
    param <- new_samples[, i]
    output_matrix[i, ] <- c(mean(param), sd(param), unname(quantile(param, probs = c(0.05, 0.5, 0.95))))
  }
  output_df <- as.data.frame(output_matrix)
  colnames(output_df) <- c("mean", "sd", "lower", "median", "upper")
  output_df$name <- repeat_strings(calpar.info$names, calpar.info$plength) # problem here
  # output_df$name <- append_los_with_index(names)
  output_df <- output_df[c(6,1,2,3,4,5)] # reorder
  return(output_df)
}



#' Run model with prior/posterior samples
#'
#' @param randsp prior/posterior samples
#' @param calpar.info info of vector of parameters to be calibrated
#' @param vparameters 
#' @param vlist 
#' @param calib_target calibration targets
#' @param valid.target validation targets
#' @param mode either `ll`, which gives likelihood, or `target`, which gives the
#' model predicted and observed target
#' @param num_targets size of the output target
#'
#' @return
#' @export
#'
#' @examples
run_model_with_sample <- function(randsp, calpar.info, vparameters, vlist, calib_target,
                                  valid.target, mode, num_targets = 144, verbal = TRUE){
	type <- match.arg(mode, choices = c("ll", "target"))
#browser()	
	output_ll <- rep(NA, nrow(randsp))
	output_model <- matrix(nrow = nrow(randsp), ncol = num_targets)
	output_obs <- matrix(nrow = nrow(randsp), ncol = num_targets)
	output_outa <- array(NA, dim = c(nrow(randsp), 348, 28, 30))
	if(verbal){
		pb <- progress_bar$new(format = "  Running ODE [:bar] :percent, eta: :eta",
													 total = nrow(randsp))
	}
	for(i in 1:nrow(randsp)){
		if(verbal){pb$tick()}
#browser()		
		obj_output <- obj(calib.par = randsp[i, ], calpar.info = calpar.info, 
											fixed = vparameters, fixed.list =vlist, calib.target = calib.target, 
											valid.target = valid.target)
#browser()	
		output_ll[i] <- obj_output$loglik
		output_model[i,] <- obj_output$model
		output_obs[i,] <- obj_output$obs
		output_outa[i,,,] <- obj_output$outa #Output for diagnostic plot, optional for calibration
	}
		
	result <- switch(type,
				 ll = output_ll,
				 target = list(model = output_model, obs = output_obs[1,], outa = output_outa))# 

	return(result)
}



#' Store posterior table for prediction
#' 
#' Store the posterior table in a similar format as `calpar`, so that the code
#' for incidence projection could be used directly
#'
#' @param posterior posterior distributions
#' @param calpar list of free parameters
#' @param loc location to store the posterior table
#'
#' @return
#' @export
#'
#' @examples
store_posterior <- function(posterior, calpar, loc){
  post_table <- calpar
  post_table <- calpar %>%
    select(c(par, gender, ethnicity, risk, sexual.intensity, cd4, pe, lower, upper,
             dist, no.par)) %>%
    mutate(lower = posterior$lower, upper = posterior$upper, pe = posterior$median)
  write.csv(post_table, loc)
  return(post_table)
}


#' Add parameter to the list of free parameters
#'
#' @param names list of parameters that we wanted to add to calibrate
#' @param calpar list of free parameters
#' @param par_info info of the list of free parameters
#' @param WB location of excel file containing information of all parameters
#' @param CITY city in which parameters are calibrated
#'
#' @return updated `calpar` and `calpar.info`
#' @export
#'
#' @examples
add_param <- function(names, par_info, WB, CITY = "SEA"){
#  aggr_params <- data.frame(matrix(nrow = 0, ncol = 15))
  aggr_params <- data.frame(matrix(nrow = 0, ncol = 16))
  colnames(aggr_params) <- c("par", "gender", "ethnicity", "risk", "sex.risk.age.group","cd4",
  													 "year", "pe", "lower", "upper", "dist", 
  													 "no.par", "dist_fit", "par1", "par2", "par3")
 
  for(name in names){
    tab_loc <- par_info$tab[match(name, par_info$parameter)]
    tab <- read.xlsx(WB, sheet=tab_loc)
#browser()     
    param_to_add <- tab %>%
      filter(city == CITY | city == "na") %>%
      mutate(par = name) 
    
#browser()     
    add_col_na_to_df <- function(df, cols){
      for(col in cols){
        if(!col %in% colnames(df)){
          df[col] = "na"
        }
      }
      return(df)
    }
     
    param_to_add <- add_col_na_to_df(param_to_add, c("cd4", "sex.risk.age.group", "year",
                                                     "risk", "gender", "ethnicity"))
    
    param_to_add <- select(param_to_add, c(par, gender, ethnicity, risk, cd4, sex.risk.age.group, year, pe, lower, upper, dist, 
                                           no.par, dist_fit, par1, par2, par3))
    # param_to_add <- tab %>%
    #   filter(city == CITY | city == "na") %>%
    #   mutate(cd4 = ifelse("cd4" %in% names(.), .data$cd4, "na")) %>%
    #   mutate(sexual.intensity = ifelse("sexual.intensity" %in% names(.), .data$sexual.intensity, "na")) %>%
    # 	mutate(risk = ifelse("risk" %in% names(.), .data$risk, "na")) %>%
    # 	mutate(gender = ifelse("gender" %in% names(.), .data$gender, "na")) %>%
    # 	mutate(ethnicity = ifelse("ethnicity" %in% names(.), .data$ethnicity, "na")) %>%
    #   mutate(par = name) %>%
    #   select(c(par, gender, ethnicity, risk, cd4, sexual.intensity, pe, lower, upper, dist, 
    #            no.par, dist_fit, par1, par2, par3))
    aggr_params <- distinct(rbind(aggr_params, param_to_add))
  }
  
  calpar.info <- create_calpar_info(aggr_params)
#browser() #Check "sex.risk.age.group", "year" is working
  list(calpar = aggr_params,
       calpar.info = calpar.info)
}



#' Return system time
#'
#' @return
#' @export
#'
#' @examples
add_timestamp <- function(){
	format(Sys.time(), "%Y-%m-%d_%H-%M")
}


#' Descriptive stat for samples
#'
#' @return
#' @export
#'
#' @examples
descriptive_stats <- function(samples) {
	samples <- as.data.frame(samples)
	
	# Keep original column order
	col_order <- names(samples)
	
	stats_long <- samples %>%
		pivot_longer(everything()) %>%
		group_by(name) %>%
		summarise(
			mean  = mean(value),
			sd    = sd(value),
			median = quantile(value, 0.5),
			lower = quantile(value, 0.025),
			upper = quantile(value, 0.975),
			.groups = "drop"
		)%>%
		mutate(name = factor(name, levels = col_order)) %>%
		arrange(name)
	# Print all rows
	print(stats_long, n = nrow(stats_long))
	
	# Return the table
	return(stats_long)
}

#' Incoporate burn-in period
#'
#' @return
#' @export
#'
#' @examples
create_burnin <- function(all_samples,
													all_lls,
													burn_frac = 0.5) {
	
	stopifnot(length(all_samples) == length(all_lls))
	
	n_chains <- length(all_lls)
	
	burned_lls <- vector("list", n_chains)
	burned_samples <- vector("list", n_chains)
	
	for (ch in seq_len(n_chains)) {
		
		chain_ll  <- all_lls[[ch]]
		chain_smp <- all_samples[[ch]]
		
		n_iter <- length(chain_ll)
		burn_in <- floor(n_iter * burn_frac)
		
		keep_idx <- (burn_in + 1):n_iter
		
		burned_lls[[ch]]     <- chain_ll[keep_idx]
		burned_samples[[ch]] <- chain_smp[keep_idx]
	}
	
	list(
		all_lls = burned_lls,
		all_samples = burned_samples,
		n_iter = length(burned_lls[[1]])
	)
}

#' Create Sigma - Correlation between target
#'
#' @return
#' @export
#'
#' @examples
make_sigma_blockwise <- function(var_vec,
																 rho_counts = 0.3,
																 rho_prop   = 0.5) {
	
	# Indices (based on your structure)
	idx_diag_race  <- 1:3
	idx_ndiag_cal  <- 4:8
	idx_ndiag_race <- 9:11
	idx_p_late     <- 12:14
	idx_p_viral    <- 15:17
	
	K <- length(var_vec)
	Sigma <- matrix(0, K, K)
	
	# --- Counts blocks ---
	Sigma[idx_diag_race, idx_diag_race] <-
		make_block_sigma(var_vec[idx_diag_race], rho_counts)
	
	Sigma[idx_ndiag_cal, idx_ndiag_cal] <-
		make_block_sigma(var_vec[idx_ndiag_cal], rho_counts)
	
	Sigma[idx_ndiag_race, idx_ndiag_race] <-
		make_block_sigma(var_vec[idx_ndiag_race], rho_counts)
	
	# --- Proportion blocks ---
	Sigma[idx_p_late, idx_p_late] <-
		make_block_sigma(var_vec[idx_p_late], rho_prop)
	
	Sigma[idx_p_viral, idx_p_viral] <-
		make_block_sigma(var_vec[idx_p_viral], rho_prop)
	
	return(Sigma)
}

#' compute split_rhat or scale reduction factor 
#'
#' @return
#' @export
#'
#' @examples
compute_split_rhat <- function(chains_array) {
	# chains_array: m x n x p
	
	m <- dim(chains_array)[1]
	n <- dim(chains_array)[2]
	p <- dim(chains_array)[3]
	
	# split chains
	half_n <- floor(n/2)
	split_array <- array(NA, dim = c(2*m, half_n, p))
	
	for (j in 1:m) {
		split_array[j, , ]       <- chains_array[j, 1:half_n, ]
		split_array[j+m, , ]     <- chains_array[j, (half_n+1):(2*half_n), ]
	}
	
	m2 <- dim(split_array)[1]
	n2 <- dim(split_array)[2]
	
	rhat <- numeric(p)
	
	for (k in 1:p) {
		chain_matrix <- split_array[, , k]
		
		chain_means <- rowMeans(chain_matrix)  # Chain means (m x p)
		B <- n2 * var(chain_means)             # Between-chain covariance B
		W <- mean(apply(chain_matrix, 1, var)) # Within-chain covariance W
		var_hat <- ((n2 - 1)/n2)*W + B/n2      # Function 11.3
		
		rhat[k] <- sqrt(var_hat / W)           # Function 11.4
	}
	
	return(rhat)
}

#' Autocorrelation-Based ESS (Stan-style)
#'
#' @return
#' @export
#'
#' @examples
# compute_ess <- function(x) {
# 	x <- as.numeric(x)
# 	N <- length(x)
# 	x <- x - mean(x)
# 	
# 	acov <- acf(x, type = "covariance", plot = FALSE)$acf
# 	acov <- as.numeric(acov)
# 	
# 	tau <- 0
# 	for (k in seq(1, length(acov)-1, by=2)) {
# 		if (k+1 > length(acov)) break
# 		pair_sum <- acov[k+1] + acov[k+2]
# 		if (pair_sum < 0) break
# 		tau <- tau + pair_sum
# 	}
# 	
# 	var_hat <- acov[1] + 2*tau
# 	ess <- N * acov[1] / var_hat
# 	
# 	return(ess)
# }
compute_ess <- function(x) {
	
	x <- as.numeric(x)
	N <- length(x)
	
	# Remove NA if any
	x <- x[!is.na(x)]
	N <- length(x)
	
	if (N < 3) return(NA_real_)
	
	# If chain is constant → ESS = 1
	if (var(x) == 0) return(1)
	
	x_centered <- x - mean(x)
	
	acov <- acf(x_centered,
							type = "covariance",
							plot = FALSE,
							lag.max = N - 1)$acf
	
	acov <- as.numeric(acov)
	
	# If variance estimate invalid
	if (is.na(acov[1]) || acov[1] <= 0)
		return(1)
	
	tau <- 0
	
	for (k in seq(1, length(acov) - 1, by = 2)) {
		
		if ((k + 1) > length(acov)) break
		
		pair_sum <- acov[k + 1] + acov[k + 2]
		
		if (is.na(pair_sum) || pair_sum < 0)
			break
		
		tau <- tau + pair_sum
	}
	
	var_hat <- acov[1] + 2 * tau
	
	if (var_hat <= 0 || is.na(var_hat))
		return(1)
	
	ess <- N * acov[1] / var_hat
	
	# ESS cannot exceed total draws
	ess <- min(ess, N)
	
	return(ess)
}

#' ESS for all parameters across chains:
#'
#' @return
#' @export
#'
#' @examples
compute_ess_all <- function(chains_array) {
	m <- dim(chains_array)[1]
	n <- dim(chains_array)[2]
	p <- dim(chains_array)[3]
	
	ess_vals <- numeric(p)
	
	for (k in 1:p) {
		combined_chain <- as.vector(chains_array[, , k])
		ess_vals[k] <- compute_ess(combined_chain)
	}
	
	return(ess_vals)
}

#' Acceptance Rate
#'
#' @return
#' @export
#'
#' @examples
compute_acceptance_rate <- function(chain_ll) {
	mean(diff(chain_ll) != 0)
}

#' Acceptance Rate - multiple chains
#'
#' @return
#' @export
#'
#' @examples
compute_acceptance_all <- function(multichain_ll_list) {
	sapply(multichain_ll_list, function(chain_ll) {
		mean(diff(chain_ll) != 0)
	})
}

#' Run MCMC Diagnostics After Burn-in
#'
#' @param burned_samples List of chains (each chain = list of iteration draws)
#' @param burned_lls List of per-chain log-likelihood vectors
#' @return diagnostics list
#' @export
run_diagnostics <- function(burned_samples,
														burned_lls) {
	
	stopifnot(length(burned_samples) == length(burned_lls))
#browser()	
	m <- length(burned_samples)                # number of chains
	n <- length(burned_lls[[1]])               # iterations per chain - post burn
	
	# ensure equal lengths
	stopifnot(all(sapply(burned_lls, length) == n))
	stopifnot(all(sapply(burned_samples, length) == n))
	
	# Number of parameters ( accounts for stratification)
	p <- ncol(burned_samples[[1]][[1]])
	
	# Create array: iterations x chains x parameters
	chains_array <- array(NA_real_, dim = c(n, m, p))
	
	for (j in seq_len(m)) {
		chain_matrix <- do.call(rbind, burned_samples[[j]])
		chains_array[, j, ] <- chain_matrix
	}
	
	# Diagnostics
	rhat_vals_classic      <- compute_split_rhat(chains_array)
#	ess_vals_classic <- compute_ess_all(chains_array)
	
	# Modern diagnostics in Posterior library
	# Convert to posterior draws object 
	library(posterior)
	draws_obj <- as_draws_array(chains_array)
	
	rhat_vals    <- rhat(draws_obj)        # modern split Rhat
	ess_vals     <- ess_bulk(draws_obj)    # bulk ESS
	ess_tail_vals <- ess_tail(draws_obj)   # optional but recommended
	
	accept_rates <- compute_acceptance_all(burned_lls)
	
	cat("\n==== MCMC Diagnostics ====\n")
	cat("Max R-hat classic:", max(rhat_vals_classic, na.rm = TRUE), "\n")
	cat("Min R-hat:", min(rhat_vals, na.rm = TRUE), "\n")
	cat("Max R-hat:", max(rhat_vals, na.rm = TRUE), "\n")
	cat("Min ESS:", min(ess_vals, na.rm = TRUE), "\n")
	cat("Mean Acceptance Rate:", mean(accept_rates, na.rm = TRUE), "\n")
	cat("==========================\n\n")
	
	list(
		rhat = rhat_vals,
		ess_bulk = ess_vals,
		ess_tail = ess_tail_vals,
		acceptance_rate = accept_rates,
		summary = list(
			min_rhat = min(rhat_vals, na.rm = TRUE),
			max_rhat = max(rhat_vals, na.rm = TRUE),
			min_bulk_ess = min(ess_vals, na.rm = TRUE),
			min_tail_ess = min(ess_tail_vals, na.rm = TRUE),
			mean_acceptance = mean(accept_rates, na.rm = TRUE)
		)
	)
}