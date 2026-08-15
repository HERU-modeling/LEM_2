require(deSolve)

obj <- function(calib.par, calpar.info, fixed, fixed.list, calib.target, valid.target, num_targets = 144) { 
	# Attain inputs for objective function
	names.gp   = as.vector(fixed$names.gp)
	names18    = as.vector(fixed$names18)
	n.gp       = length(fixed$names.gp)
	state.name = fixed$state.name
	leng   = calpar.info$plength
	names  = calpar.info$names

  # Load target
	diag.obs  = calib.target$diag18.obs [-1, ] # Remove 2016 data
	ndiag.obs = calib.target$ndiag18.obs[-1, ]
	death.obs = calib.target$death18.obs[-1, ]
	inc.all   = valid.target$obs.inc.all
	inc.msm   = valid.target$obs.inc.msm
	
################### Replace original value with new sampled value #####	
	vlist = fixed.list
	vparameters = fixed
		
	par2   = as.list(leng)
#browser()	#Check calib.par = new_prior_samples
	for(l in 1:length(leng)){
		if (l==1) {
			par2[[l]]   = calib.par[1:leng[l]]
		}
		else {
			par2[[l]]   = calib.par[(sum(leng[1:(l-1)])+1):sum(leng[1:l])]
		}
	}
	
	names(par2)   = names
	
	# update parameter list with new free parameter values and save in par3
	for (l in seq_along(names)) {
		if (par_info[par_info$parameter == names[l], ]$stratification == 0) {
			vlist[[names[l]]]$pe <- par2[[names[l]]]
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 1) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$cd4       == calpar[calpar$par == names[l], ][j, ]$cd4
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 1.1) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$gender    == calpar[calpar$par == names[l], ][j, ]$gender
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 1.2) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$ethnicity == calpar[calpar$par == names[l], ][j, ]$ethnicity
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 1.3) { 
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$risk      == calpar[calpar$par == names[l], ][j, ]$risk
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 3) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$gender  == calpar[calpar$par == names[l], ][j, ]$gender &
					vlist[[names[l]]]$ethnicity          == calpar[calpar$par == names[l], ][j, ]$ethnicity &
					vlist[[names[l]]]$risk               == calpar[calpar$par == names[l], ][j, ]$risk
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 3.5) {
			for (j in 1:leng[l]) {
				if (calpar[calpar$par == names[l], ][j, ]$component != "reach") {
					vlist[[names[l]]][vlist[[names[l]]]$component == calpar[calpar$par == names[l], ][j, ]$component, ]$pe <- par2[names[l]][[1]][j]
				} else {
					vlist[[names[l]]][vlist[[names[l]]]$gender    == calpar[calpar$par == names[l], ][j, ]$gender & 
														vlist[[names[l]]]$ethnicity == calpar[calpar$par == names[l], ][j, ]$ethnicity & 
														vlist[[names[l]]]$component == calpar[calpar$par == names[l], ][j, ]$component, ]$pe <- par2[names[l]][[1]][j]
				}
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 4) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$gender  == calpar[calpar$par == names[l], ][j, ]$gender &
					vlist[[names[l]]]$ethnicity          == calpar[calpar$par == names[l], ][j, ]$ethnicity &
					vlist[[names[l]]]$risk               == calpar[calpar$par == names[l], ][j, ]$risk &
					vlist[[names[l]]]$sex.risk.age.group == calpar[calpar$par == names[l], ][j, ]$sex.risk.age.group
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 4.5) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$gender  == calpar[calpar$par == names[l], ][j, ]$gender &
					vlist[[names[l]]]$ethnicity          == calpar[calpar$par == names[l], ][j, ]$ethnicity &
					vlist[[names[l]]]$risk               == calpar[calpar$par == names[l], ][j, ]$risk &
					vlist[[names[l]]]$cd4                == calpar[calpar$par == names[l], ][j, ]$cd4
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
		if (par_info[par_info$parameter == names[l], ]$stratification == 5) {
			for (j in 1:leng[l]) {
				match_idx <- vlist[[names[l]]]$gender  == calpar[calpar$par == names[l], ][j, ]$gender &
					vlist[[names[l]]]$ethnicity          == calpar[calpar$par == names[l], ][j, ]$ethnicity &
					vlist[[names[l]]]$risk               == calpar[calpar$par == names[l], ][j, ]$risk &
					vlist[[names[l]]]$sex.risk.age.group == calpar[calpar$par == names[l], ][j, ]$sex.risk.age.group &
					vlist[[names[l]]]$year               == calpar[calpar$par == names[l], ][j, ]$year
				vlist[[names[l]]]$pe[match_idx] <- par2[[names[l]]][j]
			}
		}
	}

	# Replace Model Parameters with new values that was sampled
	source("ModelCoreModules/Parameter_input/Update_vlist_to_vparameters.R")
	vparameters <- par_fix(vparameters, vlist)
	result      <- par_pop(vparameters, vlist)  #Function part of Update_vlist_to_vparameters
	vparameters <- result$vparameters
	x           <- result$x
	vparameters <- add_intervention_to_param(vparameters, "No interventions")
	par3 = vparameters
	
  # Run the model with new sampled values
	out_euler <- euler(x, vt, ode_model, par3)[ ,-1]
	outa = array(out_euler[-1, ], dim = c(n, length(state.name), n.gp)) 
	dimnames(outa)[[2]] = state.name
	dimnames(outa)[[3]] = names.gp
	outn = outa[ ,1:19, ] 

################### Extract Model Outputs & Observed data #####	
	# Create Equation to Modify outputs to different strata
	# group11 = function(group18){
	#   g11 = matrix(0,length(end_yr_ind_calibration),11)
	#   g11[ ,4] = rowSums(group18[ ,4:6])
	#   g11[ ,5] = rowSums(group18[ ,7:12])
	#   g11[ ,c(1:3, 6:11)] = group18[ ,c(1:3, 13:18)]
	#   return(g11)
	# }
	group18 = function(group30){
	  g18 = matrix(0, length(end_yr_ind_calibration), 18)
	  for (i in 1:18){
	    ind = which(rname %in% names18[i])
	    g18[ ,i] = rowSums(group30[ ,ind, drop = FALSE]) #must include drop= FALSE
	    colnames(g18) <- names18 
	  }
	  return(g18)
	}
	group.race = function(group30){
	  g.race = matrix(0,length(end_yr_ind_calibration),3)
	  g.race[ ,1] = rowSums(group30[ ,gp.gn$white])  #White
	  g.race[ ,2] = rowSums(group30[ ,gp.gn$black])  #Black
	  g.race[ ,3] = rowSums(group30[ ,gp.gn$hisp])   #Hispanic
	  colnames(g.race) <- c("White","Black", "Hisp")
	  return(g.race)
	}
	group.calibration = function(group30){
	  g.cal = matrix(0,length(end_yr_ind_calibration),5)
	  g.cal[ ,1] = rowSums(group30[ ,gp.gn$het.m])  #M/HET
	  g.cal[ ,2] = rowSums(group30[ ,gp.gn$het.f])  #F/HET
	  g.cal[ ,3] = rowSums(group30[ ,gp.gn$young])   #15-34/MSM,now include MWID
	  g.cal[ ,4] = rowSums(group30[ ,gp.gn$old])  #35-34/MSM,now include MWID
	  g.cal[ ,5] = rowSums(group30[ ,gp.gn$pwid])   #PWID only
	  colnames(g.cal) <- c("het.m","het.f","young","old","pwid")
	  return(g.cal)
	}
	
	### Output all model for Diagnostic plots
	# Target 1: Total Diagnosed PLHIV
	diag         = apply(outa[end_yr_ind_calibration, 7:16, ], c(1,3), sum) # Prevalence (Ok to be cumulative) 12:Dec2017, 24:Dec2018, 36:Dec2019, 48:Dec2020, 60:Dec2021, 48:Dec2022
	colnames(diag)  = names.gp
	diag.all     = rowSums(diag)
	# 18 Strata
	diag18.model    <- group18(diag)
	#Calibration strata
	diag.race.model <- round(group.race(diag)) 
	diag.cal.model  <- round(group.calibration(diag)) 
	# Obs- Calibration target
	diag.obs     = calib.target$diag.obs 
	
	# Target 2: New diagnosed among infected
	ndiag        = rbind(outa[12, "diag", ],  apply(outa[end_yr_ind_calibration, "diag", ], 2, diff)) # New Diag: cumulative for up to t=12 for 2017 ndiag value, after, difference between t=24-12 give 2018 value,..., t72-60=2022 value 
	colnames(ndiag)  = names.gp
	ndiag.all    = rowSums(ndiag) 
	# 18 Strata
	ndiag18.model    <- group18(ndiag)
	#Calibration strata
	ndiag.race.model <- round(group.race(ndiag)) 
	ndiag.cal.model  <- round(group.calibration(ndiag)) 
	# Obs- Calibration target
	ndiag.obs     = calib.target$ndiag.obs 
	
	# Target 4: (Proportion) Yearly late diagnoses: diag_late <- psi[6] * I3 / diagm: psi[3] * Ia + psi[4] * I1 + psi[5] * I2 + psi[6] * I3 (Target)
	ldiag  <- rbind(outa[12, "diag_late", ],  apply(outa[end_yr_ind_calibration, "diag_late", ], 2, diff)) # Exactly same as New diag (only that diag late is for I3 only), then divided by ndiag
	p.diag_late.all = ldiag/ndiag

	p.diag_late.w   = rowSums(ldiag[, gp.gn$white])/rowSums(ndiag[, gp.gn$white])
	p.diag_late.b   = rowSums(ldiag[, gp.gn$black])/rowSums(ndiag[, gp.gn$black])
	p.diag_late.h   = rowSums(ldiag[, gp.gn$hisp ])/rowSums(ndiag[, gp.gn$hisp ])
	p.diag_late.cal.model   = cbind(p.diag_late.w,p.diag_late.b,p.diag_late.h)
	# Obs- Calibration target
	p.diag_late.cal.obs     = calib.target$late_diag.obs  
	
	# Target 5: (Proportion) Yearly Viral suppression 2017 to 2021
	p.viral_sup.all = apply(outa[end_yr_ind_calibration, c(11,12), ], c(1,3), sum)*mean(vparameters$prop_vs_T) / apply(outa[end_yr_ind_calibration, c(7:16), ], c(1,3), sum) #T1, T2
	p.viral_sup.w   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$white], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$white][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$white], c(1,3), sum))
	p.viral_sup.b   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$black], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$black][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$black], c(1,3), sum))
	p.viral_sup.h   = rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$hisp ], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$hisp ][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$hisp ], c(1,3), sum))
	p.viral_sup.cal.model   = cbind(p.viral_sup.w,p.viral_sup.b,p.viral_sup.h)
	# Obs- Calibration target
	p.viral_sup.cal.obs     = calib.target$viral_sup.obs 
	
	
  n.model= round(c(diag.race.model, ndiag.cal.model, ndiag.race.model))
  p.model= c(p.diag_late.cal.model, p.viral_sup.cal.model)
  n.obs  = round(c(diag.obs  , ndiag.obs))
  p.obs  = c(p.diag_late.cal.obs  , p.viral_sup.cal.obs)
  model = c(n.model, p.model)
  obs   = c(n.obs  , p.obs)


################### Obtained likelihood of the model outputs from the new sampled input when compared to observed targets #####	
	#Convert p.obs to pseudo-count
  model = c(n.model, p.model*100)
  obs   = c(n.obs  , p.obs*100)
  mean_norm <- obs
  variance_norm <- obs*10

  ## 1. Define Year wt
  # Number of years
  yr_calib <- 6
  year_weight <- numeric(yr_calib)
  year_weight[1:2] <- (1/3) / 2   # years 1–2
  year_weight[3:4] <- 0           # years 3–4
  year_weight[5:6] <- (2/3) / 2   # years 5–6
  
  ## 2. Organize your model / obs by year
  K <- length(model)/yr_calib   # targets per year
  model_mat <- matrix(model        , nrow = yr_calib, ncol = K, byrow = FALSE)
  obs_mat   <- matrix(obs          , nrow = yr_calib, ncol = K, byrow = FALSE)
  var_mat   <- matrix(variance_norm, nrow = yr_calib, ncol = K, byrow = FALSE)

  make_sigma_blockwise <- function(var_vec,
  																	rho_target = rho_target,
  																	rho_race  = rho_race) {
  	
  	K <- length(var_vec)
  	stopifnot(K == 17)
  	
  	sd_vec <- sqrt(var_vec)
  	
  	sigma_block <- matrix(0, K, K)
  	
  	# ---- Define blocks ----
  	idx_target <- list(
  		diag_race  = 1:3,
  		ndiag_cal  = 4:8,
  		ndiag_race = 9:11,
  		p_late     = 12:14,
  		p_viral    = 15:17
  	)
  	
  	# ---- Define race groups ----
  	idx_race <- list(
  		W = c(1,9,12,15),
  		B = c(2,10,13,16),
  		H = c(3,11,14,17)
  	)
	
  	# ---- Add target covariance ----
  	for (target in idx_target) {
# browser() 		
  		for (i in target) {
  			for (j in target) {
  				if (i != j) {
  					sigma_block[i,j] <- sigma_block[i,j] + rho_target * sd_vec[i] * sd_vec[j]
  				}
  			}
  		}
  	}
# browser() 	
  	# ---- Add race covariance ----
  	for (race in idx_race) {
# browser() 		
  		for (i in race) {
  			for (j in race) {
  				if (i != j) {
  					sigma_block[i,j] <- sigma_block[i,j] + rho_race * sd_vec[i] * sd_vec[j]
  				}
  			}
  		}
  	}
# browser() 	
  	# ---- Set diagonal ----
  	diag(sigma_block) <- var_vec
  	
  	# ---- Ensure symmetry ----
  	sigma_block<- (sigma_block + t(sigma_block)) / 2 # This should be unneeded
# browser()   	
  	return(sigma_block)
  }
  
  
  ## 5. Per-year MVN likelihood with block-wise
  loglik_year <- function(model_row,
  												obs_row,
  												var_row,
  												rho_target = 0,
  												rho_race   = 0) {
  	sigma_block <- make_sigma_blockwise(var_row,
  																rho_target = rho_target,
  																rho_race   = rho_race)
 	
  	loglik_year <- mvtnorm::dmvnorm(model_row,
  									 mean  = obs_row,
  									 sigma = sigma_block,
  									 log   = TRUE)
#browser()  	#Check the same as loglik_year = correspond sum of values in dnorm(model, mean = mean_norm, sd = sd_norm, log = TRUE)  when all rho = 0
  	return(loglik_year)
  }
  # A note regarding target weights, used for weighted average: Cannot simply weighted average likelihood
  # The correct way is to Inflate the variance of targets you want to downweight
  # If goal is to calibrate only to 3 targets (O_T to p.viral_sup.W/B/H) and completely ignore the other 14
  # The only need to run dmvnorm w subset the likelihood to only those 3 targets.
  
  ## 6. Full weighted likelihood
  loglik_vec <- numeric(yr_calib)
  
  for (t in 1:yr_calib) {
  	loglik_vec[t] <- loglik_year(model_mat[t, ],
  															 obs_mat[t, ],
  															 var_mat[t, ],
  															 rho_target = 0.6,
  															 rho_race   = 0.1) # rho_target, rho_race: from correlation.R
  	
  }
#browser()  # Check sum(loglik_vec) = sum(dnorm(model, mean = mean_norm, sd = sd_norm, log = TRUE)) when rho =0
  # Weighted log-likelihood
  loglik_total <- sum(year_weight * loglik_vec) 
  exp_ll <- exp(loglik_total)
  if (any(is.nan(exp_ll) | is.infinite(exp_ll) | exp_ll == 0)) {   
  	cat("exp(loglik_total) is NaN, Inf, or numerically underflowed to 0 — entering debug mode\n")  
  	browser()
 } 
return(list(	model  = model,	loglik = loglik_total,	obs    = obs, outa=outa))# 

}

