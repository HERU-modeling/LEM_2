#####################################################################################################################################
## Cascade CEA Model - get number pf HIV infections and susceptibles for different time periods
## Latest update: Oct 2022
####################################################################################################################################
## Output: a size 33 list
## 1-6: "Infections.total-20Y", "SuscPY-over20Y", "Infections.total-10Y", 
##      "SuscPY-over10Y", "Infections.total-5Y", "SuscPY-over5Y"
## 7-32: infections from 2015 -2040
###################################################################################################################################

#' Get number of HIV infections and susceptible for different time periods
#'
#' @param simul.in Simulation data containing information about infections.
#' @param start.period Starting period for evaluating infections.
#' @param end.period Ending period for evaluating infections.
#' @param sust.period Sustainment period.
#'
#' @return a list (out.inf.ls) containing matrices for the yearly number of infections for different racial or ethnic groups
#' @export
#'
#' @examples
time.period.infections <- function(simul.in,
                                   start.period,
                                   end.period,
                                   sust.period){
	with(as.list(vparameters),{
  #simul.in: Include value at t=0 at 1, t=12 at 13,..., t=348 at 349 (months time step saved in vt)
	
	# Target 1: Total Diagnosed PLHIV
		diag         = apply(simul.in[end_yr_ind, 7:16, ], c(1,3), sum)
		diag.all     = rowSums(diag)
		# 3 strata
		diag.w       = rowSums(diag[ , gp.gn$white])
		diag.b       = rowSums(diag[ , gp.gn$black])
		diag.h       = rowSums(diag[ , gp.gn$hisp])
		# 18 Strata
		diag18.model = matrix(0, length(end_yr_ind), 18)
		for (i in 1:18){
			ind = which(rname %in% names18[i])
			diag18.model[ ,i] = rowSums(diag[ ,ind, drop = FALSE]) #must include drop= FALSE
		}
		
	# Target 2: New diagnosed among infected
		ndiag        = rbind(simul.in[12, "diag", ],  apply(simul.in[end_yr_ind, "diag", ], 2, diff))
		ndiag.all    = rowSums(ndiag) # total
		# 3 strata
		ndiag.w       = rowSums(ndiag[ , gp.gn$white])
		ndiag.b       = rowSums(ndiag[ , gp.gn$black])
		ndiag.h       = rowSums(ndiag[ , gp.gn$hisp])
		# 18 Strata
		ndiag18.model = matrix(0, length(end_yr_ind), 18)
		for (i in 1:18){
			ind = which(rname %in% names18[i])
			ndiag18.model[ ,i] = rowSums(ndiag[ ,ind, drop = FALSE]) #must include drop= FALSE
		}
		
	# Target 3: (Number) Yearly all-cause mortality PLHIV	- deathm
		death        = rbind(simul.in[12,"death", ],  apply(simul.in[end_yr_ind, "death", ], 2, diff)) 
		death.all    = rowSums(death) # total: #Should be the same as sum(rowSums(mo*y[,1:16]))
		# 3 strata
		death.w       = rowSums(death[ , gp.gn$white])
		death.b       = rowSums(death[ , gp.gn$black])
		death.h       = rowSums(death[ , gp.gn$hisp])
		# 18 Strata
		death18.model = matrix(0, length(end_yr_ind), 18)
		for (i in 1:18){
			ind = which(rname %in% names18[i])
			death18.model[ ,i] = rowSums(death[ ,ind, drop = FALSE]) #must include drop= FALSE
		}

		# Target 4: (Proportion) Yearly late diagnoses	- late_diag	
		mon.late_diag        <- rbind(simul.in[1, "diag_late", ],  apply(simul.in[vt, "diag_late", ], 2, diff))
		late_diag.all        <- mapply(function(start, end) mean(mon.late_diag[start:end]), end_yr_ind-11, end_yr_ind)
		# 3 strata
		late_diag.w          <- mapply(function(start, end) mean(mon.late_diag[start:end,gp.gn$white]), end_yr_ind-11, end_yr_ind)
		late_diag.b          <- mapply(function(start, end) mean(mon.late_diag[start:end,gp.gn$black]), end_yr_ind-11, end_yr_ind)
		late_diag.h          <- mapply(function(start, end) mean(mon.late_diag[start:end,gp.gn$hisp]), end_yr_ind-11, end_yr_ind)
		# 18 Strata
		late_diag18.model = matrix(0, length(end_yr_ind), 18)
		for (i in 1:18){
			ind = which(rname %in% names18[i])
			late_diag18.model[ ,i] = mapply(function(start, end) mean(mon.late_diag[start:end,ind]), end_yr_ind-11, end_yr_ind)
		}
	
		#revisit		
		# Target 5: (Proportion) Yearly Viral suppression 2017 to 2021	- viral_sup
		mon.viral_sup        <- rbind(simul.in[1, "viral_sup", ],  apply(simul.in[vt, "viral_sup", ], 2, diff))
		viral_sup.all        <- mapply(function(start, end) mean(mon.viral_sup[start:end]), end_yr_ind-11, end_yr_ind)
		# 3 strata
		viral_sup.w          <- mapply(function(start, end) mean(mon.viral_sup[start:end,gp.gn$white]), end_yr_ind-11, end_yr_ind)
		viral_sup.b          <- mapply(function(start, end) mean(mon.viral_sup[start:end,gp.gn$black]), end_yr_ind-11, end_yr_ind)
		viral_sup.h          <- mapply(function(start, end) mean(mon.viral_sup[start:end,gp.gn$hisp]), end_yr_ind-11, end_yr_ind)
		# 18 Strata
		viral_sup18.model = matrix(0, length(end_yr_ind), 18)
		for (i in 1:18){
			ind = which(rname %in% names18[i])
			viral_sup18.model[ ,i] = mapply(function(start, end) mean(mon.viral_sup[start:end,ind]), end_yr_ind-11, end_yr_ind)
		}
		# Reformat for output
	
		out.target    <- cbind(diag.all, ndiag.all, death.all, late_diag.all, viral_sup.all)
			
		
#		simul.in <- simul.in[-1,,] # Remove initial t=0 value
  # Get number of infections (incidence): 
  # (1) over 20-year evaluation period: 2025-2045
  t.int.infs.eval <- sum(simul.in[end.period, c("inc_bo", "inc_bs", "inc_g"), ]) - 
  	sum(simul.in[start.period, c("inc_bo", "inc_bs", "inc_g"), ])
  
  out.inf.eval    <- as.matrix(t.int.infs.eval)
  colnames(out.inf.eval) <- paste0("Infections.total-20Y")
  
  # (2) over 10-year implementation period: 2025-2035 (Dec, 2035)
  t.int.infs.impl <- sum(simul.in[(start.period + (sust.period * 12)), c("inc_bo", "inc_bs", "inc_g"), ]) - 
    sum(simul.in[start.period, c("inc_bo", "inc_bs", "inc_g"), ])
  
  out.inf.impl    <- as.matrix(t.int.infs.impl)
  colnames(out.inf.impl) <- paste0("Infections.total-10Y")
  
  # (3) over 5-year period: 2025-2030 (Dec, 2030)
  t.int.infs.5y   <- sum(simul.in[(start.period + (5 * 12)), c("inc_bo", "inc_bs", "inc_g"), ]) - 
                              sum(simul.in[start.period, c("inc_bo", "inc_bs", "inc_g"), ])
  
  out.inf.5y      <- as.matrix(t.int.infs.5y)
  colnames(out.inf.5y)   <- paste0("Infections.total-5Y")

  # Get person-years among susceptibles over same time periods ?
  # (1) over entire evaluation period
  denom.sum.eval  <- as.matrix(sum(simul.in[start.period:end.period, 1:2, names.gp]) / (12 * 100000))
  colnames(denom.sum.eval) <- paste0("SuscPY-over20Y")

  # (2) over 10-year implementation period
  denom.sum.impl  <- as.matrix(sum(simul.in[start.period:(start.period + sust.period * 12), 1:2, names.gp])/ (12 * 100000))
  colnames(denom.sum.impl) <- paste0("SuscPY-over10Y")

  # (3) over 5-year period
  denom.sum.5y  <- as.matrix(sum(simul.in[start.period:(start.period + 5 * 12), 1:2, names.gp])/ (12 * 100000))
  colnames(denom.sum.5y)   <- paste0("SuscPY-over5Y")
  
  # Get yearly number of infections (2020 - 2045)
  yr_end               <- seq(36, end.period + 12, by=12) # Vector 27:36=1/2015;48=1/2016;60=1/2017;...;336=1/2040;348=1/2041 (Note: 0=1/2012;12=1/2013;24=1/2014)
  t.int.infs.eval      <- apply((simul.in[yr_end, c("inc_bo", "inc_bs", "inc_g"), ]), c(1), sum) # Still Vector 37
  out.inf.yr           <- t(as.matrix(diff(t.int.infs.eval))) # Vector 26, because take the difference: Yr_2015 = Value 1/2016 (t=48) - Value 1/2015 (t=36); Yr_2040 = Value 1/2041 (t=348) - Value 1/2015 (t=336)
  colnames(out.inf.yr) <- paste0(c(2020:lyr))
  
  out.inf    <- cbind(out.inf.eval, denom.sum.eval, out.inf.impl, denom.sum.impl, out.inf.5y, denom.sum.5y)
  out.inf.ls <- list (out.inf = out.inf, out.inf.yr = out.inf.yr)
browser()
  return(out.inf.ls)
	})
}