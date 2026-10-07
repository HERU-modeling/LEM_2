rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

source("Calibration/6_Others/Utility_functions.R") 	
#source("Model_development/7.Red_strata/7.Red_strata-LoadBaselineWorkspace-Stra.R")
source("ModelCoreModules/Parameter_input/LoadBaselineWorkspace.R")

CITY <- "MAR" # MIA1, NYC1 for original, KIN, PHI, BAL

# source("Model_development/7.Red_strata/Calibration/7.load_targets.R")
# source("ModelCoreModules/Parameter_input/_targets.R") #Update Target as needed - Part of par_pop.R 
source("ModelCoreModules/Parameter_input/Module-Data.input.R")
source("Calibration/4_Visualization/plots.R")

interv <- c(17) # Pick 0 for no intervention for NYC and MIA
total.int <- length(interv) 

library(ggpubr)

plot_line <- function(data, name, TITLE = NA, CITY = CITY, start_year = 2017, y_min = NULL) {
	df <- as.data.frame(data); colnames(df) <- name
	df$Year <- seq(from = start_year, by = 1, length.out = nrow(df))
	
	df_long <- df |>
		tidyr::pivot_longer(cols = name, names_to = "Legend", values_to = "Value")
	
	label_data <- dplyr::filter(df_long, Year == min(Year) | Year == max(Year))
	
	g <- ggplot(df_long, aes(x = Year, y = Value, color = Legend)) +
		geom_line(size = 1.2) +
		geom_text(data = label_data, aes(label = round(Value, 1)),
							vjust = -0.5, size = 3.5, show.legend = FALSE) +
		theme_minimal(base_size = 14) +
		labs(title = paste(TITLE, CITY), x = "Year", y = TITLE, color = "Legend") +
		scale_color_brewer(palette = "Set1") +
		theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
	
	if (!is.null(y_min)) g <- g + expand_limits(y = y_min)
	g  # return, don't print
}

#Remove these - streamline them to be same as the calibration obj.R, as now we no longer needed the first year (which were conceive to show the prevalence match at t=0, which cannot be index normally)
#end_yr_ind   = c(c(12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,300,312,324,336,348)+1)
#end_yr_ind_plus   = c(1,end_yr_ind)

##########################################################
# Write function
extract_outputs <- function(outa) {

	diag  <- apply(outa[end_yr_ind_calibration, 7:16, ], c(1,3), sum) # Prevalence (Ok to be cumulative) 12:Dec2017, 24:Dec2018, 36:Dec2019, 48:Dec2020, 60:Dec2021, 48:Dec2022
	# Total inc
	base_w <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$white],dim = c(1, 3, 10),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$white) )
	diff_w <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$white],		c(2,3),		diff	)
	inc_w <- abind(base_w, diff_w, along = 1)
	base_b <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$black],dim = c(1, 3, 10),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$black) )
	diff_b <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$black],		c(2,3),		diff	)
	inc_b <- abind(base_b, diff_b, along = 1)
	base_h <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$hisp ],dim = c(1, 3, 10),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$hisp ) )
	diff_h <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$hisp ],		c(2,3),		diff	)
	inc_h <- abind(base_h, diff_h, along = 1)
	
	base_young <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$young],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$young) ) # 1: t=12, 3: inc_bs/bo/g, 6: young strata
	diff_young <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$young],		c(2,3),		diff	)
	inc_young  <- abind(base_young, diff_young, along = 1)
	base_old <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$old],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$old) )
	diff_old <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$old],		c(2,3),		diff	)
	inc_old  <- abind(base_old, diff_old, along = 1)
	base_pwid <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$pwid],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$pwid ) )
	diff_pwid <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$pwid ],		c(2,3),		diff	)
	inc_pwid <- abind(base_pwid, diff_pwid, along = 1)
	base_het.m <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$het.m],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$het.m) ) # 1: t=12, 3: inc_bs/bo/g, 6: young strata
	diff_het.m <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$het.m],		c(2,3),		diff	)
	inc_het.m  <- abind(base_het.m, diff_het.m, along = 1)
	base_het.f <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$het.f],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$het.f) )
	diff_het.f <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$het.f],		c(2,3),		diff	)
	inc_het.f  <- abind(base_het.f, diff_het.f, along = 1)
	base_low <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$low],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$low) ) # 1: t=12, 3: inc_bs/bo/g, 6: young strata
	diff_low <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$low],		c(2,3),		diff	)
	inc_low  <- abind(base_low , diff_low , along = 1)
	base_high <- array(outa[12, c("inc_bo","inc_bs","inc_g"), gp.gn$high],dim = c(1, 3, 6),dimnames = list("12",	c("inc_bo","inc_bs","inc_g"),	gp.gn$high) )
	diff_high <- apply(		outa[end_yr_ind, c("inc_bo","inc_bs","inc_g"), gp.gn$high],		c(2,3),		diff	)
	inc_high  <- abind(base_high, diff_high, along = 1)
# browser()
	list(
		diag_race  = cbind(rowSums(diag[, gp.gn$white]),
											 rowSums(diag[, gp.gn$black]),
											 rowSums(diag[, gp.gn$hisp])),
	ndiag  <- rbind(outa[12, "diag", ],  apply(outa[end_yr_ind_calibration, "diag", ], 2, diff)), # New Diag: cumulative for up to t=12 for 2017 ndiag value, after, difference between t=24-12 give 2018 value,..., t72-60=2022 value 
		ndiag_race = cbind(
#			rowSums(apply(outa[end_yr_ind_plus,"diag",],2, diff)[, gp.gn$white]),
#			rowSums(apply(outa[end_yr_ind_plus,"diag",],2, diff)[, gp.gn$black]),
#			rowSums(apply(outa[end_yr_ind_plus,"diag",],2, diff)[, gp.gn$hisp])
			rowSums(ndiag[, gp.gn$white]),
			rowSums(ndiag[, gp.gn$black]),
			rowSums(ndiag[, gp.gn$hisp])
			),
		ndiag_risk = cbind(
			rowSums(ndiag[, gp.gn$young]), #Before gp.gn$Y.msm
			rowSums(ndiag[, gp.gn$old]), #Before O.msm
			rowSums(ndiag[, gp.gn$pwid]), #Before all.pwid
			rowSums(ndiag[, gp.gn$het.m]),
			rowSums(ndiag[, gp.gn$het.f]),
			rowSums(ndiag[, gp.gn$het])
		),
	ldiag  <- rbind(outa[12, "diag_late", ],  apply(outa[end_yr_ind_calibration, "diag_late", ], 2, diff)), # Exactly same as New diag (only that diag late is for I3 only), then divided by ndiag
			diag_late = cbind(
			rowSums(ldiag[, gp.gn$white])/rowSums(ndiag[, gp.gn$white]),
			rowSums(ldiag[, gp.gn$black])/rowSums(ndiag[, gp.gn$black]),
			rowSums(ldiag[, gp.gn$hisp ])/rowSums(ndiag[, gp.gn$hisp ])
		),
		viral_sup = cbind(
			rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$white], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$white][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$white], c(1,3), sum)),
			rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$black], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$black][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$black], c(1,3), sum)),
			rowSums(apply(outa[end_yr_ind_calibration, c(11,12), gp.gn$hisp ], c(1,3), sum))*vparameters$prop_vs_T[gp18.gn$hisp ][1]/rowSums(apply(outa[end_yr_ind_calibration, c(7:16), gp.gn$hisp ], c(1,3), sum))
		),
		inc = cbind(
			rowSums(rbind(outa[12, c("inc_bo"), ],apply((outa[end_yr_ind, c("inc_bo"), ]), 2, diff))),
			rowSums(rbind(outa[12, c("inc_bs"), ],apply((outa[end_yr_ind, c("inc_bs"), ]), 2, diff))),
			rowSums(rbind(outa[12, c("inc_g"),  ],apply((outa[end_yr_ind, c("inc_g"),  ]), 2, diff)))
		),
		# ir_source = cbind(
		# 	rowSums(rbind(outa[12, c("inc_bo"), ],apply((outa[end_yr_ind, c("inc_bo"), ]), 2, diff)))/rowSums(apply(outa[end_yr_ind, c(1:16),               ], c(1,3), sum)),
		# 	rowSums(rbind(outa[12, c("inc_bs"), ],apply((outa[end_yr_ind, c("inc_bs"), ]), 2, diff)))/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$all.msm ], c(1,3), sum)),
		# 	rowSums(rbind(outa[12, c("inc_g"),  ],apply((outa[end_yr_ind, c("inc_g"),  ]), 2, diff)))/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$all.pwid], c(1,3), sum))
		# ),
      ir_risk = cbind(
      	apply(inc_young,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$young ], c(1,3), sum)), #Note the denomiantor
      	apply(inc_old  ,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$old   ], c(1,3), sum)),
      	apply(inc_pwid ,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$pwid  ], c(1,3), sum)),
      	apply(inc_het.m,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$het.m ], c(1,3), sum)), #Note the denomiantor
      	apply(inc_het.f,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$het.f ], c(1,3), sum)),
      	apply(inc_low  ,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$low  ], c(1,3), sum)), #Note the denomiantor
      	apply(inc_high ,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$high ], c(1,3), sum))
     ),
		ir_race = cbind(
			apply(inc_w,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$white ], c(1,3), sum)),
			apply(inc_b,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$black ], c(1,3), sum)),
			apply(inc_h,1,sum)/rowSums(apply(outa[end_yr_ind, c(1:16), gp.gn$hisp  ], c(1,3), sum))
		),
     # prop_T = rowSums(apply(outa[, 11:13, ], c(1,3), sum))/rowSums(apply(outa[, 7:16, ], c(1,3), sum)) #Prop T (11:13) over sum D+T+O (7:16) = rowSums(diag) at specific end_yr_ind_calib
		 prop_high.inf.among.het = 		rowSums(apply(outa[end_yr_ind, 3:16, gp.gn$high ], c(1,3), sum))/rowSums(apply(outa[end_yr_ind, 3:16, gp.gn$het ], c(1,3), sum)),
     prop_T_race = cbind(
     	rowSums(apply(outa[end_yr_ind, 11:13, gp.gn$white], c(1,3), sum))/rowSums(apply(outa[end_yr_ind, 7:16, gp.gn$white], c(1,3), sum)),
     	rowSums(apply(outa[end_yr_ind, 11:13, gp.gn$black], c(1,3), sum))/rowSums(apply(outa[end_yr_ind, 7:16, gp.gn$black], c(1,3), sum)),
     	rowSums(apply(outa[end_yr_ind, 11:13, gp.gn$hisp ], c(1,3), sum))/rowSums(apply(outa[end_yr_ind, 7:16, gp.gn$hisp ], c(1,3), sum))
		)
	)
	  # browser()	# Check apply(inc_w,1,sum)+apply(inc_b,1,sum)+apply(inc_h,1,sum) =  rowSums(rbind(outa[12, c("inc_bo"), ],apply((outa[end_yr_ind, c("inc_bo"), ]), 2, diff)))+rowSums(rbind(outa[12, c("inc_bs"), ],apply((outa[end_yr_ind, c("inc_bs"), ]), 2, diff)))+rowSums(rbind(outa[12, c("inc_g"),  ],apply((outa[end_yr_ind, c("inc_g"),  ]), 2, diff))) 
}

##########################################################
# Target
vparameters$pop.diag # First year 2016
calib.target$diag.obs # Remaining years 2016
new_row <- data.frame(
	W = sum(vparameters$pop.diag[gp18.gn$white]),
	B = sum(vparameters$pop.diag[gp18.gn$black]),
	H = sum(vparameters$pop.diag[gp18.gn$hisp])
)
diag.obs <- rbind(new_row, calib.target$diag.obs)

##########################################################
# Run baseline
cc <- interv[1]
out_base <- intervention.model.combination(input.parameters = vparameters,
																					 current.int = interventions[cc])
base <- extract_outputs(out_base)

##########################################################
# Set whether to use calibrated parameters from Excel
use.cali_par.xlsx <- "N"  
#browser()
if (use.cali_par.xlsx == "Y") {
	## READ in calibrated free parameters ##
	calpar <- read.xlsx("Model_development/7.Red_strata/7.Red_strata-Par_adj.xlsx", sheet = CITY) # Update the excel! Later on this is will be where we compile All the MANUAL adjustment
#	calpar <- read.xlsx("Model_development/7.Red_strata/7.Red_strata-Par_adj-Azzam.xlsx", sheet = CITY)
	names <- unique(calpar$par) # names of free.par
	leng <- as.numeric(table(factor(calpar$par, levels = names))) # length of each free.par
	
	# Assigning parameter values for selected free parameters
	par2 <- as.list(leng)
	for (l in seq_along(leng)) {
		if (l == 1) {
			par2[[l]] <- calpar$pe[1:leng[l]]
		} else {
			par2[[l]] <- calpar$pe[(sum(leng[1:(l - 1)]) + 1):sum(leng[1:l])]
		}
	}
	names(par2) <- names
	
	for (l in seq_along(names)) {
		#browser(expr= {l = 11})
		#		if (length(unlist(vparameters[names[l]])) == 1) { #Don't really want to update vparameters as this should be done through vlist_to_vparameters
		#			vparameters[[names[l]]] <- par2[[names[l]]]
		#		}
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
browser()  # check vlist, par2 for calibrated free parameter
	
} else {
	# Equivalent to the current use.cali_par.xlsx
 # browser()
	# vlist$mu_adj$pe[vlist$mu_adj$ethnicity == "w"] <- 0.0025
	# vlist$mu_adj$pe[vlist$mu_adj$ethnicity == "b"] <- 0.0025
	# vlist$mu_adj$pe[vlist$mu_adj$ethnicity == "h"] <- 0.001# -0.0003138994
	# vlist$mor_I2$pe[vlist$mor_I2$ethnicity == "b"] <- vlist$mor_I2$pe[vlist$mor_I2$ethnicity == "b"]*10
	# vlist$mor_I3$pe[vlist$mor_I3$ethnicity == "b"] <- vlist$mor_I3$pe[vlist$mor_I3$ethnicity == "b"]*10
	# vlist$mor_I2$pe[vlist$mor_I2$ethnicity == "h"] <- vlist$mor_I2$pe[vlist$mor_I2$ethnicity == "h"]*10
	# vlist$mor_I3$pe[vlist$mor_I3$ethnicity == "h"] <- vlist$mor_I3$pe[vlist$mor_I3$ethnicity == "h"]*10
	
	## Prop.hiv.aware - Differ By city
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "m" & vlist$prop.HIV.aware$risk == "pwid"] <- 0.95
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "m" & vlist$prop.HIV.aware$risk == "msm"] <- 0.95
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "m" & vlist$prop.HIV.aware$risk == "mwid"]<-0.95
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "m" & vlist$prop.HIV.aware$risk == "het"] <- 0.95
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "f" & vlist$prop.HIV.aware$risk == "pwid"]<-0.95
# vlist$prop.HIV.aware$pe[vlist$prop.HIV.aware$gender == "f" & vlist$prop.HIV.aware$risk == "het"] <-0.95

	## noG - Mostly the same across - except F/B/HEt/High
	# vlist$noG$pe[vlist$noG$risk == "msm"] <- 0.0643                                        #0.0643
	# vlist$noG$pe[vlist$noG$risk == "het" & vlist$noG$sex.risk.age.group == "low"] <- 0.0896 #0.0896
	# vlist$noG$pe[vlist$noG$ethnicity == "w" & vlist$noG$gender == "m"] <- 0.2              #0.670
	# vlist$noG$pe[vlist$noG$ethnicity == "b" & vlist$noG$gender == "m"] <- 0.7               #0.633
	# vlist$noG$pe[vlist$noG$ethnicity == "h" & vlist$noG$gender == "m"] <- 0.5               #0.579
	# vlist$noG$pe[vlist$noG$ethnicity == "w" & vlist$noG$gender == "f"] <- 0.2               #0.656
	# vlist$noG$pe[vlist$noG$ethnicity == "b" & vlist$noG$gender == "f"] <- 0.8               #0.545
	# vlist$noG$pe[vlist$noG$ethnicity == "h" & vlist$noG$gender == "f"] <- 0.5               #0.474
	
	# ## nsG - Same across all cities, same between Old & Y, W >H >B in term of partner
	# vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "old" & vlist$nsG$ethnicity == "w"] <- 0.12 #0.436
	# vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "old" & vlist$nsG$ethnicity == "b"] <- 0.4 #0.25
	# vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "old" & vlist$nsG$ethnicity == "h"] <- 0.4 #0.35
	vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "young" & vlist$nsG$ethnicity == "w"] <- 0.3 #0.436
	# vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "young" & vlist$nsG$ethnicity == "b"] <- 0.9 #0.25
	# vlist$nsG$pe[vlist$nsG$sex.risk.age.group == "young" & vlist$nsG$ethnicity == "h"] <- 0.8 #0.35

	## uio - Same across all cities
	# vlist$uio$pe[vlist$uio$gender  == "f"]<- 0.25     #0.207
	# vlist$uio$pe[vlist$uio$gender == "m"] <- 0.45     #0.286
	# vlist$uio$pe[vlist$uio$risk   == "msm"] <- 0.1 #0.291
   
	## uis - Same across all cities, even between old young
	# vlist$uis$pe[vlist$uis$sex.risk.age.group == "old"]  <- 0.8 #0.668
	# vlist$uis$pe[vlist$uis$sex.risk.age.group == "young"]<- 0.2 #0.668

	## d - Same across all cities -Pick over s as deemed more realistially to be different
	# vlist$d$pe <- 40 #75
	## s - same across all cities: 0.028

	#v2,v3 -Same across all cities
	# vlist$v2$pe <- 0.05  # 0.00878
	# vlist$v3$pe <- 0.05 # 0.018595296
	
	# #psi - same across cities: decrease to increase late diagnoses
	# vlist$psi$pe[vlist$psi$ethnicity == "w"] <- 0.035 # 0.02
	# vlist$psi$pe[vlist$psi$ethnicity == "b"] <- 0.035 # 0.02
	# vlist$psi$pe[vlist$psi$ethnicity == "h"] <- 0.035 # 0.02

	#3.3.1-ARTDropout500 T1_O: Similar by region 
	# vlist$T1_O$pe[vlist$T1_O$ethnicity == "w"] <- 0.007 #0.009195213
	# vlist$T1_O$pe[vlist$T1_O$ethnicity == "b"] <- 0.007 #0.008708514
	# vlist$T1_O$pe[vlist$T1_O$ethnicity == "h"] <- 0.007 #0.00563768

	#3.3.2-ARTReInitiation O_T: Similar by region - Pick over T1_O for less involvement
	# vlist$O_T$pe[vlist$O_T$ethnicity == "w"] <- 0.04 #0.01766174
	# vlist$O_T$pe[vlist$O_T$ethnicity == "b"] <- 0.04 #0.02072576
	# vlist$O_T$pe[vlist$O_T$ethnicity == "h"] <- 0.04 #0.01934913

	# vlist$I3_multi$pe[vlist$I3_multi$ethnicity == "w"] <- 0.8
	# vlist$I3_multi$pe[vlist$I3_multi$ethnicity == "b"] <- 0.8
	# vlist$I3_multi$pe[vlist$I3_multi$ethnicity == "h"] <- 0.8
	# vlist$prop_vs_T$pe[vlist$prop_vs_T$ethnicity == "w"] <- 0.95
	# vlist$prop_vs_T$pe[vlist$prop_vs_T$ethnicity == "b"] <- 0.8
	# vlist$prop_vs_T$pe[vlist$prop_vs_T$ethnicity == "h"] <- 0.9
	
	# vlist$sigmaO.FM$pe <- vlist$sigmaO.MF$pe
	# vlist$sigmaO.MF$pe <- vlist$sigmaO.FM$pe
	# vlist$sigmaS$pe <- vlist$sigmaS$pe *2
	
	# vlist$p1$pe     <- 0.9 #0.7-0.8 # increase help reduce IR trend, but already close to the ceiling
	# vlist$phi1_$pe  <- 0.001#0.1456954 - increase increase IR, not much impact, bring things up parallel no impact on trends
	# vlist$alpha1_$pe <- 0.9 #0.066163 # increase help reduce IR trend, but ROI is less and less
	# vlist$ass.eO$pe[vlist$ass.eO$sex.risk.age.group == "low"] <- 0 # ~ 0.6
	# vlist$ass.eO$pe[vlist$ass.eO$sex.risk.age.group == "high"] <- 0 # ~ 0.6
	# vlist$mix.prob.HL$pe <- 0 # 0.01
	
	# vlist$ass.eS$pe <- 0.2 # Same with other ass.eO, mixed effect, decrease -> reduce for B/H but increase for W ( i)
	# vlist$epsilonS$pe <- 0.95 # 0.68 - 95% reduction helps, but then again the trends for ndaig drop while ldiag increase
	# vlist$kappa_H$pe  <- 0.5 #  0.8
	# vlist$delta_sex$pe <- 0.99 # 0.95 - increase to 99% reduction helps
	# vlist$T2_T1$pe <- vlist$T2_T1$pe*10
	# vlist$T1_T2$pe <- vlist$T1_T2$pe/10
	# vlist$mor_T3$pe <- vlist$mor_T3$pe *10
	# 							 
	# vlist$prop.prep$pe[vlist$prop.prep$year == "2021"] <- vlist$prop.prep$pe[vlist$prop.prep$year == "2017"]
	# vlist$prop.ssp$pe[vlist$prop.ssp$year == "2021"] <- vlist$prop.ssp$pe[vlist$prop.ssp$year == "2017"]
	# vlist$prop.oat$pe[vlist$prop.oat$year == "2021"] <- vlist$prop.oat$pe[vlist$prop.oat$year == "2017"]
	# vlist$p1$pe[vlist$p1$year == "2022"]             <- vlist$p1$pe[vlist$p1$year == "2017"]
	
 # browser()	
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "m" & vlist$prop.high.sus$ethnicity == "w" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.0456 - slight lower
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "m" & vlist$prop.high.sus$ethnicity == "b" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.1067 - slight higher
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "m" & vlist$prop.high.sus$ethnicity == "h" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.053 - Higher (x2 times)
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "f" & vlist$prop.high.sus$ethnicity == "w" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.0167 - Lower for W (x1/2)
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "f" & vlist$prop.high.sus$ethnicity == "b" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.0206 - Lower for B (x1/2)
 # vlist$prop.high.sus$pe[vlist$prop.high.sus$gender == "f" & vlist$prop.high.sus$ethnicity == "h" & vlist$prop.high.sus$risk == "het"] <- 0.1 # 0.0109 - Higher for H (x8 times)
 # 
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "m" & vlist$prop.high.inf$ethnicity == "w" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.0912 - slight lower
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "m" & vlist$prop.high.inf$ethnicity == "b" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.2134 - slight higher
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "m" & vlist$prop.high.inf$ethnicity == "h" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.106 - Higher (x2 times)
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "f" & vlist$prop.high.inf$ethnicity == "w" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.0334 - Lower for W (x1/2)
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "f" & vlist$prop.high.inf$ethnicity == "b" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.0412 - Lower for B (x1/2)
 # vlist$prop.high.inf$pe[vlist$prop.high.inf$gender == "f" & vlist$prop.high.inf$ethnicity == "h" & vlist$prop.high.inf$risk == "het"] <- 0.5 # 0.0218 - slight lower
	
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "m" & vlist$T1_T2$ethnicity == "w"] <-c(0.0114,0.0077,0.0097,0.0091)  #0.01568938 0.01104667 0.01335415 0.01299676
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "m" & vlist$T1_T2$ethnicity == "b"] <-c(0.0144,0.0097,0.0122,0.0115)  #0.01589330 0.01119741 0.01353252 0.01317089
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "m" & vlist$T1_T2$ethnicity == "h"] <-c(0.0142,0.0096,0.0121,0.0114)  #0.01821785 0.01279905 0.01549150 0.01506880
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "f" & vlist$T1_T2$ethnicity == "w"] <-c(0.0114,0.0091)  # 0.01568938 0.01299676
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "f" & vlist$T1_T2$ethnicity == "b"] <-c(0.0144,0.0115) # 0.01589330 0.01317089
	# vlist$T1_T2$pe[vlist$T1_T2$gender == "f" & vlist$T1_T2$ethnicity == "h"] <-c(0.0142,0.0114)  # 0.01821785 0.01506880
	# 
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "m" & vlist$T1_T3$ethnicity == "w"] <-c(5.431057e-05,2.425084e-05,3.363399e-05,4.071277e-05)  #1.634494e-04 9.199502e-05 1.297629e-04 1.157771e-04
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "m" & vlist$T1_T3$ethnicity == "b"] <-c(8.369407e-05,3.749081e-05,5.160176e-05,6.330265e-05)  #0.0001872988 0.0001065733 0.0001498670 0.0001330940
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "m" & vlist$T1_T3$ethnicity == "h"] <-c(7.664574e-05,3.451266e-05,4.690573e-05,5.881245e-05)  #1.515281e-04 8.458902e-05 1.197753e-04 1.068987e-040
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "f" & vlist$T1_T3$ethnicity == "w"] <-c(5.431057e-05,4.071277e-05)  # 0.0001634494 0.0001157771
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "f" & vlist$T1_T3$ethnicity == "b"] <-c(8.369407e-05,6.330265e-05)  # 0.0001872988 0.0001330940
	# vlist$T1_T3$pe[vlist$T1_T3$gender == "f" & vlist$T1_T3$ethnicity == "h"] <-c(7.664574e-05,5.881245e-05)  # 0.0001515281 0.0001068987
	# 
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "m" & vlist$T2_T3$ethnicity == "w"] <-c(0.006927903,0.004291840,0.005444633,0.005463597)  #0.008136512 0.005901296 0.007025376 0.006836824
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "m" & vlist$T2_T3$ethnicity == "b"] <-c(0.008296383,0.005136189,0.006517093,0.006543455)  #0.007806070 0.005654422 0.006733779 0.006556683
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "m" & vlist$T2_T3$ethnicity == "h"] <-c(0.007445199,0.004608214,0.005845962,0.005873486)  #0.007024237 0.005061862 0.006055548 0.005871417
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "f" & vlist$T2_T3$ethnicity == "w"] <-c(0.006927903,0.005463597)  # 0.008136512 0.006836824
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "f" & vlist$T2_T3$ethnicity == "b"] <-c(0.008296383,0.006543455)  # 0.007806070 0.006556683
	# vlist$T2_T3$pe[vlist$T2_T3$gender == "f" & vlist$T2_T3$ethnicity == "h"] <-c(0.007445199,0.005873486)  # 0.007024237 0.005871417
	
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "m" & vlist$T2_T1$ethnicity == "w"] <-c(0.01688293,0.01974120,0.01910241,0.01746286)*10  #0.01427164 0.01708957 0.01552887 0.01572851
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "m" & vlist$T2_T1$ethnicity == "b"] <-c(0.01691035,0.01978452,0.01913763,0.01749784)*10  #0.01480496 0.01773416 0.01611310 0.01631958
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "m" & vlist$T2_T1$ethnicity == "h"] <-c(0.01771010,0.02070959,0.02003449,0.01832210)*10  #0.01270797 0.01518633 0.01382489 0.01397766
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "f" & vlist$T2_T1$ethnicity == "w"] <-c(0.01688293,0.01746286)*10  # 0.01427164 0.01572851
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "f" & vlist$T2_T1$ethnicity == "b"] <-c(0.01691035,0.01749784)*10  # 0.01480496 0.01631958
	# vlist$T2_T1$pe[vlist$T2_T1$gender == "f" & vlist$T2_T1$ethnicity == "h"] <-c(0.0177101,0.0183221)*10  # 0.01270797 0.01397766
	# 
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "m" & vlist$T3_T1$ethnicity == "w"] <-c(0.01688293,0.01974120,0.01910241,0.01746286)*10  #0.01427164 0.01708957 0.01552887 0.01572851
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "m" & vlist$T3_T1$ethnicity == "b"] <-c(0.01691035,0.01978452,0.01913763,0.01749784)*10  #0.01480496 0.01773416 0.01611310 0.01631958
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "m" & vlist$T3_T1$ethnicity == "h"] <-c(0.01771010,0.02070959,0.02003449,0.01832210)*10  #0.01270797 0.01518633 0.01382489 0.01397766
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "f" & vlist$T3_T1$ethnicity == "w"] <-c(0.01688293,0.01746286)*10  # 0.01427164 0.01572851
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "f" & vlist$T3_T1$ethnicity == "b"] <-c(0.01691035,0.01749784)*10  # 0.01480496 0.01631958
	# vlist$T3_T1$pe[vlist$T3_T1$gender == "f" & vlist$T3_T1$ethnicity == "h"] <-c(0.0177101,0.0183221)*10  # 0.01270797 0.01397766
	# 
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "m" & vlist$T3_T2$ethnicity == "w"] <-c(0.01688293,0.01974120,0.01910241,0.01746286)*10  #0.01427164 0.01708957 0.01552887 0.01572851
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "m" & vlist$T3_T2$ethnicity == "b"] <-c(0.01691035,0.01978452,0.01913763,0.01749784)*10  #0.01480496 0.01773416 0.01611310 0.01631958
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "m" & vlist$T3_T2$ethnicity == "h"] <-c(0.01771010,0.02070959,0.02003449,0.01832210)*10  #0.01270797 0.01518633 0.01382489 0.01397766
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "f" & vlist$T3_T2$ethnicity == "w"] <-c(0.01688293,0.01746286)*10  # 0.01427164 0.01572851
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "f" & vlist$T3_T2$ethnicity == "b"] <-c(0.01691035,0.01749784)*10  # 0.01480496 0.01631958
	# vlist$T3_T2$pe[vlist$T3_T2$gender == "f" & vlist$T3_T2$ethnicity == "h"] <-c(0.0177101,0.0183221)*10  # 0.01270797 0.01397766
	
	# vlist$ass.eO$pe[vlist$ass.eO$ethnicity == "w"]<- c(0.717,0.481) #0.692 0.540
	# vlist$ass.eO$pe[vlist$ass.eO$ethnicity == "b"]<- c(0.635,0.740) #0.763 0.668
	# vlist$ass.eO$pe[vlist$ass.eO$ethnicity == "h"]<- c(0.505,0.358) #0.615 0.491
	
	# 3.2.2-DiagOnART500	
	# vlist$phi1_$pe[vlist$phi1_$ethnicity == "w"] <- 0.7 #0.1456954
	# vlist$phi1_$pe[vlist$phi1_$ethnicity == "b"] <- 0.5 #0.1456954
	# vlist$phi1_$pe[vlist$phi1_$ethnicity == "h"] <- 0.5 #0.1456954
	# vlist$phi12$pe <- 1.219122# 1.738
	# vlist$phi13$pe <- 1.891889# 2.866239
	
	# vlist$mor_T2$pe <- c(0.001402457, 0.000627895, 0.000627895, 0.000984972, 0.000984972, 0.000895693, 0.001131864, 0.000506915, 0.000506915, 0.000795637, 0.000795637, 0.000722606, 0.000768049,0.000343906, 0.000343906, 0.000540286, 0.000540286, 0.000489929, 0.001422027, 0.000907699, 0.001147956, 0.000732457, 0.000779927, 0.000497106) #0.001323565 0.000594733 0.000594733 0.000935709 0.000935709 0.000843974 0.001070396 0.000481267 0.000481267 0.000757995 0.000757995 0.000682048 0.000729865 0.000328039 0.000328039 0.000517844 0.000517844 0.000464147 0.001346843 0.000857871 0.001089742 0.000693567 0.000744824 0.000472884
	# vlist$mor_T3$pe <- c(0.004256345, 0.002953857, 0.002953857, 0.004749318, 0.004749318, 0.002647791, 0.003836675, 0.002663693, 0.002663693, 0.004280957, 0.004280957, 0.002387884, 0.003281019,0.002275964, 0.002275964, 0.003661118, 0.003661118, 0.002040071, 0.004570418, 0.002843129, 0.004118037, 0.002562942, 0.003521186, 0.002189316) #0.007982172 0.005556425 0.005556425 0.008908776 0.008908776 0.004980047 0.007182611 0.005004484 0.005004484 0.008017165 0.008017165 0.004485448 0.006162769 0.004285613 0.004285613 0.006877311 0.006877311 0.003841511 0.008559039 0.005339918 0.007698117 0.004807281 0.006605118 0.004117095
}

# Incoporate the NEW values above into vparameters
source("ModelCoreModules/Parameter_input/Update_vlist_to_vparameters.R")
vparameters <- par_fix(vparameters, vlist)

# read in initial population
#	source("Model_development/7.Red_strata/7.Red_strata-Par_pop.R")	#Source do not use update vlist value so have to convert it into function
result      <- par_pop(vparameters, vlist)  #Function part of Update_vlist_to_vparameters
vparameters <- result$vparameters
x           <- result$x
# browser()
##########################################################
#Run sensitivities
vparams_sens <- vparameters
out_sens <- intervention.model.combination(input.parameters = vparams_sens,
																					 current.int = interventions[cc])
sens <- extract_outputs(out_sens)

# ##########################################################
# #Build plots (baseline vs sens vs target), 3×2 grid
# W_diag <- plot_line(cbind(base$diag_race[1:6,1], sens$diag_race[1:6,1], calib.target$diag.obs[,"W"]),
# 										c("Baseline","Sens","Target"),
# 										TITLE="Total diagnoses W", CITY=CITY, y_min=0)
# B_diag <- plot_line(cbind(base$diag_race[1:6,2], sens$diag_race[1:6,2], calib.target$diag.obs[,"B"]),
# 										c("Baseline","Sens","Target"),
# 										TITLE="Total diagnoses B", CITY=CITY, y_min=0)
# H_diag <- plot_line(cbind(base$diag_race[1:6,3], sens$diag_race[1:6,3], calib.target$diag.obs[,"H"]),
# 										c("Baseline","Sens","Target"),
# 										TITLE="Total diagnoses H", CITY=CITY, y_min=0)
# 
# W_ndiag <- plot_line(cbind(base$ndiag_race[1:6,1], sens$ndiag_race[1:6,1], calib.target$ndiag.obs[,"W"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="New diagnoses W", CITY=CITY, y_min=0) + theme(legend.position = "none")
# B_ndiag <- plot_line(cbind(base$ndiag_race[1:6,2], sens$ndiag_race[1:6,2], calib.target$ndiag.obs[,"B"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="New diagnoses B", CITY=CITY, y_min=0) + theme(legend.position = "none")
# H_ndiag <- plot_line(cbind(base$ndiag_race[1:6,3], sens$ndiag_race[1:6,3], calib.target$ndiag.obs[,"H"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="New diagnoses H", CITY=CITY, y_min=0) + theme(legend.position = "none")
# 
# Young_ndiag <- plot_line(cbind(base$ndiag_risk[1:6,1], sens$ndiag_risk[1:6,1], calib.target$ndiag.obs[,"Young"]),
# 											c("Baseline","Sens","Target"),
# 											TITLE="New diagnoses Young", CITY=CITY, y_min=0) + theme(legend.position = "none")
# Old_ndiag <- plot_line(cbind(base$ndiag_risk[1:6,2], sens$ndiag_risk[1:6,2],  calib.target$ndiag.obs[,"Old"]),
# 											c("Baseline","Sens","Target"),
# 											TITLE="New diagnoses Old", CITY=CITY, y_min=0) + theme(legend.position = "none")
# PWID_ndiag <- plot_line(cbind(base$ndiag_risk[1:6,3], sens$ndiag_risk[1:6,3], calib.target$ndiag.obs[,"PWID"]),
# 											c("Baseline","Sens","Target"),
# 											TITLE="New diagnoses PWID", CITY=CITY, y_min=0) + theme(legend.position = "none")
# het.m_ndiag <- plot_line(cbind(base$ndiag_risk[1:6,4], sens$ndiag_risk[1:6,4], calib.target$ndiag.obs[,"M.HET"]),
# 											c("Baseline","Sens","Target"),
# 											TITLE="New diagnoses het.m", CITY=CITY, y_min=0) + theme(legend.position = "none")
# het.f_ndiag <- plot_line(cbind(base$ndiag_risk[1:6,5], sens$ndiag_risk[1:6,5], calib.target$ndiag.obs[,"F.HET"]),
# 											c("Baseline","Sens","Target"),
# 											TITLE="New diagnoses het.f", CITY=CITY, y_min=0) + theme(legend.position = "none")
# 
# W_diag_late <- plot_line(cbind(base$diag_late[1:6,1], sens$diag_late[1:6,1], calib.target$late_diag.obs[,"W"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="late diagnoses W", CITY=CITY, y_min=0) + theme(legend.position = "none")
# B_diag_late <- plot_line(cbind(base$diag_late[1:6,2], sens$diag_late[1:6,2], calib.target$late_diag.obs[,"B"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="late diagnoses B", CITY=CITY, y_min=0) + theme(legend.position = "none")
# H_diag_late <- plot_line(cbind(base$diag_late[1:6,3], sens$diag_late[1:6,3], calib.target$late_diag.obs[,"H"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="late diagnoses H", CITY=CITY, y_min=0) + theme(legend.position = "none")
# 
# W_viral_sup <- plot_line(cbind(base$viral_sup[1:6,1], sens$viral_sup[1:6,1], calib.target$viral_sup.obs[,"W"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="viral sup W", CITY=CITY, y_min=0) + theme(legend.position = "none")
# B_viral_sup <- plot_line(cbind(base$viral_sup[1:6,2], sens$viral_sup[1:6,2], calib.target$viral_sup.obs[,"B"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="viral sup B", CITY=CITY, y_min=0) + theme(legend.position = "none")
# H_viral_sup <- plot_line(cbind(base$viral_sup[1:6,3], sens$viral_sup[1:6,3], calib.target$viral_sup.obs[,"H"]),
# 										 c("Baseline","Sens","Target"),
# 										 TITLE="viral sup H", CITY=CITY, y_min=0) + theme(legend.position = "none")
# 
# inc_bo <- plot_line(cbind(base$inc[,1], sens$inc[,1]),
# 												 c("Baseline","Sens"),
# 												 TITLE="Inc bo", CITY=CITY, y_min=0) + theme(legend.position = "none")
# inc_bs <- plot_line(cbind(base$inc[,2], sens$inc[,2]),
# 												 c("Baseline","Sens"),
# 												 TITLE="Inc bs", CITY=CITY, y_min=0) + theme(legend.position = "none")
# inc_g  <- plot_line(cbind(base$inc[,3], sens$inc[,3]),
# 												 c("Baseline","Sens"),
# 												 TITLE="Inc g", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_bo <- plot_line(cbind(base$ir_source[,1], sens$ir_source[,1])*100000,
# 										c("Baseline","Sens"),
# 										TITLE="IR bo per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_bs <- plot_line(cbind(base$ir_source[,2], sens$ir_source[,2])*100000,
# 										c("Baseline","Sens"),
# 										TITLE="IR bs per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_g  <- plot_line(cbind(base$ir_source[,3], sens$ir_source[,3])*100000,
# 										c("Baseline","Sens"),
# 										TITLE="IR g per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_w <- plot_line(cbind(base$ir_race[,1], sens$ir_race[,1])*100000,
# 									 c("Baseline","Sens"),
# 									 TITLE="IR W per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_b <- plot_line(cbind(base$ir_race[,2], sens$ir_race[,2])*100000,
# 									 c("Baseline","Sens"),
# 									 TITLE="IR B per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")
# ir_h  <- plot_line(cbind(base$ir_race[,3], sens$ir_race[,3])*100000,
# 									 c("Baseline","Sens"),
# 									 TITLE="IR H per 100k", CITY=CITY, y_min=0) + theme(legend.position = "none")

# g_all <- ggpubr::ggarrange(
# 	W_diag, B_diag, H_diag,
# 	W_ndiag, B_ndiag, H_ndiag,
# 	Young_ndiag, Old_ndiag,
# 	PWID_ndiag, het.m_ndiag, het.f_ndiag,
# 	W_diag_late, B_diag_late, H_diag_late,
# 	W_viral_sup, B_viral_sup, H_viral_sup,
# 	inc_bo, inc_bs, inc_g,
# 	ir_bo, ir_bs, ir_g,
# 	ir_w, ir_b, ir_h,
# 	ncol=5, nrow=6, common.legend=TRUE, legend="right"
# #	,labels=c("A","B","C","D","E","F","G","H")
# )

make_calibration_plots <- function(base, sens = NULL, CITY) {

	safe_cbind3 <- function(b, s, t) {
		if (is.null(s)) return(cbind(b, t))
		cbind(b, s, t)
	}

	safe_cbind2 <- function(b, s) {
		if (is.null(s)) return(cbind(b))
		cbind(b, s)
	}

	labs3 <- function() if (is.null(sens)) c("Baseline","Target") else c("Baseline","Sens","Target")
	labs2 <- function() if (is.null(sens)) c("Baseline") else c("Baseline","Sens")

	# ================= Diagnoses =================
	W_diag <- plot_line(safe_cbind3(base$diag_race[1:6,1], sens$diag_race[1:6,1], calib.target$diag.obs[,"W"]),
											labs3(), TITLE="Total diagnoses W", CITY=CITY, y_min=0)

	B_diag <- plot_line(safe_cbind3(base$diag_race[1:6,2], sens$diag_race[1:6,2], calib.target$diag.obs[,"B"]),
											labs3(), TITLE="Total diagnoses B", CITY=CITY, y_min=0)

	H_diag <- plot_line(safe_cbind3(base$diag_race[1:6,3], sens$diag_race[1:6,3], calib.target$diag.obs[,"H"]),
											labs3(), TITLE="Total diagnoses H", CITY=CITY, y_min=0)

	# ================= New diagnoses (race) =================
	W_ndiag <- plot_line(safe_cbind3(base$ndiag_race[1:6,1], sens$ndiag_race[1:6,1], calib.target$ndiag.obs[,"W"]),
											 labs3(), TITLE="New diagnoses W", CITY=CITY, y_min=0) + theme(legend.position="none")

	B_ndiag <- plot_line(safe_cbind3(base$ndiag_race[1:6,2], sens$ndiag_race[1:6,2], calib.target$ndiag.obs[,"B"]),
											 labs3(), TITLE="New diagnoses B", CITY=CITY, y_min=0) + theme(legend.position="none")

	H_ndiag <- plot_line(safe_cbind3(base$ndiag_race[1:6,3], sens$ndiag_race[1:6,3], calib.target$ndiag.obs[,"H"]),
											 labs3(), TITLE="New diagnoses H", CITY=CITY, y_min=0) + theme(legend.position="none")

	# ================= New diagnoses (risk) =================
	Young_ndiag <- plot_line(safe_cbind3(base$ndiag_risk[1:6,1], sens$ndiag_risk[1:6,1], calib.target$ndiag.obs[,"Young"]),
											 labs3(), TITLE="New diagnoses Young", CITY=CITY, y_min=0) + theme(legend.position="none")
	
	Old_ndiag   <- plot_line(safe_cbind3(base$ndiag_risk[1:6,2], sens$ndiag_risk[1:6,2], calib.target$ndiag.obs[,"Old"]),
											 labs3(), TITLE="New diagnoses Old"  , CITY=CITY, y_min=0) + theme(legend.position="none")
	
	PWID_ndiag  <- plot_line(safe_cbind3(base$ndiag_risk[1:6,3], sens$ndiag_risk[1:6,3], calib.target$ndiag.obs[,"PWID"]),
											 labs3(), TITLE="New diagnoses PWID" , CITY=CITY, y_min=0) + theme(legend.position="none")
	
	het.m_ndiag <- plot_line(safe_cbind3(base$ndiag_risk[1:6,4], sens$ndiag_risk[1:6,4], calib.target$ndiag.obs[,"M.HET"]),
													 labs3(), TITLE="New diagnoses het.m", CITY=CITY, y_min=0) + theme(legend.position="none")
	
	het.f_ndiag <- plot_line(safe_cbind3(base$ndiag_risk[1:6,5], sens$ndiag_risk[1:6,5], calib.target$ndiag.obs[,"F.HET"]),
													 labs3(), TITLE="New diagnoses het.f", CITY=CITY, y_min=0) + theme(legend.position="none")

	# ================= Late diagnoses =================
	W_diag_late <- plot_line(safe_cbind3(base$diag_late[1:6,1], sens$diag_late[1:6,1], calib.target$late_diag.obs[,"W"]),
													 labs3(), TITLE="late diagnoses W", CITY=CITY, y_min=0) + theme(legend.position="none")

	B_diag_late <- plot_line(safe_cbind3(base$diag_late[1:6,2], sens$diag_late[1:6,2], calib.target$late_diag.obs[,"B"]),
													 labs3(), TITLE="late diagnoses B", CITY=CITY, y_min=0) + theme(legend.position="none")

	H_diag_late <- plot_line(safe_cbind3(base$diag_late[1:6,3], sens$diag_late[1:6,3], calib.target$late_diag.obs[,"H"]),
													 labs3(), TITLE="late diagnoses H", CITY=CITY, y_min=0) + theme(legend.position="none")

	# ================= Viral suppression =================
	W_viral_sup <- plot_line(safe_cbind3(base$viral_sup[1:6,1], sens$viral_sup[1:6,1], calib.target$viral_sup.obs[,"W"]),
													 labs3(), TITLE="viral sup W", CITY=CITY, y_min=0) + theme(legend.position="none")

	B_viral_sup <- plot_line(safe_cbind3(base$viral_sup[1:6,2], sens$viral_sup[1:6,2], calib.target$viral_sup.obs[,"B"]),
													 labs3(), TITLE="viral sup B", CITY=CITY, y_min=0) + theme(legend.position="none")

	H_viral_sup <- plot_line(safe_cbind3(base$viral_sup[1:6,3], sens$viral_sup[1:6,3], calib.target$viral_sup.obs[,"H"]),
													 labs3(), TITLE="viral sup H", CITY=CITY, y_min=0) + theme(legend.position="none")

	# ================= Incidence =================
	# inc_bo <- plot_line(safe_cbind2(base$inc[,1], sens$inc[,1]), labs2(), TITLE="Inc bo", CITY=CITY, y_min=0) + theme(legend.position="none")
	# inc_bs <- plot_line(safe_cbind2(base$inc[,2], sens$inc[,2]), labs2(), TITLE="Inc bs", CITY=CITY, y_min=0) + theme(legend.position="none")
	# inc_g  <- plot_line(safe_cbind2(base$inc[,3], sens$inc[,3]), labs2(), TITLE="Inc g",  CITY=CITY, y_min=0) + theme(legend.position="none")
	# 
	# ir_bo <- plot_line(safe_cbind2(base$ir_source[,1], sens$ir_source[,1])*100000, labs2(), TITLE="IR bo per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	# ir_bs <- plot_line(safe_cbind2(base$ir_source[,2], sens$ir_source[,2])*100000, labs2(), TITLE="IR bs per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	# ir_g  <- plot_line(safe_cbind2(base$ir_source[,3], sens$ir_source[,3])*100000, labs2(), TITLE="IR g per 100k",  CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_young <- plot_line(safe_cbind2(base$ir_risk[,1], sens$ir_risk[,1])*100000, labs2(), TITLE="IR young per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_old   <- plot_line(safe_cbind2(base$ir_risk[,2], sens$ir_risk[,2])*100000, labs2(), TITLE="IR old per 100k"  , CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_pwid  <- plot_line(safe_cbind2(base$ir_risk[,3], sens$ir_risk[,3])*100000, labs2(), TITLE="IR pwid per 100k" , CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_het.m <- plot_line(safe_cbind2(base$ir_risk[,4], sens$ir_risk[,4])*100000, labs2(), TITLE="IR het.m per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_het.f <- plot_line(safe_cbind2(base$ir_risk[,5], sens$ir_risk[,5])*100000, labs2(), TITLE="IR het.f per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_low   <- plot_line(safe_cbind2(base$ir_risk[,6], sens$ir_risk[,6])*100000, labs2(), TITLE="IR het.low per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_high  <- plot_line(safe_cbind2(base$ir_risk[,7], sens$ir_risk[,7])*100000, labs2(), TITLE="IR het.high per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	
	ir_w <- plot_line(safe_cbind2(base$ir_race[,1], sens$ir_race[,1])*100000, labs2(), TITLE="IR W per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_b <- plot_line(safe_cbind2(base$ir_race[,2], sens$ir_race[,2])*100000, labs2(), TITLE="IR B per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	ir_h <- plot_line(safe_cbind2(base$ir_race[,3], sens$ir_race[,3])*100000, labs2(), TITLE="IR H per 100k", CITY=CITY, y_min=0) + theme(legend.position="none")
	# ================= Proportion =================
	prop_high.inf.among.het <- plot_line(safe_cbind2(base$prop_high.inf.among.het, sens$prop_high.inf.among.het),
																			 labs2(), TITLE="prop_high.inf.among.het", CITY=CITY, y_min=0)
	prop_T_W <- plot_line(safe_cbind2(base$prop_T_race[,1], sens$prop_T_race[,1]),
											labs2(), TITLE="Prop T Diag_PLHIV W", CITY=CITY, y_min=0)
	prop_T_B <- plot_line(safe_cbind2(base$prop_T_race[,2], sens$prop_T_race[,2]),
												labs2(), TITLE="Prop T Diag_PLHIV B", CITY=CITY, y_min=0)
	prop_T_H <- plot_line(safe_cbind2(base$prop_T_race[,3], sens$prop_T_race[,3]),
												labs2(), TITLE="Prop T Diag_PLHIV H", CITY=CITY, y_min=0)
	
	ggpubr::ggarrange(
		W_diag, B_diag, H_diag,
		W_ndiag, B_ndiag, H_ndiag,
		Young_ndiag, Old_ndiag, PWID_ndiag, het.m_ndiag, het.f_ndiag,
		W_diag_late, B_diag_late, H_diag_late,
		W_viral_sup, B_viral_sup, H_viral_sup,
#		inc_bo, inc_bs, inc_g,
		# ir_bo, ir_bs, ir_g,
# ir_young,ir_old,ir_pwid,ir_het.m,ir_het.f,ir_low ,ir_high,
		ir_w, ir_b, ir_h,
# prop_high.inf.among.het ,
    # prop_T_W,prop_T_B,prop_T_H,
#		ncol = 5, nrow = 5, 
    common.legend = TRUE, legend = "right"
	)
}

g_all <- make_calibration_plots(base, sens, CITY)


if (interactive()) {
	if (Sys.info()[["sysname"]] == "Windows") windows(width=20, height=12) 
	else if (Sys.info()[["sysname"]] == "Darwin") quartz(width=20, height=12)
	else x11(width=20, height=12, type="cairo")
}
print(g_all)
# ggplot2::ggsave("DSA_total_and_new_diag_3x2.png", g_all, width=16, height=10, dpi=300)