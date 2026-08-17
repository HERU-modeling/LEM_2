######################
## READ in calibrated free parameters ##
## Derive free parameter information and values ##
calpar <- read.xlsx("Data_files/calib_par/cali_par.xlsx", sheet = CITY) # Load free parameters
names <- unique(calpar$par) # names of free.par
leng <- as.numeric(table(factor(calpar$par, levels = names))) # length of each free.par

# Assigning parameter values for selected free parameters #
par2 <- as.list(leng)
for (l in seq_along(leng)) {
	if (l == 1) {
		par2[[l]] <- calpar$pe[1:leng[l]]
	} else {
		par2[[l]] <- calpar$pe[(sum(leng[1:(l - 1)]) + 1):sum(leng[1:l])]
	}
}
names(par2) <- names


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
#browser()
source("ModelCoreModules/Parameter_input/Update_vlist_to_vparameters.R")
vparameters <- par_fix(vparameters, vlist)

# read in initial population
result      <- par_pop(vparameters, vlist)  #Function part of Update_vlist_to_vparameters
vparameters <- result$vparameters
x           <- result$x
#browser()