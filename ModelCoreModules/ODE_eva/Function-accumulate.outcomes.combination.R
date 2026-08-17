## Cascade CEA Model - Costs & QALYs calculations
## Latest update: Oct 2022
#############################################################################
## Output: 
## out.sum: list of size 12 contains interesting quantities for CEA
#############################################################################

#' Calculate Costs & QALYs of interventions
#'
#' @param input.parameters Provide in RunModelAndAccum in load_run, often vparameter
#' @param pop.in A tensor (3D matrix) with 1st dimension being time, 2nd dimension being the 19 HIV states, and 3rd dimension being the 42 subgroups.
#' @param current.int Provide in RunModelAndAccum in load_run, often a vector list of intervention (or no intervention)
#' @param sums.only The function has a conditional check on sums.only to determine whether to compute and return detailed results (sums.only == FALSE) or only aggregate sums (sums.only == TRUE).
#'
#' @return out.sum A list of size 12 contains total QALYs and total costs
#' @export
#'
#' @examples
accum.outcomes.combination <- function(input.parameters, 
																			 pop.in, 
																			 current.int,
																			 sums.only = TRUE){

	tf.ls   <- true.false.interventions.combination(current.int)
	
	if (!all(current.int == c("No interventions"))){
		# Set intervention parameters for scale*effect 
		int.pop.eff <- scale.interventions.combination(current.int, tf.ls)
	}	
	

	
 with(as.list(c(Int.Scale.ls, Int.Eff.ls, Int.Baseline.ls, All.Costs.ls, input.parameters, gp.gn)),{
	
	# Holding arrays
	QALYs.period            <- vector("list", dim(pop.in)[1])
	hru.costs.period        <- vector("list", dim(pop.in)[1])
	art.costs.period        <- vector("list", dim(pop.in)[1])
	art.ini.costs.period    <- vector("list", dim(pop.in)[1])
	oat.costs.period        <- vector("list", dim(pop.in)[1])
	prep.costs.period       <- vector("list", dim(pop.in)[1])
	prep.tests.costs.period <- vector("list", dim(pop.in)[1])
	tests.costs.period      <- vector("list", dim(pop.in)[1])
	
	# Intervention costs
	int.costs.period        <- vector("list", dim(pop.in)[1])
	int.impl.costs.period   <- vector("list", dim(pop.in)[1])
	int.sust.costs.period   <- vector("list", dim(pop.in)[1])
	
	for(i in 1:(dim(pop.in)[1])){
		
		CurrentPeriod <- pop.in[i, , ]
		rownames(CurrentPeriod) <- state.name[c(1:16)]
		
		# Only consider evaluation period for accumulating outcomes 2020-2040
		if(i > int.start & i <= (int.start + (20 * 12))){
			
			# QALYs
			QALYs.period[[i]]            <- ((CurrentPeriod * StateQALYs)     / ((1 + Discounting / 12) ^ (i - int.start)))
			
			# STATE COSTS
			hru.costs.period[[i]]        <- ((CurrentPeriod * StateCostsHRU)  / ((1 + Discounting / 12) ^ (i - int.start)))
			art.costs.period[[i]]        <- ((CurrentPeriod * StateCostsART)  / ((1 + Discounting / 12) ^ (i - int.start)))
			oat.costs.period[[i]]        <- ((CurrentPeriod * StateCostsOAT)  / ((1 + Discounting / 12) ^ (i - int.start)))
			prep.costs.period[[i]]       <- ((CurrentPeriod * StateCostsPrEP) / ((1 + Discounting / 12) ^ (i - int.start)))
			prep.tests.costs.period[[i]] <- ((CurrentPeriod * psi.p * StateCostsTests.psi.p) / 
																			 	((1 + Discounting / 12) ^ (i - int.start)))
			
			####Initialize intervention cost, int, impl and sust cost = 0 outside impl&sust period
			int.costs.period[[i]]      <- (CurrentPeriod * 0)
			int.impl.costs.period[[i]] <- (CurrentPeriod * 0)
			int.sust.costs.period[[i]] <- (CurrentPeriod * 0)

			
			############################    
			# TESTING COST: time-varying testing rates

			## TESTING
			if (i <= 12) psi.t <- psi.m[1, ]          # Jan 2017 - Dec 2017 
			if (i > 12 & i <= 24) psi.t <- psi.m[2, ] # Jan 2018 - Dec 2018 
			if (i > 24 & i <= 36) psi.t <- psi.m[3, ] # Jan 2019 - Dec 2019 
			if (i > 36) psi.t <- psi.m[4, ]           # Dec 2020 onward
				
			psi <- matrix(0, nrow = 30, ncol = 16)
			colnames(psi) <- vparameters$state.name[1:16]		
			psi[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t	
#			psi[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi[ , c("S1", "Ia", "I1", "I2", "I3")] * (1-prop.prep.t) + psi.p * prop.prep.t # Compared to ODE, skip impact of interventions & PrEP tesing rate, need revisit
			psi[ , c("I2")] <- psi[ , c("I2")] + v2
			psi[ , c("I3")] <- psi[ , c("I3")] + v3
			psi <- t(psi)
		
			# Testing costs
#			tests.costs.period[[i]] <- ((CurrentPeriod * psi.matrix * StateCostsTests.psi) /
#																		((1 + Discounting / 12) ^ (i - int.start)))
			tests.costs.period[[i]] <- ((CurrentPeriod * psi * StateCostsTests.psi) /
																		((1 + Discounting / 12) ^ (i - int.start)))			
			
	
			## INTERVENTION 2020 -2030
			if(tf.ls$intervention.testing.modif==TRUE & i > int.start & i <= int.end){
					if(is.element(c("Opt-out testing (ER)"), current.int)){
						IntCosts.Tests <- IntCosts.test.optout.ER
						IntCosts.impl.Tests <- IntCosts.impl.test.optout.ER
						IntCosts.sust.Tests <- IntCosts.sust.test.optout.ER
						IntCosts.impl.lump.Tests <- IntCosts.impl.lump.test.optout.ER
						IntCosts.sust.lump.Tests <- IntCosts.sust.lump.test.optout.ER
						
						temp.int     <- IntEff.test.optout.ER # final scale-up multiplier
						scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
															rep(temp.int, (int.end - int.start - scale.up.period))) # scale-up multiplier vector 
						#      scale.v      <- IntScale.test.optout.ER["S1", ] # original level 42 group
						scale.m      <- IntScale.test.optout.ER
						scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.optout.ER) * IntSetting.test.optout.ER / 4,
																												nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
						scale.m["S1", ]                   <-        colSums(IntScale.test.optout.ER) * (1 - IntSetting.test.optout.ER)
						#      scale.m.temp <- t(replicate((int.sus * 12), scale.v)) # total scale-up month as rows, 42 groups in col      
						scale.m.temp <- replicate((int.sus * 12), scale.m)
						scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))

						int.test.eff <- scale.m.temp * (scale.up.v - 1) + 1
						
						psi.temp <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.temp) <- vparameters$state.name[1:16]						
						psi.temp[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.temp <- psi.temp * t(int.test.eff[(i - int.start), ,])
						psi.temp[ ,c("I2")] <- psi.temp[ , c("I2")] + v2
						psi.temp[ ,c("I3")] <- psi.temp[ , c("I3")] + v3
						psi.temp <- t(psi.temp)
						
#						psi.matrix.full <- matrix(0, nrow = 30, ncol = 16)
#						colnames(psi.matrix.full) <- state.name[-c(17:21)]
#						psi.matrix.full[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi
#						psi.matrix.full <- t(psi.matrix.full)
#						psi.matrix.full <- psi.matrix.full * int.test.eff[scale.up.period, ]
#						psi.matrix.full["I2", ] <- psi.matrix.full["I2", ] + v2
#						psi.matrix.full["I3", ] <- psi.matrix.full["I3", ] + v3
						psi.full <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.full) <- vparameters$state.name[1:16]						
						psi.full[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.full <- psi.full * t(int.test.eff[(i - int.start), ,])
						psi.full[ ,c("I2")] <- psi.full[ , c("I2")] + v2
						psi.full[ ,c("I3")] <- psi.full[ , c("I3")] + v3
						psi.full <- t(psi.full)
						
						# Incremental testing costs
#						testing.costs.mod <- ((CurrentPeriod * (psi.matrix.temp - psi.matrix) * StateCostsTests.psi) 
#																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						testing.costs.mod <- ((CurrentPeriod * (psi.temp - psi) * StateCostsTests.psi) 
																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
#						int.costs.mod <- ((CurrentPeriod * psi.matrix.temp * IntCosts.Tests) 
#															/ ((1 + Discounting / 12) ^ (i - int.start)))
						int.costs.mod <- ((CurrentPeriod * psi.temp * IntCosts.Tests) 
															/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if(i <= (scale.up.period + int.start)){
#							int.impl.costs.temp <- (((CurrentPeriod * psi.matrix.full * IntCosts.impl.Tests) + IntCosts.impl.lump.Tests) /
#																				((1 + Discounting / 12) ^ (i - int.start)))
							int.impl.costs.temp <- (((CurrentPeriod * psi.full * IntCosts.impl.Tests) + IntCosts.impl.lump.Tests) /
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
#							int.sust.costs.temp <- (((CurrentPeriod * psi.matrix.temp * IntCosts.sust.Tests) + IntCosts.sust.lump.Tests)/ 
#																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (((CurrentPeriod * psi.temp * IntCosts.sust.Tests) + IntCosts.sust.lump.Tests)/ 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- (int.costs.mod + testing.costs.mod) + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
					}
					if(is.element(c("Opt-out testing (PC)"), current.int)){
						IntCosts.Tests <- IntCosts.test.optout.PC
						IntCosts.impl.Tests <- IntCosts.impl.test.optout.PC
						IntCosts.sust.Tests <- IntCosts.sust.test.optout.PC
						IntCosts.impl.lump.Tests <- IntCosts.impl.lump.test.optout.PC
						IntCosts.sust.lump.Tests <- IntCosts.sust.lump.test.optout.PC
						
						temp.int     <- IntEff.test.optout.PC
						scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
															rep(temp.int, (int.end - int.start - scale.up.period)))
						#      scale.v      <- IntScale.test.optout.PC["S1", ]
						scale.m      <- IntScale.test.optout.PC
						scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.optout.PC) * IntSetting.test.optout.PC / 4,
																												nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
						scale.m["S1", ]                   <-        colSums(IntScale.test.optout.PC) * (1 - IntSetting.test.optout.PC)      
						#      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
						scale.m.temp <- replicate((int.sus * 12), scale.m)
						scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
						int.test.eff <- scale.m.temp * (scale.up.v - 1) + 1
						
						psi.temp <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.temp) <- vparameters$state.name[1:16]						
						psi.temp[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.temp <- psi.temp * t(int.test.eff[(i - int.start), ,])
						psi.temp[ ,c("I2")] <- psi.temp[ , c("I2")] + v2
						psi.temp[ ,c("I3")] <- psi.temp[ , c("I3")] + v3
						psi.temp <- t(psi.temp)
						
						psi.full <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.full) <- vparameters$state.name[1:16]						
						psi.full[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.full <- psi.full * t(int.test.eff[(i - int.start), ,])
						psi.full[ ,c("I2")] <- psi.full[ , c("I2")] + v2
						psi.full[ ,c("I3")] <- psi.full[ , c("I3")] + v3
						psi.full <- t(psi.full)
						
						# Incremental testing costs
						testing.costs.mod <- ((CurrentPeriod * (psi.temp - psi) * StateCostsTests.psi) 
																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
						int.costs.mod <- ((CurrentPeriod * psi.temp * IntCosts.Tests) 
															/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if(i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- (((CurrentPeriod * psi.full * IntCosts.impl.Tests) + IntCosts.impl.lump.Tests)/
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- (((CurrentPeriod * psi.temp * IntCosts.sust.Tests) + IntCosts.sust.lump.Tests)/ 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- (int.costs.mod + testing.costs.mod) + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
					} 
					if(is.element(c("EMR testing reminder"), current.int)){
						IntCosts.Tests <- IntCosts.test.EMR
						IntCosts.impl.Tests <- IntCosts.impl.test.EMR
						IntCosts.sust.Tests <- IntCosts.sust.test.EMR
						IntCosts.impl.lump.Tests <- IntCosts.impl.lump.test.EMR
						IntCosts.sust.lump.Tests <- IntCosts.sust.lump.test.EMR
						
						temp.int     <- IntEff.test.EMR
						scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
															rep(temp.int, (int.end - int.start - scale.up.period)))
						scale.m      <- IntScale.test.EMR
						reweight    <- sum(inits[,c("Ia","I1","I2","I3")])/sum(inits[,c("S1")])*((1-IntPosRate.test.EMR)/IntPosRate.test.EMR)
						scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.EMR) * 1/(1 + reweight) / 4,
																												nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
						scale.m["S1", ]                   <-        colSums(IntScale.test.EMR) * reweight/(1 + reweight)
						scale.m.temp <- replicate((int.sus * 12), scale.m)
						scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
						int.test.eff <- scale.m.temp * (scale.up.v - 1) + 1
#browser() #Compare this to 	int.pop.eff$int.test.eff 					
						psi.temp <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.temp) <- vparameters$state.name[1:16]						
						psi.temp[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.temp <- psi.temp * t(int.test.eff[(i - int.start), ,])
						psi.temp[ ,c("I2")] <- psi.temp[ , c("I2")] + v2
						psi.temp[ ,c("I3")] <- psi.temp[ , c("I3")] + v3
						psi.temp <- t(psi.temp)
						
						psi.full <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.full) <- vparameters$state.name[1:16]						
						psi.full[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.full <- psi.full * t(int.test.eff[(i - int.start), ,])
						psi.full[ ,c("I2")] <- psi.full[ , c("I2")] + v2
						psi.full[ ,c("I3")] <- psi.full[ , c("I3")] + v3
						psi.full <- t(psi.full)
						
						# Incremental testing costs
						testing.costs.mod <- ((CurrentPeriod * (psi.temp - psi) * StateCostsTests.psi) 
																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
						int.costs.mod <- ((CurrentPeriod * psi.temp * IntCosts.Tests) 
															/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if(i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- (((CurrentPeriod * psi.full * IntCosts.impl.Tests) + IntCosts.impl.lump.Tests)/
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- (((CurrentPeriod * psi.temp * IntCosts.sust.Tests) + IntCosts.sust.lump.Tests)/ 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- (int.costs.mod + testing.costs.mod) + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
					}
					if(is.element(c("Nurse-initiated testing"), current.int)){
						IntCosts.Tests <- IntCosts.test.nurse
						IntCosts.impl.Tests <- IntCosts.impl.test.nurse
						IntCosts.sust.Tests <- IntCosts.sust.test.nurse
						IntCosts.impl.lump.Tests <- IntCosts.impl.lump.test.nurse
						IntCosts.sust.lump.Tests <- IntCosts.sust.lump.test.nurse
						
						temp.int     <- IntEff.test.nurse
						scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
															rep(temp.int, (int.end - int.start - scale.up.period)))
						scale.m      <- IntScale.test.nurse
						reweight     <- sum(inits[,c("Ia","I1","I2","I3")])/sum(inits[,c("S1")])*((1-IntPosRate.test.nurse)/IntPosRate.test.nurse)
						scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.nurse) * 1/(1 + reweight) / 4,
																												nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
						scale.m["S1", ]                   <-        colSums(IntScale.test.nurse) * reweight/(1 + reweight)      
						scale.m.temp <- replicate((int.sus * 12), scale.m)
						scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
						int.test.eff <- scale.m.temp * (scale.up.v - 1) + 1
#browser()	#Compare 		int.pop.eff$int.test.eff			
						psi.temp <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.temp) <- vparameters$state.name[1:16]						
						psi.temp[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.temp <- psi.temp * t(int.test.eff[(i - int.start), ,])
						psi.temp[ ,c("I2")] <- psi.temp[ , c("I2")] + v2
						psi.temp[ ,c("I3")] <- psi.temp[ , c("I3")] + v3
						psi.temp <- t(psi.temp)
						
						psi.full <- matrix(0, nrow = 30, ncol = 16)
						colnames(psi.full) <- vparameters$state.name[1:16]						
						psi.full[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t
						psi.full <- psi.full * t(int.test.eff[(i - int.start), ,])
						psi.full[ ,c("I2")] <- psi.full[ , c("I2")] + v2
						psi.full[ ,c("I3")] <- psi.full[ , c("I3")] + v3
						psi.full <- t(psi.full)
						
						# Incremental testing costs
						testing.costs.mod <- ((CurrentPeriod * (psi.temp - psi) * StateCostsTests.psi) 
																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
						int.costs.mod <- ((CurrentPeriod * psi.temp * IntCosts.Tests) 
															/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if(i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- (((CurrentPeriod * psi.full * IntCosts.impl.Tests) + IntCosts.impl.lump.Tests)/
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- (((CurrentPeriod * psi.temp * IntCosts.sust.Tests) + IntCosts.sust.lump.Tests)/ 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- (int.costs.mod + testing.costs.mod) + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
					}
					if(is.element(c("OAT integrated testing"), current.int)){
						IntCosts.Tests <- IntCosts.test.OAT
						IntCosts.impl.Tests <- IntCosts.impl.test.OAT
						IntCosts.sust.Tests <- IntCosts.sust.test.OAT
						
						temp.int     <- IntEff.test.OAT
						scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
															rep(temp.int, (int.end - int.start - scale.up.period)))
						#      scale.v      <- IntScale.test.OAT["S1", ]
						scale.m      <- IntScale.test.OAT
						reweight     <- sum(inits[,c("Ia","I1","I2","I3")])/sum(inits[,c("S1")])*((1-IntPosRate.test.OAT)/IntPosRate.test.OAT)
						scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.OAT) * 1/(1 + reweight) / 4,
																												nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
						scale.m["S1", ]                   <-        colSums(IntScale.test.OAT) * reweight/(1 + reweight)     
						#      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
						scale.m.temp <- replicate((int.sus * 12), scale.m)
						scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
						
						int.test.eff <- scale.m.temp * (scale.up.v - 1) + 1
#browser()	#Compare 		int.pop.eff$int.test.eff
            # int.test.eff only apply to prop.oat
#						psi.matrix.temp <- psi.matrix
#						psi.matrix.temp["S1", all.pwid] <- (psi[all.pwid] * int.test.eff[(i - int.start), all.pwid])     * prop.oat.12 +  psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.temp["Ia", all.pwid] <- (psi[all.pwid] * int.test.eff[(i - int.start), all.pwid])     * prop.oat.12 +  psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.temp["I1", all.pwid] <- (psi[all.pwid] * int.test.eff[(i - int.start), all.pwid])     * prop.oat.12 +  psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.temp["I2", all.pwid] <- (psi[all.pwid] * int.test.eff[(i - int.start), all.pwid] + v2)* prop.oat.12 + (psi[all.pwid] + v2)* (1-prop.oat.12)
#						psi.matrix.temp["I3", all.pwid] <- (psi[all.pwid] * int.test.eff[(i - int.start), all.pwid] + v3)* prop.oat.12 + (psi[all.pwid] + v3)* (1-prop.oat.12)
						psi.temp <- psi
						psi.temp["S1", all.pwid] <- (psi.t[all.pwid] * int.test.eff[(i - int.start),"S1",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.temp["Ia", all.pwid] <- (psi.t[all.pwid] * int.test.eff[(i - int.start),"Ia",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.temp["I1", all.pwid] <- (psi.t[all.pwid] * int.test.eff[(i - int.start),"I1",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.temp["I2", all.pwid] <- (psi.t[all.pwid] * int.test.eff[(i - int.start),"I2",all.pwid] + v2)* prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] + 
																																											(psi.t[all.pwid] + v2)* (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.temp["I3", all.pwid] <- (psi.t[all.pwid] * int.test.eff[(i - int.start),"I3",all.pwid] + v3)* prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] + 
																																											(psi.t[all.pwid] + v3)* (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						
#						psi.matrix.full <- psi.matrix
#						psi.matrix.full["S1", all.pwid] <- (psi[all.pwid] * int.test.eff[scale.up.period, all.pwid])     * prop.oat.12 + psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.full["Ia", all.pwid] <- (psi[all.pwid] * int.test.eff[scale.up.period, all.pwid])     * prop.oat.12 + psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.full["I1", all.pwid] <- (psi[all.pwid] * int.test.eff[scale.up.period, all.pwid])     * prop.oat.12 + psi[all.pwid]     * (1-prop.oat.12)
#						psi.matrix.full["I2", all.pwid] <- (psi[all.pwid] * int.test.eff[scale.up.period, all.pwid] + v2)* prop.oat.12 + (psi[all.pwid] + v2)* (1-prop.oat.12)
#						psi.matrix.full["I3", all.pwid] <- (psi[all.pwid] * int.test.eff[scale.up.period, all.pwid] + v3)* prop.oat.12 + (psi[all.pwid] + v3)* (1-prop.oat.12)
						psi.full <- psi
						psi.full["S1", all.pwid] <- (psi.t[all.pwid] * int.test.eff[scale.up.period,"S1",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.full["Ia", all.pwid] <- (psi.t[all.pwid] * int.test.eff[scale.up.period,"Ia",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.full["I1", all.pwid] <- (psi.t[all.pwid] * int.test.eff[scale.up.period,"I1",all.pwid])     * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] +  
																																												psi.t[all.pwid]     * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.full["I2", all.pwid] <- (psi.t[all.pwid] * int.test.eff[scale.up.period,"I2",all.pwid] + v2)* prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] + 
																																											(psi.t[all.pwid] + v2)* (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
						psi.full["I3", all.pwid] <- (psi.t[all.pwid] * int.test.eff[scale.up.period,"I3",all.pwid] + v3)* prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] + 
																																											(psi.t[all.pwid] + v3)* (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid])
#browser()	#Change due to OAT strata removal					
						# Incremental testing costs
#						testing.costs.mod <- ((CurrentPeriod * (psi.matrix.temp - psi.matrix) * StateCostsTests.psi) 
#																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						testing.costs.mod <- ((CurrentPeriod * (psi.temp - psi) * StateCostsTests.psi) 
																	/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
#						int.costs.mod <- ((CurrentPeriod * psi.matrix.temp * IntCosts.Tests) 
#															/ ((1 + Discounting / 12) ^ (i - int.start)))
						int.costs.mod <- ((CurrentPeriod * psi.temp * IntCosts.Tests) 
															/ ((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if(i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- ((CurrentPeriod * psi.full * IntCosts.impl.Tests) /
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- ((CurrentPeriod * psi.temp * IntCosts.sust.Tests) / 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- (int.costs.mod + testing.costs.mod) + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
					}
				} # END OF TESTING LOOP
			
			
			############################ 
			# ART ENGAGEMENT COST
			## INCLUDES ART-INITIATION RELATED COSTS, i.e. CD4 & vL tests, FOR EVERY NEW ART EPISODE
			## Same as STATE COSTS above. 
			## Placed here as psi is time-varying
			art.ini.matrix <- matrix(0, nrow = 30, ncol = 16)
			colnames(art.ini.matrix) <- state.name[c(1:16)]
			phi   = matrix(0, length(names.gp), 3)
			alpha = matrix(0, length(names.gp), 3)
			alpha.re = numeric(length(names.gp))
			theta.t1O = numeric(length(names.gp))
			theta.t2O = numeric(length(names.gp))
			theta.t3O = numeric(length(names.gp))
			
			for (q in 1:18){
				ind = which(rname %in% names18[q])
				phi[ind, ]    = c(phi1[q],   phi2[q],   phi3[q])
				alpha[ind, ]  = c(alpha1[q], alpha2[q], alpha3[q])
				alpha.re[ind] = O_T[q]
				theta.t1O[ind]= T1_O[q]  
				theta.t2O[ind]= T2_O[q] 
				theta.t3O[ind] = T3_O[q]
			}
			
			if (tf.ls$intervention.testing.modif==TRUE & i > int.start & i <= int.end){
#				art.ini.matrix[ , "I1"] <- psi.matrix.temp["I1", ] * phi[ , 1] 
#				art.ini.matrix[ , "I2"] <- psi.matrix.temp["I2", ] * phi[ , 2]
#				art.ini.matrix[ , "I3"] <- psi.matrix.temp["I3", ] * phi[ , 3]
				art.ini.matrix[ , "I1"] <- psi.temp["I1", ] * phi[ , 1] 
				art.ini.matrix[ , "I2"] <- psi.temp["I2", ] * phi[ , 2]
				art.ini.matrix[ , "I3"] <- psi.temp["I3", ] * phi[ , 3]
			} else {
#				art.ini.matrix[ , "I1"] <- psi.matrix["I1", ] * phi[ , 1]
#				art.ini.matrix[ , "I2"] <- psi.matrix["I2", ] * phi[ , 2]
#				art.ini.matrix[ , "I3"] <- psi.matrix["I3", ] * phi[ , 3]
				art.ini.matrix[ , "I1"] <- psi["I1", ] * phi[ , 1]
				art.ini.matrix[ , "I2"] <- psi["I2", ] * phi[ , 2]
				art.ini.matrix[ , "I3"] <- psi["I3", ] * phi[ , 3]
			}
			art.ini.matrix[ , c("D1", "D2", "D3")]   <- alpha[ , 1:3]
			art.ini.matrix[ , c("T1", "T2", "T3")]   <- -cbind(theta.t1O, theta.t2O, theta.t3O) # MINUS?
			art.ini.matrix[ , c("O1", "O2", "O3") ]  <- alpha.re
			art.ini.matrix                           <- t(art.ini.matrix)
			art.ini.matrix[c("T1", "T2", "T3"), all.pwid] <- art.ini.matrix[c("T1", "T2", "T3"), all.pwid] * theta.o.oat * prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid] + 
				                                                             art.ini.matrix[c("T1", "T2", "T3"), all.pwid] * (1-prop.oat.m[nrow(prop.oat.m),gp.gn$all.pwid]) 
			#theta.o.oat only apply to prop.oat.12 proportion
		
			# ART INITIATION COSTS
			art.ini.costs.period[[i]] <- ((CurrentPeriod * art.ini.matrix * StateCostsARTinit) /
																			((1 + Discounting / 12) ^ (i - int.start)))
#browser()		
			## ART INITIATION INTERVENTION
			if (tf.ls$intervention.initiation.modif==TRUE & i > int.start & i <= int.end){
					art.ini.matrix.temp <- art.ini.matrix
					art.ini.matrix.temp["D1", ] = alpha[ , 1] * int.pop.eff$int.art.ini.eff[(i - int.start), ]
					art.ini.matrix.temp["D2", ] = alpha[ , 2] * int.pop.eff$int.art.ini.eff[(i - int.start), ]
					art.ini.matrix.temp["D3", ] = alpha[ , 3] * int.pop.eff$int.art.ini.eff[(i - int.start), ]
#browser()	#Check int.pop.eff$int.art.ini.eff, tf.ls$intervention.initiation.modif==TRUE Default?
					# Incremental ART initiation costs
					art.ini.costs.mod <- ((CurrentPeriod * (art.ini.matrix.temp - art.ini.matrix) * StateCostsARTinit) /
																	((1 + Discounting / 12) ^ (i - int.start))) 
					
					# Intervention costs
					int.costs.mod <- (((CurrentPeriod * prop.link) * IntScale.ARTinit["D1", ] * IntCosts.ARTinit) / 
															((1 + Discounting / 12) ^ (i - int.start)))
#browser()					
					# Implementation & sustainment costs
					if (i <= (scale.up.period + int.start)){
						int.impl.costs.temp <- (((CurrentPeriod * prop.link) * IntScale.ARTinit["D1", ] * IntCosts.impl.ARTinit) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- ((CurrentPeriod * prop.link) * 0)
						# END OF SCALE-UP PERIOD
					} else { 
						int.impl.costs.temp <- ((CurrentPeriod * prop.link) * 0)
						int.sust.costs.temp <- (((CurrentPeriod * prop.link) * IntScale.ARTinit["D1", ] * IntCosts.sust.ARTinit) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
					}
					int.costs.period[[i]] <- int.costs.mod + art.ini.costs.mod + int.costs.period[[i]]
					int.impl.costs.period[[i]] <-  int.impl.costs.temp + int.impl.costs.period[[i]]   
					int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]  
#browser()				
					# END OF SUSTAINMENT PERIOD
			} # END OF ART INITIATION LOOP
			
			## ART RETENTION INTERVENTION
			if (tf.ls$intervention.retention.modif==TRUE & i > int.start & i <= int.end){
					if(is.element(c("ART retention"), current.int)){
						this.scale <- IntScale.ARTreten["T1", ]
					} else if (is.element(c("ART retention, targeted"), current.int)){
						IntCosts.ARTreten = IntCosts.ARTreten.target
						IntCosts.impl.ARTreten = IntCosts.impl.ARTreten.target
						IntCosts.sust.ARTreten = IntCosts.sust.ARTreten.target
						this.scale <- IntScale.ARTreten.target["T3", ]
					}
#browser()	# check int.pop.eff$int.art.ret.eff	#This part seems clearly wrong, cannot stop at this.scale, have to get to int.art.ret.eff before * CurrentPeriod 
					# Intervention costs
					int.costs.mod <- ((CurrentPeriod * this.scale * IntCosts.ARTreten ) / 
															((1 + Discounting / 12) ^ (i - int.start)))
					
					# Implementation & sustainment costs
					if (i <= (scale.up.period + int.start)){
						int.impl.costs.temp <- ((CurrentPeriod * this.scale * IntCosts.impl.ARTreten) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- (CurrentPeriod * 0)
						# END OF SCALE-UP PERIOD
					} else {
						int.impl.costs.temp <- (CurrentPeriod * 0)
						int.sust.costs.temp <- ((CurrentPeriod * this.scale * IntCosts.sust.ARTreten) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
					}
					int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
					int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
					int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
					# END OF SUSTAINMENT PERIOD
			} # END OF ART RETENTION LOOP 
			
			## ART DROPOUT EMR-PROMPT INTERVENTION
			if(tf.ls$intervention.dropout.modif==TRUE & i > int.start & i <= int.end){
#browser()	# check int.pop.eff$int.art.drop.eff - Same here CurrentPeriod * IntScale.ARTEMRprompt["T1", ] does not make sense 
				    # Unless we only look at the additional cost, 
					# Intervention costs, thus we look only look at value above 1, but then should this apply to other interventions too?
					int.costs.mod <- ((CurrentPeriod * IntScale.ARTEMRprompt["T1", ] * IntCosts.ARTEMRprompt ) / 
															((1 + Discounting / 12) ^ (i - int.start)))
					
					# Implementation & sustainment costs
					if (i <= (scale.up.period + int.start)){
						int.impl.costs.temp <- ((CurrentPeriod * IntScale.ARTEMRprompt["T1", ] * IntCosts.impl.ARTEMRprompt) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- (CurrentPeriod * 0)
						# END OF SCALE-UP PERIOD
					} else {
						int.impl.costs.temp <- (CurrentPeriod * 0)
						int.sust.costs.temp <- ((CurrentPeriod * IntScale.ARTEMRprompt["T1", ] * IntCosts.sust.ARTEMRprompt) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
					}
					int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
					int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]  
					int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
					# END OF SUSTAINMENT PERIOD
			} # END OF ART ART DROPOUT EMR-PROMPT INTERVENTION LOOP
			
			## IMMEDIATE ART INTERVENTION
			if(tf.ls$intervention.immediateART.modif==TRUE & i > int.start & i <= int.end){
					art.ini.matrix.temp <- art.ini.matrix
					if(is.element(c("RAPID ART"), current.int)){
						
#						art.ini.matrix.temp[c("I1", "I2", "I3"), ] = psi.matrix[c("I1", "I2", "I3") , ] * t(phi[ , 1:3]) * int.pop.eff$int.art.imm.eff[(i - int.start), ]
						art.ini.matrix.temp[c("I1", "I2", "I3"), ] = psi[c("I1", "I2", "I3") , ] * t(phi[ , 1:3]) * int.pop.eff$int.art.imm.eff[(i - int.start), ]
					} else if (is.element(c("RAPID ART, targeted"), current.int)){
#						art.ini.matrix.temp["I3", ] = psi.matrix["I3", ] * phi[ , 3] * int.pop.eff$int.art.imm.eff[(i - int.start), ]
						art.ini.matrix.temp["I3", ] = psi["I3", ] * phi[ , 3] * int.pop.eff$int.art.imm.eff[(i - int.start), ]
browser()	# check int.pop.eff$int.art.imm.eff					
						IntCosts.immART = IntCosts.immART.target
						IntCosts.impl.immART = IntCosts.impl.immART.target
						IntCosts.sust.immART = IntCosts.sust.immART.target
					}
					
					# Incremental ART initiation costs
					art.ini.costs.mod <- ((CurrentPeriod * (art.ini.matrix.temp - art.ini.matrix) * StateCostsARTinit) /
																	((1 + Discounting / 12) ^ (i - int.start))) 
					
					# Intervention costs
#					int.costs.mod <- (((CurrentPeriod * psi.matrix) * IntScale.immART["I1", ] * IntCosts.immART ) / 
#															((1 + Discounting / 12) ^ (i - int.start)))
					int.costs.mod <- (((CurrentPeriod * psi) * IntScale.immART["I1", ] * IntCosts.immART ) / 
															((1 + Discounting / 12) ^ (i - int.start)))
					
					# Implementation & sustainment costs
					if (i <= (scale.up.period + int.start)){
						
#						int.impl.costs.temp <- (((CurrentPeriod * psi.matrix) * IntScale.immART["I1", ] * IntCosts.impl.immART) / 
#																			((1 + Discounting / 12) ^ (i - int.start)))
#						int.sust.costs.temp <- ((CurrentPeriod * psi.matrix) * 0)
						int.impl.costs.temp <- (((CurrentPeriod * psi) * IntScale.immART["I1", ] * IntCosts.impl.immART) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- ((CurrentPeriod * psi) * 0)
						# END OF SCALE-UP PERIOD
					} else {
#						int.impl.costs.temp <- ((CurrentPeriod * psi.matrix) * 0)
#						int.sust.costs.temp <- (((CurrentPeriod * psi.matrix) * IntScale.immART["I1", ] * IntCosts.sust.immART) / 
#																			((1 + Discounting / 12) ^ (i - int.start)))
						int.impl.costs.temp <- ((CurrentPeriod * psi) * 0)
						int.sust.costs.temp <- (((CurrentPeriod * psi) * IntScale.immART["I1", ] * IntCosts.sust.immART) / 
																			((1 + Discounting / 12) ^ (i - int.start)))
					}
					int.costs.period[[i]] <- int.costs.mod + art.ini.costs.mod + int.costs.period[[i]]
					int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
					int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
					# END OF SUSTAINMENT PERIOD
			} # END OF IMMEDIATE ART INTERVENTION LOOP
			
			## ART RE-INITIATION INTERVENTION
			if(tf.ls$intervention.reinitiation.modif==TRUE & i > int.start & i <= int.end){
					if(is.element(c("ART re-initiation"), current.int)){
						temp.int <- IntEff.reART
						scale.v <- IntScale.reART["O1", ]
						scale.m.temp <- t(replicate((int.sus * 12), scale.v))
						scale.up.v <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
														rep(temp.int, (int.end - int.start - scale.up.period)))
						int.art.reini.eff <- scale.m.temp * (scale.up.v - 1) + 1
						
						art.ini.matrix.temp <- art.ini.matrix
						art.ini.matrix.temp["O1", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
						art.ini.matrix.temp["O2", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
						art.ini.matrix.temp["O3", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
#browser()	# check int.pop.eff$int.art.reini.eff							
						# Incremental ART initiation costs
						art.ini.costs.mod <- ((CurrentPeriod * (art.ini.matrix.temp - art.ini.matrix) * StateCostsARTinit) /
																		((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
						int.costs.mod <- ((CurrentPeriod * IntScale.reART["O1", ] * IntCosts.reART ) / 
																((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if (i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- ((CurrentPeriod * IntScale.reART["O1", ] * IntCosts.impl.reART) / 
																				((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- ((CurrentPeriod * IntScale.reART["O1", ] * IntCosts.sust.reART) / 
																				((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- int.costs.mod + art.ini.costs.mod + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]    
						# END OF RE-ART LOOP
						
					} else if(is.element(c("ART re-linkage"), current.int)){
						temp.int <- IntEff.relink
						scale.v <- IntScale.relink["O1", ]
						scale.m.temp <- t(replicate((int.sus * 12), scale.v))
						scale.up.v <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
														rep(temp.int, (int.end - int.start - scale.up.period)))
						int.art.reini.eff <- scale.m.temp * (scale.up.v - 1) + 1
						
						art.ini.matrix.temp <- art.ini.matrix
						art.ini.matrix.temp["O1", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
						art.ini.matrix.temp["O2", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
						art.ini.matrix.temp["O3", ] = alpha.re * int.art.reini.eff[(i - int.start), ]
#browser()	# check int.pop.eff$int.art.reini.eff						
						# Incremental ART initiation costs
						art.ini.costs.mod <- ((CurrentPeriod * (art.ini.matrix.temp - art.ini.matrix) * StateCostsARTinit) /
																		((1 + Discounting / 12) ^ (i - int.start)))
						
						# Intervention costs
						int.costs.mod <- ((CurrentPeriod * IntScale.relink["O1", ] * IntCosts.relink ) / 
																((1 + Discounting / 12) ^ (i - int.start)))
						
						# Implementation & sustainment costs
						if (i <= (scale.up.period + int.start)){
							int.impl.costs.temp <- ((IntCosts.impl.lump.relink) / ((1 + Discounting / 12) ^ (i - int.start)))
							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
						} else {
							int.impl.costs.temp <- (CurrentPeriod * 0)
							int.sust.costs.temp <- ((IntCosts.sust.relink) / ((1 + Discounting / 12) ^ (i - int.start)))
						}
						int.costs.period[[i]] <- int.costs.mod + art.ini.costs.mod + int.costs.period[[i]]
						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF RE-LINKAGE LOOP
				} # END OF SUSTAINMENT PERIOD
			} # END OF ART RE-ENGAGEMENT LOOP
			
			############################   
			# OAT INTERVENTION
			if(tf.ls$intervention.oat.modif==TRUE & i > int.start & i <= int.end){
#					pop.pwid = pop.pwid
#					nOAT.m   = nOAT.m
#					nOAT     = round(c(nOAT.m[nrow(nOAT.m), 1:3] * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[nrow(nOAT.m), 4:6]))
#					if(is.element(c("OAT with BUP"), current.int)){
#						IntCosts.OAT = IntCosts.BUP
#						IntCosts.impl.OAT = IntCosts.impl.BUP
#						IntCosts.sust.OAT = IntCosts.sust.BUP
#browser() #Check while this int impact costs.art.ini.sum						
#						scale.growth <- IntScale.BUP
#						scale.v <- IntBaseline.BUP
#						scale.m.temp <- t(replicate((int.sus * 12), scale.v))
#						scale.up.v <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
#														rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))  # Growth rate on OAT due to intervention (1 is SQ)
#						int.oat.eff <- scale.m.temp * (scale.up.v - 1) + t(replicate(int.sus * 12, nOAT.m[nrow(nOAT.m),]))   # Extrapolate value from nOAT.m due to interventions
						
#						nOAT.int = round(c(as.vector(t(int.oat.eff[(i - int.start), 1:3])) * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), as.vector(t(int.oat.eff[(i - int.start), 4:6]))))
#browser()						
						# revisit
#            prop.oat     =    nOAT /pop.pwid
#						prop.oat.int = nOAT.int/pop.pwid
#						if(any(prop.oat.int > 0.95 * 0.727)) {
#							while (any(prop.oat.int > 0.95 * 0.727 + 0.00001)){
#								sat            <- which(prop.oat.int >= 0.95 * 0.727)
#								nOAT.int[-sat] <- nOAT.int[-sat] + sum(nOAT.int[sat] - pop.pwid[sat] * 0.95 * 0.727) * (pop.pwid[-sat] / sum(pop.pwid[-sat]))
#								nOAT.int[sat]  <- pop.pwid[sat] * 0.95 * 0.727
#								prop.oat.int   <- nOAT.int/pop.pwid
#							}
#						}
#browser()						
#						oat.e     = prop.oat     / (1 - prop.oat    ) * oat.q 
#						oat.e.int = prop.oat.int / (1 - prop.oat.int) * oat.q 
						## to address entry rate > 1 or <0
#						oat.e[which(oat.e > 0.9 | oat.e <= 0)] <- 0.9
#						oat.e.int[which(oat.e.int > 0.9 | oat.e.int <= 0)] <- 0.9
						##
#						rname = gsub(paste(c("/low", "/high","/young","/old"), collapse="|"), "", names.gp)
#						pwid.name = names18[grep("PWID", names18)]
#						oat.e.all.m.int <- oat.e.all.m <- matrix(0, nrow = 19, ncol = 42)
#						for (q in 1:length(pwid.name)){
#							ind = which(rname %in% pwid.name[q])
#							oat.e.all.m[ , intersect(ind, off.oat)]     = oat.e[q]
#							oat.e.all.m.int[ , intersect(ind, off.oat)] = oat.e.int[q]
#						}
						## Replace oat.e with prop.oat						
#						prop.oat.all.m.int <- prop.oat.all.m <- matrix(0, nrow = 16, ncol = 30)
#						for (q in 1:length(pwid.name)){
#							ind = which(rname %in% pwid.name[q])
#							prop.oat.all.m[ , intersect(ind, pwid)]     = prop.oat[q]
#							prop.oat.all.m.int[ , intersect(ind, pwid)] = prop.oat.int[q]
#						}						
#						int.costs.mod <- ((CurrentPeriod  * (oat.e.all.m.int - oat.e.all.m) * IntCosts.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#						int.costs.mod <- ((CurrentPeriod  * (prop.oat.all.m.int - prop.oat.all.m) * IntCosts.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#browser()						
						# Implementation & sustainment costs
#						if (i <= (scale.up.period + int.start)){
#							nOAT.int.full <- round(c(as.vector(t(int.oat.eff[scale.up.period, 1:3])) * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), as.vector(t(int.oat.eff[scale.up.period, 4:6]))))
							
#							prop.oat.full = nOAT.int.full/pop.pwid
#							if(any(prop.oat.full > 0.95 * 0.727)) {
#								while (any(prop.oat.full > 0.95 * 0.727 + 0.00001)){
#									sat <- which(prop.oat.full >= 0.95 * 0.727)
#									nOAT.int.full[-sat] <- nOAT.int.full[-sat] + sum(nOAT.int.full[sat] - pop.pwid[sat] * 0.95 * 0.727) * (pop.pwid[-sat] / sum(pop.pwid[-sat]))
#									nOAT.int.full[sat] <- pop.pwid[sat] * 0.95 * 0.727
#									prop.oat.full = nOAT.int.full/pop.pwid
#								}
#							}
#browser()							
#							oat.e.int.full = prop.oat.full / (1 - prop.oat.full) * oat.q
							## to address entry rate > 1 or <0
#							oat.e.int.full[which(oat.e.int.full > 0.9 | oat.e.int.full <= 0)] <- 0.9
							##
#							oat.e.all.m.int.full = oat.e.all.m
#							for (q in 1:length(pwid.name)){
#								ind = which(rname %in% pwid.name[q])
#								oat.e.all.m.int.full[ , intersect(ind, off.oat)] = oat.e.int.full[q]
#							}
#             prop.oat.all.m.full <- prop.oat.all.m 
#						  for (q in 1:length(pwid.name)){
#								ind = which(rname %in% pwid.name[q])
#								prop.oat.all.m.full[ , intersect(ind, pwid)] = prop.oat.full[q]
#							}
#             int.impl.costs.temp <- ((CurrentPeriod * oat.e.all.m.int.full * IntCosts.impl.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))							
#              int.impl.costs.temp <- ((CurrentPeriod * prop.oat.all.m.full * IntCosts.impl.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
#						} else {
#							int.impl.costs.temp <- (CurrentPeriod * 0)
#							int.sust.costs.temp <- ((CurrentPeriod * oat.e.all.m.int * IntCosts.sust.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))							
#							int.sust.costs.temp <- ((CurrentPeriod * prop.oat.all.m.int  * IntCosts.sust.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#						}
#						int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
#						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]    
#						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
#browser() # Done with BUP intervention 						
#					}
#					if(is.element(c("OAT with methadone"), current.int)){
#						IntCosts.OAT = IntCosts.MET
#						IntCosts.impl.OAT = IntCosts.impl.MET
#						IntCosts.sust.OAT = IntCosts.sust.MET
						
#						scale.growth <- IntScale.MET
#						scale.v <- IntBaseline.MET
#						scale.m.temp <- t(replicate((int.sus * 12), scale.v))
#						scale.up.v <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
#														rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))  # Growth rate on OAT due to intervention (1 is SQ)
#						int.oat.eff <- scale.m.temp * (scale.up.v - 1) + t(replicate(int.sus * 12, nOAT.m[nrow(nOAT.m),]))   # Extrapolate value from nOAT.m due to interventions

#						nOAT.int = round(c(as.vector(t(int.oat.eff[(i - int.start), 1:3])) * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), as.vector(t(int.oat.eff[(i - int.start), 4:6]))))
						
#						prop.oat     =    nOAT /pop.pwid
#						prop.oat.int = nOAT.int/pop.pwid
#						if(any(prop.oat.int > 0.95 * 0.727)) {
#							while (any(prop.oat.int > 0.95 * 0.727 + 0.00001)){
#								sat <- which(prop.oat.int >= 0.95 * 0.727)
#								nOAT.int[-sat] <- nOAT.int[-sat] + sum(nOAT.int[sat] - pop.pwid[sat] * 0.95 * 0.727) * (pop.pwid[-sat] / sum(pop.pwid[-sat]))
#								nOAT.int[sat] <- pop.pwid[sat] * 0.95 * 0.727
#								prop.oat.int = nOAT.int/pop.pwid
#							}
#						}
#browser()						
#						oat.e     = (nOAT / pop.pwid) / (1 - (nOAT / pop.pwid)) * oat.q 
#						oat.e.int = prop.oat.int / (1 - prop.oat.int) * oat.q 
						## to address entry rate > 1 or <0
#						oat.e[which(oat.e > 0.9 | oat.e <= 0)] <- 0.9
#						oat.e.int[which(oat.e.int > 0.9 | oat.e.int <= 0)] <- 0.9
						##
#						rname = gsub(paste(c("/low", "/high", "/young", "/old"), collapse="|"), "", names.gp)
#						pwid.name = names18[grep("PWID", names18)]
#						oat.e.all.m.int <- oat.e.all.m <- matrix(0, nrow = 19, ncol = 42)
#						for (q in 1:length(pwid.name)){
#							ind = which(rname %in% pwid.name[q])
#							oat.e.all.m[ , intersect(ind, off.oat)]     = oat.e[q]
#							oat.e.all.m.int[ , intersect(ind, off.oat)] = oat.e.int[q]
#						prop.oat.all.m.int <- prop.oat.all.m <- matrix(0, nrow = 16, ncol = 30)
#						for (q in 1:length(pwid.name)){
#							ind = which(rname %in% pwid.name[q])
#							prop.oat.all.m[ , intersect(ind, pwid)]     = prop.oat[q]
#							prop.oat.all.m.int[ , intersect(ind, pwid)] = prop.oat.int[q]
#						}		
						
#						int.costs.mod <- ((CurrentPeriod  * (oat.e.all.m.int - oat.e.all.m) * IntCosts.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#						int.costs.mod <- ((CurrentPeriod  * (prop.oat.all.m.int - prop.oat.all.m) * IntCosts.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#browser()						
						# Implementation & sustainment costs
#						if (i <= (scale.up.period + int.start)){
#							nOAT.int.full <- round(c(as.vector(t(int.oat.eff[scale.up.period, 1:3])) * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), as.vector(t(int.oat.eff[scale.up.period, 4:6]))))
							
#							prop.oat.full = nOAT.int.full/pop.pwid
#							if(any(prop.oat.full > 0.95 * 0.727)) {
#								while (any(prop.oat.full > 0.95 * 0.727 + 0.00001)){
#									sat <- which(prop.oat.full >= 0.95 * 0.727)
#									nOAT.int.full[-sat] <- nOAT.int.full[-sat] + sum(nOAT.int.full[sat] - pop.pwid[sat] * 0.95 * 0.727) * (pop.pwid[-sat] / sum(pop.pwid[-sat]))
#									nOAT.int.full[sat] <- pop.pwid[sat] * 0.95 * 0.727
#									prop.oat.full = nOAT.int.full/pop.pwid
#								}
#							}
#							oat.e.int.full = prop.oat.full / (1 - prop.oat.full) * oat.q
							## to address entry rate > 1 or <0
#							oat.e.int.full[which(oat.e.int.full > 0.9 | oat.e.int.full <= 0)] <- 0.9
							##
#							oat.e.all.m.int.full = oat.e.all.m
#							for (q in 1:length(pwid.name)){
#								ind = which(rname %in% pwid.name[q])
#								oat.e.all.m.int.full[ , intersect(ind, off.oat)] = oat.e.int.full[q]
#							}
#              prop.oat.all.m.full <- prop.oat.all.m 
#              for (q in 1:length(pwid.name)){
#	                ind = which(rname %in% pwid.name[q])
#	                prop.oat.all.m.full[ , intersect(ind, pwid)] = prop.oat.full[q]
#             }
#							int.impl.costs.temp <- ((CurrentPeriod  * oat.e.all.m.int.full * IntCosts.impl.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#              int.impl.costs.temp <- ((CurrentPeriod * prop.oat.all.m.full * IntCosts.impl.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#							int.sust.costs.temp <- (CurrentPeriod * 0)
							# END OF SCALE-UP PERIOD
#						} else {
#							int.impl.costs.temp <- (CurrentPeriod * 0)
#							int.sust.costs.temp <- ((CurrentPeriod * oat.e.all.m.int * IntCosts.sust.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#							int.sust.costs.temp <- ((CurrentPeriod * prop.oat.all.m.int  * IntCosts.sust.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
#						}
#						int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
#						int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]    
#						int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
						# END OF SUSTAINMENT PERIOD
#					}
#			} # END OF OAT LOOP
			prop.oat.matrix = prop.oat.matrix.temp = prop.oat.matrix.full = matrix(0, nrow=16, ncol=30)
			prop.oat.temp <- int.pop.eff$int.oat.eff[i - int.start,]
			
			prop.oat.matrix[ , ]       =  matrix(IntBaseline.OAT/16, nrow = 16, ncol = 30, byrow = TRUE)
			prop.oat.matrix.temp[ , ]  = rep(prop.oat.temp/16 , each=16)
			
			int.costs.mod <- ((CurrentPeriod * (prop.oat.matrix.temp - prop.oat.matrix) * IntCosts.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
			
			if (i <= (scale.up.period + int.start)){
				prop.oat.full <- int.pop.eff$int.oat.eff[scale.up.period,]
				prop.oat.matrix.full[ , ] <- rep(prop.oat.full/16, each = 16)
				
				int.impl.costs.temp <- ((CurrentPeriod * prop.oat.matrix.full * IntCosts.impl.OAT) / 
																	((1 + Discounting / 12) ^ (i - int.start)))
				int.sust.costs.temp <- (CurrentPeriod * 0)
				
				# END OF SCALE-UP PERIOD
			} else {
				int.impl.costs.temp <- (CurrentPeriod * 0)
				int.sust.costs.temp <- ((CurrentPeriod * prop.oat.matrix.temp * IntCosts.sust.OAT) / ((1 + Discounting / 12) ^ (i - int.start)))
			}
			
			int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
			int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
			int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
			# END OF SUSTAINMENT PERIOD
	 # END OF OAT LOOP
			}
			
			############################   
			# SSP INTERVENTION
			if(tf.ls$intervention.ssp.modif==TRUE & i > int.start & i <= int.end){
#					v.ssp.int <- matrix(0, nrow=16, ncol=30)
					prop.ssp.matrix    = prop.ssp.matrix.temp = prop.ssp.matrix.full =  matrix(0, nrow=16, ncol=30)
#					v.ssp.int[ , all.pwid] <- (int.pop.eff$int.ssp.eff[i - int.start] - IntBaseline.SSP) / 12 / sum(CurrentPeriod[ , all.pwid]) * CurrentPeriod[ , all.pwid]
					prop.ssp.temp <- int.pop.eff$int.ssp.eff[i - int.start,]
					prop.ssp.matrix[ , ]       = rep(IntBaseline.SSP/16, each=16) #colSums(prop.ssp.matrix.temp) = prop.ssp.temp
					prop.ssp.matrix.temp[ , ]  = rep(prop.ssp.temp  /16 , each=16)  # All 16 states can be on SSP, index for all.pwid no longer required
					
#					int.costs.mod <- ((v.ssp.int * IntCosts.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
					int.costs.mod <- ((CurrentPeriod * (prop.ssp.matrix.temp - prop.ssp.matrix) * IntCosts.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
					
					# Implementation & sustainment costs
					if (i <= (scale.up.period + int.start)){
#						v.ssp.int.full <- matrix(0, nrow=16, ncol=30)
#						v.ssp.int.full[ , all.pwid] <- (int.pop.eff$int.ssp.eff[scale.up.period] - IntBaseline.SSP) / 12 / sum(CurrentPeriod[ , all.pwid]) * CurrentPeriod[ , all.pwid]
						prop.ssp.full            <- int.pop.eff$int.ssp.eff[scale.up.period,] - IntBaseline.SSP
						prop.ssp.matrix.full[ , ]<- rep(prop.ssp.full/16, each=16)
						
#						int.impl.costs.temp <- ((v.ssp.int.full * IntCosts.impl.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
						int.impl.costs.temp <- ((CurrentPeriod * prop.ssp.matrix.full * IntCosts.impl.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- (CurrentPeriod * 0)
						# END OF SCALE-UP PERIOD
					} else {
						int.impl.costs.temp <- (CurrentPeriod * 0)
#						int.sust.costs.temp <- ((v.ssp.int * IntCosts.sust.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
						int.sust.costs.temp <- ((CurrentPeriod * prop.ssp.matrix.temp * IntCosts.sust.SSP) / ((1 + Discounting / 12) ^ (i - int.start)))
					}
					int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
					int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
					int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
				  # END OF SUSTAINMENT PERIOD
			} # END OF SSP LOOP
			
			############################   
			# PrEP INTERVENTION
			if(tf.ls$intervention.prep.modif==TRUE & i > int.start & i <= int.end){
				#				  eta.matrix      = eta.matrix.temp = matrix(0, nrow=16, ncol=30)
				#					prep.proportion <- prep.proportion
				#					msm.h.scep      <- msm.h.scep
				#					prep.vol  <- int.pop.eff$int.prep.eff[i - int.start] * prep.proportion
				#					prop.prep <- int.pop.eff$int.prep.eff[i - int.start] * prep.proportion / msm.h.scep
				prop.prep.matrix = prop.prep.matrix.temp = prop.prep.matrix.full = matrix(0, nrow=16, ncol=30)
				prop.prep.temp <- int.pop.eff$int.prep.eff[i - int.start,]
				
				#					while (any(prop.prep > 0.92001)){
				#						sat            <- which(prop.prep >= 0.92)
				#						prep.vol[-sat] <- prep.vol[-sat] + sum(prep.vol[sat] - msm.h.scep[sat] * 0.92) * (msm.h.scep[-sat] / sum(msm.h.scep[-sat]))
				#						prep.vol[sat]  <- msm.h.scep[sat] * 0.92
				#						prop.prep      <- prep.vol/msm.h.scep
				#					}
				
				#          eta.matrix[c(1:2) , Y.allmsm] = rep(-log(1 - IntBaseline.prep * prep.proportion / msm.h.scep)/12 ,each=2)
				#          eta.matrix.temp[c(1:2) , Y.allmsm] = rep(prop.prep/(1-prop.prep) * wp, each=2)
				
				prop.prep.matrix[c(1:2) , ]       = rep(IntBaseline.prep/2, each=2) # Divided by 2 so that colSums(prop.prep.matrix) = IntBaseline.prep
				prop.prep.matrix.temp[c(1:2) , ]  = rep(prop.prep.temp/2 , each=2) # colSums(prop.prep.matrix.temp) = prop.prep.temp; prop.prep.temp, IntBaseline.prep are all those on PrEP divided among 2 S states
				
				## to address entry rate > 1 (causing NA)
				#					eta.matrix.temp[which(is.na(eta.matrix.temp))] <- 0.9
				#					eta.matrix.temp[eta.matrix.temp > 0.9] <- 0.9
				##
				#browser() 			
				#					int.costs.mod <- ((CurrentPeriod * (eta.matrix.temp - eta.matrix) * IntCosts.prep) / ((1 + Discounting / 12) ^ (i - int.start)))
				int.costs.mod <- ((CurrentPeriod * (prop.prep.matrix.temp - prop.prep.matrix) * IntCosts.prep) / ((1 + Discounting / 12) ^ (i - int.start)))
				#browser()	#Check int.costs.mod					
				# Implementation & sustainment costs
				if (i <= (scale.up.period + int.start)){
					#						eta.matrix.full <- eta.matrix
					#						prop.prep.full  <- int.pop.eff$int.prep.eff[scale.up.period] * prep.proportion / msm.h.scep
					#						prep.vol.full   <- int.pop.eff$int.prep.eff[scale.up.period] * prep.proportion
					#browser() 	
					#							while (any(prop.prep.full > 0.92001)){
					#							sat                 <- which(prop.prep.full >= 0.92)
					#							prep.vol.full[-sat] <- prep.vol.full[-sat] + sum(prep.vol.full[sat] - msm.h.scep[sat] * 0.92) * (msm.h.scep[-sat] / sum(msm.h.scep[-sat]))
					#							prep.vol.full[sat]  <- msm.h.scep[sat] * 0.92
					#							prop.prep.full      <- prep.vol.full/msm.h.scep
					#						}
					
					#						eta.matrix.full[c(1:2) , Y.allmsm] = rep(prop.prep.full/(1-prop.prep.full) * wp, each = 2)
					
					prop.prep.full                   <- int.pop.eff$int.prep.eff[scale.up.period,]
					prop.prep.matrix.full[c(1:2) , ] <- rep(prop.prep.full/2, each = 2)
					
					#browser()	# Check	eta.matrix.full for use later in 	int.impl.costs.temp
					## to address entry rate > 1 (causing NA)
					#						eta.matrix.full[which(is.na(eta.matrix.full))] <- 0.9
					#						eta.matrix.full[eta.matrix.full > 0.9] <- 0.9
					##
					#						int.impl.costs.temp <- ((CurrentPeriod * eta.matrix.full * IntCosts.impl.prep) / 
					#																			((1 + Discounting / 12) ^ (i - int.start)))
					int.impl.costs.temp <- ((CurrentPeriod * prop.prep.matrix.full * IntCosts.impl.prep) / 
																		((1 + Discounting / 12) ^ (i - int.start)))
					int.sust.costs.temp <- (CurrentPeriod * 0)
					#browser() # Check	int.impl.costs.temp, int.sust.costs.temp (=0 for earlier i period)							
					# END OF SCALE-UP PERIOD
				} else {
					int.impl.costs.temp <- (CurrentPeriod * 0)
					#						int.sust.costs.temp <- ((CurrentPeriod * eta.matrix.temp * IntCosts.sust.prep) / 
					#																			((1 + Discounting / 12) ^ (i - int.start)))
					int.sust.costs.temp <- ((CurrentPeriod * prop.prep.matrix.temp * IntCosts.sust.prep) / 
																		((1 + Discounting / 12) ^ (i - int.start)))
				}
				int.costs.period[[i]] <- int.costs.mod + int.costs.period[[i]]
				int.impl.costs.period[[i]] <- int.impl.costs.temp + int.impl.costs.period[[i]]
				int.sust.costs.period[[i]] <- int.sust.costs.temp + int.sust.costs.period[[i]]
				# END OF SUSTAINMENT PERIOD
				#browser(expr = {i > (scale.up.period + int.start)}) #For PrEP, int.sust.costs.temp = 0 due to IntCosts.sust.prep =0
			} # END OF PrEP LOOP
				
			
			############################
		} else { # All outcomes = 0 outside of evaluation period
			# QALYs
			QALYs.period[[i]]            <- (CurrentPeriod * 0)
			# STATE COSTS
			hru.costs.period[[i]]        <- (CurrentPeriod * 0)
			art.costs.period[[i]]        <- (CurrentPeriod * 0)
			art.ini.costs.period[[i]]    <- (CurrentPeriod * 0)
			oat.costs.period[[i]]        <- (CurrentPeriod * 0)
			tests.costs.period[[i]]      <- (CurrentPeriod * 0)
			prep.costs.period[[i]]       <- (CurrentPeriod * 0)
			prep.tests.costs.period[[i]] <- (CurrentPeriod * 0)
			# INTERVENTION COSTS
			int.costs.period[[i]]        <- (CurrentPeriod * 0)
			int.impl.costs.period[[i]]   <- (CurrentPeriod * 0)
			int.sust.costs.period[[i]]   <- (CurrentPeriod * 0)
		}
		
		if(tf.ls$intervention.testing.modif      == FALSE &
			 tf.ls$intervention.initiation.modif   == FALSE &
			 tf.ls$intervention.retention.modif    == FALSE &
			 tf.ls$intervention.dropout.modif      == FALSE &
			 tf.ls$intervention.immediateART.modif == FALSE &
			 tf.ls$intervention.reinitiation.modif == FALSE &
			 tf.ls$intervention.ssp.modif          == FALSE &
			 tf.ls$intervention.oat.modif          == FALSE &
			 tf.ls$intervention.prep.modif         == FALSE){
			# INTERVENTION COSTS
			int.costs.period[[i]]      <- (CurrentPeriod * 0)
			int.impl.costs.period[[i]] <- (CurrentPeriod * 0)
			int.sust.costs.period[[i]] <- (CurrentPeriod * 0)
		}
	}
	#### END OF LOOP OVER PERIODS
	
	if(sums.only==FALSE){
		
		QALYs            <- abind(QALYs.period, along = 0)
		# Stratified costs
		costs.hru        <- abind(hru.costs.period, along = 0)
		costs.art        <- abind(art.costs.period, along = 0)
		costs.art.ini    <- abind(art.ini.costs.period, along = 0)
		costs.oat        <- abind(oat.costs.period, along = 0)
		costs.prep       <- abind(prep.costs.period, along = 0)
		costs.prep.tests <- abind(prep.tests.costs.period, along = 0)
		costs.test       <- abind(tests.costs.period, along = 0)
		# Stratified intervention costs
		int.costs        <- abind(int.costs.period, along = 0)
		int.impl.costs   <- abind(int.impl.costs.period, along = 0)
		int.sust.costs   <- abind(int.sust.costs.period, along = 0)
		
		# Sums
		QALYs.sum            <- sum(unlist(QALYs.period))
		costs.hru.sum        <- sum(unlist(hru.costs.period))
		costs.art.sum        <- sum(unlist(art.costs.period))
		costs.art.ini.sum    <- sum(unlist(art.ini.costs.period))
		costs.oat.sum        <- sum(unlist(oat.costs.period))
		costs.prep.sum       <- sum(unlist(prep.costs.period))
		costs.prep.tests.sum <- sum(unlist(prep.tests.costs.period))
		costs.test.sum       <- sum(unlist(tests.costs.period))
		int.costs.sum        <- sum(unlist(int.costs.period))
		int.impl.costs.sum   <- sum(unlist(int.impl.costs.period))
		int.sust.costs.sum   <- sum(unlist(int.sust.costs.period))
		
		# Total costs
		costs.total.sum <- (costs.hru.sum + costs.art.sum + costs.art.ini.sum +
													costs.oat.sum + costs.prep.sum + costs.prep.tests.sum +
													costs.test.sum + int.costs.sum +
													int.impl.costs.sum + int.sust.costs.sum)
		
		# In list
		out.sum <- list(QALYs                = QALYs, 
										costs.hru            = costs.hru,
										costs.art            = costs.art,
										costs.art.ini        = costs.art.ini,
										costs.oat            = costs.oat,
										costs.prep           = costs.prep,
										costs.prep.tests     = costs.prep.tests,
										costs.test           = costs.test,
										int.costs            = int.costs,
										int.impl.costs       = int.impl.costs,
										int.sust.costs       = int.sust.costs,
										QALYs.sum            = QALYs.sum, 
										costs.total.sum      = costs.total.sum,
										costs.hru.sum        = costs.hru.sum,
										costs.art.sum        = costs.art.sum,
										costs.art.ini.sum    = costs.art.ini.sum,
										costs.oat.sum        = costs.oat.sum,
										costs.prep.sum       = costs.prep.sum,
										costs.prep.tests.sum = costs.prep.tests.sum,
										costs.test.sum       = costs.test.sum,
										int.costs.sum        = int.costs.sum,
										int.impl.costs.sum   = int.impl.costs.sum,
										int.sust.costs.sum   = int.sust.costs.sum)
		
	} else if (sums.only==TRUE){
		
		# Sums
		QALYs.sum            <- sum(unlist(QALYs.period))
		costs.hru.sum        <- sum(unlist(hru.costs.period))
		costs.art.sum        <- sum(unlist(art.costs.period))
		costs.art.ini.sum    <- sum(unlist(art.ini.costs.period))
		costs.oat.sum        <- sum(unlist(oat.costs.period))
		costs.prep.sum       <- sum(unlist(prep.costs.period))
		costs.prep.tests.sum <- sum(unlist(prep.tests.costs.period))
		costs.test.sum       <- sum(unlist(tests.costs.period))
		int.costs.sum        <- sum(unlist(int.costs.period))
		int.impl.costs.sum   <- sum(unlist(int.impl.costs.period))
		int.sust.costs.sum   <- sum(unlist(int.sust.costs.period))
		
		# Total costs
		costs.total.sum <- (costs.hru.sum + costs.art.sum + costs.art.ini.sum +
													costs.oat.sum + costs.prep.sum + costs.prep.tests.sum +
													costs.test.sum + int.costs.sum +
													int.impl.costs.sum + int.sust.costs.sum)
		
		# In list
		out.sum <- list(QALYs.sum            = QALYs.sum, 
										costs.total.sum      = costs.total.sum,
										costs.hru.sum        = costs.hru.sum,
										costs.art.sum        = costs.art.sum,
										costs.art.ini.sum    = costs.art.ini.sum,
										costs.oat.sum        = costs.oat.sum,
										costs.prep.sum       = costs.prep.sum,
										costs.prep.tests.sum = costs.prep.tests.sum,
										costs.test.sum       = costs.test.sum,
										int.costs.sum        = int.costs.sum,
										int.impl.costs.sum   = int.impl.costs.sum,
										int.sust.costs.sum   = int.sust.costs.sum)
	}
	
	
	return(out.sum) 
})
}