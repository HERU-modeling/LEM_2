##############################################################################
## Cascade CEA Model - Interventions Scale-up Coefficient
## Latest update: Oct 2022
##############################################################################
## Output: a list int.pop.eff containing parameter scale change multiplier due
## to the intervention
##############################################################################
# interventions <- c("Opt-out testing (ER)", "Opt-out testing (PC)",
#                    "EMR testing reminder", "Nurse-initiated testing", "OAT integrated testing",
#                    "ART initiation", "ART retention", "ART retention, targeted",
#                    "ART EMR reminder", "RAPID ART",
#                    "ART re-initiation", "ART re-linkage",
#                    "SSP", "OAT with BUP", "OAT with methadone", "PrEP",
#                    "No interventions")
#
# Interventions, suffix
# "Opt-out testing (ER)", test.optout.ER
# "Opt-out testing (PC)", test.optout.PC
# "EMR testing reminder", test.EMR
# "Nurse-initiated testing", test.nurse
# "OAT integrated testing", test.OAT
# 
# "ART initiation", ARTinit
# "ART retention", ARTreten
# "ART retention, targeted", ARTreten.target
# "ART EMR reminder", ARTEMRprompt
# "RAPID ART", immART
# 
# "ART re-initiation", reART
# "ART re-linkage", relink
# 
# "SSP", SSP
# "OAT with BUP", BUP
# "OAT with methadone", MET
# "PrEP", prep
#
#############################################################################

#' Interventions Scale-up Coefficient
#'
#' @param current.int a vector of all the interventions investigated
#' @param tf.ls a list Value of logical vectors for intervention modeling
#' @param scale.up.period scale up increasing period, after fixed level
#'
#' @return a list int.pop.eff containing parameter scale change multiplier due to the intervention
#' @export
#'
#' @examples
scale.interventions.combination <- function(current.int,  tf.ls,
                                            scale.up.period = 18){
  
  # Set scale inputs, 18 scale up increasing period, after fixed level
	# Int.Scale.ls, Int.Eff.ls, Int.Baseline.ls are loaded from data.input
	# int.end int.start are set in loadBaselineWorkspace
#browser() #Check current.int  
	with(as.list(c(Int.Scale.ls, Int.Eff.ls, Int.Baseline.ls, vparameters)),{
  if (tf.ls$intervention.testing.modif == TRUE){
#    int.test.eff   <- matrix(0, nrow=(int.end-int.start), ncol = length(names.gp))
  	#Need ways to figure out int.test.eff: 3dim time 120, state 16, strata 30
    int.test.eff   <- array(0, dim = c(
    	n=(int.end-int.start), length(vparameters$state.name[1:16]),
    	length(vparameters$names.gp)
    ))
    if (is.element(c("Opt-out testing (ER)"), current.int)){ 
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
      
      int.test.eff <- scale.m.temp * (scale.up.v - 1) + int.test.eff 
    } 
    if (is.element(c("Opt-out testing (PC)"), current.int)){ 
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
      
      int.test.eff <- scale.m.temp * (scale.up.v - 1) + int.test.eff
    } 
    if (is.element(c("EMR testing reminder"), current.int)){ 
      # Intervention effectiveness
			temp.int     <- IntEff.test.EMR
			scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
												rep(temp.int, (int.end - int.start - scale.up.period)))
			# Scale of delivery: Change v to m since this is a matrix
#     scale.v      <- IntScale.test.EMR["S1", ]
      scale.m      <- IntScale.test.EMR
#browser()      
      reweight    <- sum(inits[,c("Ia","I1","I2","I3")])/sum(inits[,c("S1")])*((1-IntPosRate.test.EMR)/IntPosRate.test.EMR)
			scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.EMR) * 1/(1 + reweight) / 4,
																									nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
			scale.m["S1", ]                   <-        colSums(IntScale.test.EMR) * reweight/(1 + reweight)
#browser()	#Check colSums(IntScale.test.EMR) vs colSums(scale.m)
			scale.m.temp <- replicate((int.sus * 12), scale.m)
			scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
			
			int.test.eff <- scale.m.temp * (scale.up.v - 1) + int.test.eff
		
    }
    if (is.element(c("Nurse-initiated testing"), current.int)){ 
      temp.int     <- IntEff.test.nurse
      scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
      									rep(temp.int, (int.end - int.start - scale.up.period)))
#      scale.v      <- IntScale.test.nurse["S1", ]
      scale.m      <- IntScale.test.nurse
      reweight     <- sum(inits[,c("Ia","I1","I2","I3")])/sum(inits[,c("S1")])*((1-IntPosRate.test.nurse)/IntPosRate.test.nurse)
      scale.m[c("Ia","I1","I2","I3"), ] <- matrix(colSums(IntScale.test.nurse) * 1/(1 + reweight) / 4,
      																								 nrow = 4,	ncol = ncol(scale.m),	byrow = TRUE)
      scale.m["S1", ]                   <-        colSums(IntScale.test.nurse) * reweight/(1 + reweight)      
#      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
      scale.m.temp <- replicate((int.sus * 12), scale.m)
      scale.m.temp <- aperm(scale.m.temp, perm = c(3, 1, 2))
     
      int.test.eff <- scale.m.temp * (scale.up.v - 1) + int.test.eff
    }
    if (is.element(c("OAT integrated testing"), current.int)){ 
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
    	
      int.test.eff <- scale.m.temp * (scale.up.v - 1) + int.test.eff
    }
    int.test.eff   <- int.test.eff + 1
  } else {int.test.eff = NULL}
 
  if (tf.ls$intervention.initiation.modif == TRUE){
    temp.int       <- IntEff.ARTinit
    scale.v        <- IntScale.ARTinit["D1", ]
    scale.m.temp   <- t(replicate((int.sus * 12), scale.v))
    scale.up.v     <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                        rep(temp.int, (int.end - int.start - scale.up.period)))
    int.art.ini.eff <- scale.m.temp * (scale.up.v - 1) + 1
  } else {int.art.ini.eff = NULL}
  
  if (tf.ls$intervention.retention.modif == TRUE){
    if (is.element(c("ART retention"), current.int)) {
      temp.int    <- IntEff.ARTreten
      scale.v     <- IntScale.ARTreten["T1", ]
    } else if (is.element(c("ART retention, targeted"), current.int)) {
      temp.int    <- IntEff.ARTreten.target
      scale.v     <- IntScale.ARTreten.target["T3", ]
    }
    scale.m.temp  <- t(replicate((int.sus * 12), scale.v))
    scale.up.v    <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                       rep(temp.int, (int.end - int.start - scale.up.period)))
    int.art.ret.eff <- scale.m.temp * (scale.up.v - 1) + 1
  } else {int.art.ret.eff = NULL}
  
  if (tf.ls$intervention.dropout.modif == TRUE) {
    temp.int      <- IntEff.ARTEMRprompt
    scale.v       <- IntScale.ARTEMRprompt["T1", ]
    scale.m.temp  <- t(replicate((int.sus * 12), scale.v))
    scale.up.v    <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                       rep(temp.int, (int.end - int.start - scale.up.period)))
    int.art.drop.eff <- scale.m.temp * (scale.up.v - 1) + 1
  } else {int.art.drop.eff = NULL}
  
  if (tf.ls$intervention.immediateART.modif == TRUE) {
    temp.int      <- IntEff.immART
    scale.v       <- IntScale.immART["I1", ]
    scale.m.temp  <- t(replicate((int.sus * 12), scale.v))
    scale.up.v    <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                       rep(temp.int, (int.end - int.start - scale.up.period)))
    int.art.imm.eff <- scale.m.temp * (scale.up.v - 1) + 1
  } else {int.art.imm.eff = NULL}
  
  if (tf.ls$intervention.reinitiation.modif == TRUE) {
    int.art.reini.eff <- matrix(0, nrow=(int.end-int.start), ncol = length(names.gp))
    if (is.element(c("ART re-initiation"), current.int)){ 
      temp.int     <- IntEff.reART
      scale.v      <- IntScale.reART["O1", ]
      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
      scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                        rep(temp.int, (int.end - int.start - scale.up.period)))
      int.art.reini.eff <- scale.m.temp * (scale.up.v - 1) + int.art.reini.eff
    } 
    if (is.element(c("ART re-linkage"), current.int)){ 
      temp.int     <- IntEff.relink
      scale.v      <- IntScale.relink["O1", ]
      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
      scale.up.v   <- c(seq(1, temp.int, by = ((temp.int - 1) / scale.up.period))[-1], 
                        rep(temp.int, (int.end - int.start - scale.up.period)))
      int.art.reini.eff <- scale.m.temp * (scale.up.v - 1) + int.art.reini.eff
    }
    int.art.reini.eff <- int.art.reini.eff + 1 
  } else {int.art.reini.eff = NULL}
  
  if (tf.ls$intervention.ssp.modif == TRUE) {
    scale.v        <- IntBaseline.SSP
    scale.m.temp   <- t(replicate((int.sus * 12), scale.v))
    scale.growth   <- IntScale.SSP
    scale.up.v     <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
                        rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))
    int.ssp.eff    <- scale.m.temp * scale.up.v
    int.ssp.eff[int.ssp.eff > IntEff.SSP] <- IntEff.SSP 
  } else {int.ssp.eff = NULL}
    
  if (tf.ls$intervention.oat.modif == TRUE){
    int.oat.eff    <- matrix(0, nrow=(int.end-int.start), ncol = length(IntBaseline.OAT))
#    if(is.element(c("OAT with BUP"), current.int)) {
#      scale.growth <- IntScale.BUP
#      scale.v      <- IntBaseline.BUP
#      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
#      scale.up.v   <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
#                        rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))
#      int.oat.eff  <- scale.m.temp * (scale.up.v - 1) + int.oat.eff
#   } 
#    if(is.element(c("OAT with methadone"), current.int)) {
#      scale.growth <- IntScale.MET
#      scale.v      <- IntBaseline.MET
#      scale.m.temp <- t(replicate((int.sus * 12), scale.v))
#      scale.up.v   <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
#                        rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))
#      int.oat.eff  <- scale.m.temp * (scale.up.v - 1) + int.oat.eff
#    }
#    int.oat.eff    <- int.oat.eff + t(replicate(nrow(int.oat.eff), IntBaseline.BUP))  + t(replicate(nrow(int.oat.eff), IntBaseline.MET))
    
		 if(any(current.int %in% c("OAT with BUP", "OAT with methadone"))) {
    	scale.v      <- IntBaseline.OAT
    	scale.m.temp <- t(replicate((int.sus * 12), scale.v))
    	scale.growth <- IntScale.OAT
    	scale.up.v   <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
    										rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))
    	int.oat.eff  <- scale.m.temp * scale.up.v
    	int.oat.eff[int.oat.eff > IntEff.OAT] <- IntEff.OAT 
		 }
  } else {int.oat.eff = NULL}
  
  if (tf.ls$intervention.prep.modif == TRUE){
    scale.v        <- IntBaseline.prep
    scale.m.temp   <- t(replicate((int.sus * 12), scale.v))
    scale.growth   <- IntScale.prep
    scale.up.v     <- c(cumprod(c(1, rep((1 + scale.growth), scale.up.period)))[-1],
                        rep(((1 + scale.growth)^scale.up.period), (int.end - int.start - scale.up.period)))
    int.prep.eff   <- scale.m.temp * scale.up.v
    int.prep.eff[int.prep.eff > IntEff.prep] <- IntEff.prep 
    } else {int.prep.eff = NULL}
  
  if(is.element(c("No interventions"), current.int)){
    warning("There are no interventions being modeled")
  }
  
  int.pop.eff <- list(int.test.eff = int.test.eff,         int.art.ini.eff = int.art.ini.eff, int.art.ret.eff = int.art.ret.eff, 
                      int.art.drop.eff = int.art.drop.eff, int.art.imm.eff = int.art.imm.eff, int.art.reini.eff = int.art.reini.eff,
                      int.ssp.eff = int.ssp.eff,           int.oat.eff = int.oat.eff,         int.prep.eff = int.prep.eff)
  
  # when it happens, revisit
  if(any(lengths(int.pop.eff)==0)){
    int.pop.eff <- int.pop.eff[-which(lengths(int.pop.eff)==0)]
  }
  return(int.pop.eff)
	})
}
