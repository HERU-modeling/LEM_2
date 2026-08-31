#################################################################################
# An R script to define ODE's of HIV model need to be solved
# Input the parameter values of the ODE with or without interventions,
# return the rate of change of all the states
# Updated: July 2026
##################################################################################
#
# Modifications to model time-varying parameters for interventions:
#
# 1. Immediate ART (phi)
# 2. ART initiation (alpha)
# 3. ART re-initiation (alpha prime)
# 4. ART retention (theta)
# 5. OAT: Proportion on OAT (FoI)
# 6. SSP: Proportion on Safer supply (FoI)
# 7. PrEP: Proportion on PrEP (FoI)
# 8. Testing rate (psi)
##################################################################################

#' ode_model
#'
#' \code{ode_model} defining ODE equations that need to be solved
#'
#' @param x vector of length 720 containing the number of individuals in 30 groups and 24 states (S1, S2, Sp, Ia, I1, I2, I3, Iap, Ip, Da, D1, D2, D3, T1, T2, T3, O1, O2, O3, inc_bo, inc_bs, inc_g, diag, diag_late, positive_test, viral_sup, death)
#' @param t time-period (month)
#' @param vparameters a list not a vector (defined in "Model-1-Module-Data.input.R")
#' @return
#' out: a vector of health state occupancy rate of change, HIV incidence, and diagnoses
#' @export

ode_model <- function(t, x, vparameters) {
  with(as.list(c(x, vparameters, gp.gn)), {

    n.gp <- length(names.gp)
    n.msm<- length(names.msm)
    # no. of opposite sex partners for all 30 groups 16 states
    no <- matrix(0, n.gp, 16)
    noG[noG < 0] <- 0 # can it be less than zero? need to revisit
    no[, 1:6] <- noG # for susceptible/infected
    no[, 7:16] <- noG * (1 - epsilonS) # for diagnosed/on ART/off ART, reduced by a fraction epsilonS

    # no. of same sexual partners
    ns <- matrix(0, n.msm, 16) 
    nsG[nsG < 0] <- 0
    ns[, 1:6] <- nsG # for susceptible/infected
    ns[, 7:16] <- nsG * (1 - epsilonS) # for diagnosed/on ART/off ART, reduced by a fraction epsilonS

    # probability of transmission by sex for 14 HIV+ states, on ART reduction by delta, 16 by 1
    # 14 HIV+ states, Ia,I1,I2,I3,Da,D1,D2,D3,T1,T2,T3,O1,O2,O3)
    sigmaFM <- c(sigmaO.FM, sigmaO.FM, sigmaO.FM[2:4] * (1 - delta_H), sigmaO.FM[2:4])
    sigmaMF <- c(sigmaO.MF, sigmaO.MF, sigmaO.MF[2:4] * (1 - delta_H), sigmaO.MF[2:4])
    sigmaM <- c(sigmaS, sigmaS, sigmaS[2:4] * (1 - delta_M), sigmaS[2:4])

    # probability of transmission per shared injection for 14 HIV+ states, on ART reduction by delta
    tau.all <- c(tau, tau, tau[2:4] * (1 - delta_I), tau[2:4])

    # condom use adjustment term for each group on 16 states,
    # (1-probability of condom use opposite sex * probability of condom effectiveness opposite sex)
    # uio 30 * 1, uis 18 *1,  kappa a number
    uio[uio > 1] <- 1 # need to revisit by bigger than 1
    uis[uis > 1] <- 1
    uoC <- matrix(1 - uio * kappa_H, n.gp , 16) # same across different HIV states
    usC <- matrix(1 - uis * kappa_M, n.msm, 16) # same across different HIV states

    # mortality, 30*16
    mo <- cbind(
      mor_S, mor_S, mor_Ia, mor_I1, mor_I2, mor_I3,
      mor_Ia, mor_I1, mor_I2, mor_I3, mor_T1, mor_T2, mor_T3,
      mor_I1, mor_I2, mor_I3
    )

    # x: 720 (30x24) initial values for 30 groups, 22 states (S1, S2, Ia, I1, I2, I3, 
    # Da, D1, D2, D3, T1, T2, T3, O1, O2, O3, inc_bo, inc_bs, inc_g, diagm, diag_late, pos_rate, viral_sup, deathm) 
    # 1:16 for states from S1 to O3
    # convert x to a matrix with each row representing a group.
    y <- matrix(x, nrow = n.gp, byrow = TRUE) # Convert  vector x to a matrix of 42 by 24 columns
    y2 <- y[, 1:16] # selects only the first 16 columns of the matrix y and assigns them to a new matrix y2
    row.names(y2) <- names.gp

# browser(expr = {t >= 0})     
if (sum(y[, 10], na.rm = TRUE) / sum(y, na.rm = TRUE) > 0.9) {
    	cat("⚠️ sum(y[,10]) / sum(y) > 0.9 — entering debug mode\n")
    	browser() 
}  
    # change n.gp group names without OAT, low, high
    rname <- gsub(paste(c("/low", "/high", "/old", "/young"), collapse = "|"), "", names.gp)

    # pwid group names without OAT, low, high
    pwid.name <- names18[grep("PWID", names18)]

    # hiv test rate per month
    # assign the parameters into all 42 groups
    phi <- matrix(0, n.gp, 3) # immediate ART initiation, cols 3 CD4 levels
    alpha <- matrix(0, n.gp, 3) # 	>= 30-days post-diagnosis ART initiation
    alpha.re <- numeric(n.gp) # ART re-initiation
    theta.t12 <- numeric(n.gp)
    theta.t13 <- numeric(n.gp)
    theta.t21 <- numeric(n.gp)
    theta.t23 <- numeric(n.gp)
    theta.t31 <- numeric(n.gp)
    theta.t32 <- numeric(n.gp)
    theta.t1O <- numeric(n.gp)
    theta.t2O <- numeric(n.gp)
    theta.t3O <- numeric(n.gp)
    mu_adj    <- numeric(n.gp)
    mu_mat    <- numeric(n.gp) # maturation-out rate
    old.e.sus <- numeric(n.gp)
    old.e.hiv <- numeric(n.gp)


    ## TO UPDATE NEW ART ENGAGEMENT ##
    if (t <= 12) {                 # At & Before Dec 2017
      phi1.t <- phi1[1, ]
      phi2.t <- phi2[1, ]
      phi3.t <- phi3[1, ]
      alpha1.t <- alpha1[1, ]
      alpha2.t <- alpha2[1, ]
      alpha3.t <- alpha3[1, ]
    } else if (t > 12 & t <= 24) { # Jan 2018 - Dec 2018
      phi1.t <- phi1[2, ]
      phi2.t <- phi2[2, ]
      phi3.t <- phi3[2, ]
      alpha1.t <- alpha1[2, ]
      alpha2.t <- alpha2[2, ]
      alpha3.t <- alpha3[2, ]
    } else if (t > 24 & t <= 36) { # Jan 2019 - Dec 2019
      phi1.t <- phi1[3, ]
      phi2.t <- phi2[3, ]
      phi3.t <- phi3[3, ]
      alpha1.t <- alpha1[3, ]
      alpha2.t <- alpha2[3, ]
      alpha3.t <- alpha3[3, ]
    } else if (t > 36 & t <= 48) { # Jan 2020 - Dec 2020
      phi1.t <- phi1[4, ]
      phi2.t <- phi2[4, ]
      phi3.t <- phi3[4, ]
      alpha1.t <- alpha1[4, ]
      alpha2.t <- alpha2[4, ]
      alpha3.t <- alpha3[4, ]
    } else if (t > 48 & t <= 60) { # Jan 2021 - Dec 2021
      phi1.t <- phi1[5, ]
      phi2.t <- phi2[5, ]
      phi3.t <- phi3[5, ]
      alpha1.t <- alpha1[5, ]
      alpha2.t <- alpha2[5, ]
      alpha3.t <- alpha3[5, ]
    } else if (t > 60) { # Jan 2022 onward
      phi1.t <- phi1[6, ]
      phi2.t <- phi2[6, ]
      phi3.t <- phi3[6, ]
      alpha1.t <- alpha1[6, ]
      alpha2.t <- alpha2[6, ]
      alpha3.t <- alpha3[6, ]
    }
    ## END UPDATING NEW ART ENGAGEMENT ##
    for (i in 1:18) {
      ind <- which(rname %in% names18[i])
      phi[ind, 1] <- phi1.t[i]
      phi[ind, 2] <- phi2.t[i]
      phi[ind, 3] <- phi3.t[i]
      alpha[ind, 1] <- alpha1.t[i]
      alpha[ind, 2] <- alpha2.t[i]
      alpha[ind, 3] <- alpha3.t[i]
      alpha.re[ind] <- O_T[i]
      theta.t12[ind] <- T1_T2[i]
      theta.t13[ind] <- T1_T3[i]
      theta.t21[ind] <- T2_T1[i]
      theta.t23[ind] <- T2_T3[i]
      theta.t31[ind] <- T3_T1[i]
      theta.t32[ind] <- T3_T2[i]
      theta.t1O[ind] <- T1_O[i]
      theta.t2O[ind] <- T2_O[i]
      theta.t3O[ind] <- T3_O[i]
      mu_adj[ind] <- vparameters$mu_adj[i]
      mu_mat[ind] <- vparameters$mu_mat[i]
      old.e.sus[ind] <- vparameters$old.e.sus[i]
      old.e.hiv[ind] <- vparameters$old.e.hiv[i]      
      prop_vs_T[ind] <- vparameters$prop_vs_T[i]
    }
# browser(expr = {t >= 60})

    ### Time-varying number of OAT initiators ###
   if (t <= 12)          prop.oat.t <- prop.oat.m[1,] # 2017 Dec and before
   if (t > 12 & t <= 24) prop.oat.t <- prop.oat.m[2,] # 2018 Jan - Dec
   if (t > 24 & t <= 36) prop.oat.t <- prop.oat.m[3,] # 2019 Jan - Dec
   if (t > 36 & t <= 48) prop.oat.t <- prop.oat.m[4,] # 2020 Jan - Dec
   if (t > 48)           prop.oat.t <- prop.oat.m[5,] # 2021 Jan and later


   if (intervention.oat.modif == TRUE & t > int.start) {
	  prop.oat.hold <- prop.oat.t
	  if (t <= int.end) {
	  	 prop.oat.t <- int.oat.eff[t - int.start,]
	   } else {
	     prop.oat.t  <- prop.oat.hold
  	 }
    }

    ### Time-varying number of SSP ###
		if (t <= 12)          prop.ssp.t <- prop.ssp.m[1,] # 2017 Dec and before
		if (t > 12 & t <= 24) prop.ssp.t <- prop.ssp.m[2,] # 2018 Jan - Dec
		if (t > 24 & t <= 36) prop.ssp.t <- prop.ssp.m[3,] # 2019 Jan - Dec
		if (t > 36 & t <= 48) prop.ssp.t <- prop.ssp.m[4,] # 2020 Jan - Dec
		if (t > 48)           prop.ssp.t <- prop.ssp.m[5,] # 2021 Jan and later
   
    ############################
    # SSP, intervention
    if (intervention.ssp.modif == TRUE & t > int.start) {
    	prop.ssp.hold <- prop.ssp.t
    	if (t <= int.end) {
    		prop.ssp.t <- int.ssp.eff[t - int.start,]
    	} else {
    		prop.ssp.t <- prop.ssp.hold
    	}
    }
    ############################
    # assign shared infjections for all n.gp groups
    s.all <- numeric(n.gp)

    for (i in seq_along(names.pwid)) {
      ind <- which(rname %in% names.pwid[i])
      s.all[ind] <- s[i]
    }
    
    # Weighted average impact of OAT: Reduce mortality
     mo <- mo * vparameters$oat_mor * prop.oat.t + mo * (1-prop.oat.t)

  
    # sufficient contact rate (y2=y[,1:19], excluding incidence/new diagnosis)
    # output: matrix of c(lambda, gamma, beta_o,  beta_s),
    # each of them with dim (n.gp,2) with row of each group, second column for PrEP
    # lambda- over all transition, gamma - needle sharing, beta_o - heterosexual, beta_s - same sex
    foi <- FoI(
      y = y2, no = no, uoC = uoC, ns = ns, usC = usC, eO = ass.eO, eS = ass.eS, sigmaFM = sigmaFM, sigmaMF = sigmaMF, sigmaM = sigmaM, tau = tau.all,
      eff.prep = eff.prep, d = d, s = s.all, eff.oat = eff.oat, 
      prop.ssp.t = prop.ssp.t, eff.ssp = eff.ssp, s.multi = s.multi, bal, prop.oat.t = prop.oat.t,
      adherence.prep = adherence.prep, mix.prob.HL = mix.prob.HL, mix.prob.YO = mix.prob.YO
    )
# browser(expr = {t >= 0})
    ############################

    ## TO UPDATE NEW PrEP ##
    if (t <= 12)          prop.prep.t <- prop.prep.m[1,] # 2017 Dec and before
    if (t > 12 & t <= 24) prop.prep.t <- prop.prep.m[2,] # 2018 Jan - Dec
    if (t > 24 & t <= 38) prop.prep.t <- prop.prep.m[3,] # 2019 Jan - Dec & 2 first months 2020
    if (t > 38 & t <= 50) prop.prep.t <- prop.prep.m[4,] # 2020 Mar - 2021 Feb - Covid
    if (t > 50 & t <= 60) prop.prep.t <- prop.prep.m[5,] # 2021 Mar - Dec
    if (t > 60)           prop.prep.t <- prop.prep.m[6,] # 2022 Jan and later
    
    ## END UPDATING NEW PrEP ##

    # PrEP, intervention - Move BEFORE to ensure prep.total[8], thus prep.total.t value is update with PrEP intervention
    # 2020 onward (t > int.start=96)
    if (intervention.prep.modif == TRUE & t > int.start) {
    	#    	prep.total.hold <- prep.total[8]
    	prop.prep.hold <- prop.prep.t
    	if (t <= int.end) {
    		#    		prep.total[8] <- int.prep.eff[t - int.start]
    		prop.prep.t <- int.prep.eff[t - int.start,]
    	} else {
    		#    		prep.total[8] <- prep.total.hold
    		prop.prep.t  <- prop.prep.hold
    	}
    }
    ############################
    # psi.t: time-varying screening rate
    if (t <= 38) psi.t <- psi.m[1, ]
    if (t > 38 & t <= 39) psi.t <- psi.m[2, ]
    if (t > 39 & t <= 40) psi.t <- psi.m[3, ]
    if (t > 40 & t <= 41) psi.t <- psi.m[4, ]
    if (t > 41 & t <= 42) psi.t <- psi.m[5, ]
    if (t > 42 & t <= 43) psi.t <- psi.m[6, ]
    if (t > 43 & t <= 44) psi.t <- psi.m[7, ]
    if (t > 44) psi.t <- psi.m[8, ]
    
    # Testing, 30*16
     psi <- matrix(0, nrow = 30, ncol = 16)
     colnames(psi) <- vparameters$state.name[1:16]		
     psi[ , c("S1", "Ia", "I1", "I2", "I3")] <- psi.t

    ############################
    # psi (screening rate), intervention, increase baseline psi
    if (intervention.testing.modif == TRUE & t > int.start) {
      psi.hold <- psi
      if (t <= int.end) {
        psi <- psi * t(int.test.eff[(t - int.start), ,])
      } else {
        psi <- psi.hold
      }
    }  
 
    # Include psi.p into psi
		psi[ , c("S1","S2")] <- psi[ , c("S1","S2")] * (1-prop.prep.t) + psi.p * prop.prep.t #correspond with Sp, no longer consider Iap, Ip per Xiao 10-21
		# Include v2,v3 into psi
    psi[ , c("I2")] <- psi[ , c("I2")] + v2
    psi[ , c("I3")] <- psi[ , c("I3")] + v3

    ############################
    # alpha (>= 30-days post-diagnosis ART initiation), intervention
    if (intervention.initiation.modif == TRUE & t > int.start) {
      alpha.hold <- alpha
      if (t <= int.end) {
        alpha <- alpha * int.art.ini.eff[(t - int.start), ]
      } else {
        alpha <- alpha.hold
      }
    }

    ############################
    # INCREASED RETENTION, intervention
    if (intervention.retention.modif == TRUE & t > int.start) {
      theta1.hold <- theta.t1O
      theta2.hold <- theta.t2O
      theta3.hold <- theta.t3O

      if (t <= int.end) {
        if (any(current.int == "ART retention")) {
          theta.t1O <- theta.t1O * (1 / int.art.ret.eff[(t - int.start), ])
          theta.t2O <- theta.t2O * (1 / int.art.ret.eff[(t - int.start), ])
          theta.t3O <- theta.t3O * (1 / int.art.ret.eff[(t - int.start), ])
        } else if (any(current.int == "ART retention, targeted")) {
          theta.t3O <- theta.t3O * (1 / int.art.ret.eff[(t - int.start), ])
        }
      } else {
        theta.t1O <- theta1.hold
        theta.t2O <- theta2.hold
        theta.t3O <- theta3.hold
      }
    }

    # DECREASED DROPOUT
    if (intervention.dropout.modif == TRUE & t > int.start) {
      theta1.hold <- theta.t1O
      theta2.hold <- theta.t2O
      theta3.hold <- theta.t3O
      if (t <= int.end) {
        theta.t1O <- theta.t1O * int.art.drop.eff[(t - int.start), ]
        theta.t2O <- theta.t2O * int.art.drop.eff[(t - int.start), ]
        theta.t3O <- theta.t3O * int.art.drop.eff[(t - int.start), ]
      } else {
        theta.t1O <- theta1.hold
        theta.t2O <- theta2.hold
        theta.t3O <- theta3.hold
      }
    }
    
    # Weighted average impact of OAT: Reduce drop out
    theta.t1O <- theta.t1O * theta.o.oat * prop.oat.t + theta.t1O * (1 - prop.oat.t) 
    theta.t2O <- theta.t2O * theta.o.oat * prop.oat.t + theta.t2O * (1 - prop.oat.t) 
    theta.t3O <- theta.t3O * theta.o.oat * prop.oat.t + theta.t3O * (1 - prop.oat.t) 

    ############################
    # phi (immediate ART initiation), intervention
    if (intervention.immediateART.modif == TRUE & t > int.start) {
      phi.hold <- phi
      if (t <= int.end) {
        phi <- phi * int.art.imm.eff[(t - int.start), ]
      } else {
        phi <- phi.hold
      }
    }


    ############################
    # alpha.re (ART re-initiation rate), intervention
    if (intervention.reinitiation.modif == TRUE & t > int.start) {
      alpha.re.hold <- alpha.re
      if (t <= int.end) {
        alpha.re <- alpha.re * int.art.reini.eff[(t - int.start), ]
      } else {
        alpha.re <- alpha.re.hold
      }
    }
   
    ############################
    #Goal: Adjust pop change so that 1) respect the input in regards to races , 2) but
    # maintain the proportion of each risk group as close as the initial population as possible
if (prop.adj == TRUE) {    
# browser(expr = {t >= 0})    
    target.pop <- model.pop <- numeric(length(vparameters$names18))
    tot.group <- rowSums(y2)
   #target.pop: Initial 18 strata population, based on init.collap
    target.pop[intersect(gp18.gn$white, gp18.gn$msm)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$white, gp18.gn$msm)]/sum(init.collap)
    target.pop[intersect(gp18.gn$black, gp18.gn$msm)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$black, gp18.gn$msm)]/sum(init.collap)
    target.pop[intersect(gp18.gn$hisp , gp18.gn$msm)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$hisp , gp18.gn$msm)]/sum(init.collap)
    target.pop[intersect(gp18.gn$white, gp18.gn$mwid)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$white, gp18.gn$mwid)]/sum(init.collap)
    target.pop[intersect(gp18.gn$black, gp18.gn$mwid)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$black, gp18.gn$mwid)]/sum(init.collap)
    target.pop[intersect(gp18.gn$hisp , gp18.gn$mwid)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$hisp , gp18.gn$mwid)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$male)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$male)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$male)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$female)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$female)]/sum(init.collap)
    target.pop[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group) * init.collap[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$female)]/sum(init.collap)
    target.pop[intersect(gp18.gn$white, gp18.gn$het.m)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$white, gp18.gn$het.m)]/sum(init.collap)
    target.pop[intersect(gp18.gn$black, gp18.gn$het.m)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$black, gp18.gn$het.m)]/sum(init.collap)
    target.pop[intersect(gp18.gn$hisp , gp18.gn$het.m)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$hisp , gp18.gn$het.m)]/sum(init.collap)
    target.pop[intersect(gp18.gn$white, gp18.gn$het.f)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$white, gp18.gn$het.f)]/sum(init.collap)
    target.pop[intersect(gp18.gn$black, gp18.gn$het.f)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$black, gp18.gn$het.f)]/sum(init.collap)
    target.pop[intersect(gp18.gn$hisp , gp18.gn$het.f)]  <- sum(tot.group) * init.collap[intersect(gp18.gn$hisp , gp18.gn$het.f)]/sum(init.collap)
   #model.pop:  Current 18 strata popluation, resulted from ODE
    model.pop[intersect(gp18.gn$white, gp18.gn$msm)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$msm)])
    model.pop[intersect(gp18.gn$black, gp18.gn$msm)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$msm)])
    model.pop[intersect(gp18.gn$hisp , gp18.gn$msm)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$msm)])
    model.pop[intersect(gp18.gn$white, gp18.gn$mwid)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$mwid)])
    model.pop[intersect(gp18.gn$black, gp18.gn$mwid)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$mwid)])
    model.pop[intersect(gp18.gn$hisp , gp18.gn$mwid)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$mwid)])
    model.pop[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$pwid.m)])
    model.pop[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$pwid.m)])
    model.pop[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$male)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$pwid.m)])
    model.pop[intersect(intersect(gp18.gn$white, gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$pwid.f)])
    model.pop[intersect(intersect(gp18.gn$black, gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$pwid.f)])
    model.pop[intersect(intersect(gp18.gn$hisp , gp18.gn$pwid), gp18.gn$female)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$pwid.f)])
    model.pop[intersect(gp18.gn$white, gp18.gn$het.m)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$het.m)])
    model.pop[intersect(gp18.gn$black, gp18.gn$het.m)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$het.m)])
    model.pop[intersect(gp18.gn$hisp , gp18.gn$het.m)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$het.m)])
    model.pop[intersect(gp18.gn$white, gp18.gn$het.f)]  <- sum(tot.group[intersect(gp.gn$white,gp.gn$het.f)])
    model.pop[intersect(gp18.gn$black, gp18.gn$het.f)]  <- sum(tot.group[intersect(gp.gn$black,gp.gn$het.f)])
    model.pop[intersect(gp18.gn$hisp , gp18.gn$het.f)]  <- sum(tot.group[intersect(gp.gn$hisp ,gp.gn$het.f)])
    
    # ---- Define lists ----
    race_list <- list(
    	white = gp18.gn$white,
    	black = gp18.gn$black,
    	hisp  = gp18.gn$hisp
    )
    risk_list <- list(
    	msm  = gp18.gn$all.msm,
    	pwid = gp18.gn$pwid,
    	het  = gp18.gn$het
    )
    
    pop.change.18.new <- numeric(length(pop.change.18))
    #Goal: Adjust pop change so that 1) respect the input in regards to races , 2) but
    # maintain the proportion of each risk group as close as the initial population as possible

    # Loop by race 
    
    for (race in names(race_list)) {
    	
    	idx_race <- race_list[[race]]
    	
    	# total change for this race (must be preserved)
    	race_total_change <- sum(pop.change.18[idx_race])   #Step 1: Compute total change for that race
   	
    	# current and target totals by risk within this race
    	current_risk_tot <- sapply(risk_list, function(risk_idx) {
    		sum(model.pop[intersect(idx_race, risk_idx)])    	})
    	
    	target_risk_tot <- sapply(risk_list, function(risk_idx) {
    		sum(target.pop[intersect(idx_race, risk_idx)])    	})
    	
    	# proportions
    	current_prop <- current_risk_tot / sum(current_risk_tot) #Step 2: Compute current risk proportions within that race
    	target_prop  <- target_risk_tot / sum(target_risk_tot)   #Step 3: Compute target risk proportions within that race
    	
    	# gap at risk level (inside race only)
    	risk_gap <- target_prop - current_prop
    	
    	# allocate change depending on direction - If gap = 0 (the case for t=0), assign pop change to each risk group equally; otherwise allocate change based on gap/difference (model vs target pop) size
    	if (race_total_change > 0) {
    		pos_gap <- pmax(risk_gap, 0)
    		if (sum(pos_gap) < 1e-12) { #Ignore noise
    			weights <- rep(1/length(pos_gap), length(pos_gap))
    		} else {
    			weights <- pos_gap / sum(pos_gap)
    		}
    	} else {
    		neg_gap <- pmax(-risk_gap, 0)
    		if (sum(neg_gap) < 1e-12) {
    			weights <- rep(1/length(neg_gap), length(neg_gap))
    		} else {
    			weights <- neg_gap / sum(neg_gap)
    		}
    	}
    	
    	risk_change <- race_total_change * weights  # Step 4: Allocate the race's total change across its risks according to risk gap (i.e., possible for W/Het to make up for W/PWID, but not B/Het or B/PWID)
# Note when race_total_change > 0; We are adding people to this race, thus We should only add to underrepresented risks, vice versa. If we remove people for this race race_total_change < 0, We should only remove from overrepresented risks.
    	# distribute inside each (race, risk)
    	for (i in seq_along(risk_list)) {
    		
    		idx_cell <- intersect(idx_race, risk_list[[i]])
    		if (length(idx_cell) == 0) next
  		
    		inside_prop <- target.pop[idx_cell] / sum(target.pop[idx_cell]) # This step indirectly ensure that proportion for msm.Y & msm.O (in MSM), pwid.f & pwid.m (in pwid), het.f & het.m (in Het) is also mantain
    		pop.change.18.new[idx_cell] <- risk_change[i] * inside_prop
    	}
    }
    
    # Checks 
    # total preserved
    stopifnot(abs(sum(pop.change.18.new) - sum(pop.change.18)) < 1e-8)
    
    # race totals preserved
    for (race in race_list) {
    	stopifnot(
    		abs(sum(pop.change.18.new[race]) - 
    					sum(pop.change.18[race])) < 1e-8
    	)
    }
}
else {pop.change.18.new = pop.change.18}    
    # Create placeholder & new parameter
    pop.change.30 <- numeric(n.gp)
    prop.high.sus.all <- numeric(n.gp)
    prop.young.sus.all <- numeric(n.gp)
    
    # Go from 18 to 30 
    for (i in 1:18) {
    	ind <- which(rname %in% names18[i])
    	#Set up required proportion, based on init conversion in model initial, try to re-use same value & keep it simple, pick prop.high.sus rather than prop.high.inf by random
    	prop.high.sus.all[ind] = prop.high.sus[i]
    	prop.young.sus.all[ind] = prop.young.sus[i]
    	#Populate pop.change.30 itself using value from pop.change.18.new	
    	pop.change.30[ind] <- pop.change.18.new[i]
    }
    # Distribute pop.change across high-young vs low-old
    pop.change.30[high] = pop.change.30[high] *as.numeric(prop.high.sus.all[high])
    pop.change.30[low]  = pop.change.30[low]   *(1-as.numeric(prop.high.sus.all[low]))
    pop.change.30[old]  = 0 # All msm pop.change got to young
    
    # Adjustment to make sure high-young stay const over time
    tot.group <- rowSums(y2)

    # Adjustment to make sure high-young stay const over time
    dif.het.high  <- (tot.group[low] + tot.group[high]) * as.numeric(prop.high.sus.all[high]) - tot.group[high]
    pop.change.30[high]  <- pop.change.30[high]  + dif.het.high
    pop.change.30[low]   <- pop.change.30[low] - dif.het.high

############################
# This code session is backbone for the "out" 30*21 result matrix
    # out is the population difference between time t and t+1
    out <- matrix(0, n.gp, 24)
		# Code below, populate the 30*28 matrix "out", 30 rows are indexed by [i, ] representing risk groups
    # PWID or HET
    for (i in c(pwid, het)) {
      out[i, ] <- ode.list(
        y[i, ], foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
				mu_adj[i], ws, psi[i,],
        theta.ai, theta.ad,
        phi[i, ], alpha[i, ], alpha.re[i],
        theta.1, theta.2,
        theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
        theta.t31[i], theta.t32[i],
        theta.t1O[i], theta.t2O[i], theta.t3O[i],
        mo[i, ], pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
      )
    }
    for (j in seq_along(young)) {
		i <- young[j]	
  	  out[i, ] <- ode.list.MSMyoung(
		    y[i, ], y[old[j], ],
		    foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
		    mu_mat[i], mu_adj[i], ws, psi[i,], # mu_adj not needed for young
        theta.ai, theta.ad,
		    phi[i, ], alpha[i, ], alpha.re[i],
		    theta.1, theta.2,
		    theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
		    theta.t31[i], theta.t32[i],
        theta.t1O[i], theta.t2O[i], theta.t3O[i], 
        mo[i, ], mo[old[j], ], old.e.sus[i], old.e.hiv[i],
		    pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
      )
     }
    for (j in seq_along(old)) {
	   i <- old[j]
	    out[i, ] <- ode.list.MSMold(
		    y[i, ], y[young[j], ], 
		    foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
		    mu_mat[i], mu_adj[i], ws, psi[i,],
        theta.ai, theta.ad,
		    phi[i, ], alpha[i, ], alpha.re[i],
		    theta.1, theta.2,
		    theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
		    theta.t31[i], theta.t32[i],
		    theta.t1O[i], theta.t2O[i], theta.t3O[i], 
        mo[i, ], old.e.sus[i], old.e.hiv[i],
		    pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
	     ) 
      }
    list(as.vector(t(out)))
# browser()      
  })
}
