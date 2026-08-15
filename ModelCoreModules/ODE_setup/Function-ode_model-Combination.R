#################################################################################
# An R script to define ODE's of HIV model need to be solved
# Input the parameter values of the ODE with or without interventions,
# return the rate of change of all the states
# Updated: Oct 2022
##################################################################################
#
# Modifications to model time-varying parameters for interventions:
#
# 1. Testing rate (psi)
# 2. ART initiation (alpha)
# 3. ART retention (theta)
# 4. Immediate ART (phi)
# 5. ART re-initiation (alpha prime)
# 6. SSP: Number of syringes available (FoI)
# 7. OAT: Number of individuals receiving OAT (FoI)
# 8. PrEP: Rate of PrEP initiation (FoI)
##################################################################################

#' ode_model
#'
#' \code{ode_model} defining ODE equations that need to be solved
#'
#' @param x vector of length 1008 containing the number of individuals in 42 groups and 24 states (S1, S2, Sp, Ia, I1, I2, I3, Iap, Ip, Da, D1, D2, D3, T1, T2, T3, O1, O2, O3, inc_bo, inc_bs, inc_g, diagm, deathm)
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

#browser(expr = {t >= 347}) #rownames(mo) <- vparameters$names.gp, round(mo,4), check y at t=1 and t =349 using browser    
if (sum(y[, 10], na.rm = TRUE) / sum(y, na.rm = TRUE) > 0.9) {
    	cat("⚠️ sum(y[,10]) / sum(y) > 0.9 — entering debug mode\n")
    	browser() 
}  
    # change n.gp group names without OAT, low, high
    rname <- gsub(paste(c("/low", "/high", "/old", "/young"), collapse = "|"), "", names.gp)

    # pwid group names without OAT, low, high
    pwid.name <- names18[grep("PWID", names18)]
#browser()
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
#    rho.m.all <- numeric(n.gp) # entry rates for in-migration
#    rho.all <- numeric(n.gp) # entry rates


    # Set rho.m to 0 for 2016 and onwards
#    if (t > 48) rho.m <- rep(0, 18)


    ## TO UPDATE NEW ART ENGAGEMENT ##
    if (t <= 48) {
      phi1.t <- phi1[1, ]
      phi2.t <- phi2[1, ]
      phi3.t <- phi3[1, ]
      alpha1.t <- alpha1[1, ]
      alpha2.t <- alpha2[1, ]
      alpha3.t <- alpha3[1, ]
    } else if (t > 48 & t <= 60) {
      phi1.t <- phi1[2, ]
      phi2.t <- phi2[2, ]
      phi3.t <- phi3[2, ]
      alpha1.t <- alpha1[2, ]
      alpha2.t <- alpha2[2, ]
      alpha3.t <- alpha3[2, ]
    } else if (t > 60 & t <= 72) {
      phi1.t <- phi1[3, ]
      phi2.t <- phi2[3, ]
      phi3.t <- phi3[3, ]
      alpha1.t <- alpha1[3, ]
      alpha2.t <- alpha2[3, ]
      alpha3.t <- alpha3[3, ]
    } else if (t > 72) {
      phi1.t <- phi1[4, ]
      phi2.t <- phi2[4, ]
      phi3.t <- phi3[4, ]
      alpha1.t <- alpha1[4, ]
      alpha2.t <- alpha2[4, ]
      alpha3.t <- alpha3[4, ]
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
#      rho.m.all[ind] <- rho.m[i]
#      rho.all[ind] <- rho[i]
    }


    ### Time-varying number of OAT initiators ###
    # MWID (pop.pwid 1-3), PWID + M (pop.pwid 4-6), PWID + F
#    if (t <= 12) nOAT <- round(c(nOAT.m[1, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[1, 4:6]))
#    if (t > 12 & t <= 24) nOAT <- round(c(nOAT.m[2, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[2, 4:6]))
#    if (t > 24 & t <= 36) nOAT <- round(c(nOAT.m[3, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[3, 4:6]))
#    if (t > 36 & t <= 48) nOAT <- round(c(nOAT.m[4, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[4, 4:6]))
#    if (t > 48 & t <= 60) nOAT <- round(c(nOAT.m[5, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[5, 4:6]))
#    if (t > 60 & t <= 72) nOAT <- round(c(nOAT.m[6, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[6, 4:6]))
#    if (t > 72 & t <= 84) nOAT <- round(c(nOAT.m[7, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[7, 4:6]))
#    if (t > 84) nOAT <- round(c(nOAT.m[8, 1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[8, 4:6]))

    ############################
    # oat, intervention
#   if (intervention.oat.modif == TRUE & t > int.start) {
#      nOAT.hold <- nOAT
#      if (t <= int.end) {
        # intervention change the value of n.OAT
#        nOAT.v <- int.oat.eff[(t - int.start), ]
#        nOAT <- round(c(nOAT.v[1:3] * c(pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6]), 1 - pop.pwid[1:3] / (pop.pwid[1:3] + pop.pwid[4:6])), nOAT.v[4:6]))
#      } else {
#        nOAT <- nOAT.hold
#      }
#    }

    ############################
    # proportion of oat: number of OAT clients/PWID population
#    prop.oat <- nOAT / pop.pwid

    # oat.e: OAT entry rate for 9 groups (collapsing onOAT/offOAT, low/high)
#    if ((t >= int.start) & (t < int.end) & (intervention.oat.modif == TRUE)) {
      # adjustment to ensure (yearly) proportion on OAT among OUD (opioid use disorder) never exceed 95% (caped),
      # 72.7% OUD among PWID, applied in Ideal scenario
      # need to revisit this part
 #     if (any(prop.oat > 0.95 * 0.727)) {
 #        while (any(prop.oat > 0.95 * 0.727 + 0.00001)) {
 #         sat <- which(prop.oat >= 0.95 * 0.727)
 #         nOAT[-sat] <- nOAT[-sat] + sum(nOAT[sat] - pop.pwid[sat] * 0.95 * 0.727) * (pop.pwid[-sat] / sum(pop.pwid[-sat]))
 #         nOAT[sat] <- pop.pwid[sat] * 0.95 * 0.727
 #         prop.oat <- nOAT / pop.pwid
 #       }
 #     }
 #   }
#    oat.e <- prop.oat / (1 - prop.oat) * oat.q
    # entry rate cap 90%, revisit
#    oat.e[which(oat.e > 0.9 | oat.e <= 0)] <- 0.9
#browser(expr = {t >= 97}) # 97 is when intervention start in accumulate.outcomes.combination, needed to get the comparable nOAT (& prop.oat)
   if (t <= 12)          prop.oat.t <- prop.oat.m[1,] # 2017 Dec and before
   if (t > 12 & t <= 24) prop.oat.t <- prop.oat.m[2,] # 2018 Jan - Dec
   if (t > 24 & t <= 36) prop.oat.t <- prop.oat.m[3,] # 2019 Jan - Dec
   if (t > 36 & t <= 48) prop.oat.t <- prop.oat.m[4,] # 2020 Jan - Dec
   if (t > 48)           prop.oat.t <- prop.oat.m[5,] # 2021 Jan and later

## END UPDATING NEW PrEP ##

# PrEP, intervention - Move BEFORE to ensure prep.total[8], thus prep.total.t value is update with PrEP intervention
# 2020 onward (t > int.start=96)
    
if (intervention.oat.modif == TRUE & t > int.start) {
	prop.oat.hold <- prop.oat.t
	if (t <= int.end) {
		prop.oat.t <- int.oat.eff[t - int.start,]
	} else {
		prop.oat.t  <- prop.oat.hold
	}
}

    ### Time-varying number of SSP ###
#   if (length(v.ssp) != 1) {
#      if (t <= 48) v.ssp.t <- v.ssp[1]
#      if (t > 48 & t <= 60) v.ssp.t <- v.ssp[2]
#      if (t > 60 & t <= 72) v.ssp.t <- v.ssp[3]
#      if (t > 72) v.ssp.t <- v.ssp[4]
#    } else {
#      v.ssp.t <- v.ssp
#    }
     
    # NEW: Use Proportion on SSP instead of syringe through SSP
		if (t <= 12)          prop.ssp.t <- prop.ssp.m[1,] # 2017 Dec and before
		if (t > 12 & t <= 24) prop.ssp.t <- prop.ssp.m[2,] # 2018 Jan - Dec
		if (t > 24 & t <= 36) prop.ssp.t <- prop.ssp.m[3,] # 2019 Jan - Dec
		if (t > 36 & t <= 48) prop.ssp.t <- prop.ssp.m[4,] # 2020 Jan - Dec
		if (t > 48)           prop.ssp.t <- prop.ssp.m[5,] # 2021 Jan and later
   
    ############################
    # SSP, intervention
#    if (intervention.ssp.modif == TRUE & t > int.start) {
#      v.ssp.hold <- v.ssp.t
#      if (t <= int.end) {
#        v.ssp.t <- int.ssp.eff[t - int.start]
#      } else {
#        v.ssp.t <- v.ssp.hold
#      }
#    }
 
    if (intervention.ssp.modif == TRUE & t > int.start) {
    	prop.ssp.hold <- prop.ssp.t
    	if (t <= int.end) {
    		prop.ssp.t <- int.ssp.eff[t - int.start,]
    	} else {
    		prop.ssp.t <- prop.ssp.hold
    	}
    }
    ############################
    # ssp coverage: #syringe distributed/(PWID population*d)
#    cov.ssp <- v.ssp.t / (sum(y2[all.pwid, ]) * d * 12)

    # set the maximum coverage as 1, revisit
 #   if (cov.ssp > 1) {
 #     cov.ssp <- 1
 #   }

    # assign oat entry rate and coverage of ssp for all n.gp groups
#   oat.e.all <- numeric(n.gp)
#   cov.ssp.all <- numeric(n.gp)
    s.all <- numeric(n.gp)
#    prop.oat.all <- numeric(n.gp)

    # names.pwid: pwid group names without OAT, low, high (among names18), 9vp
    for (i in seq_along(names.pwid)) {
      ind <- which(rname %in% names.pwid[i])
#     oat.e.all[ind] <- oat.e[i]
#     cov.ssp.all[ind] <- cov.ssp    
      s.all[ind] <- s[i]
#      prop.oat.all[ind] <- prop.oat[i]
    }
    
    # Weighted average impact of OAT: Reduce mortality
#    mo <- mo * vparameters$oat_mor * prop.oat.all + mo * (1-prop.oat.all)
     mo <- mo * vparameters$oat_mor * prop.oat.t + mo * (1-prop.oat.t)

if (sum(y[, 10], na.rm = TRUE) / sum(y, na.rm = TRUE) > 0.9) {
  cat("⚠️ sum(y[,10]) / sum(y) > 0.9 — entering debug mode\n")
  browser() 
}     
#browser()    
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

    # eta: time-varying PrEP entry rate
#    eta <- numeric(30)
#    if (t <= 6) {
#      eta[young] <- 0
#    } else if (t > 6 & t <= 12) {
#      eta[young] <- rep(-log(1 - prep.total.t * prep.proportion / msm.h.scep) / 6, 2)
#    } else if (t > 12 & t <= 48) {
#      eta[young] <- rep(-log(1 - prep.total.t * prep.proportion / msm.h.scep) / 12, 2)
#    } else if (t > 48) { ## Since 2016, modeling PrEP by proportions
#    	prep.vol  <- prep.total.t * prep.proportion               # no. of people on prep by race
#    	prop.prep <- prep.total.t * prep.proportion / msm.h.scep # proportion of msm high risk susceptible on PrEP by race
      # adjustment to ensure proportion on PrEP never exceed 92% of high-risk MSM in each pop group, applied in Ideal scenario
#      while (any(prop.prep > 0.92001)) {
#        sat <- which(prop.prep >= 0.92)
#        prep.vol[-sat] <- prep.vol[-sat] + sum(prep.vol[sat] - msm.h.scep[sat] * 0.92) * (msm.h.scep[-sat] / sum(msm.h.scep[-sat]))
#        prep.vol[sat] <- msm.h.scep[sat] * 0.92
#        prop.prep <- prep.vol / msm.h.scep
#      }
#      eta[young] <- rep(prop.prep / (1 - prop.prep) * wp, 2)
#browser() #get prop.prep
      # NA check and print out when it is happened.
#      if (anyNA(eta)) warning(sprintf("Current time period = \"%s\". \n", t))
#      if (anyNA(eta)) {
#        warning(sprintf(
#          "PrEP coverage by race/ethnicity = \"%s\". \n",
#          prep.total[8] * prep.proportion / msm.h.scep
#        ))
#      }
#      if (anyNA(eta)) stop(sprintf(">100% PrEP coverage. Current time period = \"%s\". \n", i))
      # eta cap at 0.9. Why?
#      eta[which(is.na(eta))] <- 0.9
#      eta[eta > 0.9] <- 0.9
#    }
    

    
#browser(expr = {t > 0})
    # psi: time-varying screening rate
#    if (t <= 12) psi <- psi.m[1, ]
#    if (t > 12 & t <= 24) psi <- psi.m[2, ]
#    if (t > 24 & t <= 36) psi <- psi.m[3, ]
#    if (t > 36) psi <- psi.m[4, ]
    ############################
    # psi (screening rate), intervention
#    if (intervention.testing.modif == TRUE & t > int.start) {
#   	psi.hold <- psi
#    	if (t <= int.end) {
#    		psi <- psi * int.test.eff[(t - int.start), ]
#    		psi[psi > 0.99] <- 0.99
#    	} else {
#    		psi <- psi.hold
#    	}
#    }
    
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


#browser(expr = {t >= int.start})       
    ############################
    # psi (screening rate), intervention, increase baseline psi
    if (intervention.testing.modif == TRUE & t > int.start) {
      psi.hold <- psi
      if (t <= int.end) {
        psi <- psi * t(int.test.eff[(t - int.start), ,])
#       psi[psi > 0.99] <- 0.99
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
#    theta.t1O <- theta.t1O * theta.o.oat * prop.oat.all + theta.t1O * (1 - prop.oat.all) 
#    theta.t2O <- theta.t2O * theta.o.oat * prop.oat.all + theta.t2O * (1 - prop.oat.all) 
#    theta.t3O <- theta.t3O * theta.o.oat * prop.oat.all + theta.t3O * (1 - prop.oat.all) 
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
#browser(expr = {t >= 210})    	
    	# allocate change depending on direction - If gap = 0 (the case for t=0), assign pop change to each risk group equally; otherwise allocate change based on gap/difference (model vs target pop) size
    	if (race_total_change > 0) {
    		pos_gap <- pmax(risk_gap, 0)
    		if (sum(pos_gap) < 1e-12) { #Ignore noise
    			weights <- rep(1/length(pos_gap), length(pos_gap))
    		} else {
#browser()    			
    			weights <- pos_gap / sum(pos_gap)
    		}
    	} else {
    		neg_gap <- pmax(-risk_gap, 0)
    		if (sum(neg_gap) < 1e-12) {
    			weights <- rep(1/length(neg_gap), length(neg_gap))
    		} else {
#browser()    			
    			weights <- neg_gap / sum(neg_gap)
    		}
    	}
    	
    	risk_change <- race_total_change * weights  # Step 4: Allocate the race's total change across its risks according to risk gap (i.e., possible for W/Het to make up for W/PWID, but not B/Het or B/PWID)
#browser(expr = {t >= 210})    	# Note when race_total_change > 0; We are adding people to this race, thus We should only add to underrepresented risks, vice versa. If we remove people for this race race_total_change < 0, We should only remove from overrepresented risks.
    	# distribute inside each (race, risk)
    	for (i in seq_along(risk_list)) {
    		
    		idx_cell <- intersect(idx_race, risk_list[[i]])
    		if (length(idx_cell) == 0) next
  		
    		inside_prop <- target.pop[idx_cell] / sum(target.pop[idx_cell]) # This step indirectly ensure that proportion for msm.Y & msm.O (in MSM), pwid.f & pwid.m (in pwid), het.f & het.m (in Het) is also mantain
    		pop.change.18.new[idx_cell] <- risk_change[i] * inside_prop
#browser(expr = {t >= 210})      		
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
 #browser()   
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
# browser(expr = {t == 200})
    # Distribute pop.change across high-young vs low-old
    pop.change.30[high] = pop.change.30[high] *as.numeric(prop.high.sus.all[high])
    pop.change.30[low]  = pop.change.30[low]   *(1-as.numeric(prop.high.sus.all[low]))
# browser(expr = {t == 200})
    pop.change.30[old]  = 0 # All msm pop.change got to young
    
    # Adjustment to make sure high-young stay const over time
    tot.group <- rowSums(y2)

    # Adjustment to make sure high-young stay const over time
    dif.het.high  <- (tot.group[low] + tot.group[high]) * as.numeric(prop.high.sus.all[high]) - tot.group[high]
    pop.change.30[high]  <- pop.change.30[high]  + dif.het.high
    pop.change.30[low]   <- pop.change.30[low] - dif.het.high
    
    
#browser(expr = {t > 200}) #Check sum(pop.change.30) = sum(pop.change.18)

############################
# This code session is backbone for the "out" 30*21 result matrix
    # out is the population difference between time t and t+1
    out <- matrix(0, n.gp, 28)
#browser()
		# Code below, populate the 30*24 matrix "out", 42 rows are indexed by [i, ] representing risk groups, determine if patient on OAT and which ODE to use
    # PWID or HET
    for (i in c(pwid, het)) {
#cat("Start Not MSM: i =",i, "\n")
#browser(expr = {i = 30})     	#Test index pwid, het
      out[i, ] <- ode.list(
        y[i, ], foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
#       rho.all[i],mu_mat[i], rho.m.all[i],  
#				ws, wp, eta[i], # eta for PrEP entry
				mu_adj[i], ws, psi[i,],
#       psi[i], psi.p, v2, v3,
        theta.ai, theta.ad,
        phi[i, ], alpha[i, ], alpha.re[i],
        theta.1, theta.2,
        theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
        theta.t31[i], theta.t32[i],
        theta.t1O[i], theta.t2O[i], theta.t3O[i],
        mo[i, ], pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
      )
#cat("End Not MSM: i =",i, "\n")
    }
#browser(expr = {t >= 210})
#for (i in c(young)) {
for (j in seq_along(young)) {
		i <- young[j]	
#browser(expr = {i = 4}) 
  	out[i, ] <- ode.list.MSMyoung(
		y[i, ], y[old[j], ],
		foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
#		rho.all[i], mu_mat[i], rho.m.all[i], 
#		ws, wp, eta[i], # eta for PrEP entry
		mu_mat[i], mu_adj[i], ws, psi[i,], # mu_afj not needed for young
#   psi[i], psi.p, v2, v3, 
    theta.ai, theta.ad,
		phi[i, ], alpha[i, ], alpha.re[i],
		theta.1, theta.2,
		theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
		theta.t31[i], theta.t32[i],
    theta.t1O[i], theta.t2O[i], theta.t3O[i], 
    mo[i, ], mo[old[j], ], old.e.sus[i], old.e.hiv[i],
		pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
	)
#	cat("End MSM/young: i =",i, "\n")
}
#browser(expr = {t >= 210})
for (j in seq_along(old)) {
	  i <- old[j]
#browser()
#browser(expr = {i = 1})	#check old_mor
	  out[i, ] <- ode.list.MSMold(
		y[i, ], y[young[j], ], 
		foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
#		rho.all[i], mu_mat[i], rho.m.all[i],
#   ws, wp, eta[i], # eta for PrEP entry
		mu_mat[i], mu_adj[i], ws, psi[i,],
#   psi[i], psi.p, v2, v3, 
    theta.ai, theta.ad,
		phi[i, ], alpha[i, ], alpha.re[i],
		theta.1, theta.2,
		theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
		theta.t31[i], theta.t32[i],
		theta.t1O[i], theta.t2O[i], theta.t3O[i], 
    mo[i, ], old.e.sus[i], old.e.hiv[i],
		pop.change.30[i], prop.prep.t[i], prop_vs_T[i]
	)
#	cat("End MSM/old: i =",i, "\n")
}
    # MSM/PWID NOT on OAT (to add OAT entry/dropout)
    # Why above use for_in, whereas with OAT, we use for_in_seq_along? Need j to index both oat and off.oat at the same time y[i, ], y[oat[j], ] where i <- off.oat[j]
#    for (j in seq_along(oat)) {
      # indices for off.oat group, off.oat
#      i <- off.oat[j]
#      out[i, ] <- ode.list.offOAT(
#        y[i, ], y[oat[j], ],
#        foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
#        rho.all[i], ws, wp, mu_mat[i], mo[i, ], eta[i], # eta for PrEP entry
#        psi[i], psi.p, theta.ai, theta.ad,
#        phi[i, ], v2, v3, alpha[i, ], alpha.re[i],
#        theta.1, theta.2,
#        theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
#        theta.t31[i], theta.t32[i],
#        theta.t1O[i], theta.t2O[i], theta.t3O[i],
#        oat.e.all[i], oat.q, rho.m.all[i]
#      )
#cat("Off OAT: i =",i,"j =",j, "\n")
#    }
#browser()
    # MSM/PWID on OAT (to add OAT entry/dropout)
#    for (j in seq_along(oat)) {
      # indices for on oat group, oat
#      i <- oat[j]
#      out[i, ] <- ode.list.OAT(
#        y[i, ], y[off.oat[j], ],
#        foi$beta_o[i, ], foi$beta_s[i, ], foi$gamma[i, ],
#        rho.all[i], ws, wp, mu_mat[i], mo[i, ], eta[i], # eta for PrEP entry
#        psi[i], psi.p, theta.ai, theta.ad,
#        phi[i, ], v2, v3, alpha[i, ], alpha.re[i],
#        theta.1, theta.2,
#       theta.t12[i], theta.t13[i], theta.t21[i], theta.t23[i],
#        theta.t31[i], theta.t32[i],
#        theta.t1O[i] * theta.o.oat, theta.t2O[i] * theta.o.oat, theta.t3O[i] * theta.o.oat,
#        oat.e.all[i], oat.q, rho.m.all[i]
#      )
#cat("On OAT: i =",i,"j =",j, "\n")
#    }
#browser(expr = {t >= 210}) # Fill Pop_entry and Save_matrix to check output

#############################################
    # revisit the population adjustment part
    if (prop.adj == TRUE) {
# adjust susceptible population S1 to maintain constant proportion of risk groups among male and female. Why only S1? cause only new people into pop can start here?
      # total population for each group
      tot.group <- rowSums(y2)

      # col S1,S2,Sp
      dif <- matrix(0, nrow = 42, ncol = 3)

      # adjust for risk group proportion among male population, using the initial pop
      # Calculate the initial proportion (to be kept constant) of each M risk group within total M population (init.group.prop <- as.vector(rowSums(init)/sum(init))
      msm.prop <- sum(init.group.prop[msm]) / sum(init.group.prop[male]) 
      mwid.prop <- sum(init.group.prop[mwid]) / sum(init.group.prop[male]) # MSM and PWID
      m_pwid.prop <- sum(init.group.prop[pwid.m]) / sum(init.group.prop[male]) # HET and PWID
      m_het.prop <- sum(init.group.prop[het.m]) / sum(init.group.prop[male]) # Het male

      # Difference between the population if CONSTANT proportion was applied (sum(tot.group[male]) * ....prop) 
      # VS current calculated value (sum(tot.group[...]), obtained from y2: tot.group <-rowsum(y2)) for each risk group
      dif.msm.total <- sum(tot.group[male]) * msm.prop - sum(tot.group[msm])
      dif.mwid.total <- sum(tot.group[male]) * mwid.prop - sum(tot.group[mwid])
      dif.m_pwid.total <- sum(tot.group[male]) * m_pwid.prop - sum(tot.group[pwid.m])
      dif.m_het.total <- sum(tot.group[male]) * m_het.prop - sum(tot.group[het.m])

      # midu_off = mwid_off idu.m_off = pwid.m_off - Getting the index value
      mwid_off <- intersect(mwid, off.oat)
      pwid.m_off <- intersect(pwid.m, off.oat)
#browser() # Check t, dif[msm, 1], ... diff
      # Index y2 by risk group stratification, divided it by sum of all stratification to get what I assume to be the DISTRIBUTION of dif...total among all states (columns) and risk stratification (rows) (y2[msm, ] / sum(tot.group[msm])) - Dimension: risk stratifications (6 for MSM) * 19 states ;Check: Sum of that =1
			# Using this DISTRIBUTION to distribute the difference in population among all states (columns) and risk stratification (rows): dif.msm.total * (y2[msm, ] / sum(tot.group[msm])) -> 1* matrix = matrix (with same dimension).
			# Using rowSums to sum values across 19 states into vector with length of the risk stratification (MSM is 6) --> Required to fit the length of the corresponding subsect dif[,]
			# Using the value of the vector to update the corresponded value of the dif matrix index by risk group, dif[msm, 1] refer to row S1, index by MSM
      dif[msm, 1] <- dif[msm, 1] + (rowSums(dif.msm.total * (y2[msm, ] / sum(tot.group[msm])))) # Populate dif matrix (row index by MSM)
      dif[mwid_off, 1] <- dif[mwid_off, 1] + (rowSums(dif.mwid.total * (y2[mwid_off, ] / sum(tot.group[mwid_off]))))
      dif[pwid.m_off, 1] <- dif[pwid.m_off, 1] + (rowSums(dif.m_pwid.total * (y2[pwid.m_off, ] / sum(tot.group[pwid.m_off])))) # when t=0, only this is differ from 0
      dif[het.m, 1] <- dif[het.m, 1] + (rowSums(dif.m_het.total * (y2[het.m, ] / sum(tot.group[het.m]))))
#browser() # Check dif[msm, 1], ... diff
      # adjust for risk group proportion among female population#
      f_pwid.prop <- sum(init.group.prop[pwid.f]) / sum(init.group.prop[female])
      f_het.prop <- sum(init.group.prop[het.f]) / sum(init.group.prop[female])

      dif.f_pwid.total <- sum(tot.group[female]) * f_pwid.prop - sum(tot.group[pwid.f])
      dif.f_het.total <- sum(tot.group[female]) * f_het.prop - sum(tot.group[het.f])

      # idu.f_off = pwid.f_off
      pwid.f_off <- intersect(pwid.f, off.oat)
#browser() # Check dif[pwid.f_off, 1], ... diff     
			#Why index y2[pwid.f_off, 1:3]? Not add up to 1 like y2[pwid.f_off, ]
			# y2[, 1:3] would mean we only care about # people in 3 susceptible states, but why here an not above. After all, pwid.f_off ~ pwid.m_off,het.f~het.m
      dif[pwid.f_off, 1] <- dif[pwid.f_off, 1] + (rowSums(dif.f_pwid.total * (y2[pwid.f_off, 1:3] / sum(tot.group[pwid.f_off]))))
      dif[het.f, 1] <- dif[het.f, 1] + (rowSums(dif.f_het.total * (y2[het.f, 1:3] / sum(tot.group[het.f]))))

# When making these adjustments, it seems that it not that we "balance" by taking from group to add to the other 
# (i.e, to make the total population before and after adjustment the same, i.e., 1000 people before adj, will still be 1000 after adj)
# Rather, we "top up" the group (using dif) in order to keep the proportion constant, but at the expense of keeping the population constant   
# (i.e., 1000 bef adj will bow be 1005 after adj, with 5 divided among dif groups to ensure const prop)
# Though rather it can be the other way, since we are "topping of" 'out' representing the "change" population (i.e., dS1, rather than S1)
# Since the change here is often a reduction, we actually "remove" people (1000 before adj, now 995 after adj), and vice versa
# This is corroborated by the fact that I can see predicted # of people in S1 predicted from ODE alone (w/o adj) is always larger than 
# the number after adj (54502.788 vs 54499.3405), as a consequent of dS1 after adjust is always larger.(note this is for [1,1] other can be different)
      

# adjust to have constant high and low risk proportion among susceptible, adjustment made to S1
      # adjust for high-risk % among MSM population
      tot.group.sus <- rowSums(y2[, 1:3])

      # high risk 25%. Why? Assume that sum of these 2 groups give # of allmsm, why multiply this sum with 25% and substract the H.allmsms?
      dif.msm <- (tot.group.sus[old] + tot.group.sus[young]) * 0.25 - tot.group.sus[young]

      ## CHANGED ADDING MSM.H / REMOVING MSM.L ONLY TO S1 ##
      # So increase the change ( remember we modifying 'out') in High risk, and decrease the chance in low risk by the same amount
      # Also, some of the index between H.allmsm, L.allmsm overlap with previous groups, so dif[] will be overwritten?
      dif[young, 1] <- dif[young, 1] + dif.msm
      dif[old, 1] <- dif[old, 1] - dif.msm

      # adjust for high-risk % among HET population het.h = H.het
      H.het <- setdiff(het, L.het)
      dif.het <- (tot.group.sus[L.het] + tot.group.sus[H.het]) * prop.high.sus[13:18] - tot.group.sus[H.het]
      
      ## CHANGED ADDING H.het / REMOVING L.het ONLY TO S1 ##
      dif[H.het, 1] <- dif[H.het, 1] + dif.het
      dif[L.het, 1] <- dif[L.het, 1] - dif.het

      # change only made to S1.  can just define a vector size 42 by 1, revisit.
#browser() #Check this browser (to the prev one) can see the value diff being overwritten, this is done without account for previous balance between const prop, sig change with L. & H.allmsms after adj
      out[, 1:3] <- out[, 1:3] + dif
#browser() #check t, dif (42x first 3 columns), check out (and out[1,1]), check y (and y[1,1])
    }

    if (inf.prop.adj == TRUE) {
# adjust infected population to keep high-risk % among MSM and HET constant among PLHIV
      # adjustment all made within gender and ethnicity

      tot.group <- rowSums(y2[, 4:19])

      # 16 states for PLHIV
      dif <- matrix(0, nrow = 42, ncol = 16)

      # adjust for high-risk % among MSM population
      dif.msm <- (tot.group[old] + tot.group[young]) * c(prop.high.inf[1:6], prop.high.inf[4:6]) - tot.group[young]
      dif[young, ] <- dif[young, ] + (dif.msm * (y2[, 4:19][young, ] / rowSums(y2[, 4:19][young, ]))) # Difference between here and previous code is that, we spread the chance across 15 states (4:19) instead of just 1
      dif[old, ] <- dif[old, ] - (dif.msm * (y2[, 4:19][old, ] / rowSums(y2[, 4:19][old, ]))) # Similar to first adjustment, distribute *by multiply) dif.msm over a 9x16 matrix, which itself results from divided a 9x16 matrix by a vector of 9 (i.e., rowsums) 
#browser()
      # adjust for high-risk % among HET population
      H.het <- setdiff(het, L.het)
      dif.het <- (tot.group[L.het] + tot.group[H.het]) * prop.high.inf[13:18] - tot.group[H.het]
      dif[H.het, ] <- dif[H.het, ] + (dif.het * (y2[H.het, 4:19] / rowSums(y2[H.het, 4:19])))
      dif[L.het, ] <- dif[L.het, ] - (dif.het * (y2[L.het, 4:19] / rowSums(y2[L.het, 4:19])))

      out[, 4:19] <- out[, 4:19] + dif
    }

#cat("End Function-ode_model-Combination", t ,"\n") #Delete
# This ODE function will be called by euler for each time point t (also the index that was passed through), i.e., the number of vt
    
    list(as.vector(t(out)))
   
  })
}
