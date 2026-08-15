##################################################################################
# An R script to specify the model ode equations for non PWID (HET and MSM)
# Details see supplement of P12
# Updated: Oct 14, 2022
##################################################################################

## ODE parameters ##
# x: a vector containing the number of individuals in each HIV state
#   x=c(S1,S2,Sp,Ia,I1,I2,I3,Iap,Ip,Da,D1,D2,D3,T1,T2,T3,O1,O2,O3)
# lambda: the sufficient contact rate
# lambda.p: the sufficient contact rate for those on PrEP.
# rho: the entry rate (assume the same for all groups)
# 	rho.m: population in-migration rate
# 1/ws: average duration uninfected individuals in compartment S remain
#   identified after screening
# 1/wp: average duration uninfected individuals in compartment Sp remain on PrEP
# mu_mat: maturation rates (age>65)
# mo: mortality rates for 19 states
# psi: screening rates for susceptible or infected
# psi.p: screening rates for Iap (acute HIV on PrEP)
# theta.ai: transition rate from acute states (Ia) to chronic state (CD4>=500; I1)
# theta.ad: transition rate from acute states (Da) to chronic state (CD4>=500; D1)
# phi: % infected of people receiving ART once diagnosed
# v2,v3: symptom-based case finding rate for I2 and I3, respectively
# alpha: ART initiation rate for D1,D2, and D3
# alpha.re: ART re-initiation rate for O1,O2, and O3
# theta.1: HIV disease progression rate for those not on ART, from I1/D1 to I2/D2
# theta.2: HIV disease progression rate for those not on ART, from I2/D2 to I3/D3
# theta.t: transition probabilities for those on ART, theta.t12,theta.t12, theta.t21,theta.t23, theta.t31, theta.t32
# theta.o: ART dropout probability from states T1,T2,T3, theta.t1O, theta.t2O, theta.t3O
# eta: PrEP entry rate

#' Function-ode.list
#'
#' \code{Function-ode.list} Define ODEs for MSM and Heterosexual
#'
#' @param x vector containing the number of individuals in each HIV state = c(S1,S2,Sp,Ia,I1,I2,I3,Iap,Ip,Da,D1,D2,D3,T1,T2,T3,O1,O2,O3)
#' @param beta_o Sufficient contact rate for individuals involved in opposite sex
#' @param beta_s Sufficient contact rate for individuals involved in same sex
#' @param gamma Sufficient contact rate for individuals involved in needle sharing
#' @param rho population entry rate
#' @param ws average duration (1/ws) uninfected individuals in compartment S remain identified after screening
#' @param wp average duration (1/wp) uninfected individuals in compartment Sp remain identified on PrEP
#' @param mu_mat 
#' @param mo mortality rates for 19 states
#' @param eta PrEP entry rate
#' @param psi screening (testing?) rates for susceptible or infected
#' @param psi.p screening (testing?) rates for infected-acute on PrEP
#' @param theta.ai transition rate from acute states (Ia) to chronic undiagnosed (I1)
#' @param theta.ad transition rate from acute states (Id) to chronic diagnosed (D1)
#' @param phi proportion of infected receiving ART within first month of diagnosis
#' @param v2 symptom-based case finding rate for I2
#' @param v3 symptom-based case finding rate for I3
#' @param alpha ART initiation rate (post 30 days) for D1, D2, and D3
#' @param alpha.re ART re-initiation rate for O1, O2, and O3
#' @param theta.1 HIV disease progression rate for those not on ART, from I1/D1 to I2/D2
#' @param theta.2 HIV disease progression rate for those not on ART, from I2/D2 to I3/D3
#' @param theta.t12 transition rate(?) for those on ART T1 -> T2
#' @param theta.t13 transition rate(?) for those on ART T1 -> T3
#' @param theta.t21 transition rate(?) for those on ART T2 -> T1
#' @param theta.t23 transition rate(?) for those on ART T2 -> T3
#' @param theta.t31 transition rate(?) for those on ART T3 -> T1
#' @param theta.t32 transition rate(?) for those on ART T3 -> T2
#' @param theta.t1O transition rate(?) for those on ART T1 -> off
#' @param theta.t2O transition rate(?) for those on ART T2 -> off
#' @param theta.t3O transition rate(?) for those on ART T3 -> off
#' @param rho population entry rate
#' @param rho.m population in-migration rate (HIV+??)
#' @return
#' out: a vector of health state occupancy rate of change, HIV incidence, and diagnoses
#' @export

ode.list.MSMyoung <- function(x, x.old, beta_o, beta_s, gamma, mu_mat.old, mu_adj, ws, psi,
															#psi.p, v2, v3, 
															theta.ai, theta.ad, phi,
															alpha, alpha.re, theta.1, theta.2, theta.t12, theta.t13,
															theta.t21, theta.t23, theta.t31, theta.t32,
															theta.t1O, theta.t2O, theta.t3O, 
															mo, mo.old ,old.e.sus, old.e.hiv, pop.change.30, prop.prep.t,prop_vs_T) {
  # Should consider indexing with names instead of numbers
  # x is a vector containing the number in each health states at time t-1 (correct ?)
  # Extract individual numbers for each health state (later apply transition rates to these numbers)
  S1 <- x[1]
  S2 <- x[2]
  Ia <- x[3]
  I1 <- x[4]
  I2 <- x[5]
  I3 <- x[6]
  Da <- x[7]
  D1 <- x[8]
  D2 <- x[9]
  D3 <- x[10]
  T1 <- x[11]
  T2 <- x[12]
  T3 <- x[13]
  O1 <- x[14]
  O2 <- x[15]
  O3 <- x[16]

  # ART initiation and re-initiation
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  alpha.re1 <- alpha.re2 <- alpha.re3 <- alpha.re

  # defined on page 4 of supplement of p12, o - opposite sex, s- same sex
  # Sufficient contact rate defined in 'CascadeCEA-Model-0-Function-FoI'
  # Combination of transmission via hetero- and same-sex contact and injection
  # Matrix of dim (n.gp, 2) - 2 for PrEP and non-PrEP
#  lambda <- sum(beta_o[1], beta_s[1], gamma[1]) # non-PrEP
#  lambda.p <- sum(beta_o[2], beta_s[2], gamma[2]) # PrEP
  lambda <- sum(beta_o[1], beta_s[1], gamma[1]) * (1-prop.prep.t) + sum(beta_o[2], beta_s[2], gamma[2]) * prop.prep.t

 
  ### 
 #dS1 <- pop.change.30/12 - psi[1] * S1  + ws * S2 - lambda * S1 - mo[1] * S1 - old.e[1] * S1
  dS1 <- pop.change.30/12 - psi[1] * S1  + ws * S2 - lambda * S1              - old.e.sus * S1 + (mu_mat.old + mo.old[1]) * x.old[1] + (mu_mat.old + mo.old[2]) * x.old[2] # mo.old is needed same as x.old, as it is diffrent, mu_mat are the same across young and old so it is nor necessary
#browser()
  # 
 #dS2 <- psi[1] * S1 - ws * S2 - lambda * S2 - mo[2] * S2 - old.e[1] * S2
  dS2 <- psi[1] * S1 - ws * S2 - lambda * S2              - old.e.sus * S2  
  # 
  #dSp <- -wp * Sp - lambda.p * Sp - mo[3] * Sp + eta * (S2 + S1) - old.e[1] * Sp
#browser()  
  # 
  dIa <- lambda * (S1 + S2) - psi[3] * Ia - theta.ai * Ia - mo[3] * Ia - old.e.hiv * Ia
  # 
  dI1 <- theta.ai * Ia - psi[4] * I1 - theta.1 * I1 - mo[4] * I1 - old.e.hiv * I1 
  # 
  dI2 <- theta.1 * I1 - psi[5] * I2 - theta.2 * I2 - mo[5] * I2 - old.e.hiv * I2 
  # 
  dI3 <- theta.2 * I2 - psi[6] * I3 - mo[6] * I3 - old.e.hiv * I3
  # 
  #dIap <- lambda.p * Sp - psi.p * Iap - theta.ai * Iap - mo[8] * Iap - old.e[2] * Iap
  # 
  #dIp <- theta.ai * Iap - psi.p * Ip - mo[9] * Ip - old.e[2] * Ip
  
  # 
  dDa <- psi[3] * Ia - theta.ad * Da - mo[7] * Da - old.e.hiv * Da
  # 
  dD1 <- theta.ad * Da + psi[4] * (1 - phi[1]) * I1 - theta.1 * D1 - alpha1 * D1 - mo[8] * D1 - old.e.hiv * D1 
  # 
  dD2 <- theta.1 * D1 + psi[5] * (1 - phi[2]) * I2 - theta.2 * D2 - alpha2 * D2 - mo[9] * D2 - old.e.hiv * D2 
  # 
  dD3 <- theta.2 * D2 + psi[6] * (1 - phi[3]) * I3 - alpha3 * D3 - mo[10] * D3 - old.e.hiv * D3 
  
  # 
  dT1 <- theta.t21 * T2 + theta.t31 * T3 + alpha1 * D1 + alpha.re1 * O1 + psi[4] * phi[1] * I1 -
  					(theta.t12 + theta.t13 + theta.t1O) * T1 - mo[11] * T1 - old.e.hiv * T1 
  # 
  dT2 <- theta.t12 * T1 + theta.t32 * T3 + alpha2 * D2 + alpha.re2 * O2 + psi[5] * phi[2] * I2 -
  					(theta.t21 + theta.t23 + theta.t2O) * T2 - mo[12] * T2 - old.e.hiv * T2 
  # 
  dT3 <- theta.t13 * T1 + theta.t23 * T2 + alpha3 * D3 + alpha.re3 * O3 + psi[6] * phi[3] * I3 -
  					(theta.t31 + theta.t32 + theta.t3O) * T3 - mo[13] * T3 - old.e.hiv * T3 
  
  # Transition to O1 (off ART CD4 > 500)
  dO1 <- theta.t1O * T1 - theta.1 * O1 - alpha.re1 * O1 - mo[14] * O1 - old.e.hiv * O1 
  # Transition to O2 (off ART CD4 200-500)
  dO2 <- theta.t2O * T2 + theta.1 * O1 - (alpha.re2 + theta.2) * O2 - mo[15] * O2 - old.e.hiv * O2 
  # Transition to O3 (off ART CD4 < 200)
  dO3 <- theta.t3O * T3 + theta.2 * O2 - alpha.re3 * O3 - mo[16] * O3 - old.e.hiv * O3 

  # Overall incidence from heterosexual contact (susceptible moving to infected)
  # inc_bo <- beta_o[1] * (S1 + S2) + beta_o[2] * Sp  
  inc_bo <- beta_o[1] * (S1 + S2) * (1-prop.prep.t) + beta_o[2] * (S1 + S2) * prop.prep.t
  # Overall incidence from homosexual contact (susceptible moving to infected)
  # inc_bs <- beta_s[1] * (S1 + S2) + beta_s[2] * Sp 
  inc_bs <- beta_s[1] * (S1 + S2) * (1-prop.prep.t) + beta_s[2] * (S1 + S2) * prop.prep.t
  # Overall incidence from injection (susceptible moving to infected)
  # inc_g <- gamma[1] * (S1 + S2) + gamma[2] * Sp  
  inc_g <- gamma[1] * (S1 + S2) * (1-prop.prep.t)+ gamma[2] * (S1 + S2) * prop.prep.t
  
  # (Number) Yearly Incidence/ new diagnoses among infected
  #diagm <- psi * (Ia + I1 + I2 + I3) + psi.p * (Iap + Ip) + v2 * I2 + v3 * I3
  #diagm <- psi * (Ia + I1 + I2 + I3) + v2 * I2 + v3 * I3
  diag <- psi[3] * x[3] + psi[4] * x[4] + psi[5] * x[5] + psi[6] * x[6] 
  # mortality among PLHIV (D,T,O not I)
  death <- sum(mo[7:16] * x[7:16])
  diag_late <- psi[6] * x[6]
  D_T <- alpha1 * D1 + alpha2 * D2 + alpha3 * D3
  I_T <- psi[4] * phi[1] * I1 + psi[5] * phi[2] * I2 + psi[6] * phi[3] * I3
  I_D <- psi[3] * Ia + psi[4] * (1 - phi[1]) * I1 +  psi[5] * (1 - phi[2]) * I2 + psi[6] * (1 - phi[3]) * I3
  death_D	<- sum(mo[7:10] * x[7:10])
  
  # (Proportion)
  positive_test <- diag/ (psi[1] * S1 + diag) #Assume no false positive among S and no false negative among I
  viral_sup <- sum(x[11:12])*prop_vs_T/sum(x[7:16]) #T1 / (Da + D1 + D2 + D3 + T1 + T2 + T3 + O1 + O2 + O3)
  
  # Vector of overall net monthly rates of change for each health state at time t, HIV incidence, HIV diagnoses and deaths
  out <- c(
  	dS1, dS2, dIa, dI1, dI2, dI3, dDa, dD1, dD2, dD3, dT1, dT2, dT3, dO1, dO2, dO3,
  	inc_bo, inc_bs, inc_g, diag, diag_late, positive_test, viral_sup, death, D_T, I_T,I_D, death_D
  )
  # the output should be the same order as the initial values
  # TO-DO: consider indexing by name instead of number
#browser() # Check out (lot of negative), Get diagm, get total infected (Ia + I1 + I2 + I3) - result is for a stratification (42 of them)
  ## adjustment for negative cell(?) - out is the change right, no negative means that number in the health states cannot increase?
  ## need to revisit later
  out[1:16][which((out[1:16] + x[1:16]) < 0)] <- -x[1:16][which((out[1:16] + x[1:16]) < 0)]

  return(out)
}
