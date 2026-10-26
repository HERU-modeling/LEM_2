## load in fixed parameters 1d or high-d

city_sp <- grep("Y", par_info$city_specific)
non_city_sp <- grep("N", par_info$city_specific)


## Parametrization for single-dimension parameters, i.e. dimension = 1

# city-scpecific parameters
for (i in city_sp) {
  vlist[[i]] <- subset(read.xlsx(WB, sheet = par_info$tab[i]), city == CITY)[, -1]
  if (par_info[i, ]$dimension == 1) {
    vparameters[[i]] <- vlist[[i]]$pe
  }
}

# non-city-scpecific parameters
for (i in non_city_sp) {
  vlist[[i]] <- read.xlsx(WB, sheet = par_info$tab[i])
  if (par_info[i, ]$dimension == 1) {
    vparameters[[i]] <- vlist[[i]]$pe
  }
}

## String variables for variable/group/target names
namelist <- read.xlsx(WB, sheet = "common")


vparameters$names.gp <- namelist$names.gp # group names, 30 groups: fully stratified 30 groups
vparameters$names18 <- namelist$names18[!is.na(namelist$names18)] # group names, 18 groups: gender*risk group*ethnicity (aggregate)
vparameters$names.pwid <- namelist$names.pwid[!is.na(namelist$names.pwid)] # group names, 9 PWID groups: gender*risk group*ethnicity (aggregate)
vparameters$names.msm <- namelist$names.msm[!is.na(namelist$names.msm)] # group names, 12 MSM and MWID groups: gender*risk group*ethnicity * risk level * OAT (fully stratified)
vparameters$names.e <- namelist$names.e[!is.na(namelist$names.e)] # group names, 6 groups: high/low*ethnicity (highly aggregate)
vparameters$names.e2 <- namelist$names.e2[!is.na(namelist$names.e2)] # group names, 6 groups: young/old*ethnicity (highly aggregate)
vparameters$state.name <- c(
  "S1", "S2", "Ia", "I1", "I2", "I3",
  "Da", "D1", "D2", "D3", "T1", "T2", "T3", "O1", "O2", "O3",
  "inc_bo", "inc_bs", "inc_g", "diag", "diag_late", "positive_test", "viral_sup", "death"
)
# names for the 16 health states plus 5 states (incidence opposite sex, same sex, needle sharing, diagnose, death)


# return indices for the subgroups we are interested among all groups
source("ModelCoreModules/Parameter_input/Group.number.R")

## Parametrization for parameters with multiple dimensions

# ART effectiveness for transmission via sex (heterosexual & homosexual)
vparameters$delta_H <- vparameters$delta_M <- vlist$delta_sex$pe

# transition rate from acute to chronic states for infected and diagnosed
vparameters$theta.ad <- vparameters$theta.ai

## Parameter in pwid domain ##
pwid.gn <- pwid.gn.fun(vparameters$names.pwid)

# s: 6 to 9, mwid = pwid  (proportion of shared injections)
vparameters$s <- numeric(length(vparameters$names.pwid))
vparameters$s <- rep(vlist$s$pe, length(vparameters$names.pwid))

## Parameter in MSM domain ##
msm.gn <- msm.gn.fun(vparameters$names.msm)

# nsG: 6 to 12, number of sexual partners msm: assume mwid = msm 
vparameters$nsG <- numeric(length(vparameters$names.msm))
 
vparameters$nsG[intersect(msm.gn$old  , msm.gn$white)] <- with(vlist$nsG, pe[sex.risk.age.group == "old" & ethnicity == "w"]) 
vparameters$nsG[intersect(msm.gn$old  , msm.gn$black)] <- with(vlist$nsG, pe[sex.risk.age.group == "old" & ethnicity == "b"]) 
vparameters$nsG[intersect(msm.gn$old  , msm.gn$hisp )] <- with(vlist$nsG, pe[sex.risk.age.group == "old" & ethnicity == "h"]) 
vparameters$nsG[intersect(msm.gn$young, msm.gn$white)] <- with(vlist$nsG, pe[sex.risk.age.group == "young" & ethnicity == "w"]) 
vparameters$nsG[intersect(msm.gn$young, msm.gn$black)] <- with(vlist$nsG, pe[sex.risk.age.group == "young" & ethnicity == "b"]) 
vparameters$nsG[intersect(msm.gn$young, msm.gn$hisp )] <- with(vlist$nsG, pe[sex.risk.age.group == "young" & ethnicity == "h"]) 

# uis: 6 to 12, probability of condom use msm: Before mwid = msm 
vparameters$uis <- numeric(length(vparameters$names.msm))
vparameters$uis[msm.gn$old  ] <- with(vlist$uis, pe[sex.risk.age.group == "old"]) 
vparameters$uis[msm.gn$young] <- with(vlist$uis, pe[sex.risk.age.group == "young"]) 

# proportion young infected undiagnosed (I only):
vparameters$prop.young.I <- numeric(length(vparameters$names.msm))
vparameters$prop.young.I[intersect(msm.gn$young, msm.gn$msm)]  <- with(vlist$prop.young.I, pe[sex.risk.age.group == "young" & risk == "msm"]) 
vparameters$prop.young.I[intersect(msm.gn$young, msm.gn$mwid)] <- with(vlist$prop.young.I, pe[sex.risk.age.group == "young" & risk == "mwid"])


## Parameter in pop18 (race*gender*risk group) domain ##
gp18.gn <- gp18.gn.fun(vparameters$names18)
	
# pop.diag: 18 to 18 (Diagnosed PLHIV might not be aware)
vparameters$pop.diag <- numeric(length(vparameters$names18))
vparameters$pop.diag[intersect(gp18.gn$pwid, gp18.gn$male)]   <- with(vlist$pop.diag, pe[risk == "pwid" & gender == "m"])
vparameters$pop.diag[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$pop.diag, pe[risk == "pwid" & gender == "f"])
vparameters$pop.diag[gp18.gn$msm]   <- with(vlist$pop.diag, pe[risk == "msm"])
vparameters$pop.diag[gp18.gn$mwid]  <- with(vlist$pop.diag, pe[risk == "mwid"])
vparameters$pop.diag[gp18.gn$het.f] <- with(vlist$pop.diag, pe[risk == "het" & gender == "f"])
vparameters$pop.diag[gp18.gn$het.m] <- with(vlist$pop.diag, pe[risk == "het" & gender == "m"])

# pop.change.18.prop
vparameters$pop.change.18.prop <- numeric(length(vparameters$names18))
vparameters$pop.change.18.prop[intersect(gp18.gn$pwid, gp18.gn$male)]   <- with(vlist$pop.change.18.prop, pe[risk == "pwid" & gender == "m"])
vparameters$pop.change.18.prop[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$pop.change.18.prop, pe[risk == "pwid" & gender == "f"])
vparameters$pop.change.18.prop[gp18.gn$msm]   <- with(vlist$pop.change.18.prop, pe[risk == "msm"])
vparameters$pop.change.18.prop[gp18.gn$mwid]  <- with(vlist$pop.change.18.prop, pe[risk == "mwid"])
vparameters$pop.change.18.prop[gp18.gn$het.f] <- with(vlist$pop.change.18.prop, pe[risk == "het" & gender == "f"])
vparameters$pop.change.18.prop[gp18.gn$het.m] <- with(vlist$pop.change.18.prop, pe[risk == "het" & gender == "m"])

# T1_O: 12 to 18, m=f for pwid, het (ART dropout rate CD4 >= 500)
vparameters$T1_O <- numeric(length(vparameters$names18))
vparameters$T1_O[gp18.gn$white] <- with(vlist$T1_O, pe[ethnicity == "w"])
vparameters$T1_O[gp18.gn$black] <- with(vlist$T1_O, pe[ethnicity == "b"])
vparameters$T1_O[gp18.gn$hisp]  <- with(vlist$T1_O, pe[ethnicity == "h"])

# T2_O: 12 to 18, m=f for pwid, het (ART dropout rate CD4:200-499)
vparameters$T2_O <- numeric(length(vparameters$names18))
#vparameters$T2_O[gp18.gn$white] <- with(vlist$T2_O, pe[ethnicity == "w"])
#vparameters$T2_O[gp18.gn$black] <- with(vlist$T2_O, pe[ethnicity == "b"])
#vparameters$T2_O[gp18.gn$hisp]  <- with(vlist$T2_O, pe[ethnicity == "h"])
vparameters$T2_O <- vparameters$T1_O * vparameters$T12_O
for (i in seq_len(length(vparameters$T2_O))) {
	if (vparameters$T2_O[i] < 0) {
		warning(paste0(
			"vparameters$T2_O = ", vparameters$T2_O[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		vparameters$T2_O[i] <- 0.05
	}
	
	if (vparameters$T2_O[i] > 1) {
		warning(paste0(
			"vparameters$T2_O = ", vparameters$T2_O[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		vparameters$T2_O[i] <- 0.95
	}
}

# T3_O: 12 to 18, m=f for pwid, het (ART dropout rate CD4 < 200)
vparameters$T3_O <- numeric(length(vparameters$names18))
vparameters$T3_O <- vparameters$T1_O * vparameters$T13_O
for (i in seq_len(length(vparameters$T3_O))) {
	if (vparameters$T3_O[i] < 0) {
		warning(paste0(
			"vparameters$T3_O = ", vparameters$T3_O[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		vparameters$T3_O[i] <- 0.05
	}
	
	if (vparameters$T3_O[i] > 1) {
		warning(paste0(
			"vparameters$T3_O = ", vparameters$T3_O[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		vparameters$T3_O[i] <- 0.95
	}
}

# T1_T2: 18 to 18 (disease progression rate from T1 to T2)
vparameters$T1_T2 <- numeric(length(vparameters$names18))
vparameters$T1_T2[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T1_T2, pe[risk == "pwid" & gender == "m"])
vparameters$T1_T2[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T1_T2, pe[risk == "pwid" & gender == "f"])
vparameters$T1_T2[gp18.gn$msm] <- with(vlist$T1_T2, pe[risk == "msm"])
vparameters$T1_T2[gp18.gn$mwid] <- with(vlist$T1_T2, pe[risk == "mwid"])
vparameters$T1_T2[gp18.gn$het.f] <- with(vlist$T1_T2, pe[risk == "het" & gender == "f"])
vparameters$T1_T2[gp18.gn$het.m] <- with(vlist$T1_T2, pe[risk == "het" & gender == "m"])


# T1_T3: 18 to 18 (disease progression rate from T1 to T3)
vparameters$T1_T3 <- numeric(length(vparameters$names18))
vparameters$T1_T3[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T1_T3, pe[risk == "pwid" & gender == "m"])
vparameters$T1_T3[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T1_T3, pe[risk == "pwid" & gender == "f"])
vparameters$T1_T3[gp18.gn$msm] <- with(vlist$T1_T3, pe[risk == "msm"])
vparameters$T1_T3[gp18.gn$mwid] <- with(vlist$T1_T3, pe[risk == "mwid"])
vparameters$T1_T3[gp18.gn$het.f] <- with(vlist$T1_T3, pe[risk == "het" & gender == "f"])
vparameters$T1_T3[gp18.gn$het.m] <- with(vlist$T1_T3, pe[risk == "het" & gender == "m"])


# T2_T3: 18 to 18 (disease progression rate from T2 to T3)
vparameters$T2_T3 <- numeric(length(vparameters$names18))
vparameters$T2_T3[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T2_T3, pe[risk == "pwid" & gender == "m"])
vparameters$T2_T3[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T2_T3, pe[risk == "pwid" & gender == "f"])
vparameters$T2_T3[gp18.gn$msm] <- with(vlist$T2_T3, pe[risk == "msm"])
vparameters$T2_T3[gp18.gn$mwid] <- with(vlist$T2_T3, pe[risk == "mwid"])
vparameters$T2_T3[gp18.gn$het.f] <- with(vlist$T2_T3, pe[risk == "het" & gender == "f"])
vparameters$T2_T3[gp18.gn$het.m] <- with(vlist$T2_T3, pe[risk == "het" & gender == "m"])


# T2_T1: 18 to 18 (disease progression rate from T2 to T1)
vparameters$T2_T1 <- numeric(length(vparameters$names18))
vparameters$T2_T1[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T2_T1, pe[risk == "pwid" & gender == "m"])
vparameters$T2_T1[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T2_T1, pe[risk == "pwid" & gender == "f"])
vparameters$T2_T1[gp18.gn$msm] <- with(vlist$T2_T1, pe[risk == "msm"])
vparameters$T2_T1[gp18.gn$mwid] <- with(vlist$T2_T1, pe[risk == "mwid"])
vparameters$T2_T1[gp18.gn$het.f] <- with(vlist$T2_T1, pe[risk == "het" & gender == "f"])
vparameters$T2_T1[gp18.gn$het.m] <- with(vlist$T2_T1, pe[risk == "het" & gender == "m"])


# T3_T1: 18 to 18 (disease progression rate from T3 to T1)
vparameters$T3_T1 <- numeric(length(vparameters$names18))
vparameters$T3_T1[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T3_T1, pe[risk == "pwid" & gender == "m"])
vparameters$T3_T1[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T3_T1, pe[risk == "pwid" & gender == "f"])
vparameters$T3_T1[gp18.gn$msm] <- with(vlist$T3_T1, pe[risk == "msm"])
vparameters$T3_T1[gp18.gn$mwid] <- with(vlist$T3_T1, pe[risk == "mwid"])
vparameters$T3_T1[gp18.gn$het.f] <- with(vlist$T3_T1, pe[risk == "het" & gender == "f"])
vparameters$T3_T1[gp18.gn$het.m] <- with(vlist$T3_T1, pe[risk == "het" & gender == "m"])


# T3_T2: 18 to 18 (disease progression rate from T3 to T2)
vparameters$T3_T2 <- numeric(length(vparameters$names18))
vparameters$T3_T2[intersect(gp18.gn$pwid, gp18.gn$male)] <- with(vlist$T3_T2, pe[risk == "pwid" & gender == "m"])
vparameters$T3_T2[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$T3_T2, pe[risk == "pwid" & gender == "f"])
vparameters$T3_T2[gp18.gn$msm] <- with(vlist$T3_T2, pe[risk == "msm"])
vparameters$T3_T2[gp18.gn$mwid] <- with(vlist$T3_T2, pe[risk == "mwid"])
vparameters$T3_T2[gp18.gn$het.f] <- with(vlist$T3_T2, pe[risk == "het" & gender == "f"])
vparameters$T3_T2[gp18.gn$het.m] <- with(vlist$T3_T2, pe[risk == "het" & gender == "m"])

# O_T: 12 to 18, m=f for pwid, het (ART re-initiation)
vparameters$O_T <- numeric(length(vparameters$names18))
vparameters$O_T[gp18.gn$white] <- with(vlist$O_T, pe[ethnicity == "w"])
vparameters$O_T[gp18.gn$black] <- with(vlist$O_T, pe[ethnicity == "b"])
vparameters$O_T[gp18.gn$hisp]  <- with(vlist$O_T, pe[ethnicity == "h"])

# mu_adj: 6 to 18, (by ethnicity, account for maturation + migration rate among PLHV)
vparameters$mu_adj <- numeric(length(vparameters$names18))
vparameters$mu_adj[gp18.gn$white] <- as.numeric(with(vlist$mu_adj, pe[ethnicity == "w"]))
vparameters$mu_adj[gp18.gn$black] <- as.numeric(with(vlist$mu_adj, pe[ethnicity == "b"]))
vparameters$mu_adj[gp18.gn$hisp]  <- as.numeric(with(vlist$mu_adj, pe[ethnicity == "h"]))

# immediate ART initiation (I -> T)
# phi:immediate ART initiators
# p1: proportion of treated PLHIV linked to HIV care following diagnosis
# phi_: proportion of immediate ART initiation among newly diagnosed PLHIV on HIV care

# not immediate ART initiation (D -> T)
# alpha: no immediate ART initiation rate
# alpha_: ART initiation rate among diagnosed PLHIV on HIV care >= 30 days
# p2: proportion of diagnosed but untreated PLHIV that are linked to care

# phi1, phi2, phi3, alpha1, alpha2, alpha3: derived from phi1_, phi2_, phi3_, p1, alpha1_, alpha2_, alpha3_
# phi1_, phi2_, phi3_: 7 to 18, msm*ethnicity, pwid, mwid, het*gender
# alpha1_, alpha2_, alpha3_: 8 to 18, msm*ethnicity, pwid, het.f, mwid. het.m.w=het.m.h, het.m.b
# phi1 =phi1_*p1; phi2 =phi2_*p1; phi3 =phi3_*p1;
# alpha1 =alpha1_*(p1*(1-phi1_))/(1-p1*phi1_); alpha2 =alpha2_*(p1*(1-phi1_))/(1-p1*phi2_);
# alpha3 =alpha3_*(p1*(1-phi3_))/(1-p1*phi3_)

phi1_ <- numeric(length(vparameters$names18))
phi2_ <- numeric(length(vparameters$names18))
phi3_ <- numeric(length(vparameters$names18))
p1 <- numeric(length(vparameters$names18))
alpha1_ <- numeric(length(vparameters$names18))
alpha2_ <- numeric(length(vparameters$names18))
alpha3_ <- numeric(length(vparameters$names18))

# Populate phi1_
phi1_[gp18.gn$white]   <- with(vlist$phi1_, pe[ethnicity == "w"])
phi1_[gp18.gn$black]   <- with(vlist$phi1_, pe[ethnicity == "b"])
phi1_[gp18.gn$hisp]    <- with(vlist$phi1_, pe[ethnicity == "h"])


# Populate phi2_
phi2_                  <- phi1_ * vparameters$phi12
for (i in seq_len(length(phi2_))) {
	if (phi2_[i] < 0) {
		warning(paste0(
			"phi2_ = ", phi2_[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		phi2_[i] <- 0.05
	}
	
	if (phi2_[i] > 1) {
		warning(paste0(
			"phi2_ = ", phi2_[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		phi2_[i] <- 0.95
	}
}

# Populate phi3_
phi3_                  <- phi1_ * vparameters$phi13
for (i in seq_len(length(phi3_))) {
	if (phi3_[i] < 0) {
		warning(paste0(
			"phi3_ = ", phi3_[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		phi3_[i] <- 0.05
	}
	
	if (phi3_[i] > 1) {
		warning(paste0(
			"phi3_ = ", phi3_[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		phi3_[i] <- 0.95
	}
}

# Populate p1
p1[gp18.gn$white]                                  <- with(vlist$p1, pe[ethnicity == "w" & year == "2017"])
p1[gp18.gn$black]                                  <- with(vlist$p1, pe[ethnicity == "b" & year == "2017"])
p1[gp18.gn$hisp ]                                  <- with(vlist$p1, pe[ethnicity == "h" & year == "2017"])

p1.update        <- p1 # Original data
p1.update[gp18.gn$white]                           <- with(vlist$p1, pe[ethnicity == "w" & year == "2022"])
p1.update[gp18.gn$black]                           <- with(vlist$p1, pe[ethnicity == "b" & year == "2022"])
p1.update[gp18.gn$hisp ]                           <- with(vlist$p1, pe[ethnicity == "h" & year == "2022"])

# Populate alpha1_
alpha1_[gp18.gn$white] <- with(vlist$alpha1_, pe[ethnicity == "w"])
alpha1_[gp18.gn$black] <- with(vlist$alpha1_, pe[ethnicity == "b"])
alpha1_[gp18.gn$hisp]  <- with(vlist$alpha1_, pe[ethnicity == "h"])


# Populate alpha2_
alpha2_                 <- alpha1_ * vparameters$alpha12
for (i in seq_len(length(alpha2_))) {
	if (alpha2_[i] < 0) {
		warning(paste0(
			"alpha2_ = ", alpha2_[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		alpha2_[i] <- 0.05
	}
	
	if (alpha2_[i] > 1) {
		warning(paste0(
			"alpha2_ = ", alpha2_[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		alpha2_[i] <- 0.95
	}
}

# Populate alpha3_
alpha3_                 <- alpha1_ * vparameters$alpha13
for (i in seq_len(length(alpha3_))) {
	if (alpha3_[i] < 0) {
		warning(paste0(
			"alpha3_ = ", alpha3_[i],	" < 0 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.05."))
		alpha3_[i] <- 0.05
	}
	
	if (alpha3_[i] > 1) {
		warning(paste0(
			"alpha3_ = ", alpha3_[i],	" > 1 in ", CITY,
			" for ", vparameters$names18[i],". Setting it to 0.95."))
		alpha3_[i] <- 0.95
	}
}

# Generate first p1 yr data vector
phi1 <- phi1_ * p1
phi2 <- phi2_ * p1
phi3 <- phi3_ * p1
alpha1 <- alpha1_ * (p1 * (1 - phi1_)) / (1 - p1 * phi1_)
alpha2 <- alpha2_ * (p1 * (1 - phi2_)) / (1 - p1 * phi2_)
alpha3 <- alpha3_ * (p1 * (1 - phi3_)) / (1 - p1 * phi3_) # 1(year) * 18 (gend*ethn*risk)

# Generate first p1 yr data vector
phi1.update = phi1_ * p1.update
phi2.update = phi2_ * p1.update
phi3.update = phi3_ * p1.update
alpha1.update = alpha1_ * (p1.update * (1 - phi1_)) / (1 - p1.update * phi1_)
alpha2.update = alpha2_ * (p1.update * (1 - phi2_)) / (1 - p1.update * phi2_)
alpha3.update = alpha3_ * (p1.update * (1 - phi3_)) / (1 - p1.update * phi3_)


# Convert data vector to matrix
length     <- length(vlist$p1$year)
first.year <- vlist$p1$year[1]
last.year  <- vlist$p1$year[length/2+1]
time = as.numeric(last.year) - as.numeric(first.year) + 1
phi1_upmx <- phi2_upmx <- phi3_upmx <- alpha1_upmx <- alpha2_upmx <- alpha3_upmx <- matrix(0, nrow = time, ncol = length(phi1))
for (uu in 1:18){
	phi1_upmx[ , uu]   <- seq(phi1[uu], phi1.update[uu], length.out = time)
	phi2_upmx[ , uu]   <- seq(phi2[uu], phi2.update[uu], length.out = time)
	phi3_upmx[ , uu]   <- seq(phi3[uu], phi3.update[uu], length.out = time)
	alpha1_upmx[ , uu] <- seq(alpha1[uu], alpha1.update[uu], length.out = time)
	alpha2_upmx[ , uu] <- seq(alpha2[uu], alpha2.update[uu], length.out = time)
	alpha3_upmx[ , uu] <- seq(alpha3[uu], alpha3.update[uu], length.out = time)
}

# Save matrix
vparameters$phi1 <- phi1_upmx
vparameters$phi2 <- phi2_upmx
vparameters$phi3 <- phi3_upmx
vparameters$alpha1 <- alpha1_upmx
vparameters$alpha2 <- alpha2_upmx
vparameters$alpha3 <- alpha3_upmx # 4 (years 2015-18) * 18 (gend*ethn*risk)

# prop.HIV.aware: 6 to 18, % aware of HIV status
vparameters$prop.HIV.aware <- numeric(length(vparameters$names18))
vparameters$prop.HIV.aware[gp18.gn$msm ] <- with(vlist$prop.HIV.aware, pe[risk == "msm"]) 
vparameters$prop.HIV.aware[gp18.gn$mwid] <- with(vlist$prop.HIV.aware, pe[risk == "mwid"]) 
vparameters$prop.HIV.aware[intersect(gp18.gn$pwid, gp18.gn$male  )] <- with(vlist$prop.HIV.aware, pe[gender == "m" & risk == "pwid"]) 
vparameters$prop.HIV.aware[intersect(gp18.gn$pwid, gp18.gn$female)] <- with(vlist$prop.HIV.aware, pe[gender == "f" & risk == "pwid"])
vparameters$prop.HIV.aware[gp18.gn$het.m] <- with(vlist$prop.HIV.aware, pe[gender == "m" & risk == "het"]) 
vparameters$prop.HIV.aware[gp18.gn$het.f] <- with(vlist$prop.HIV.aware, pe[gender == "f" & risk == "het"])

# prop.ever.art: 9 to 18, 1-D data to 3-D data, use male white het as the reference level = p1
# initial value for proportion ever on ART
# p1: PLHIV linked to HIV care following diagnosis, derived from city-/state-level surveillance database
vparameters$prop.ever.art[gp18.gn$white] <- with(vlist$prop.ever.art, pe[ethnicity == "w"])
vparameters$prop.ever.art[gp18.gn$black] <- with(vlist$prop.ever.art, pe[ethnicity == "b"])
vparameters$prop.ever.art[gp18.gn$hisp]  <- with(vlist$prop.ever.art, pe[ethnicity == "h"])
#vparameters$prop.ever.art <- p1


# prop.curr.art: 3 to  18,
# initial value for proportion currently on ART, by ethnicity
vparameters$prop.curr.art <- numeric(length(vparameters$names18))
vparameters$prop.curr.art[gp18.gn$white] <- with(vlist$prop.curr.art, pe[ethnicity == "w"])
vparameters$prop.curr.art[gp18.gn$black] <- with(vlist$prop.curr.art, pe[ethnicity == "b"])
vparameters$prop.curr.art[gp18.gn$hisp] <- with(vlist$prop.curr.art, pe[ethnicity == "h"])

# prop.high.sus: 7 to  18, 0.25 certain for msm, mpwid, NA for pwid
# proportion high risk susceptible
vparameters$prop.high.sus <- numeric(length(vparameters$names18))
vparameters$prop.high.sus[intersect(gp18.gn$het.m,gp18.gn$white)] <- with(vlist$prop.high.sus, pe[gender == "m" & ethnicity == "w" & risk == "het"])
vparameters$prop.high.sus[intersect(gp18.gn$het.m,gp18.gn$black)] <- with(vlist$prop.high.sus, pe[gender == "m" & ethnicity == "b" & risk == "het"])
vparameters$prop.high.sus[intersect(gp18.gn$het.m,gp18.gn$hisp)]  <- with(vlist$prop.high.sus, pe[gender == "m" & ethnicity == "h" & risk == "het"])
vparameters$prop.high.sus[intersect(gp18.gn$het.f,gp18.gn$white)] <- with(vlist$prop.high.sus, pe[gender == "f" & ethnicity == "w" & risk == "het"])
vparameters$prop.high.sus[intersect(gp18.gn$het.f,gp18.gn$black)] <- with(vlist$prop.high.sus, pe[gender == "f" & ethnicity == "b" & risk == "het"])
vparameters$prop.high.sus[intersect(gp18.gn$het.f,gp18.gn$hisp)]  <- with(vlist$prop.high.sus, pe[gender == "f" & ethnicity == "h" & risk == "het"])

# proportion high risk infected
vparameters$prop.high.inf <- numeric(length(vparameters$names18))
vparameters$prop.high.inf[intersect(gp18.gn$het.m,gp18.gn$white)] <- with(vlist$prop.high.inf, pe[gender == "m" & ethnicity == "w" & risk == "het"])
vparameters$prop.high.inf[intersect(gp18.gn$het.m,gp18.gn$black)] <- with(vlist$prop.high.inf, pe[gender == "m" & ethnicity == "b" & risk == "het"])
vparameters$prop.high.inf[intersect(gp18.gn$het.m,gp18.gn$hisp)]  <- with(vlist$prop.high.inf, pe[gender == "m" & ethnicity == "h" & risk == "het"])
vparameters$prop.high.inf[intersect(gp18.gn$het.f,gp18.gn$white)] <- with(vlist$prop.high.inf, pe[gender == "f" & ethnicity == "w" & risk == "het"])
vparameters$prop.high.inf[intersect(gp18.gn$het.f,gp18.gn$black)] <- with(vlist$prop.high.inf, pe[gender == "f" & ethnicity == "b" & risk == "het"])
vparameters$prop.high.inf[intersect(gp18.gn$het.f,gp18.gn$hisp)]  <- with(vlist$prop.high.inf, pe[gender == "f" & ethnicity == "h" & risk == "het"])

# proportion young susceptible:
vparameters$prop.young.sus <- numeric(length(vparameters$names18))
vparameters$prop.young.sus[intersect(gp18.gn$msm,gp18.gn$white)] <- with(vlist$prop.young.sus, pe[ethnicity == "w" & risk == "msm"])
vparameters$prop.young.sus[intersect(gp18.gn$msm,gp18.gn$black)] <- with(vlist$prop.young.sus, pe[ethnicity == "b" & risk == "msm"])
vparameters$prop.young.sus[intersect(gp18.gn$msm,gp18.gn$hisp)]  <- with(vlist$prop.young.sus, pe[ethnicity == "h" & risk == "msm"])
vparameters$prop.young.sus[intersect(gp18.gn$mwid,gp18.gn$white)] <- with(vlist$prop.young.sus, pe[ethnicity == "w" & risk == "mwid"])
vparameters$prop.young.sus[intersect(gp18.gn$mwid,gp18.gn$black)] <- with(vlist$prop.young.sus, pe[ethnicity == "b" & risk == "mwid"])
vparameters$prop.young.sus[intersect(gp18.gn$mwid,gp18.gn$hisp)]  <- with(vlist$prop.young.sus, pe[ethnicity == "h" & risk == "mwid"])

# proportion young infected:
vparameters$prop.young.inf <- numeric(length(vparameters$names18))
vparameters$prop.young.inf[intersect(gp18.gn$msm,gp18.gn$white)] <- with(vlist$prop.young.inf, pe[ethnicity == "w" & risk == "msm"])
vparameters$prop.young.inf[intersect(gp18.gn$msm,gp18.gn$black)] <- with(vlist$prop.young.inf, pe[ethnicity == "b" & risk == "msm"])
vparameters$prop.young.inf[intersect(gp18.gn$msm,gp18.gn$hisp)]  <- with(vlist$prop.young.inf, pe[ethnicity == "h" & risk == "msm"])
vparameters$prop.young.inf[intersect(gp18.gn$mwid,gp18.gn$white)] <- with(vlist$prop.young.inf, pe[ethnicity == "w" & risk == "mwid"])
vparameters$prop.young.inf[intersect(gp18.gn$mwid,gp18.gn$black)] <- with(vlist$prop.young.inf, pe[ethnicity == "b" & risk == "mwid"])
vparameters$prop.young.inf[intersect(gp18.gn$mwid,gp18.gn$hisp)]  <- with(vlist$prop.young.inf, pe[ethnicity == "h" & risk == "mwid"])

# mu_mat: 3 to 18,
# Maturation among S MSM/old
vparameters$mu_mat <- numeric(length(vparameters$names18))
vparameters$mu_mat[gp18.gn$white] <- with(vlist$mu_mat, pe[ethnicity == "w"])
vparameters$mu_mat[gp18.gn$black] <- with(vlist$mu_mat, pe[ethnicity == "b"])
vparameters$mu_mat[gp18.gn$hisp]  <- with(vlist$mu_mat, pe[ethnicity == "h"])

# old.e.sus: 3 to 18,
# Maturation among S MSM/old
vparameters$old.e.sus <- numeric(length(vparameters$names18))
vparameters$old.e.sus[gp18.gn$white] <- with(vlist$old.e.sus, pe[ethnicity == "w"])
vparameters$old.e.sus[gp18.gn$black] <- with(vlist$old.e.sus, pe[ethnicity == "b"])
vparameters$old.e.sus[gp18.gn$hisp]  <- with(vlist$old.e.sus, pe[ethnicity == "h"])

# old.e.plhiv: 3 to 18,
# Maturation among S MSM/old
vparameters$old.e.hiv <- numeric(length(vparameters$names18))
vparameters$old.e.hiv[gp18.gn$white] <- with(vlist$old.e.hiv, pe[ethnicity == "w"])
vparameters$old.e.hiv[gp18.gn$black] <- with(vlist$old.e.hiv, pe[ethnicity == "b"])
vparameters$old.e.hiv[gp18.gn$hisp]  <- with(vlist$old.e.hiv, pe[ethnicity == "h"])

# 3.2.2-PropViralSupOnART: 3 to 18,
# For viral suppression target calculation
vparameters$prop_vs_T <- numeric(length(vparameters$names18))
vparameters$prop_vs_T[gp18.gn$white] <- with(vlist$prop_vs_T, pe[ethnicity == "w"])
vparameters$prop_vs_T[gp18.gn$black] <- with(vlist$prop_vs_T, pe[ethnicity == "b"])
vparameters$prop_vs_T[gp18.gn$hisp]  <- with(vlist$prop_vs_T, pe[ethnicity == "h"])

# 1.2.2-I3PropMultiplier: 3 to 18,
# For late diagnoses target calculation
vparameters$I3_multi <- numeric(length(vparameters$names18))
vparameters$I3_multi[gp18.gn$white] <- with(vlist$I3_multi, pe[ethnicity == "w"])
vparameters$I3_multi[gp18.gn$black] <- with(vlist$I3_multi, pe[ethnicity == "b"])
vparameters$I3_multi[gp18.gn$hisp]  <- with(vlist$I3_multi, pe[ethnicity == "h"])

## Parameter in 30 full gp domain ##
gp.gn <- gp.gn.fun(vparameters$names.gp)
# prop.link: % of diag linked to care (42 groups), no OAT and risk level differences
vparameters$prop.link <- numeric(length(vparameters$names.gp))
rname <- gsub(paste(c("/low", "/high","/old", "/young"), collapse = "|"), "", vparameters$names.gp)
for (i in 1:18) {
  ind <- which(rname %in% vparameters$names18[i])
  vparameters$prop.link[ind] <- p1[i]
}

# noG: 18 to 30, msm=mwid, pwid=high.het (number of sexual partners opposite sex)
vparameters$noG <- numeric(length(vparameters$names.gp))
vparameters$noG[gp.gn$all.msm] <- with(vlist$noG, pe[risk == "msm"])
vparameters$noG[gp.gn$low]     <- with(vlist$noG, pe[risk == "het" & sex.risk.age.group == "low"])
vparameters$noG[gp.gn$H.het.m] <- with(vlist$noG, pe[gender == "m" & risk == "het" & sex.risk.age.group == "high"])
vparameters$noG[gp.gn$H.het.f] <- with(vlist$noG, pe[gender == "f" & risk == "het" & sex.risk.age.group == "high"])
vparameters$noG[gp.gn$pwid] <- with(vlist$noG, pe[risk == "het" & sex.risk.age.group == "low"]) #* vparameters$noG.pwid


# uio: 12 to 30, no risk diff in low,  pwid=high het, msm=mwid, H.het.f=L.het.f (prob condom use opposite sex)
vparameters$uio <- numeric(length(vparameters$names.gp))
vparameters$uio[gp.gn$all.msm] <- with(vlist$uio, pe[risk == "msm"])
vparameters$uio[union(gp.gn$het.m, gp.gn$pwid.m)]   <- with(vlist$uio, pe[risk == "het" & gender == "m"])
vparameters$uio[union(gp.gn$het.f, gp.gn$pwid.f)]   <- with(vlist$uio, pe[risk == "het" & gender == "f"])

# mor_S: 13 to 30, all pwid equal, low=high (mortality baseline), on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_S <- numeric(length(vparameters$names.gp))
vparameters$mor_S[gp.gn$pwid] <- with(vlist$mor_S, pe[risk == "pwid"])[1]
vparameters$mor_S[gp.gn$Y.msm] <- with(vlist$mor_S, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_S[gp.gn$O.msm] <- with(vlist$mor_S, pe[risk == "msm" & sex.risk.age.group == "old"])
vparameters$mor_S[gp.gn$Y.mwid] <- with(vlist$mor_S, pe[risk == "mwid" & sex.risk.age.group == "young"])
vparameters$mor_S[gp.gn$O.mwid] <- with(vlist$mor_S, pe[risk == "mwid" & sex.risk.age.group == "old"])
vparameters$mor_S[gp.gn$het.m] <- with(vlist$mor_S, pe[risk == "het" & gender == "m"])
vparameters$mor_S[gp.gn$het.f] <- with(vlist$mor_S, pe[risk == "het" & gender == "f"])


# Mortality rates for PLHIV were differentiated by HIV progression stages (acute and three CD4 levels)
# and whether the individual is on ART

# mor_Ia: = mor_S (mor acute infected)
vparameters$mor_Ia <- vparameters$mor_S

# mor_I1: 13 to 30, all pwid=, (mor infected CD4 >= 500)
vparameters$mor_I1 <- vparameters$mor_S

# whether is pwid and whether on OAT
# mor_I2 (mor infected CD4: 200-499): 1 to 42, only diff by risk, on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_I2 <- numeric(length(vparameters$names.gp))
vparameters$mor_I2[gp.gn$pwid]   <- with(vlist$mor_I2, pe[risk == "pwid"]) * vparameters$pwid_mor2
vparameters$mor_I2[gp.gn$Y.msm]  <- with(vlist$mor_I2, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_I2[gp.gn$O.msm]  <- with(vlist$mor_I2, pe[risk == "msm" & sex.risk.age.group == "old"])
vparameters$mor_I2[gp.gn$Y.mwid] <- with(vlist$mor_I2, pe[risk == "mwid" & sex.risk.age.group == "young"]) * vparameters$pwid_mor2
vparameters$mor_I2[gp.gn$O.mwid] <- with(vlist$mor_I2, pe[risk == "mwid" & sex.risk.age.group == "old"])   * vparameters$pwid_mor2
vparameters$mor_I2[gp.gn$het]    <- with(vlist$mor_I2, pe[risk == "het"])

# mor_I3 (mor infected CD4 < 200): 1 to 30, only diff by risk, on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_I3 <- numeric(length(vparameters$names.gp))
vparameters$mor_I3[gp.gn$pwid]   <- with(vlist$mor_I3, pe[risk == "pwid"]) * vparameters$pwid_mor3
vparameters$mor_I3[gp.gn$Y.msm]  <- with(vlist$mor_I3, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_I3[gp.gn$O.msm]  <- with(vlist$mor_I3, pe[risk == "msm" & sex.risk.age.group == "old"])
vparameters$mor_I3[gp.gn$Y.mwid] <- with(vlist$mor_I3, pe[risk == "mwid" & sex.risk.age.group == "young"]) * vparameters$pwid_mor3
vparameters$mor_I3[gp.gn$O.mwid] <- with(vlist$mor_I3, pe[risk == "mwid" & sex.risk.age.group == "old"])   * vparameters$pwid_mor3
vparameters$mor_I3[gp.gn$het]    <- with(vlist$mor_I3, pe[risk == "het"])


# mor_T1 (mor infected CD4 >=500 on ART): 18 to 30, on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_T1 <- numeric(length(vparameters$names.gp))
vparameters$mor_T1[gp.gn$pwid.m] <- with(vlist$mor_T1, pe[risk == "pwid" & gender == "m"])
vparameters$mor_T1[gp.gn$pwid.f] <- with(vlist$mor_T1, pe[risk == "pwid" & gender == "f"])
vparameters$mor_T1[gp.gn$Y.msm] <- with(vlist$mor_T1, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_T1[gp.gn$O.msm] <- with(vlist$mor_T1, pe[risk == "msm" & sex.risk.age.group == "old"])

vparameters$mor_T1[gp.gn$Y.mwid] <- with(vlist$mor_T1, pe[risk == "mwid" & sex.risk.age.group == "young"])
vparameters$mor_T1[gp.gn$O.mwid] <- with(vlist$mor_T1, pe[risk == "mwid" & sex.risk.age.group == "old"])
vparameters$mor_T1[gp.gn$het.m] <- with(vlist$mor_T1, pe[risk == "het" & gender == "m"])
vparameters$mor_T1[gp.gn$het.f] <- with(vlist$mor_T1, pe[risk == "het" & gender == "f"])


# mor_T2 (mor infected CD4: 200-499 on ART): 18 to 30, on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_T2 <- numeric(length(vparameters$names.gp))
vparameters$mor_T2[gp.gn$pwid.m] <- with(vlist$mor_T2, pe[risk == "pwid" & gender == "m"])
vparameters$mor_T2[gp.gn$pwid.f] <- with(vlist$mor_T2, pe[risk == "pwid" & gender == "f"])
vparameters$mor_T2[gp.gn$Y.msm] <- with(vlist$mor_T2, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_T2[gp.gn$O.msm] <- with(vlist$mor_T2, pe[risk == "msm" & sex.risk.age.group == "old"])
vparameters$mor_T2[gp.gn$Y.mwid] <- with(vlist$mor_T2, pe[risk == "mwid" & sex.risk.age.group == "young"])
vparameters$mor_T2[gp.gn$O.mwid] <- with(vlist$mor_T2, pe[risk == "mwid" & sex.risk.age.group == "old"])
vparameters$mor_T2[gp.gn$het.m] <- with(vlist$mor_T2, pe[risk == "het" & gender == "m"])
vparameters$mor_T2[gp.gn$het.f] <- with(vlist$mor_T2, pe[risk == "het" & gender == "f"])

# mor_T3 (mor infected CD4<200 on ART): 18 to 30, on OAT is lower (implement in ODE as part of prop.oat)
vparameters$mor_T3 <- numeric(length(vparameters$names.gp))
vparameters$mor_T3[gp.gn$pwid.m] <- with(vlist$mor_T3, pe[risk == "pwid" & gender == "m"])
vparameters$mor_T3[gp.gn$pwid.f] <- with(vlist$mor_T3, pe[risk == "pwid" & gender == "f"])
vparameters$mor_T3[gp.gn$Y.msm] <- with(vlist$mor_T3, pe[risk == "msm" & sex.risk.age.group == "young"])
vparameters$mor_T3[gp.gn$O.msm] <- with(vlist$mor_T3, pe[risk == "msm" & sex.risk.age.group == "old"])
vparameters$mor_T3[gp.gn$Y.mwid] <- with(vlist$mor_T3, pe[risk == "mwid" & sex.risk.age.group == "young"])
vparameters$mor_T3[gp.gn$O.mwid] <- with(vlist$mor_T3, pe[risk == "mwid" & sex.risk.age.group == "old"])
vparameters$mor_T3[gp.gn$het.m] <- with(vlist$mor_T3, pe[risk == "het" & gender == "m"])
vparameters$mor_T3[gp.gn$het.f] <- with(vlist$mor_T3, pe[risk == "het" & gender == "f"])

# prop.S2: 21 to 30, LOW: mwid=msm=het, HIGH: msm=mwid (proportion of screened susceptible)
vparameters$prop.S2 <- numeric(length(vparameters$names.gp))
vparameters$prop.S2[intersect(gp.gn$msm, gp.gn$white)] <- with(vlist$prop.S2, pe[risk == "msm" & ethnicity == "w"]) 
vparameters$prop.S2[intersect(gp.gn$msm, gp.gn$black)] <- with(vlist$prop.S2, pe[risk == "msm" & ethnicity == "b"]) 
vparameters$prop.S2[intersect(gp.gn$msm, gp.gn$hisp )] <- with(vlist$prop.S2, pe[risk == "msm" & ethnicity == "h"]) 
vparameters$prop.S2[intersect(gp.gn$mwid, gp.gn$white)] <- with(vlist$prop.S2, pe[risk == "mwid" & ethnicity == "w"]) 
vparameters$prop.S2[intersect(gp.gn$mwid, gp.gn$black)] <- with(vlist$prop.S2, pe[risk == "mwid" & ethnicity == "b"]) 
vparameters$prop.S2[intersect(gp.gn$mwid, gp.gn$hisp )] <- with(vlist$prop.S2, pe[risk == "mwid" & ethnicity == "h"]) 
vparameters$prop.S2[intersect(gp.gn$pwid, gp.gn$white)]    <- with(vlist$prop.S2, pe[risk == "pwid" & ethnicity == "w"]) #m & f
vparameters$prop.S2[intersect(gp.gn$pwid, gp.gn$black)]    <- with(vlist$prop.S2, pe[risk == "pwid" & ethnicity == "b"]) #m & f
vparameters$prop.S2[intersect(gp.gn$pwid, gp.gn$hisp )]    <- with(vlist$prop.S2, pe[risk == "pwid" & ethnicity == "h"]) #m & f
vparameters$prop.S2[intersect(gp.gn$L.het.m, gp.gn$white)]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "w" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$L.het.m, gp.gn$black)]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "b" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$L.het.m, gp.gn$hisp )]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "h" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$H.het.m, gp.gn$white)]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "w" & sex.risk.age.group == "high"])
vparameters$prop.S2[intersect(gp.gn$H.het.m, gp.gn$black)]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "b" & sex.risk.age.group == "high"])
vparameters$prop.S2[intersect(gp.gn$H.het.m, gp.gn$hisp )]   <- with(vlist$prop.S2, pe[gender == "m" & risk == "het" & ethnicity == "h" & sex.risk.age.group == "high"])
vparameters$prop.S2[intersect(gp.gn$L.het.f, gp.gn$white)]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "w" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$L.het.f, gp.gn$black)]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "b" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$L.het.f, gp.gn$hisp )]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "h" & sex.risk.age.group == "low"])
vparameters$prop.S2[intersect(gp.gn$H.het.f, gp.gn$white)]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "w" & sex.risk.age.group == "high"])
vparameters$prop.S2[intersect(gp.gn$H.het.f, gp.gn$black)]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "b" & sex.risk.age.group == "high"])
vparameters$prop.S2[intersect(gp.gn$H.het.f, gp.gn$hisp )]   <- with(vlist$prop.S2, pe[gender == "f" & risk == "het" & ethnicity == "h" & sex.risk.age.group == "high"])


# prop.Ia, prop.Da: 1 to 30 (proportion in acute stage IA, DA)
vparameters$prop.Ia <- vparameters$prop.Da <- rep(vlist$prop.a$pe, 30)


# prop.I1, prop.I2, prop.I3:  30*1 each (initial proportion of different I stages)
vparameters$prop.I1 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.I2 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.I3 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.I1 <- as.numeric(with(vlist$prop.I, pe[cd4 == 1]))
vparameters$prop.I2 <- as.numeric(with(vlist$prop.I, pe[cd4 == 2]))
vparameters$prop.I3 <- as.numeric(with(vlist$prop.I, pe[cd4 == 3]))


# prop.D1, prop.D2, prop.D3: 30*1 each (initial proportion of different D stages)
vparameters$prop.D1 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.D2 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.D3 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)

vparameters$prop.D1[gp.gn$pwid] <- with(vlist$prop.D, pe[cd4 == 1 & risk == "pwid"])[1:3]
vparameters$prop.D1[gp.gn$msm] <- with(vlist$prop.D, pe[cd4 == 1 & risk == "msm"])
vparameters$prop.D1[gp.gn$mwid] <- with(vlist$prop.D, pe[cd4 == 1 & risk == "mwid"])
vparameters$prop.D1[gp.gn$het.m] <- with(vlist$prop.D, pe[cd4 == 1 & risk == "het" & gender == "m"])
vparameters$prop.D1[gp.gn$het.f] <- with(vlist$prop.D, pe[cd4 == 1 & risk == "het" & gender == "f"])

vparameters$prop.D2[gp.gn$pwid] <- with(vlist$prop.D, pe[cd4 == 2 & risk == "pwid"])[1:3]
vparameters$prop.D2[gp.gn$msm] <- with(vlist$prop.D, pe[cd4 == 2 & risk == "msm"])
vparameters$prop.D2[gp.gn$mwid] <- with(vlist$prop.D, pe[cd4 == 2 & risk == "mwid"])
vparameters$prop.D2[gp.gn$het.m] <- with(vlist$prop.D, pe[cd4 == 2 & risk == "het" & gender == "m"])
vparameters$prop.D2[gp.gn$het.f] <- with(vlist$prop.D, pe[cd4 == 2 & risk == "het" & gender == "f"])


vparameters$prop.D3[gp.gn$pwid] <- with(vlist$prop.D, pe[cd4 == 3 & risk == "pwid"])[1:3]
vparameters$prop.D3[gp.gn$msm] <- with(vlist$prop.D, pe[cd4 == 3 & risk == "msm"])
vparameters$prop.D3[gp.gn$mwid] <- with(vlist$prop.D, pe[cd4 == 3 & risk == "mwid"])
vparameters$prop.D3[gp.gn$het.m] <- with(vlist$prop.D, pe[cd4 == 3 & risk == "het" & gender == "m"])
vparameters$prop.D3[gp.gn$het.f] <- with(vlist$prop.D, pe[cd4 == 3 & risk == "het" & gender == "f"])


# viral_sup: 3 to  18,
# viral suppression
vparameters$viral_sup <- numeric(length(vparameters$names18))
vparameters$viral_sup[gp18.gn$white] <- with(vlist$viral_sup, pe[ethnicity == "w"])
vparameters$viral_sup[gp18.gn$black] <- with(vlist$viral_sup, pe[ethnicity == "b"])
vparameters$viral_sup[gp18.gn$hisp]  <- with(vlist$viral_sup, pe[ethnicity == "h"])

# prop.T1, prop.T2, prop.T3: 30*1 each (initial proportion of different T stages)
vparameters$prop.T1 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.T2 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.T3 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)

vparameters$prop.T1[gp.gn$msm] <- with(vlist$prop.T, pe[cd4 == 1 & risk == "msm"])
vparameters$prop.T1[gp.gn$mwid] <- with(vlist$prop.T, pe[cd4 == 1 & risk == "mwid"])
vparameters$prop.T1[gp.gn$pwid] <- with(vlist$prop.T, pe[cd4 == 1 & risk == "pwid"])
vparameters$prop.T1[gp.gn$het.m] <- with(vlist$prop.T, pe[cd4 == 1 & risk == "het" & gender == "m"])
vparameters$prop.T1[gp.gn$het.f] <- with(vlist$prop.T, pe[cd4 == 1 & risk == "het" & gender == "f"])

vparameters$prop.T2[gp.gn$msm] <- with(vlist$prop.T, pe[cd4 == 2 & risk == "msm"])
vparameters$prop.T2[gp.gn$mwid] <- with(vlist$prop.T, pe[cd4 == 2 & risk == "mwid"])
vparameters$prop.T2[gp.gn$pwid] <- with(vlist$prop.T, pe[cd4 == 2 & risk == "pwid"])
vparameters$prop.T2[gp.gn$het.m] <- with(vlist$prop.T, pe[cd4 == 2 & risk == "het" & gender == "m"])
vparameters$prop.T2[gp.gn$het.f] <- with(vlist$prop.T, pe[cd4 == 2 & risk == "het" & gender == "f"])

vparameters$prop.T3 <- 1-vparameters$prop.T1-vparameters$prop.T2


# prop.O: 12*3 to 30*3 (initial proportion of different O stages)
vparameters$prop.O1 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.O2 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)
vparameters$prop.O3 <- matrix(0, nrow = length(vparameters$names.gp), ncol = 1)

vparameters$prop.O1[gp.gn$msm] <- with(vlist$prop.O, pe[cd4 == 1 & risk == "msm"])
vparameters$prop.O1[gp.gn$mwid] <- with(vlist$prop.O, pe[cd4 == 1 & risk == "mwid"])
vparameters$prop.O1[gp.gn$pwid] <- with(vlist$prop.O, pe[cd4 == 1 & risk == "pwid"])
vparameters$prop.O1[gp.gn$het.m] <- with(vlist$prop.O, pe[cd4 == 1 & risk == "het" & gender == "m"])
vparameters$prop.O1[gp.gn$het.f] <- with(vlist$prop.O, pe[cd4 == 1 & risk == "het" & gender == "f"])

vparameters$prop.O2[gp.gn$msm] <- with(vlist$prop.O, pe[cd4 == 2 & risk == "msm"])
vparameters$prop.O2[gp.gn$mwid] <- with(vlist$prop.O, pe[cd4 == 2 & risk == "mwid"])
vparameters$prop.O2[gp.gn$pwid] <- with(vlist$prop.O, pe[cd4 == 2 & risk == "pwid"])
vparameters$prop.O2[gp.gn$het.m] <- with(vlist$prop.O, pe[cd4 == 2 & risk == "het" & gender == "m"])
vparameters$prop.O2[gp.gn$het.f] <- with(vlist$prop.O, pe[cd4 == 2 & risk == "het" & gender == "f"])

vparameters$prop.O3[gp.gn$msm] <- with(vlist$prop.O, pe[cd4 == 3 & risk == "msm"])
vparameters$prop.O3[gp.gn$mwid] <- with(vlist$prop.O, pe[cd4 == 3 & risk == "mwid"])
vparameters$prop.O3[gp.gn$pwid] <- with(vlist$prop.O, pe[cd4 == 3 & risk == "pwid"])
vparameters$prop.O3[gp.gn$het.m] <- with(vlist$prop.O, pe[cd4 == 3 & risk == "het" & gender == "m"])
vparameters$prop.O3[gp.gn$het.f] <- with(vlist$prop.O, pe[cd4 == 3 & risk == "het" & gender == "f"])

## Parameter in transmissibility domain##

# transmission probability through sex (FM) by acute and 3 CD4 stages
vparameters$sigmaO.FM <- numeric(4)
vparameters$sigmaO.FM[2:4] <- vlist$sigmaO.FM$pe
vparameters$sigmaO.FM[1] <- vparameters$trans.acute * vparameters$sigmaO.FM[2]

# transmission probability through sex (MF) by acute and 3 CD4 stages
vparameters$sigmaO.MF <- numeric(4)
vparameters$sigmaO.MF[2:4] <- vlist$sigmaO.MF$pe
vparameters$sigmaO.MF[1] <- vparameters$trans.acute * vparameters$sigmaO.MF[2]

# transmission probability through sex (MM) by acute and 3 CD4 stages
vparameters$sigmaS <- numeric(4)
vparameters$sigmaS[2:4] <- vlist$sigmaS$pe
vparameters$sigmaS[1] <- vparameters$trans.acute * vparameters$sigmaS[2]


# probability of transmission per shared injection by acute and 3 CD4 stages
vparameters$tau <- numeric(4)
vparameters$tau[2:4] <- vlist$tau$pe
vparameters$tau[1] <- vparameters$trans.acute * vparameters$tau[2]

## Parameter in assortativeness domain ##
# names.e groups by ethnicity
# opposite sex assortative coefficient for each ethnicity by risk level
vparameters$ass.eO <- numeric(length(vparameters$names.e))
vparameters$ass.eO[grep("low", vparameters$names.e)] <- with(vlist$ass.eO, pe[sex.risk.age.group == "low"])
vparameters$ass.eO[grep("high", vparameters$names.e)] <- with(vlist$ass.eO, pe[sex.risk.age.group == "high"])

# homosexual assortative coefficient by ethnicity and risk level
vparameters$ass.eS <- numeric(length(vparameters$names.e2))
vparameters$ass.eS[grep("old", vparameters$names.e2)] <- with(vlist$ass.eS, pe[sex.risk.age.group == "old"])
vparameters$ass.eS[grep("young", vparameters$names.e2)] <- with(vlist$ass.eS, pe[sex.risk.age.group == "young"])


# Set up time to avoid hard-coded year
length     <- length(vlist$prop.ssp$year)
first.year <- vlist$prop.ssp$year[1]
last.year  <- vlist$prop.ssp$year[length/2+1]
time = as.numeric(last.year) - as.numeric(first.year) + 1 #Plus 1 to ensure the time accounted for first & last year and not just the difference
# Populate data for first and last year
# SSP only apply to those on PWID/MWID, the rest =0
vparameters$prop.ssp.m <- matrix(0, nrow = time, ncol = length(vparameters$names.gp))
vparameters$prop.ssp.m[1,gp.gn$pwid]    <- with(vlist$prop.ssp, pe[risk == "pwid" & year == first.year])
vparameters$prop.ssp.m[1,gp.gn$mwid]    <- with(vlist$prop.ssp, pe[risk == "mwid" & year == first.year])
vparameters$prop.ssp.m[time,gp.gn$pwid] <- with(vlist$prop.ssp, pe[risk == "pwid" & year == last.year])
vparameters$prop.ssp.m[time,gp.gn$mwid] <- with(vlist$prop.ssp, pe[risk == "mwid" & year == last.year])
# Populate data for years in-between
for (i in 1:length(vparameters$names.gp)) {
	vparameters$prop.ssp.m[, i] <- seq(
		from = vparameters$prop.ssp.m[1,][i],
		to = vparameters$prop.ssp.m[time,][i],
		length.out = time
	)
}

## PrEP parameter ##
# Covid disruption
covid_prep <- numeric(length(vparameters$names.gp))
covid_prep[intersect(gp.gn$white, gp.gn$young)] <- with(vlist$covid_prep, pe[ethnicity == "w" & sex.risk.age.group == "young"])
covid_prep[intersect(gp.gn$black, gp.gn$young)] <- with(vlist$covid_prep, pe[ethnicity == "b" & sex.risk.age.group == "young"])
covid_prep[intersect(gp.gn$hisp , gp.gn$young)] <- with(vlist$covid_prep, pe[ethnicity == "h" & sex.risk.age.group == "young"])
covid_prep[intersect(gp.gn$white, gp.gn$old)] <- with(vlist$covid_prep, pe[ethnicity == "w" & sex.risk.age.group == "old"])
covid_prep[intersect(gp.gn$black, gp.gn$old)] <- with(vlist$covid_prep, pe[ethnicity == "b" & sex.risk.age.group == "old"])
covid_prep[intersect(gp.gn$hisp , gp.gn$old)] <- with(vlist$covid_prep, pe[ethnicity == "h" & sex.risk.age.group == "old"])
covid_prep[intersect(gp.gn$white, union(gp.gn$pwid,gp.gn$het))] <- with(vlist$covid_prep, pe[ethnicity == "w" & risk == "het"])
covid_prep[intersect(gp.gn$black, union(gp.gn$pwid,gp.gn$het))] <- with(vlist$covid_prep, pe[ethnicity == "b" & risk == "het"])
covid_prep[intersect(gp.gn$hisp , union(gp.gn$pwid,gp.gn$het))] <- with(vlist$covid_prep, pe[ethnicity == "h" & risk == "het"])

# prop.prep: Prop on PrEP for all strata (NEW parameter)
# Set up time to avoid hard-coded year
first.year <- 2017
last.year  <- 2019
time = as.numeric(last.year) - as.numeric(first.year) + 1 #Plus 1 to ensure the time accounted for first & last year and not just the difference
# Populate data for first and last year
prop.prep.m1 <- matrix(0, nrow = time, ncol = length(vparameters$names.gp))

prop.prep.m1[1,intersect(gp.gn$msm, gp.gn$old)]    <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "old"   & year == first.year])
prop.prep.m1[1,intersect(gp.gn$msm, gp.gn$young)]  <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "young" & year == first.year])
prop.prep.m1[1,intersect(gp.gn$mwid, gp.gn$old)]   <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "old"   & year == first.year])
prop.prep.m1[1,intersect(gp.gn$mwid, gp.gn$young)] <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "young" & year == first.year])
prop.prep.m1[1,gp.gn$pwid.m]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "m" & year == first.year])
prop.prep.m1[1,gp.gn$pwid.f]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "f" & year == first.year])
prop.prep.m1[1,gp.gn$L.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "low"    & year == first.year])
prop.prep.m1[1,gp.gn$H.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "high"   & year == first.year])
prop.prep.m1[1,gp.gn$L.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "low"    & year == first.year])
prop.prep.m1[1,gp.gn$H.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "high"   & year == first.year])

prop.prep.m1[time,intersect(gp.gn$msm, gp.gn$old)]    <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "old"   & year == last.year])
prop.prep.m1[time,intersect(gp.gn$msm, gp.gn$young)]  <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "young" & year == last.year])
prop.prep.m1[time,intersect(gp.gn$mwid, gp.gn$old)]   <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "old"   & year == last.year])
prop.prep.m1[time,intersect(gp.gn$mwid, gp.gn$young)] <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "young" & year == last.year])
prop.prep.m1[time,gp.gn$pwid.m]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "m" & year == last.year])
prop.prep.m1[time,gp.gn$pwid.f]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "f" & year == last.year])
prop.prep.m1[time,gp.gn$L.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "low"    & year == last.year])
prop.prep.m1[time,gp.gn$H.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "high"   & year == last.year])
prop.prep.m1[time,gp.gn$L.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "low"    & year == last.year])
prop.prep.m1[time,gp.gn$H.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "high"   & year == last.year])

for (i in 1:length(vparameters$names.gp)) {
	prop.prep.m1[, i] <- seq(
		from = prop.prep.m1[1,][i],
		to = prop.prep.m1[time,][i],
		length.out = time
	)
}

prop.prep.covid <- prop.prep.m1[3, ]*(1+covid_prep)
# Set up time to avoid hard-coded year
first.year <- 2021
last.year  <- 2022
time = as.numeric(last.year) - as.numeric(first.year) + 1 #Plus 1 to ensure the time accounted for first & last year and not just the difference
# Populate data for first and last year
prop.prep.m2 <- matrix(0, nrow = time, ncol = length(vparameters$names.gp))

prop.prep.m2[1,intersect(gp.gn$msm, gp.gn$old)]    <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "old"   & year == first.year])
prop.prep.m2[1,intersect(gp.gn$msm, gp.gn$young)]  <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "young" & year == first.year])
prop.prep.m2[1,intersect(gp.gn$mwid, gp.gn$old)]   <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "old"   & year == first.year])
prop.prep.m2[1,intersect(gp.gn$mwid, gp.gn$young)] <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "young" & year == first.year])
prop.prep.m2[1,gp.gn$pwid.m]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "m" & year == first.year])
prop.prep.m2[1,gp.gn$pwid.f]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "f" & year == first.year])
prop.prep.m2[1,gp.gn$L.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "low"    & year == first.year])
prop.prep.m2[1,gp.gn$H.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "high"   & year == first.year])
prop.prep.m2[1,gp.gn$L.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "low"    & year == first.year])
prop.prep.m2[1,gp.gn$H.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "high"   & year == first.year])

prop.prep.m2[time,intersect(gp.gn$msm, gp.gn$old)]    <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "old"   & year == last.year])
prop.prep.m2[time,intersect(gp.gn$msm, gp.gn$young)]  <- with(vlist$prop.prep, pe[risk == "msm"  & gender == "m" & sex.risk.age.group == "young" & year == last.year])
prop.prep.m2[time,intersect(gp.gn$mwid, gp.gn$old)]   <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "old"   & year == last.year])
prop.prep.m2[time,intersect(gp.gn$mwid, gp.gn$young)] <- with(vlist$prop.prep, pe[risk == "mwid" & gender == "m" & sex.risk.age.group == "young" & year == last.year])
prop.prep.m2[time,gp.gn$pwid.m]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "m" & year == last.year])
prop.prep.m2[time,gp.gn$pwid.f]                       <- with(vlist$prop.prep, pe[risk == "pwid" & gender == "f" & year == last.year])
prop.prep.m2[time,gp.gn$L.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "low"    & year == last.year])
prop.prep.m2[time,gp.gn$H.het.m]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "m" & sex.risk.age.group == "high"   & year == last.year])
prop.prep.m2[time,gp.gn$L.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "low"    & year == last.year])
prop.prep.m2[time,gp.gn$H.het.f]                      <- with(vlist$prop.prep, pe[risk == "het" & gender == "f" & sex.risk.age.group == "high"   & year == last.year])

for (i in 1:length(vparameters$names.gp)) {
	prop.prep.m2[, i] <- seq(
		from = prop.prep.m2[1,][i],
		to = prop.prep.m2[time,][i],
		length.out = time
	)
}

vparameters$prop.prep.m <- rbind(prop.prep.m1, prop.prep.covid,prop.prep.m2)
rownames(vparameters$prop.prep.m) <- c("2017", "2018", "2019", "covid","2021","2022")


## OAT parameter ##
# Set up time to avoid hard-coded year
length     <- length(vlist$prop.oat$year)
first.year <- vlist$prop.oat$year[1]
last.year  <- vlist$prop.oat$year[length/2+1]
time = as.numeric(last.year) - as.numeric(first.year) + 1 #Plus 1 to ensure the time accounted for first & last year and not just the difference
# Populate data for first and last year
# OAT only apply to those on PWID/MWID, the rest =0
vparameters$prop.oat.m <- matrix(0, nrow = time, ncol = length(vparameters$names.gp))
vparameters$prop.oat.m[1,gp.gn$pwid]    <- with(vlist$prop.oat, pe[risk == "pwid" & year == first.year])
vparameters$prop.oat.m[1,gp.gn$mwid]    <- with(vlist$prop.oat, pe[risk == "mwid" & year == first.year])
vparameters$prop.oat.m[time,gp.gn$pwid] <- with(vlist$prop.oat, pe[risk == "pwid" & year == last.year])
vparameters$prop.oat.m[time,gp.gn$mwid] <- with(vlist$prop.oat, pe[risk == "mwid" & year == last.year])

# Populate data for years in-between
for (i in 1:length(vparameters$names.gp)) {
	vparameters$prop.oat.m[, i] <- seq(
		from = vparameters$prop.oat.m[1,][i],
		to = vparameters$prop.oat.m[time,][i],
		length.out = time
	)
}

#### Screening rate parameter ####:

# psi: directly pasted from the model calibration, no need to manipulate (baseline testing rate)
psi.initial <- numeric(length(vparameters$names.gp))
psi.initial[gp.gn$white] <- with(vlist$psi, pe[ethnicity == "w"])
psi.initial[gp.gn$black] <- with(vlist$psi, pe[ethnicity == "b"])
psi.initial[gp.gn$hisp ] <- with(vlist$psi, pe[ethnicity == "h"])

# Covid disruption
covid_psi <- numeric(length(vparameters$names.gp))
covid_psi[gp.gn$white] <- with(vlist$covid_psi, pe[ethnicity == "w"])
covid_psi[gp.gn$black] <- with(vlist$covid_psi, pe[ethnicity == "b"])
covid_psi[gp.gn$hisp ] <- with(vlist$covid_psi, pe[ethnicity == "h"])

psi.covid <- psi.initial * (1+covid_psi)

psi.covid.m <- matrix(0, nrow = 7, ncol = length(vparameters$names.gp))
for (i in 1:length(vparameters$names.gp)) {
	psi.covid.m[, i] <- seq(
		from = psi.covid[i],
		to   = psi.initial[i],
		length.out = 7
	)
}

vparameters$psi.m <- rbind(psi.initial,psi.covid.m)
rownames(vparameters$psi.m ) <- c("before covid", "Mar-20", "Apr-20", "May-20", "Jun-20", "Jul-20", "Aug-20","after covid")

rm(p1,alpha1_,alpha2_,alpha3_,phi1_,phi2_,phi3_)
