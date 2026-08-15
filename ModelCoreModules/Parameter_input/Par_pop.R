## Parameter in population demographics domain ##

# pop.male (male population by ethnicity)
vparameters$pop.male <- numeric(3)
vparameters$pop.male[1] <- with(vlist$pop.total, pe[ethnicity == "w" & gender == "m" & year == "2016"])
vparameters$pop.male[2] <- with(vlist$pop.total, pe[ethnicity == "b" & gender == "m" & year == "2016"])
vparameters$pop.male[3] <- with(vlist$pop.total, pe[ethnicity == "h" & gender == "m" & year == "2016"])

# pop.female (female population by ethnicity)
vparameters$pop.female <- numeric(3)
vparameters$pop.female[1] <- with(vlist$pop.total, pe[ethnicity == "w" & gender == "f" & year == "2016"])
vparameters$pop.female[2] <- with(vlist$pop.total, pe[ethnicity == "b" & gender == "f" & year == "2016"])
vparameters$pop.female[3] <- with(vlist$pop.total, pe[ethnicity == "h" & gender == "f" & year == "2016"])


# prop.msm.among.male: 1 to 3 (proportion of MSM for each ethnicity)
#if (CITY == "NYC") {
#  region.prop_m.msm <- numeric(5)
#  region.prop_m.msm[1:4] <- vlist$prop.msm.among.male$pe[1:4]
#  region.prop_m.msm[5] <- with(vlist$prop.msm.among.male, pe[region == "brooklyn"])
#  prop.msm.among.male <- sum(region.prop_m.msm * as.numeric(vlist$prop.msm.among.male$male.pop)) / sum(as.numeric(vlist$prop.msm.among.male$male.pop))
#  vparameters$prop.msm.among.male <- rep(prop.msm.among.male, 3)
#} else {
vparameters$prop.msm.among.male <- rep(vlist$prop.msm.among.male$pe, 3)
#}

# prop.male.among.pwid: 1 to 3 (Male proportion among all ethnic groups for PWID)
vparameters$prop.male.among.pwid <- rep(vlist$prop.male.among.pwid$pe, 3)


# prop.pwid.among.total (PWID proportion among total population by ethnicity, equally by male and female)
vparameters$prop.pwid.among.total <- numeric(3)
vparameters$prop.pwid.among.total[1] <- with(vlist$prop.pwid.among.total, pe[ethnicity == "w"])
vparameters$prop.pwid.among.total[2] <- with(vlist$prop.pwid.among.total, pe[ethnicity == "b"])
vparameters$prop.pwid.among.total[3] <- with(vlist$prop.pwid.among.total, pe[ethnicity == "h"])


# prop.mpwid.among.pwid (proportion of msm and inject drugs among all pwid population by ethnicity)
vparameters$prop.mpwid.among.pwid <- numeric(3)
vparameters$prop.mpwid.among.pwid[1] <- with(vlist$prop.mpwid.among.pwid, pe[ethnicity == "w"])
vparameters$prop.mpwid.among.pwid[2] <- with(vlist$prop.mpwid.among.pwid, pe[ethnicity == "b"])
vparameters$prop.mpwid.among.pwid[3] <- with(vlist$prop.mpwid.among.pwid, pe[ethnicity == "h"])


# prop.mpwid.among.msm (proportion of msm and inject drugs among all msm population by ethnicity)
vparameters$prop.mpwid.among.msm <- numeric(3)
vparameters$prop.mpwid.among.msm[1] <- with(vlist$prop.mpwid.among.msm, pe[ethnicity == "w"])
vparameters$prop.mpwid.among.msm[2] <- with(vlist$prop.mpwid.among.msm, pe[ethnicity == "b"])
vparameters$prop.mpwid.among.msm[3] <- with(vlist$prop.mpwid.among.msm, pe[ethnicity == "h"])


## Read in target calibration population data (2011-2015) _move to load_target?
# total diagnosed cases, 18 groups
#diag18.obs <- data.matrix(subset(data.frame(read.xlsx("Data_Files/DM/Target.xlsx", sheet = "diag18.obs", colNames = T)), city == CITY)[, -1])
#ndiag18.obs<- data.matrix(subset(data.frame(read.xlsx("Data_Files/DM/Target.xlsx", sheet = "ndiag18.obs", colNames = T)), city == CITY)[, -1]) #new diagnoses, 18 groups
#death18.obs<- data.matrix(subset(data.frame(read.xlsx("Data_Files/DM/Target.xlsx", sheet = "death18.obs", colNames =T)), city == CITY)[ ,-1]) #all-casue deaths, 18 groups
#obs.inc.all<- subset(data.frame(read.xlsx("Data_Files/DM/Target.xlsx", sheet="inc.all.obs", colNames =T)), city == CITY)[ , -1]              #new incidence, total and range
#obs.inc.msm<- subset(data.frame(read.xlsx("Data_Files/DM/Target.xlsx", sheet="inc.msm.obs", colNames =T)), city == CITY)[ , -1]              #new incidence, MSM and range

# Derive model population initials #
source("ModelCoreModules/Parameter_input/Function-model.initial.R")
#init <- model.initial(par = vparameters, diag18 = diag18.obs[1, ], gp.gn = gp.gn) # Old with above diag18.obs

init <- model.initial(par = vparameters, diag18 = vparameters$pop.diag, gp.gn = gp.gn) # 30*12 initials, initial population for 30 groups at 12 different states, inits2 within source(Function-model.initial.R)

#####
inits <- cbind(init, inc_bo = 0, inc_bs = 0, inc_g = 0, diag = 0, diag_late = 0, positive_test = 0, viral_sup = 0, death = 0, D_T=0, I_T=0 , I_D = 0, death_D = 0) # adding states to track # infections (by 3 routes of transmission), new diagnoses and deaths among PLHIV,T: this is how x becomes 42X25
x <- as.vector(t(inits)) # ode function requires init as vector (x)
vparameters$init.group.prop <- as.vector(rowSums(init) / sum(init)) #T: basically a proportion of the number of people of each group over the total people, sum() to 1

# initial population for 18 groups (collapsing young/old, low/high)
init.collap <- numeric(18)
rname <- gsub(paste(c("/low", "/high", "/old","/young"), collapse = "|"), "", vparameters$names.gp)
for (i in 1:18) {
  ind <- which(rname %in% vparameters$names18[i])
  init.collap[i] <- sum(init[ind, ])
}
# T:code aggregates the initial population data for 18 groups by collapsing groups with similar names (after removing "/OAT", "/low", or "/high" substrings) and summing up their populations
#browser()

# These two variables are required for estimating OAT and PrEP entry rates in ODE module (only use initial values) #
# number of PWID (include MWID) within the 9 collapsed groups:
# "MSM/PWID/white","MSM/PWID/black","MSM/PWID/hispanic","PWID/male/white", "PWID/male/black",
# "PWID/male/hispanic", "PWID/female/white","PWID/female/black","PWID/female/hispanic"
#vparameters$pop.pwid <- init.collap[grep("PWID", vparameters$names18)]

# Create prop.oat for use in Par_int and Par_CEA to replace parts that used to indexed by OAT, also resuse later
#nOAT.m   <- vparameters$nOAT.m
#pop.pwid <- vparameters$pop.pwid
#nOAT     <- round(c(nOAT.m[nrow(nOAT.m), 1:3] * c(pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6]), 1-pop.pwid[1:3]/(pop.pwid[1:3] + pop.pwid[4:6])), nOAT.m[nrow(nOAT.m), 4:6]))
#vparameters$prop.oat <- nOAT / pop.pwid

# Repeat first 3 values for MSM/PWID/old&young, revisit for update time-varying prop.oat
#vparameters$prop.oat.12 <- c(vparameters$prop.oat[1:3], vparameters$prop.oat)

# initial number of susceptible high-risk MSM (include MWID) within the 3 collapsed groups (race only) - Now Young/MSM
#vparameters$msm.h.scep <- rowSums(matrix(rowSums(inits[gp.gn$Y.allmsm, 1:3]), nrow = 3))

vparameters$prop.adj <- FALSE # force the risk group proportions to remain constant
vparameters$inf.prop.adj <- FALSE # force the risk level proportions to remain constant among infected
vparameters$bal <- TRUE # balance # of sexual partnership between males & females

############################
#pop.change: Pop growth is constant, convert pop.total to pop.change 
length     <- length(vlist$pop.total$year)
first.year <- vlist$pop.total$year[1]
last.year  <- vlist$pop.total$year[length/2+1]
time = as.numeric(last.year) - as.numeric(first.year)
pop.change.6    <- numeric(6)
pop.change.6[1] <- (with(vlist$pop.total, pe[ethnicity == "w" & gender == "m" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "w" & gender == "m" & year == first.year]))/time
pop.change.6[2] <- (with(vlist$pop.total, pe[ethnicity == "b" & gender == "m" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "b" & gender == "m" & year == first.year]))/time
pop.change.6[3] <- (with(vlist$pop.total, pe[ethnicity == "h" & gender == "m" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "h" & gender == "m" & year == first.year]))/time
pop.change.6[4] <- (with(vlist$pop.total, pe[ethnicity == "w" & gender == "f" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "w" & gender == "f" & year == first.year]))/time
pop.change.6[5] <- (with(vlist$pop.total, pe[ethnicity == "b" & gender == "f" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "b" & gender == "f" & year == first.year]))/time
pop.change.6[6] <- (with(vlist$pop.total, pe[ethnicity == "h" & gender == "f" & year == last.year]) - with(vlist$pop.total, pe[ethnicity == "h" & gender == "f" & year == first.year]))/time
#browser()

############################
# Convert from 6 to 18 stratification: Inspired from model.initial
vparameters$pop.change.18 <- numeric(length(vparameters$names18))

for (i in 1:3) {
	
	## Indices in pop.change.6
	male_idx   <- i
	female_idx <- i + 3
	
	## Male / Female totals
	pop_male <- pop.change.6[male_idx]
	pop_female <- pop.change.6[female_idx]
	pop_total <- pop_male + pop_female
	
	pop_msm    <- pop_male * vparameters$prop.msm.among.male[i]
	
	pop_pwid_f <- pop_total * vparameters$prop.pwid.among.total[i] *
		(1 - vparameters$prop.male.among.pwid[i])
	
	pop_pwid_m <- pop_total * vparameters$prop.pwid.among.total[i] *
		vparameters$prop.male.among.pwid[i]
	
	pop_mwid <- pop_pwid_m *
		vparameters$prop.mpwid.among.pwid[i]
	
	pop_pwid_m <- pop_pwid_m - pop_mwid
	pop_msm    <- pop_msm    - pop_mwid
	
	pop_het_m <- pop_male - pop_msm - pop_pwid_m - pop_mwid
	pop_het_f <- pop_female - pop_pwid_f
	
	vparameters$pop.change.18[gp18.gn$msm[i]]    <- pop_msm
	vparameters$pop.change.18[gp18.gn$mwid[i]]   <- pop_mwid
	vparameters$pop.change.18[intersect(gp18.gn$pwid,gp18.gn$male  )[i]]   <- pop_pwid_m
	vparameters$pop.change.18[intersect(gp18.gn$pwid,gp18.gn$female)[i]]   <- pop_pwid_f
	vparameters$pop.change.18[gp18.gn$het.m[i]]    <- pop_het_m
	vparameters$pop.change.18[gp18.gn$het.f[i]]    <- pop_het_f
#browser()	
}
#browser()  # Check sum(vparameters$pop.change.18[gp18.gn$white]), sum(vparameters$pop.change.18[gp18.gn$black]),... same as calculate from excel
############################
# # Manually inflate Pop entry/growth for MSM, MSM/PWID, PWID
# # M/W/MSM-PWID	
# vparameters$pop.change.18[gp18.gn$mwid[1]] <- vparameters$pop.change.18[gp18.gn$mwid[1]] * ifelse(vparameters$pop.inflate.rate_msm > vparameters$pop.inflate.rate_pwid, vparameters$pop.inflate.rate_msm, vparameters$pop.inflate.rate_pwid)
# # M/W/PWID	
# vparameters$pop.change.18[gp18.gn$pwid[1]] <- vparameters$pop.change.18[gp18.gn$pwid[1]] * vparameters$pop.inflate.rate_pwid
# # M/W/MSM
# vparameters$pop.change.18[gp18.gn$msm[1]]  <- vparameters$pop.change.18[gp18.gn$msm[1]]  * vparameters$pop.inflate.rate_msm
# # M/B/MSM-PWID
# vparameters$pop.change.18[gp18.gn$mwid[2]] <- vparameters$pop.change.18[gp18.gn$mwid[2]] * ifelse(vparameters$pop.inflate.rate_msm > vparameters$pop.inflate.rate_pwid, vparameters$pop.inflate.rate_msm, vparameters$pop.inflate.rate_pwid)
# # M/B/PWID
# vparameters$pop.change.18[gp18.gn$pwid[2]] <- vparameters$pop.change.18[gp18.gn$pwid[2]] * vparameters$pop.inflate.rate_pwid
# # M/B/MSM
# vparameters$pop.change.18[gp18.gn$msm[2]]  <- vparameters$pop.change.18[gp18.gn$msm[2]]  * vparameters$pop.inflate.rate_msm
# # M/H/MSM-PWID
# vparameters$pop.change.18[gp18.gn$mwid[3]] <- vparameters$pop.change.18[gp18.gn$mwid[3]] * ifelse(vparameters$pop.inflate.rate_msm > vparameters$pop.inflate.rate_pwid, vparameters$pop.inflate.rate_msm, vparameters$pop.inflate.rate_pwid)
# # M/H/PWID 
# vparameters$pop.change.18[gp18.gn$pwid[3]] <- vparameters$pop.change.18[gp18.gn$pwid[3]] * vparameters$pop.inflate.rate_pwid
# # M/H/MSM
# vparameters$pop.change.18[gp18.gn$msm[3]]  <- vparameters$pop.change.18[gp18.gn$msm[3]]  * vparameters$pop.inflate.rate_msm
# # F/W/PWID
# vparameters$pop.change.18[gp18.gn$pwid[4]] <- vparameters$pop.change.18[gp18.gn$pwid[4]] * vparameters$pop.inflate.rate_pwid
# # F/B/PWID
# vparameters$pop.change.18[gp18.gn$pwid[5]] <- vparameters$pop.change.18[gp18.gn$pwid[5]] * vparameters$pop.inflate.rate_pwid
# # F/H/PWID
# vparameters$pop.change.18[gp18.gn$pwid[6]] <- vparameters$pop.change.18[gp18.gn$pwid[6]] * vparameters$pop.inflate.rate_pwid
# #browser()
