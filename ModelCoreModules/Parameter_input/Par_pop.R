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
vparameters$prop.msm.among.male <- rep(vlist$prop.msm.among.male$pe, 3)

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

# Derive model population initials #
source("ModelCoreModules/Parameter_input/Function-model.initial.R")
init <- model.initial(par = vparameters, diag18 = vparameters$pop.diag, gp.gn = gp.gn) # 30*12 initials, initial population for 30 groups at 12 different states, inits2 within source(Function-model.initial.R)

#####
inits <- cbind(init, inc_bo = 0, inc_bs = 0, inc_g = 0, diag = 0, diag_late = 0, positive_test = 0, viral_sup = 0, death = 0) # adding states to track # infections (by 3 routes of transmission), new diagnoses and deaths among PLHIV,T: this is how x becomes 42X25
x <- as.vector(t(inits)) # ode function requires init as vector (x)
vparameters$init.group.prop <- as.vector(rowSums(init) / sum(init)) 

# initial population for 18 groups (collapsing young/old, low/high)
init.collap <- numeric(18)
rname <- gsub(paste(c("/low", "/high", "/old","/young"), collapse = "|"), "", vparameters$names.gp)
for (i in 1:18) {
  ind <- which(rname %in% vparameters$names18[i])
  init.collap[i] <- sum(init[ind, ])
}

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
}

