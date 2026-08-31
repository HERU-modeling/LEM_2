library(openxlsx)
# Prevalence
diag.obs  = data.matrix(subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="TotalDiagnosedPLHIV",  colNames =TRUE)), city == CITY & year != 2016)[ , c("W","B","H")])
# New diagnoses
ndiag.obs = data.matrix(subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="NewDiagnoses", colNames =TRUE)), city == CITY & year != 2016)[ , c("M.HET","F.HET","X15.34.MSM","X35.64.MSM","X15.34.MWID","X35.64.MWID","PWID","W","B","H")])
# (Proportion) Late diagnoses
late_diag.obs = data.matrix(subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="LateDiagnoses", colNames =TRUE)), city == CITY & year != 2016)[ , c("W","B","H")])
# (Proportion) Viral suppression
viral_sup.obs = data.matrix(subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="ViralSuppression", colNames =TRUE)), city == CITY & year != 2016)[ , c("W","B","H")])

obs.inc.all = subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="obs.inc.all", colNames =TRUE)), city == CITY)[ , c("value","low","high")]
obs.inc.msm = subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx", sheet="obs.inc.msm", colNames =TRUE)), city == CITY)[ , c("value","low","high")]
#####################################################
myFile.name <- c("Data_files/DM/Evidence-Inputs-Master")


WB <- openxlsx::loadWorkbook(paste0(myFile.name, ".xlsx"))
#From par_fix
namelist <- read.xlsx(WB, sheet = "common")
names18 <- namelist$names18[!is.na(namelist$names18)]

#####################################################
calib.target <- list(diag.obs = diag.obs, ndiag.obs = ndiag.obs, late_diag.obs = late_diag.obs, viral_sup.obs = viral_sup.obs)

valid.target <- list (obs.inc.all = obs.inc.all, obs.inc.msm = obs.inc.msm)

calib.target$ndiag.obs <- cbind(
	calib.target$ndiag.obs,
	Young = calib.target$ndiag.obs[, "X15.34.MSM"] + calib.target$ndiag.obs[, "X15.34.MWID"],
	Old   = calib.target$ndiag.obs[, "X35.64.MSM"] + calib.target$ndiag.obs[, "X35.64.MWID"]
) 	# add new cols
calib.target$ndiag.obs <- calib.target$ndiag.obs[, c("M.HET", "F.HET", "Young", "Old", "PWID", "W", "B", "H")]
#####################################################



