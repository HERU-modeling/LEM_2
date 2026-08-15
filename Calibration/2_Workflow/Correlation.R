# Input Average of all outcomes across 32 jurisdictions
Jurisdictions <- all.cities
#all.cities #c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL")

all_data_list <- list()

for (j in seq_along(Jurisdictions)) {
	
	diag.obs  <- data.matrix(
		subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx",
																sheet="TotalDiagnosedPLHIV", colNames=TRUE)),
					 city == Jurisdictions[j] & year != 2016)[ , c("W","B","H")]
	)
	
	ndiag.obs <- data.matrix(
		subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx",
																sheet="NewDiagnoses", colNames=TRUE)),
					 city == Jurisdictions[j] & year != 2016)[ , c("M.HET","F.HET","X15.34.MSM","X35.64.MSM","X15.34.MWID","X35.64.MWID","PWID","W","B","H")]
	)
	ndiag.obs <- cbind(
		ndiag.obs,
		Young = ndiag.obs[, "X15.34.MSM"] + ndiag.obs[, "X15.34.MWID"],
		Old   = ndiag.obs[, "X35.64.MSM"] + ndiag.obs[, "X35.64.MWID"]
	)
	ndiag.obs <- ndiag.obs[, c("M.HET", "F.HET", "Young", "Old", "PWID", "W", "B", "H")]
	
	late_diag.obs <- data.matrix(
		subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx",
																sheet="LateDiagnoses", colNames=TRUE)),
					 city == Jurisdictions[j] & year != 2016)[ , c("W","B","H")]
	)
	
	viral_sup.obs <- data.matrix(
		subset(data.frame(read.xlsx("Data_files/Target/calibration_targets.xlsx",
																sheet="ViralSuppression", colNames=TRUE)),
					 city == Jurisdictions[j] & year != 2016)[ , c("W","B","H")]
	)
	
	# Combine horizontally (6 × 12)
	Y_j <- cbind(
		diag.obs,
		ndiag.obs,
		late_diag.obs,
		viral_sup.obs
	)
	
	colnames(Y_j) <- c(
		"diag_W","diag_B","diag_H",
		"ndiag_M.HET","ndiag_F.HET","ndiag_Young","ndiag_Old","ndiag_PWID","ndiag_W","ndiag_B","ndiag_H",
		"late_W","late_B","late_H",
		"viral_W","viral_B","viral_H"
	)
	
	all_data_list[[j]] <- Y_j
	
}

Y_all <- do.call(rbind, all_data_list)

dim(Y_all)

#Step  — Estimate covariance and correlation
Sigma_hat <- cov(Y_all) #covariance
R_hat <- cor(scale(Y_all)) #correlation

#Step  — Average within group / within race
idx_diag  <- 1:3
idx_ndiag <- 4:11
idx_late  <- 12:14
idx_viral <- 15:17

idx_W <- c(1,9,12,15)
idx_B <- c(2,10,13,16)
idx_H <- c(3,11,14,17)

mean_offdiag <- function(M) {
	mean(M[upper.tri(M)])
}

# Within outcome correlation
rho_diag  <- mean_offdiag(R_hat[idx_diag, idx_diag])
rho_ndiag <- mean_offdiag(R_hat[idx_ndiag, idx_ndiag])
rho_late  <- mean_offdiag(R_hat[idx_late, idx_late])
rho_viral <- mean_offdiag(R_hat[idx_viral, idx_viral])

rho_target <- mean(c(rho_diag, rho_ndiag, rho_late, rho_viral)) #0.6096505 round to 0.6
# Within race correlation
rho_W <- mean_offdiag(R_hat[idx_W, idx_W])
rho_B <- mean_offdiag(R_hat[idx_B, idx_B])
rho_H <- mean_offdiag(R_hat[idx_H, idx_H])

rho_race <- mean(c(rho_W, rho_B, rho_H)) #0.0718699 round to 0.1
