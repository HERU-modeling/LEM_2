####################################################################################################
## Cascade CEA Model - setup the workspace
## 1. Load packages 
## 2. initialization values (interventions, related time values, discount rate, cities)
## 3. Load necessary functions 
####################################################################################################

## 1. load packeges 

require(openxlsx) # read/write xlsx 
require(abind)
require(deSolve) # for the function "euler" or "ode"
require(tictoc) # for timing

# packages and inputs required only when do parallel modeling
require(foreach)
require(doParallel)
cores <- detectCores()-2
registerDoParallel(cores)



## 2. Initialization values 

# GET INTERVENTIONS 
interventions <- c("Opt-out testing (ER)", "Opt-out testing (PC)",
                   "EMR testing reminder", "Nurse-initiated testing", 
                   "OAT integrated testing",
                   "ART initiation", "ART retention", "ART retention, targeted",
                   "ART EMR reminder", "RAPID ART",
                   "ART re-initiation", "ART re-linkage",
                   "SSP", "OAT with BUP", "OAT with methadone", "PrEP",
									 "No interventions")

# SET time periods
calibration.year.start <- 2017               # calibration start year 2012
calibration.year.end   <- 2022               # NEW: In the new Model it likely will be 2022, but for now 2020 to fix for year LEM 1.0
int.first.year <- 2025                       # year of intervention starts = 2020
int.sus        <- 10                         # intervention sustainment duration = 10 years
lyr            <- 2045                       # end of projection year = 2040
yr             <- calibration.year.start:lyr # year index (in years) from calibration to the end of projection

# SET month index for first year of intervention and the end of projection
start.proj     <- (int.first.year - yr[1]) * 12 # 
end.proj       <- ((lyr - yr[1]) * 12)

# Set time steps
nyr = lyr - calibration.year.start  +1  # no. of years in the model 
nyr_cal = calibration.year.end - calibration.year.start  +1  # NEW
end_yr_ind = c(12*(1:nyr))              # year index (in months)
end_yr_ind_calibration = c(12*(1:nyr_cal))              # NEW
n   = nyr*12                            # number of months between Jan 1 2012 and Jan 1, 2041
vt  = seq(0, n, 1)                      # monthly index, includes t=0


# Set month index for the intervention sustainment (2020-2030) period
int.start <- end_yr_ind[(which(int.first.year == yr, arr.ind = TRUE) - 1)]
int.end <- (int.start + (int.sus * 12)) 

# Set scale-up period length
scale.up.period <- 18

# Discount rate for cost and QALY
Discounting <- 0.03

# City list
all.cities <- c("MAR", "SFA", "LOR", "SAC", "SDI", "DUV", "TAM", "BMP", "ORF", "ATL", "COO", "MIN", "EBR", "ORL", "BAL", "SUF", "WAY", "CLA", "NYC", "NJE", "MEC", "CUY", "FRA", "HAM", "PHI", "SHE", "BEX", "DFW", "HAR", "TRV", "KIN", "DCM")

## 3. load functions 
# functions for setup ODE
source("ModelCoreModules/ODE_setup/Function-ode.list.R")        #ode equations for HET and PWID not MSM
source("ModelCoreModules/ODE_setup/Function-ode.list.MSMyoung.R")        #ode equations for young MSM
source("ModelCoreModules/ODE_setup/Function-ode.list.MSMold.R")        #ode equations for oldMSM
source("ModelCoreModules/ODE_setup/Function-FoI.R") # transition matrix 
source("ModelCoreModules/ODE_setup/Function-ode_model-Combination.R")

# functions for evaluate ODE and generate result matrix
source("ModelCoreModules/ODE_eva/Function-true.false.interventions.combination.R")
source("ModelCoreModules/ODE_eva/Function-scale.interventions.combination.R")
source("ModelCoreModules/ODE_eva/Function-intervention.model.combination.R")
source("ModelCoreModules/ODE_eva/Function-time.period.infections.R")
source("ModelCoreModules/ODE_eva/Function-accumulate.outcomes.combination.R")
source("ModelCoreModules/ODE_eva/Function-comb.eva.func.R")