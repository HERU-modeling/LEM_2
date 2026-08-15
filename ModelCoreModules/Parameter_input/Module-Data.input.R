############################################################################################
## DATA INPUT
## 1. directly load in fixed parameters (1-dim and multi-dim) and calibrated free parameters
## 2. CAE parameters
## Update: Oct 2022
###########################################################################################

## Read in data from excel file ##
myFile.name <- c("Data_files/DM/Evidence-Inputs-Master")


WB <- openxlsx::loadWorkbook(paste0(myFile.name, ".xlsx"))

par_info <- read.xlsx(WB, sheet = "parameter_info") # first tab "parameter_info" for parameter info, e.g. name, dimension
vlist <- list() # initialize intermediate vlist (as a list) to save inputs
length(vlist) <- nrow(par_info) # number of input parameters
names(vlist) <- par_info$parameter # name of input parameters
vparameters <- vlist # initialize vparameters - parameter list used in model


# read in fixed parameters
source("ModelCoreModules/Parameter_input/Par_fix.R")
# read in initial population
source("ModelCoreModules/Parameter_input/Par_pop.R")
#browser()
# read in calibrated free parameters
source("ModelCoreModules/Parameter_input/Par_cali.R")


vparameters <- vparameters[-which(sapply(vparameters, is.null))]
