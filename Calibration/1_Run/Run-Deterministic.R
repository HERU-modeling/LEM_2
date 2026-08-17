
rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
CITY <- Jurisdictions <- "MAR" # MIA1, NYC1 for original, KIN, BAL, PHI -Pick any random cities, just need to get this running
current.int = c("No interventions")
j = 1
	
## LOAD all input parameters and comparators
source("Calibration/2_Workflow/load_targets.R") #Update Target as needed - Part of par_pop.R - Need to move up here since diag18.obs is needed for initial pop
source("ModelCoreModules/Parameter_input/LoadBaselineWorkspace.R")
source("ModelCoreModules/Parameter_input/Module-Data.input.R")

intervention.model.combination <- function(input.parameters,
																					 current.int) {
	# Set TRUE/FALSE intervention parameters
	
	tf.ls <- true.false.interventions.combination(current.int)
	

		new.inputs <- c(input.parameters, tf.ls,
										current.int = list(current.int),
										int.start = 0, int.end = 0
		)


	out_euler <- deSolve::euler(y = x, times = vt, func = ode_model, parms = new.inputs)[, -1] # add parameter names to map
	#  outa <- array(out_euler[-1 , ], dim = c(n, 
	outa <- array(out_euler[   , ], dim = c(n+1, #If want t=0, to compared number of diagnosed PLHIV matched target
																					length(input.parameters$state.name),
																					length(input.parameters$names.gp)
	)) # out_euler[ , ] & n+1, if want to keep initial value 
	
	#  dimnames(outa)[[1]] <- 1:n # 1:(n+1)
	dimnames(outa)[[1]] <- 1:(n+1)
	dimnames(outa)[[2]] <- input.parameters$state.name
	dimnames(outa)[[3]] <- input.parameters$names.gp
	
	return(outa)
}

## 3. evaluate the model and generating the result matrix
outa <- intervention.model.combination(input.parameters = vparameters,
																						 current.int = interventions[16])
n.strata    <- apply(outa[ , c(1:16), ], c(1,3), sum)
# Population Plot
plot_stacked_area(time = vt, data=n.strata , name= colnames(n.strata), CITY="MAR")
# Incidence Projection plot
# Incidence Projection plot
## Add 1 dim to match existing code:
outa <- array(
	outa,
	dim = c(1, dim(outa)),
	dimnames = c(list("1"), dimnames(outa))
)
### Run compute_inc_total,plot_target_2_median
# Incidence index mapping
inc_inds <- c(
	White = gp.gn$white,
	Black = gp.gn$black,
	Hisp  = gp.gn$hisp
)
#browser()	
# Compute totals
inc_totals <- lapply(
	inc_inds,
	compute_inc_total,
	outa        = outa,
	end_yr_ind  = end_yr_ind,
	start_yr_ind = 12
)
# Combine 10 W, 10 B, 10 H strata into 1 W, 1B, 1H
inc_totals <- list(
	black = Reduce(`+`, inc_totals[grep("^Black", names(inc_totals))]),
	hisp  = Reduce(`+`, inc_totals[grep("^Hisp",  names(inc_totals))]),
	white = Reduce(`+`, inc_totals[grep("^White", names(inc_totals))])
) #inc_sum <- Reduce(`+`, inc_totals) -> will give the same value for calculate from sources & races (& risks)
inc_all = inc_totals$white +	inc_totals$black +	inc_totals$hisp

inc_array <- abind::abind(inc_totals, along = 3)

g_inc <- plot_median_quantile(inc_array, dashed_year = 2025, strata_labels = c("Black","Hispanic","White/Other"), y_label = "Number of new infections", legend_title = "Race/Ethnicity",color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	)) 
g_inc

#Population sums 
pop_sum_white   <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$white],	c(1),	sum) # Per Xiao, denominator for inc is All pop
pop_sum_black   <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$black],	c(1),	sum)
pop_sum_hisp    <- apply(	outa[ , end_yr_ind, 1:16, gp.gn$hisp ],	c(1),	sum)
pop_sum_all     <- apply(	outa[ , end_yr_ind, 1:16,  ],	c(1),	sum)


##IRs 
IR_sum_white        <- inc_totals$white        / pop_sum_white      * 100000
IR_sum_black        <- inc_totals$black        / pop_sum_black      * 100000
IR_sum_hisp         <- inc_totals$hisp         / pop_sum_hisp       * 100000
IR_sum_all          <- inc_all                 / pop_sum_all        * 100000

IR_array <- abind::abind(IR_sum_black,IR_sum_hisp,IR_sum_white, along = 3)

g_ir <- plot_median_quantile(IR_array,skip_time = 0 ,dashed_year = 2025 ,strata_labels = c("Black","Hispanic","White/Other"), y_label = "Rate of new HIV incidence per 100,000", legend_title = "Race/Ethnicity",
														 color_values = c("White/Other" = "#F8766D",	"Black" = "#619CFF","Hispanic" = "#00BA38"	))

g_ir
