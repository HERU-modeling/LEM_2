############################################################################################
## Cascade CEA Model - Evaluate the model and generate result matrix
## combination evaluation function: generate population list, incidence results, CEA results
## Latest update: Oct 2022
############################################################################################
## Output: incidence and cost-effectiveness outcomes in a vector
############################################################################################

#' Evaluate the ODE with specified interventions
#' 
#' combination evaluation function: generate population list, incidence results, CEA results
#'
#' @param input.parameters Provide in RunModelAndAccum in load_run, often vparameter
#' @param current.int Provide in RunModelAndAccum in load_run, often a vector list of intervention (or no intervention)
#' @param rname Name of the 42 population subgroupsAdded to help with RunModelAndAccum (resolve R name cannot be found)
#'
#' @return a concatenated vector containing unlisted results from the number of infections (t.ni) and the accumulated outcomes (out.ls.accum.int).
#' @export
#'
#' @examples
comb.eva.func <- function(input.parameters,
                          current.int){

  ## Evaluates the model, taking into consideration w/wo interventions
	## three-dimensional array
	## 1-dim: time (not including initial); 2-dim: different states(24); 3-dim: different groups(42)
  out.ls.int <- intervention.model.combination(input.parameters = input.parameters,
                                               current.int = current.int)
  
  ## Gather the results we want (infections, susceptible, yearly infections)
  t.ni <- time.period.infections(simul.in     = out.ls.int,
                                 start.period = start.proj,
                                 end.period   = end.proj,
                                 sust.period  = int.sus)
#browser()
  ## Calculates QALYs, costs, etc
  out.ls.accum.int <- accum.outcomes.combination(input.parameters = input.parameters,
                                                pop.in           = out.ls.int[, 1:16, ],  #Original
  																						#	pop.in           = out.ls.int[-1, 1:16, ],# Compared to original, new out.ls.int include t=0, which we now exclude to make it the same
                                                current.int      = current.int,
                                                sums.only        = TRUE)

  return(c(unlist(t.ni), unlist(out.ls.accum.int)))
}
