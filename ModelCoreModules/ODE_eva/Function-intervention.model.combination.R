#########################################################################################
## Cascade CEA Model - Solve ODE Model w or wo Interventions
## Function for running the simulation model w/ or w/o interventions
## Latest update: Oct 2022
########################################################################################
## Output: three-dimensional array
## 1-dim: time ; 2-dim: different states(24= 19 + 5[inc_bo,inc_bs,inc_g,diagm,deathm]);
## 3-dim: different groups(42)
########################################################################################

#' Cascade CEA Model - Solve ODE Model w or wo Interventions
#'
#' @param input.parameters Provide in RunModelAndAccum in load_run, often vparameter
#' @param current.int Provide in RunModelAndAccum in load_run, often a vector list of intervention (or no intervention)
#'
#' @return outa Three-dimensional array:1-dim: time ; 2-dim: different states(24= 19 + 5[inc_bo,inc_bs,inc_g,diagm,deathm]);3-dim: different groups(42)
#' @export
#'
#' @examples
intervention.model.combination <- function(input.parameters,
                                           current.int) {
  # Set TRUE/FALSE intervention parameters
  tf.ls <- true.false.interventions.combination(current.int)
#browser()
  # Run model
  if (all(current.int %in% c(
    "Opt-out testing (ER)", "Opt-out testing (PC)",
    "EMR testing reminder", "Nurse-initiated testing",
    "OAT integrated testing",
    "ART initiation",
    "ART retention", "ART retention, targeted",
    "ART EMR reminder", "RAPID ART",
    "ART re-initiation", "ART re-linkage",
    "SSP", "OAT with BUP", "OAT with methadone",
    "PrEP"
  ))) {
    # Set intervention parameters for scale*effect
    int.pop.eff <- scale.interventions.combination(current.int, tf.ls)


    new.inputs <- c(input.parameters, tf.ls,
      current.int = list(current.int),
      int.start = int.start, int.end = int.end,
      int.pop.eff
    )
  } else if (all(current.int == c("No interventions"))) {
    new.inputs <- c(input.parameters, tf.ls,
      current.int = list(current.int),
      int.start = 0, int.end = 0
    )
    warning("There are no interventions being modeled")
  }

  # Call model function and solve with Euler's method
  # Euler function part of deSolve (added prefix for R package)
  # y - the initial (state) values for the DE system, a vector
  # times - time sequence for which output is wanted; the first value of times must be the initial time, a vector
  # func - function that computes the values of the derivatives in the ODE system (the model definition) at time t
  # ode_model defined in "CascadeCEA-Model-0-Function-ode_model-Combination.R"
  # parms - parameters passed to func
  out_euler <- deSolve::euler(y = x, times = vt, func = ode_model, parms = new.inputs)[, -1] # add parameter names to map
#  outa <- array(out_euler[-1 , ], dim = c(n, 
  outa <- array(out_euler[   , ], dim = c(n+1, #If want t=0, to compared number of diagnosed PLHIV matched target
  	                                      length(input.parameters$state.name),
    																			length(input.parameters$names.gp)
  																				)) # out_euler[ , ] & n+1, if want to keep initial value 
#browser()  
#  dimnames(outa)[[1]] <- 1:n # 1:(n+1)
  dimnames(outa)[[1]] <- 1:(n+1)
  dimnames(outa)[[2]] <- input.parameters$state.name
  dimnames(outa)[[3]] <- input.parameters$names.gp
#browser()  
  return(outa)
}
