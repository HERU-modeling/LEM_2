#' Apply intervention effects to parameters
#'
#' @param input.parameters list of parameters
#' @param current.int types of interventions to be applied
#'
#' @return updated list of parameters after appyling interventions
#' @export
#'
#' @examples
add_intervention_to_param <- function(input.parameters,
																			current.int){
	# Set TRUE/FALSE intervention parameters
	tf.ls <- true.false.interventions.combination(current.int, warn = FALSE)
	
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
		int.pop.eff <- scale.interventions.combination(current.int, tf.ls, warn = FALSE)
		
		
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
	}
return(new.inputs)
}