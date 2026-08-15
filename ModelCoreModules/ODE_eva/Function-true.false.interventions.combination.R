#############################################################################
## Cascade CEA Model - Interventions True/False
## Last updated:Oct 2022
#############################################################################
## Output: Value of logical vectors for intervention modeling
## The vector indicates which interventions are applied
#############################################################################
# interventions <- c("Opt-out testing (ER)", "Opt-out testing (PC)",
#                    "EMR testing reminder", "Nurse-initiated testing", "OAT integrated testing",
#                    "ART initiation", "ART retention", "ART retention, targeted",
#                    "ART EMR reminder", "RAPID ART",
#                    "ART re-initiation", "ART re-linkage",
#                    "SSP", "OAT with BUP", "OAT with methadone", "PrEP",
#                    "No interventions")
#
# Interventions, suffix
# "Opt-out testing (ER)", test.optout.ER
# "Opt-out testing (PC)", test.optout.PC
# "EMR testing reminder", test.EMR
# "Nurse-initiated testing", test.nurse
# "OAT integrated testing", test.OAT
# 
# "ART initiation", ARTinit
# "ART retention", ARTreten
# "ART retention, targeted", ARTreten.target
# "ART EMR reminder", ARTEMRprompt
# "RAPID ART", immART
# 
# "ART re-initiation", reART
# "ART re-linkage", relink
# 
# "SSP", SSP
# "OAT with BUP", BUP
# "OAT with methadone", MET
# "PrEP", prep
#
#############################################################################

#' Interventions True/False
#'
#' @param current.int a vector of all the interventions investigated
#'
#' @return a list of Value of logical vectors for intervention modeling
#' @export
#'
#' @examples
true.false.interventions.combination <- function(current.int, warn = TRUE){
  if(any(current.int %in% c("Opt-out testing (ER)", "Opt-out testing (PC)", 
                        "EMR testing reminder", "Nurse-initiated testing", 
                        "OAT integrated testing"))){
          intervention.testing.modif <- TRUE
  } else {intervention.testing.modif <- FALSE} 
  
  if(any(current.int==c("ART initiation"))) {
          intervention.initiation.modif <- TRUE
  } else {intervention.initiation.modif <- FALSE}
  
  if(any(current.int %in% c("ART retention", "ART retention, targeted"))) {
          intervention.retention.modif <- TRUE
  } else {intervention.retention.modif <- FALSE}
  
  if(any(current.int==c("ART EMR reminder"))) {
          intervention.dropout.modif <- TRUE
  } else {intervention.dropout.modif <- FALSE}
  
  if(any(current.int %in% c("RAPID ART"))) {
          intervention.immediateART.modif <- TRUE
  } else {intervention.immediateART.modif <- FALSE}
  
  if(any(current.int %in% c("ART re-initiation", "ART re-linkage"))) {
          intervention.reinitiation.modif <- TRUE
  } else {intervention.reinitiation.modif <- FALSE}
  
  if(any(current.int==c("SSP"))) {
          intervention.ssp.modif <- TRUE
  } else {intervention.ssp.modif <- FALSE}
  
  if(any(current.int %in% c("OAT with BUP", "OAT with methadone"))) {
          intervention.oat.modif <- TRUE
  } else {intervention.oat.modif <- FALSE}
  
  if(any(current.int==c("PrEP"))) {
          intervention.prep.modif <- TRUE
  } else {intervention.prep.modif <- FALSE}
  
  if(any(current.int==c("No interventions")) & warn) {
     warning("There are no interventions being modeled")
  }
 
  tf.ls <- list(
    intervention.testing.modif      = intervention.testing.modif,
    intervention.initiation.modif   = intervention.initiation.modif,
    intervention.retention.modif    = intervention.retention.modif,
    intervention.dropout.modif      = intervention.dropout.modif,
    intervention.immediateART.modif = intervention.immediateART.modif,
    intervention.reinitiation.modif = intervention.reinitiation.modif,
    intervention.ssp.modif          = intervention.ssp.modif,
    intervention.oat.modif          = intervention.oat.modif,
    intervention.prep.modif         = intervention.prep.modif
  )
  return(tf.ls)
}