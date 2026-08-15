# LEM_2

Base version of our LEM 2 code.

## list of folders:

- Calibration

  - 1_Run:

  - 2_Workflow:

  - 3_Results:

  - 4_Visualization:

  - 5_Plots:

- Data_Files:\
  contains all the excel and rds files for parameter data inputting

- ModelCoreModules:\
  contains core scripts for the model

  - Parameter_input:\
    all the code related to data inputting

    “Module-Data.input” - main data input code

    “Group.number” - define the subgroups indices among different population groups

    “Par_fix” - load fixed parameters

    “model.initial”- function to initialize population for 42 groups and 19 states

    “Par_pop” – load population related parameters

  - ODE_setup:\
    code for computing transmission matrix and set up the ODE model

    "ode.list" - define the ODEs for non PWID (HET and MSM) groups

    "FoI" - compute the transmission matrix lambda

    "ode_model-Combination" - define the ODE systems needed to solve with or without interventions

  - ODE_eva:\
    code for evaluating the ODE model and generate the main result matrix

    "true.false.interventions.combination" - generate a true/false vector indicating which interventions are applied

    "scale.interventions.combination" - generate parameter change multiplier due to the intervention

    "intervention.model.combination" - function used to evaluate the ODE model with or without the interventions

    "time.period.infections" - get number of HIV infections and susceptibles for different time periods

    "accum.outcomes.combination" - calculate the costs and QALYs

    "comb.eva.func" - call the previous three functions and generate the result matrix

- Load_run:\
