# LEM_base_model
Base version of our LEM code.

Update (February 23, 2024): Elvis incorporates new calibration routines (ModelCoreModules/Calibration & Notebook); Tri tests model changes (Model-development). 

## list of folders: 
- Data_Files: \
contains all the excel and rds files for parameter data inputting

	For DM: "Evidence-Inputs-Master.xlsx", "Target.xlsx", "cali_par.xlsx"

	For PSA: rds files in the PSA folder

	For Ideal: "Evidence-Inputs-Master-Ideal.xlsx" in the Ideal folder

- ModelCoreModules: \
contains core scripts for the model 

	- Parameter_input: \
		all the code related to data inputting 
		
		“Module-Data.input” - main data input code 
		
		“Group.number” - define the subgroups indices among different population groups
		
		“Par_fix” - load fixed parameters 
		
		“model.initial”- function to initialize population for 42 groups and 19 states
		
		“Par_pop” – load population related parameters
		
		“Par_cali” – load calibrated free parameters
		
		“set.scale.race.gender”, “set.scale.risk.race.gender” - functions related to intervention effect parameters 
		defining
		
		“Par_int” – load intervention related parameters
		
		“Par_CEA” - load CEA related parameters 

	- Combination_list: \
	  code used to generate a complete list of combinations
	  
	  "combination.list" - construct the complete list of all the intervention combinations 
	  
	- ODE_setup: \
		code for computing transmission matrix and set up the ODE model
		
		"ode.list" - define the ODEs for non PWID (HET and MSM) groups
		
		"ode.list.OAT" - define the ODEs for PWID and MSM-PWID receiving OAT
		
		"ode.list.offOAT" - define the ODEs for PWID and MSM-PWID not receiving OAT
		
		"FoI" - compute the transmission matrix lambda 
		
		"ode_model-Combination" - define the ODE systems needed to solve with or without interventions
		
	- ODE_eva: \
		code for evaluating the ODE model and generate the main result matrix
		
		"true.false.interventions.combination" - generate a true/false vector indicating which interventions are applied
		
		"scale.interventions.combination" - generate parameter change multiplier due to the intervention
		
		"intervention.model.combination" - function used to evaluate the ODE model with or without the interventions
		
		"time.period.infections" - get number of HIV infections and susceptibles for different time periods
		
		"accum.outcomes.combination" - calculate the costs and QALYs 
		
		"comb.eva.func" - call the previous three functions and generate the result matrix
	
	- Calibration **(New 2/2024)**: \	code for implementing new calibration routines.


- Load_run: \
	load data and functions, initialization the model, and run the model to get the resulting matrix
	
	"LoadBaselineWorkspace" - load packages, initialization values, load necessary functions
	
	"Combination-Deterministic" - run parallelly and generate results with combination or single interventions for the 
	deterministic model 
	
	"NoIntervention-Determinstic" - run and generate results for the status quo (no interventions)
	
	"Combination-PSA" - run parallelly and generate results with combination or single interventions for PSA model
	
- Output: \
contains rds files for the resulting matrix after running the model 

	Outcome-All-Combination-CITY-DM.rds:Result deterministic matrix of all combination for each city
	
	Outcome-All-Single-CITY-DM.rds: Result deterministic matrix of all single intervention for each city 
	(first 16 interventions in the combination list)
	
	Outcome-No-Intervention-CITY-DM.rds: Result deterministic matrix of Status Quo for each city
	
	Outcome-NUMBERS-CITY-PSA.rds: Result PSA matrix of each combination of interventions for each city (NUMBERS indicate which 
	interventions are considered)

		
- Other_scripts: 
all scripts used to do analysis, not relevant to the main code

- Notebook **(New 2/2024)**: 
Runa and produce new calibration routines result.

- Model_development **(New 2/2024)**: 
Tests new changes to the model, implement one change at a time, all code are
independent from the main model.
