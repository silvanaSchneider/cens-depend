The codes in this archive can be used to reproduce the results in the paper 

"An approach to model clustered survival data with dependent censoring"
Silvana Schneider, Fábio Nogueira Demarqui, Enrico Antônio Colosimo, Vinícius Diniz Mayrink

Author responsible for the code: Silvana Schneider (silvana.schneider@ufrgs.br)

The code has been written using R version 3.5.3 (2019-03-11).
Packages required:
dlm_1.1-5
survival_2.43-3
rootSolve_1.7
snowfall_1.84-6.1 
snow_0.4-3
MASS_7.3-51.1
dplyr_0.7.8 
finalfit_0.9.5 


#################
## SIMULATIONS ##
#################

For the simulation study present in the paper, we considered 500 datasets for each scenario. 
As this study can take a few days, changing the number of datasets can be easily done in line 28 
(file "Table_1_fit_models.R" and file "Table_2_fit_models.R") changing replica to another value. 

The scripts to perform the simulations are in the "Simulation" folder.

Code to reproduce the simulation results presented in Secation 3 
Table 1: The key functions are contained in the R file "Table_1_fit_models.R", which load the files functions "function_data.R", "function_Weibull.r" and "function_PEM.r".
Table 2: The key functions are contained in the R file "Table_2_fit_models.R", which load the files functions "function_data.R", "function_Weibull.r" and "function_PEM.r".
Figure 1: The codes to generate Figure 1 are in the R file "Figure _1_script_correlations.R".


#################
## APPLICATION ##
#################	   

The dataset form DOPPS is confidential, so we provide artificial dataset simulated that mimics the real data and produces similar results.
The scripts to execute the analysis are in the "Application" folder.

Code to produce similar results as presented in Section 4
Table 3:  The key functions to generate summary statistics for the DOPPS data are in the R file "Table_3.R", which load the file functions "function_Weibull.R".
Table 4:  The key functions for adjusting models are in the R file "Table_4_fit_models.R", which load the files functions "function_Weibull.R" and "functions_PEM.R".


In case of questions, please contact author Silvana Schneider at silvana.schneider@ufrgs.br

