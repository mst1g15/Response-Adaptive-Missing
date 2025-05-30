# Implementing Response-Adaptive Designs when Responses are Missing: Impute or Ignore?
Welcome! This repository contains R scripts to run the simulation study in the article, "Implementing Response-Adaptive Designs when Responses are Missing: Impute or Ignore?"

This simulation study is based on the iCanQuit trial and compares the performance of a number of response-adaptive designs combined with missing data strategies. Please see the manuscript for full details. 

## Code for Simulation Study 

The repository contains the following R scripts: 

- `00_init.R` loads packages, variable names and settings.

- `01_main_simulation_function.R` contains the main wrapper function to run simulation study 

- `01_missingness_functions.R` contains functions to implement missing data strategies

- `01_simulation_functions_support.R` contains functions to support the main simulation function 

- `01_simulation_settings.R` creates .RData files which contain simulation settings. 

- `01_treatment_allocation_functions.R` contains functions to for treatment allocation to support the main simulation function 

- `02_run_simulation_example.R` demonstrates how to run the simulation for one setting 

To run the simulation with the Gittins Index, the table `GIindex_R.csv` is required but is too large to upload. 

## Output 

Results from the simulation study are provided as .RData files in the output folder. 

-  `03_forest_plots_smoking.R`  produces key figures of results 

## How to run the simulation study

You can go straight to `02_run_simulation_example.R` to select a simulation setting and run the simulation. We recommend reducing `nsim` to a low number if running on a PC. 

Please contact Mia.Tackney@mrc-bsu.cam.ac.uk for any queries/comments.
