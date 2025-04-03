# Implementing Response-Adaptive Designs when Responses are Missing: Impute or Ignore?
Welcome! This repository contains R scripts to run the simulation study in the article, "Implementing Response-Adaptive Designs when Responses are Missing: Impute or Ignore?"

This simulation study is based on the iCanQuit trial and compares the performance of a number of response-adaptive designs combined with missing data strategies. Please see the manuscript for full details. 

## Code for Simulation Study 

The repository contains the following R scripts: 

### `00_init.R`

### 01_main_simulation_function.R

### 01_missingness_functions.R contains functions to implement missing data strategies

### 01_simulation_functions_support.R contains functions to support the main simulation function 

### 01_simulation_settings.R creates .RData files which contain simulation settings. 

### 01_treatment_allocation_functions.R contains functions to for treatment allocation to support the main simulation function 

### 02_run_simulation_example.R demonstrates how to run the simulation for one setting 

### 03_forest_plots_smoking.R  produces key figures of results 

In addition, GIindex_R.csv contains Gittins Indices. 

## How to run the simulation study

You can go straight to 02_run_simulation_example.R to select a simulation setting and run the simulation. 

Please contact Mia.Tackney@mrc-bsu.cam.ac.uk for any queries/comments.
