#-------------------------------------------------------------------------------
#This script demonstrates how to run the simulation for one setting 
#-------------------------------------------------------------------------------

#load packages and settings for the simulation study 
source("00_init.R")

#load functions 
source("01_simulation_functions_support.R")
source("01_main_simulation_function.R")
source("01_missingness_functions.R")
source("01_treatment_allocation_functions.R")
GIindex <- read.csv("GIindex_R.csv", row.names=1, header=T)


all_settings = readRDS("settings/all_settings_smoke.RDS")

#number of simulation repetitions
N <- 10000

#sample size of trial 
Nt <- 1622

#number of arms
K=2

set.seed(200)



#Covariates are NOT studied in this manuscript, but code was set up to allow for 
#them - created here arbitrarily and have no impact on simulation study 
X <- c(rep(0, Nt/2), rep(1, Nt/2))
X <- sample(X)  #randomize covariates


task_id=1  #choose simulation setting 

setting_i <- all_settings[task_id,]
method_name <- as.character(setting_i$Method)
alpha <- alpha_list[as.character(setting_i$alpha)][[1]]
beta <- beta_list[as.character(setting_i$beta)][[1]]
missing_approach= setting_i$missing_approach
MI_emp <- setting_i$MI_emp
truncated = setting_i$truncated


res_sim <- run_sim(alpha, beta, N, Nt, method_name, missing_approach, setting_i, truncated, burn_in=TRUE)
