#---------------------------------------------------------
#This script creates runs the simulations on an HPC  
#Last updated 22 October 2024
#-----------------------------------------------------------


source("00_init.R")
source("01_simulation_functions.R")
source("01_missingness_functions.R")
#reticulate::source_python("indices_output.py")
reticulate::source_python("indices.py")


settings_path1 = "settings/nocov_settings.RDS" 
settings_path2 = "settings/controlled_settings.RDS" 

all_settings <- rbind(readRDS(settings_path1), readRDS(settings_path2))  

N <- 10000
Nt <- 300
set.seed(200)
K=2


X <- c(rep(0, Nt/2), rep(1, Nt/2))
#randomize
X <- sample(X)



task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(task_id)


  setting_i <- all_settings[task_id,]
  method_name <- as.character(setting_i$Method)
  alpha <- alpha_list[as.character(setting_i$alpha)][[1]]
  beta <- beta_list[as.character(setting_i$beta)][[1]]
  missing_approach= setting_i$missing_approach
  MI_emp <- setting_i$MI_emp
  controlled = setting_i$controlled
  
  
  res_sim <- run_sim(alpha, beta, N, Nt, method_name, missing_approach, setting_i, MI_emp, controlled)
  



saveRDS(res_sim[[1]], paste0("output/res_seq", task_id, ".RDS"))
saveRDS(res_sim[[2]], paste0("output/res_final", task_id, ".RDS"))
saveRDS(res_sim[[3]], paste0("output/res_pstar", task_id, ".RDS"))

