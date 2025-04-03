#Load packages------------------------------------------------------------------

#library(reticulate)
library(tidyverse)
library(matrixStats)
#library(gridExtra)
library(grid)
#library(ggpubr)
library(mice)
library(foreach)
library(gittins)
library(doParallel)
#library(detectseparation)

#load settings for simulation study --------------------------------------------
Method= c("FR","permuted_block", "neyman", "RTS", "GI", "RGI", "RGI_nt", "CB", "RBI", "RBI_nt")

beta_names <- c("Null", "Null-low", 
                "Alternative", "Alternative-low")

true_success <- data.frame(beta=beta_names, 
                           p0=c(0.2788848, 0.07242649, 0.2141650, 0.07242649), 
                           p1=c(0.2788848, 0.07242649, 0.2809003, 0.11405238))


#beta parameters for model for the response
beta_list <- list("Null" =c(-0.95, 0, 0, 0), #succcess probability of 28% in each arm 
                  "Null-low" =c(-2.55, 0, 0, 0),#succcess probability of 7% in each arm 
                  "Alternative"=c(-1.3, 0, 0.36, 0), # 21% in placebo arm while 28% in the treatment arm 
                  "Alternative-low" = c(-2.55, 0, 0.5, 0)) #success probability 7% in placebo arm and 11% in the treatment arm  


alpha_names <- c("Fully Observed", "MCAR",
                 "MAR-T1-low", "MAR-T0-high",
                 "MNAR-Y0", "MNAR-T0Y0")


#alpha parameters for the model for missingness mechanism
alpha_list <- list("Fully Observed"= "FALSE",                     
                   "MCAR" = c(-1.8, 0, 0, 0),
                   "MAR-T0-low"=c(-1.8, 0, -0.3, 0),
                   "MAR-T0-high"=c(0, 0, -4, 0),
                   "MNAR-Y0" =c(-1.5, 0, 0, -1.5), 
                   "MNAR-T0Y0" = c(-1.5, 0, -0.5, -1.5))

missing_names <- c("cc", "single_impute_current", 
                   "single_impute",  "impute_zero")
