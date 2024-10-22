#---------------------------------------------------------
#This script creates .RData files which contain simulation settings. 
#Last updated 22 October 2024
#-----------------------------------------------------------

source("00_init.R")


#beta parameters for model for the response
beta  <- c("Null", "Null_low",  "Null_high", 
           "Alternative", "Alternative_low", "Alternative_low2", "Alternative_high",
           "Coveffect", "Synergistic", "Antagonistic" )
alpha <- c("Complete", "MCAR", "MAR_X1", "MAR_X0", "MAR_T1_high", "MAR_T1_low", "MAR_T0_high",
           "MNAR_Y0", "MNAR_Y1", "MNAR_T1Y1", "MNAR_X1Y1", "MNAR_T0Y1", "MNAR_X0Y1")
missing_approach <- c("cc", 
                      "single_impute_current", "single_impute", 
                      "MI_current", "MI", 
                      "MI_covariate_current",  "MI_covariate")

#save tables of settings to explore in simulations-----------------------------------------------------------


#SIM1: no covariates. 

beta  <- c("Null", "Null_low",  "Null_high", 
           "Alternative", "Alternative_low", "Alternative_low2", "Alternative_high" )

alpha <- c("MCAR", "MAR_T1_high", "MAR_T1_low", "MAR_T0_high")

missing_approach <- c("cc", "single_impute_current", "single_impute")

nocov_complete <- expand.grid(Method=Method, beta=beta, alpha="Complete", missing_approach="cc", MI_emp=FALSE, controlled=FALSE)
nocov_mar <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)

#with mnar0
alpha <- c("MNAR_Y0", "MNAR_T1Y0", "MNAR_T0Y0")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_zero")
nocov_mnar0 <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)


#with mnar1
alpha <- c("MNAR_Y1", "MNAR_T1Y1", "MNAR_T0Y1")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_one")
nocov_mnar1 <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)


nocov_settings <- rbind(nocov_complete, nocov_mar, nocov_mnar0, nocov_mnar1)


saveRDS(nocov_settings, "settings/nocov_settings.RDS")   

#controlled settings------------------------------------------------------------------------------------------------------------------
cont_method <- c("CB", "RTS", "GI")
beta  <- c("Null", "Null_low",  "Null_high", 
           "Alternative", "Alternative_low", "Alternative_low2", "Alternative_high" )

alpha <- c("MCAR", "MAR_T1_high", "MAR_T1_low","MAR_T0_high" )

missing_approach <- c("cc", "single_impute_current", "single_impute")
controlled_complete <- expand.grid(Method=cont_method, beta=beta, alpha="Complete", missing_approach="cc", MI_emp=FALSE, controlled=TRUE)
controlled_mar <- expand.grid(Method=cont_method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=TRUE)

#with mnar - impute zero
alpha <- c("MNAR_Y0", "MNAR_T1Y0", "MNAR_T0Y0")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_zero")
controlled_mnar0 <- expand.grid(Method=cont_method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=TRUE)

#with mnar - impute one
alpha <- c("MNAR_Y1", "MNAR_T1Y1", "MNAR_T0Y1")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_one")
controlled_mnar1 <- expand.grid(Method=cont_method, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=TRUE)



controlled_settings <- rbind(controlled_complete, controlled_mar, controlled_mnar0, controlled_mnar1)


saveRDS(controlled_settings, "settings/controlled_settings.RDS") 


#nt methods for RBI and RGI, introduce UCB-----------------------------------------------------------------

beta  <- c("Null", "Null_low",  "Null_high", 
           "Alternative", "Alternative_low", "Alternative_low2", "Alternative_high" )

alpha <- c("MCAR", "MAR_T1_high", "MAR_T1_low", "MAR_T0_high")

missing_approach <- c("cc", "single_impute_current", "single_impute")

method_nt <- c("UCB", "UCB_nt", "RBI_nt", "RGI_nt")

nocov_complete <- expand.grid(Method=method_nt, beta=beta, alpha="Complete", missing_approach="cc", MI_emp=FALSE, controlled=FALSE)
nocov_mar <- expand.grid(Method=method_nt, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)

#with mnar0
alpha <- c("MNAR_Y0", "MNAR_T1Y0", "MNAR_T0Y0")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_zero")
nocov_mnar0 <- expand.grid(Method=method_nt, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)


#with mnar1
alpha <- c("MNAR_Y1", "MNAR_T1Y1", "MNAR_T0Y1")
missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_one")
nocov_mnar1 <- expand.grid(Method=method_nt, beta=beta, alpha=alpha, missing_approach=missing_approach, MI_emp=FALSE, controlled=FALSE)


nt_settings <- rbind(nocov_complete, nocov_mar, nocov_mnar0, nocov_mnar1)


saveRDS(nt_settings, "settings/nt_settings.RDS")   

#look at impute-one-----------------------------------------------------------------------------


beta  <- c("Null", "Null_low",  "Null_high", 
           "Alternative", "Alternative_low", "Alternative_low2", "Alternative_high" )

alpha <- c("MNAR_Y1", "MNAR_T1Y1")

missing_approach <- c("cc", "single_impute_current", "single_impute", "impute_one")


#SIM2: covariates
beta  <- c("Coveffect", "Synergistic","Antagonistic" )
alpha <- c("Complete", "MCAR", "MAR_T1_high", "MAR_T1_low", "MAR_X0")
missing_approach <- c("cc", "mean_impute_current", "mean_impute", "MI", "MI_covariate")
cov_mar <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach)

#with mnar 
alpha <- c("MNAR_Y0", "MNAR_T1Y0", "MNAR_X0Y0")
missing_approach <- c("cc", "mean_impute_current", "mean_impute", "MI", "MI_covariate", "impute_zero")
cov_mnar <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach)
cov_settings <- rbind(cov_mar, cov_mnar)
saveRDS(cov_settings, "settings/cov_settings.RDS")   




all_settings <- expand.grid(Method=Method, beta=beta, alpha=alpha, missing_approach=missing_approach)






saveRDS(all_settings, "settings/all_settings.RDS")   


#illustrate probabilities--------------------------------------------------------------------------------------

setting_summary(beta=c(0, 0, 0, 0)) #null case
setting_summary(beta=c(0, 0, 2, 0)) # alternative, no effect of covariate
setting_summary(beta=c(0, 1, 1, 0))  # covariate makes choice for arm 1 harder 
setting_summary(beta=c(0, -1, 1, -2)) # choice of best arm is different 

#null 
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, 0, 0, 0)) #MCAR
setting_summary(alpha=c(-1, 1, 0, 0), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-5, 0, 5, 0), beta=c(0, 0, 0, 0)) #MAR:T1
setting_summary(alpha=c(0, 0, 0, -5), beta=c(0, 0, 0, 0)) #MNAR Y0
setting_summary(alpha=c(-5, 0, 0, 5), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 0, 1, 1), beta=c(0, 0, 0, 0)) #MNAR:T1Y1


#alternative 
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, 0, 2, 0))
setting_summary(alpha=c(-1, 1, 0, 0), beta=c(0, 0, 2, 0))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, 0, 2, 0))
setting_summary(alpha=c(-5, 0, 5, 0), beta=c(0, 0, 2, 0))


#synergistic
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, 1, 1, 0))
setting_summary(alpha=c(-1, 1, 0, 0), beta=c(0, 1, 1, 0))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, 1, 1, 0))


#Antagonistic
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, -1, 1, -2))
setting_summary(alpha=c(-1, 1, 0, 0), beta=c(0, -1, 1, -2))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, -1, 1, -2))


#Coveffect
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, -2, 2, 0))
setting_summary(alpha=c(-1, 1, 0, 0), beta=c(0, -2, 2, 0))
setting_summary(alpha=c(-1, 1, 0, 1), beta=c(0, -2, 2, 0))



setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 0, 0, 0), beta=c(0, 0, 2, 0))


setting_summary(alpha=c(-1, 0, 1, 0), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 0, 1, 0), beta=c(0, 0, 2, 0))


setting_summary(alpha=c(1, 0, 0, -1), beta=c(0, 0, 0, 0))
setting_summary(alpha=c(-1, 0, 0, 1), beta=c(0, 0, 0, 0))


setting_summary(alpha=c(1, 0, 0, -1), beta=c(0, 0, 2, 0))
setting_summary(alpha=c(-1, 0, 0, 1), beta=c(0, 0, 2, 0))



setting_summary(alpha=c(1, 0, 0, -1), beta=c(-1, 0, 2, 0))
setting_summary(alpha=c(-1, 0, 0, 1), beta=c(-1, 0, 2, 0))

setting_summary(alpha=c(0, 0, 0, -5), beta=c(0, 0, 2, 0))
setting_summary(alpha=c(-5, 0, 0, 5), beta=c(0, 0.5, 1.5, 0))



setting_summary(alpha=c(0, -5, 0, 0), beta=c(0, -2, 2, 0))

