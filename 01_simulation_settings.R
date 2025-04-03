#-------------------------------------------------------------------------------
#This script creates .RData files which contain simulation settings. 
#-------------------------------------------------------------------------------

source("00_init.R")


#no-truncation settings---------------------------------------------------------
#with MCAR and MAR
core_methods= c("FR","permuted_block", "neyman", "RTS", "GI", "RGI", "CB", "RBI")
#note that "RTS" (Raw Thompson Sampling) is referred to as BRAR (Bayesian 
#response adaptive randomisation) in the manuscript

alpha_mcar_mar <- c("MCAR", "MAR-T0-low", "MAR-T0-high")
missing_approach_mcar_mar <- c("cc", "single_impute_current", "single_impute")
fully_observed <- expand.grid(Method=core_methods, beta=beta_names, alpha="Fully Observed", missing_approach="cc", truncated=FALSE)
mcar_mar <- expand.grid(Method=core_methods, beta=beta_names, alpha=alpha_mcar_mar, missing_approach=missing_approach_mcar_mar, truncated=FALSE)

#with mnar0
alpha_mnar <- c("MNAR-Y0", "MNAR-T0Y0")
#use all missing names
mnar0 <- expand.grid(Method=core_methods, beta=beta_names, alpha=alpha_mnar, missing_approach=missing_names, truncated=FALSE)

#random component_changed 
rand_comp <- expand.grid(Method=c("RGI_nt", "RBI_nt"), beta=beta_names, alpha=c("MCAR", "MAR-T0-low", "MAR-T0-high", "MNAR-Y0", "MNAR-T0Y0"), 
                         missing_approach="cc", truncated=FALSE)

non_trunc <- rbind(fully_observed, mcar_mar, mnar0, rand_comp)
#truncated settings-------------------------------------------------------------
trunc <- non_trunc %>% filter(Method %in% c("GI", "CB"))
trunc[, "truncated"] <- TRUE
all_settings_smoke <- rbind(non_trunc, trunc)
all_settings_smoke <- all_settings_smoke %>% arrange(Method)
saveRDS(all_settings_smoke, "settings/all_settings_smoke.RDS")   
#-------------------------------------------------------------------------------
