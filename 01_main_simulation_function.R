#-------------------------------------------------------------------------------
#This script contains the main wrapper function to run simulation study 
#-------------------------------------------------------------------------------

run_sim <- function(alpha, beta, N, Nt, method_name, missing_approach, setting_i, 
                    truncated=FALSE, burn_in=FALSE, ...){
  #Purpose: run simulation study with N repetitions and output results
  #Input: alpha: vector for missing data mechanism 
  #       beta: vector for data generating mechanism 
  #       N: number of repetitions for the simulation 
  #       Nt: number of participants in a trial 
  #       method_name: design used to allocate treatments 
  #       missing_approach: strategy for handling missing data 
  #       setting_i: vector summarising settings for simulation 
  #       truncated: set to TRUE if proportion of participants allocated to one arm is restricted
  #       burn-in: set to TRUE if the initial four participants are allocated with a 
  #       permuted block design and their responses are observed.
  #Output: list of three dataframes: 
  #       results_throughout: running results from sample size 1 to Nt
  #       results_mle; estimates of success probabilities at Nt
  #       results_pstar: pstar values at Nt 
  
  #set a beta(1,1) prior for both arms
  sk0=fk0=as.integer(1)
  
  #each repetition of simulation 
  results <-  foreach(i=1:N) %do% {
    
    
    #outputs from each repetition of the simulation-----------------------------
    y <- rep(NA, Nt)
    
    #record allocation probabilities for BRAR and Neyman 
    allocate_all=rep(NA, Nt)
    
    #indicator for whether observed response or not 
    #this is needed since y can sometimes include imputations 
    obs_vec <- rep(NA, Nt)
    
    #pstar (unconditional) and pstar in the observed responses
    pstar = pstar_obs = rep(NA, Nt)
    
    #record estimated probability of success and bias
    p_0 = p_1 = rep(NA, Nt)
    
    #record observed number of successes
    ONS = rep(NA, Nt)
    
    #record assigned arm 
    arm_choice = rep(NA,Nt)
    
    #initial values-------------------------------------------------------------
    s_obs0=s_obs1 = as.integer(0)
    f_obs0=f_obs1 = as.integer(0)
    
    
    #setting up treatment allocations for permuted block design-----------------
    #currently hard-coded for the sample size 1622 
    #first block is of size 4 (suitable for the for the burn-in)
    #and following blocks are of sizes 6, 4 and 2 in a random order 
    blocks <- c(4, sample(c(rep(6,100), rep(4, 249), rep(2, 11))))     
    permuted_block_list <- unlist(lapply(blocks, create_permuted_block))
    
    
    #sequential assignment of Nt participants----------------------------------
    for (t in 1:Nt){
      #  print(t)
      #find ratio for purposes of truncation  
      ratio = ifelse(t==1, 0.5, sum(arm_choice[1:(t-1)]==0)/length(arm_choice[1:(t-1)])) 
      
      
      if(burn_in==TRUE & t<=4){
        
        arm=permuted_block_list[t]
        allocation_prop=0.5
        allocate_all[t] <- allocation_prop
        
      }else{
        
        
        if(method_name=="FR"){
          
          arm=sample(c(0,1), 1)
          
        }else if(method_name=="permuted_block"){
          
          arm=permuted_block_list[t]
          
        }else if(method_name=="RTS"){
          
          allocation_prop = sum(rbeta(1000, sk0 + s_obs0, fk0 + f_obs0) > rbeta(1000, sk0 + s_obs1, fk0 + f_obs1))/1000
          
          arm = sample(c(0, 1), 1, prob=c(allocation_prop, 1-allocation_prop))
          
          allocate_all[t] <- allocation_prop
          
        }else if(method_name=="neyman"){
          
          sigma_0 = sqrt(p_0[t-1]*(1-p_0[t-1]))
          sigma_1=sqrt(p_1[t-1]*(1- p_1[t-1]))
          
          if(t==1){
            arm=sample(c(0,1), 1)
            allocate_all[t] <- 0.5
          }else if(sigma_0==0 | sigma_1==0 | is.null(sigma_0) | is.null(sigma_1) |
                   is.na(sigma_0) | is.na(sigma_1)){
            arm=sample(c(0,1), 1)
            allocate_all[t] <- 0.5
          }else{
            
            allocation_prop =sigma_0/(sigma_0+sigma_1)
            arm= sample(c(0, 1), 1, prob=c(allocation_prop, 1-allocation_prop))
            allocate_all[t] <- allocation_prop
            
          }
          
        }else if(truncated==TRUE &  ratio < 0.10){
          arm=0
        }else if(truncated==TRUE &  ratio > 0.90){
          arm=1
        }else{
          
          #all other allocation methods require computing index/prob for each arm
          
          #compute index/prob for arm 0
          num_0=ifelse(is.na(sum(arm_choice[1:(t-1)]==0)), 0, sum(arm_choice[1:(t-1)]==0))
          index0 <- calc_index(method_name, sk0=sk0, skt=s_obs0, fk0=fk0, fkt=f_obs0, K=2, num_t= num_0)
          #compute index/prob for arm1
          num_1=ifelse(is.na(sum(arm_choice[1:(t-1)]==1)), 0, sum(arm_choice[1:(t-1)]==1))
          index1 <-calc_index(method_name, sk0=sk0, skt=s_obs1, fk0=fk0, fkt=f_obs1, K=2, num_t=num_1)
          
          #select best arm, randomize in case of a tie--------------------------
          arm = which(c(index0, index1) == max(c(index0, index1)))-1
          if(length(arm)==2){
            arm=sample(c(0,1), 1)
          }
          
        }
        
      }
      
      #record arm choice
      arm_choice[t]=arm
      
      #generate response--------------------------------------------------------
      res = generate_response(x=X[t], k=arm, beta)
      
      #add response 
      y[t] = res
      
      #generate missingness-----------------------------------------------------
      
      if (length(alpha)==1){  #when there is complete data (no missingness)
        miss=0
      }else if (burn_in==TRUE & t <=4){   #initial design with participant patients - no missing data 
        miss=0
      }else{
        miss=allocate_missing(X[t], arm, res, alpha)
      }
      
      #add indicator missingess in y_obs  
      y[t] = ifelse(miss==0, res, NA)
      obs_vec[t] <- ifelse(miss==0, 1, NA)
      
      #total observed number of successes before any imputation happens
      ONS[t] <- sum(y[1:t]*obs_vec[1:t], na.rm=T)
      
      if(t==1){
        p0_t=NA
        p1_t=NA
      }else{
        p0_t=p_0[t-1]
        p1_t=p_1[t-1]
      }
      
      #handle missing data------------------------------------------------------
      
      est <- handle_missing(y_t=y[1:t], arm_t=arm_choice[1:t], 
                            p0_t, p1_t, approach=missing_approach)
      
      
      #response after missingness approach applied 
      y[1:t]  <- est$y
      p_0[t] = est$p0 
      p_1[t] =  est$p1
      
      
      s_obs0 = as.integer(est$s_obs0) 
      f_obs0= as.integer(est$f_obs0) 
      s_obs1 = as.integer(est$s_obs1) 
      f_obs1 = as.integer(est$f_obs1) 
      
      
      #calculate pstar (including those who have missing responses)
      pstar[t] <- sum(arm_choice[1:t]==1)/t
      
      
      #conditional pstar: only calculated for those observed 
      pstar_obs[t] <- sum(arm_choice[1:t]==1 & obs_vec[1:t]==1, na.rm=TRUE) / sum(obs_vec[1:t]==1, na.rm=TRUE)
      
    }
    
    
    #analysis function here
    final_est <- analysis_step(y, arm_choice, obs_vec, p_0, p_1, method, missing_approach, allocate_all)
    
    
    #return results
    c(list(arm_choice=arm_choice,y=y, obs_vec=obs_vec,
           allocate_all=allocate_all,
           p_0=p_0, p_1=p_1, 
           pstar=pstar, pstar_obs=pstar_obs, 
           ONS=ONS), #these have Nt results per iteration
      final_est)     #this has one result per iteration 
    
    
  }
  
  
  #merge all results 
  all_results <- lapply(purrr::transpose(results), function(l) do.call(rbind, l))
  
  #get summaries----------------------------------------------------------------
  
  p0_mean <- colMeans(all_results$p_0)
  p1_mean <- colMeans(all_results$p_1)
  
  pstar_mean <- colMeans(all_results$pstar, na.rm=T)
  pstar_sd <- colSds(all_results$pstar, na.rm=T)
  
  pstar_obs_mean <-  colMeans(all_results$pstar_obs, na.rm=T)
  pstar_obs_sd <- colSds(all_results$pstar_obs, na.rm=T)
  
  ONS_mean<- colMeans(all_results$ONS, na.rm=T)
  ONS_sd <- colSds(all_results$ONS, na.rm=T)
  
  #output summary results-------------------------------------------------------
  #running results from sample size up to Nt
  results_throughout <- data.frame(method=method_name, 
                                   beta_setting=as.character(setting_i$beta), 
                                   alpha_setting=as.character(setting_i$alpha),
                                   missing_approach=missing_approach, truncated=truncated,
                                   sample_size=1:Nt,  
                                   pstar_mean, pstar_sd,
                                   pstar_obs_mean, pstar_obs_sd = pstar_obs_sd,
                                   p0_mean=p0_mean, p1_mean=p1_mean,
                                   ONS_mean=ONS_mean, ONS_sd=ONS_sd)
  
  #results for final analyses, MLEs---------------------------------------------
  
  #get summaries of MLEs
  #estimates of with cc---------------------------------------------------------
  p0_cc_mean = mean(all_results$cc_p0, na.rm=TRUE)
  p0_cc_na = sum(is.na(all_results$cc_p0))   
  p0_cc_sd = sd(all_results$cc_p0, na.rm=TRUE)
  
  IPW_p0_cc_mean = mean(all_results$IPW_est0_cc, na.rm=TRUE)
  IPW_p0_cc_na = sum(is.na(all_results$IPW_est0_cc))   
  IPW_p0_cc_sd = sd(all_results$IPW_est0_cc, na.rm=TRUE)
  
  p1_cc_mean = mean(all_results$cc_p1, na.rm=TRUE)
  p1_cc_na = sum(is.na(all_results$cc_p1))   
  p1_cc_sd = sd(all_results$cc_p1, na.rm=TRUE)
  
  IPW_p1_cc_mean = mean(all_results$IPW_est1_cc, na.rm=TRUE)
  IPW_p1_cc_na = sum(is.na(all_results$IPW_est1_cc))   
  IPW_p1_cc_sd = sd(all_results$IPW_est1_cc, na.rm=TRUE)
  
  #estimates with single imputation---------------------------------------------
  p0_single_imp_mean=mean(all_results$single_imp_p0, na.rm=TRUE)
  p0_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),
                          sum(is.na(all_results$single_imp_p0)), NA)
  p0_single_imp_sd=sd(all_results$single_imp_p0, na.rm=TRUE)
  
  IPW_p0_single_imp_mean=mean(all_results$IPW_est0_imp, na.rm=TRUE)
  IPW_p0_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),
                              sum(is.na(all_results$IPW_est0_imp)), NA)
  IPW_p0_single_imp_sd=sd(all_results$IPW_est0_imp, na.rm=TRUE)
  
  #estimates of p1
  p1_single_imp_mean=mean(all_results$single_imp_p1, na.rm=TRUE)
  p1_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),
                          sum(is.na(all_results$single_imp_p1)), NA)    
  p1_single_imp_sd=sd(all_results$single_imp_p1, na.rm=TRUE)
  
  IPW_p1_single_imp_mean=mean(all_results$IPW_est1_imp, na.rm=TRUE)
  IPW_p1_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),
                              sum(is.na(all_results$IPW_est1_imp)), NA)
  IPW_p1_single_imp_sd=sd(all_results$IPW_est1_imp, na.rm=TRUE)
  
  
  
  results_mle <- data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
                            alpha_setting=as.character(setting_i$alpha),
                            missing_approach=missing_approach, truncated=truncated,
                            pstar_mean=pstar_mean[Nt], pstar_sd=pstar_sd[Nt],
                            pstar_obs_mean = pstar_obs_mean[Nt], pstar_obs_sd = pstar_obs_sd[Nt],
                            ONS_mean = ONS_mean[Nt], ONS_sd= ONS_sd[Nt],
                            p0_cc_mean, p0_cc_sd, p0_cc_na,
                            IPW_p0_cc_mean, IPW_p0_cc_sd, IPW_p0_cc_na, 
                            p1_cc_mean, p1_cc_sd, p1_cc_na,
                            IPW_p1_cc_mean, IPW_p1_cc_sd, IPW_p1_cc_na, 
                            p0_single_imp_mean, p0_single_imp_sd, p0_single_imp_na,
                            IPW_p0_single_imp_mean, IPW_p0_single_imp_sd, IPW_p0_single_imp_na, 
                            p1_single_imp_mean, p1_single_imp_sd, p1_single_imp_na, 
                            IPW_p1_single_imp_mean, IPW_p1_single_imp_sd, IPW_p1_single_imp_na)
  
  results_pstar <- data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
                              alpha_setting=as.character(setting_i$alpha),
                              missing_approach=missing_approach,  truncated=truncated, pstar=all_results$pstar[,Nt])
  
  
  return(list(results_throughout, results_mle, results_pstar))
  
}


