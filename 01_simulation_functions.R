#---------------------------------------------------------
#This script contains functions to run the simulation study
#Last updated 22 October 2024
#-----------------------------------------------------------




allocate_missing <- function(x, k, y, alpha){
#Purpose: given a missing data mechanism, values for the treatment, response and covariate, 
#         generate a missing data indicator 
#Input  : x - covariate 
#         k - arm 
#         y - response 
#         alpha - missing data mechanism
#Output:  binary indicator for missingness 
  
  des_vec <- c(1, x, k, y) #one row of design matrix 
  f <- alpha %*% des_vec  #linear predictor
  indic <- rbinom(1, 1, 1/(1+exp(-f)))
  
  return(indic)
}

generate_response <- function(x, k, beta){
  #Purpose: given a data generating mechanism, values for the treatment and covariate, 
  #         generate a response
  #Input  : x - covariate 
  #         k - arm 
  #         beta - data generating mechanism 
  #Output:  binary indicator for missingness 
  
  des_vec <- c(1, x, k, x*k) #one row of design matrix 
  f <- beta %*% des_vec  #linear predictor
  indic <- rbinom(1, 1, 1/(1+exp(-f)))
  
  return(indic)
}

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

compute_prob <- function(x, k, beta){
  
  des_vec <- c(1, x, k, x*k) #one row of design matrix 
  f <- beta %*% des_vec  #linear predictor
  prob= 1/(1+exp(-f))
  
  return(prob)
}

comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY=FALSE)
}

analysis_step <- function(y, arm_t, obs_vec, p0, p1, missing_approach){
  #Purpose: compute p0 and p1 at the end of the trial with complete case
  # and other approaches (single imputation, MI) depending on missingness approach taken 
  #Input: y - total responses, obs_vec: indicator for observed/missing, 
  #       p0 - final value of p0, p1 - final value of p1, missing approach - approach 
  #       taken for treatment assignment, i.e. cc, mean_impute 
  #Output: estimates for p0 and p1 under cc, single imputation and MI 
  #        (NA where approach not applicable)
 
  
  comp_y <- y*obs_vec
  
  s_obs0 = sum(comp_y[arm_t==0]==1, na.rm=TRUE)
  f_obs0 = sum(comp_y[arm_t==0]==0, na.rm=TRUE)
  s_obs1 = sum(comp_y[arm_t==1]==1, na.rm=TRUE)
  f_obs1 = sum(comp_y[arm_t==1]==0, na.rm=TRUE)
  
  cc_p0=s_obs0/(s_obs0+f_obs0)
  cc_p1=s_obs1/(s_obs1+f_obs1)
  
  
  if(missing_approach %in% c("single_impute_current", "single_impute", "impute_zero", "impute_one")){
    
    single_imp_p0=p0
    single_imp_p1=p1
    
  }else{
    single_imp_p0=single_imp_p1=NA
  }
  
  if(missing_approach=="MI"){
    
    MI_p0=p0
    MI_p1=p1
    
  }else{
    
    MI_p0=MI_p1=NA
  }
    
  
  
  return(list(cc_p0=cc_p0, cc_p1=cc_p1, 
              single_imp_p0=single_imp_p0, single_imp_p1=single_imp_p1, 
              MI_p0=MI_p0, MI_p1=MI_p1))
  
  
}




calc.mi.index <- function(sk0, skt, fk0, fkt, Nt, t, K, method_name){
  #purpose: compute index M times given M copies of the complete dataset (with imputations), 
  #         then take the mode of M indices (for CB and GI) or the mean of the M indices (for all others)

  M <- length(skt)
  all_index <- rep(NA, M)
  
  
  for (i in 1:M){
    all_index[i] <- match.fun(method_name)(sk0, skt[i], fk0, fkt[i], Nt, t, K) 
    
  }
  
  if(method_name %in% c("CB", "GI")){
    
    index <- Modes(all_index)
    if (length(index)>1){
      index <- sample(index, 1)
    }

  }else{
    
    
    index <- mean(all_index)
  }

  return(index)
  
}




calc.mi.index.emp <- function(sk0, skt, fk0, fkt, Nt, t, K, method_name){
  #purpose: compute index M times given M copies of the complete dataset (with imputations) for 
  #         GI and CB. Create an empirical distribution of the index (where there are q<=M) possible
  #         values. Select one of q values with weighting defined by the empirical distribution. 

  M <- length(skt)
  all_index <- rep(NA, M)
  
  
  for (i in 1:M){
    all_index[i] <- match.fun(method_name)(sk0, skt[i], fk0, fkt[i], Nt, t, K) 
    
  }
  
  lookup <- data.frame(table(all_index))
  
  lookup$Freq <- lookup$Freq/M
  #make index is numeric and not factor
  lookup$all_index <- as.numeric(levels(lookup$all_index))
  index <- sample(lookup$all_index, 1, prob=lookup$Freq)

  return(index)
  
}



show.mi.index<- function(sk0, skt, fk0, fkt, Nt, t, K, method_name){
  #put description here
  M <- length(skt)
  all_index <- rep(NA, M)
  
  
  for (i in 1:M){
    all_index[i] <- match.fun(method_name)(sk0, skt[i], fk0, fkt[i], Nt, t, K) 
    
  }
  
 return(all_index)  
  
}




run_sim <- function(alpha, beta, N, Nt, method_name, missing_approach, setting_i, 
                    MI_emp=FALSE,  controlled=FALSE, ...){
  
  
  prior0=c(1,1)
  prior1=c(1,1)
  
  #record responses 
  y <- matrix(0, N, Nt)
  
  #indicator for whether observed response or not 
  obs_vec <- rep(0, Nt)
  
  #record indices
  index0_all <- matrix(NA, N, Nt)
  index1_all <- matrix(NA, N, Nt)
  
  #record proportion of assignment to arm 1
  pi <- matrix(NA, N, Nt)
  
  #pstar (unconditional) and pstar in the observed responses
  pstar = pstar_obs = matrix(NA, N, Nt)
  
  #number of observations in each arm 
  N_obs0 = N_obs1 = matrix(NA, N, Nt)
  
  #record estimated probability of success and bias
  p_0 = p_1 = matrix(NA, N, Nt)
  
  #computing final analysis estimates and biases 
  p0_cc_Nt = bias0_cc_Nt = p1_cc_Nt = bias1_cc_Nt= rep(NA, N)
  p0_single_imp_Nt = bias0_single_imp_Nt = p1_single_imp_Nt = bias1_single_imp_Nt =  rep(NA, N)
  p0_MI_Nt = bias0_MI_Nt = p1_MI_Nt = bias1_MI_Nt = rep(NA, N)
  
  #record observed number of successes
  ONS = matrix(NA, N, Nt)
  
  #record assigned arm 
  arm_choice = matrix(NA, N, Nt)
  
  
  #record first mle 
  first_mle_p0=first_mle_p1= rep(NA, N)
  
  #each repetition of simulation 
  for (i in 1:N){
    print(i)
    
    s_obs0=s_obs1 = 0
    f_obs0=f_obs1 = 0
    
    #sequential assignment of Nt participants
    for (t in 1:Nt){
      
      #this is for GI in Python code 
      seed=sample(1:100000, 1)
      
      if(t==1){
        ratio=0.5
      }else{
        ratio=sum(arm_choice[i, 1:(t-1)]==0)/length(arm_choice[i, 1:(t-1)])
      }
      
      
      if(controlled==TRUE &  ratio < 0.10){
          arm=0
        }else if(controlled==TRUE &  ratio > 0.90){
          arm=1
        }else{
      
      
            #compute index------------------------------------------------------------
            sk0=as.integer(prior0[1])
            skt= as.integer(s_obs0)
            fk0=as.integer(prior0[2])
                fkt=as.integer(f_obs0)
            
            #if index has been computed multiple times - summarise 
            if(length(skt)>1 | length(fkt)>1){
              
              if(MI_emp==TRUE){
                index0 <- calc.mi.index.emp(sk0, skt, fk0, fkt, Nt, t, K, method_name)
                
              }else{
                index0 <- calc.mi.index(sk0, skt, fk0, fkt, Nt, t, K, method_name)
              }
              
              
            }else{
              
                num_t=sum(arm_choice[i, 1:(t-1)]==0)
              
                index0 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K, num_t) 
                
            }
    
            index0_all[i, t] =index0
            
            #compute the index for arm1
            sk0=as.integer(prior1[1])
            skt=as.integer(s_obs1)
            fk0=as.integer(prior1[2])
            fkt=as.integer(f_obs1)
    
            if(length(skt)>1 | length(fkt)>1){
              
              if(MI_emp==TRUE){
                index1 <- calc.mi.index.emp(sk0, skt, fk0, fkt, Nt, t, K, method_name)
              }else{
                 index1 <- calc.mi.index(sk0, skt, fk0, fkt, Nt, t, K, method_name)
              }
              
            }else{
              
              num_t=sum(arm_choice[i, 1:(t-1)]==1)
              
              index1 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K, num_t)
              
            }
            
            index1_all[i, t] =index1
            
            pi[i, t]=index1/(index0+index1)
            
            #select best arm, randomize in case of a tie----------------------------
            
            arm = which(c(index0, index1) == max(c(index0, index1)))-1
            if(length(arm)==2){
              arm=sample(c(0,1), 1)
            }
            
        }
            
        #record arm choice
        arm_choice[i, t]=arm
            
        
      
        #generate response------------------------------------------------------
        res = generate_response(x=X[t], k=arm, beta)
      
        #add response 
        y[i,t] = res
      
        #generate missingness--------------------------------------------------------------------
      
        if (length(alpha)==1){
          miss=0
        }else{
          miss=allocate_missing(X[t], arm, res, alpha)
        }
      
        #add indicate missingess in y_obs  
        y[i,t] = ifelse(miss==0, res, NA)
        obs_vec[t] <- ifelse(miss==0, 1, NA)
      
        #total observed number of successes before any imputation happens
        ONS[i, t] <- sum(y[i,1:t]*obs_vec[1:t], na.rm=T)
      
        if(t==1){
          p0_t=NA
          p1_t=NA
        }else{
          p0_t=p_0[i,(t-1)]
          p1_t=p_1[i,(t-1)]
        }
      
        #handle missing data------------------------------------------------------------
      
     
        est <- handle_missing(y_t=y[i,1:t], arm_t=arm_choice[i,1:t], 
                                p0_t, p1_t, approach=missing_approach)
        
      
        #response after missingness approach applied 
        y[i,1:t]  <- est$y
        p_0[i,t] = est$p0 
        p_1[i,t] =  est$p1
        
        
        s_obs0 = as.integer(est$s_obs0) 
        f_obs0= as.integer(est$f_obs0) 
        s_obs1 = as.integer(est$s_obs1) 
        f_obs1 = as.integer(est$f_obs1) 
        
        
        #calculate pstar (including those who have missing responses)
        pstar[i,t] <- sum(arm_choice[i,1:t]==1)/t
        
        
        #conditional pstar: only calculated for those observed 
        pstar_obs[i,t] <- sum(arm_choice[i,1:t]==1 & obs_vec[1:t]==1, na.rm=TRUE) / sum(obs_vec[1:t]==1, na.rm=TRUE)
  
      
      
    }
    
    
    #analysis function here
    final_est <- analysis_step(y[i,], arm_choice[i,], obs_vec, p_0[i, Nt], p_1[i, Nt], missing_approach)
    
    p0_cc_Nt[i] = final_est$cc_p0 
    bias0_cc_Nt[i] = final_est$cc_p0 - compute_prob(0, 0, beta)  #assuming no covariate 
    p1_cc_Nt[i] = final_est$cc_p1
    bias1_cc_Nt[i] = final_est$cc_p1 - compute_prob(0, 1, beta)
    p0_single_imp_Nt[i]=final_est$single_imp_p0
    bias0_single_imp_Nt[i] = final_est$single_imp_p0 - compute_prob(0, 0, beta)  #assuming no covariate 
    p1_single_imp_Nt[i]=final_est$single_imp_p1 
    bias1_single_imp_Nt[i]=final_est$single_imp_p1 - compute_prob(0, 1, beta)
    p0_MI_Nt[i]=final_est$MI_p0
    bias0_MI_Nt[i]=final_est$MI_p0 - compute_prob(0, 0, beta)  #assuming no covariate 
    p1_MI_Nt[i]=final_est$MI_p1
    bias1_MI_Nt[i]=final_est$MI_p1 - compute_prob(0, 1, beta)
    
   
    
  }
  
  
  #summarise results across repetitions 
  #indices 
  index0_mean <- colMeans(index0_all, na.rm=T)
  index0_sd <- colSds(index0_all, na.rm=T)
  index1_mean <- colMeans(index1_all, na.rm=T)
  index1_sd <- colSds(index1_all, na.rm=T)
  
  p0_mean <- colMeans(p_0, na.rm=T)
  p0_sd <- colSds(p_0, na.rm=T)
  p1_mean <- colMeans(p_1, na.rm=T)
  p1_sd <- colSds(p_1, na.rm=T)
  
  #number of patients assigned to arm 1
  N1_mean<- colMeans(t(apply(arm_choice, 1, cumsum)), na.rm=T)
  N1_sd <- colSds(t(apply(arm_choice, 1, cumsum)), na.rm=T)
  
  #pstar
  pstar_mean <- colMeans(pstar, na.rm=T)
  pstar_sd <- colSds(pstar, na.rm=T)
  pstar_obs_mean <-  colMeans(pstar_obs, na.rm=T)
  pstar_obs_sd <- colSds(pstar_obs, na.rm=T)
  
  #observed number of patients
  ONS_mean<- colMeans(ONS, na.rm=T)
  ONS_sd <- colSds(ONS, na.rm=T)
  
  
  #calculate expected bias according to Bowden and Trippa's formula 
  #expbias_p0 <- -cov(p_0[,Nt], N_obs0[,Nt])/mean(N_obs0[,Nt], na.rm=T)
  #expbias_p1 <- -cov(p_1[,Nt], N_obs1[,Nt])/mean(N_obs1[,Nt], na.rm=T)
  
  results_all <- data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
               alpha_setting=as.character(setting_i$alpha),
               missing_approach=missing_approach, MI_emp=MI_emp, controlled= controlled,
               sample_size=1:(Nt),  
               N1_mean=N1_mean, N1_sd=N1_sd,
               pstar_mean=pstar_mean, pstar_sd=pstar_sd,
               pstar_obs_mean = pstar_obs_mean, pstar_obs_sd = pstar_obs_sd,
               p0_mean=p0_mean, p0_sd=p0_sd, 
               p1_mean=p1_mean, p1_sd=p1_sd, 
               ONS_mean=ONS_mean, ONS_sd=ONS_sd, 
               index0_mean=index0_mean, index0_sd=index0_sd, 
               index1_mean=index1_mean, index1_sd=index1_sd)
  
  
  
  #results for final analysis -----------------------------------------------
  
  
  p0_cc_mean = mean(p0_cc_Nt, na.rm=TRUE)
  p0_cc_na = sum(is.na(p0_cc_Nt))   #
  p0_cc_sd = sd(p0_cc_Nt, na.rm=TRUE)
  
  bias0_cc_mean = mean(bias0_cc_Nt, na.rm=TRUE)
  bias0_cc_sd = sd(bias0_cc_Nt, na.rm=TRUE)
  
  p1_cc_mean = mean(p1_cc_Nt, na.rm=TRUE)
  p1_cc_na = sum(is.na(p1_cc_Nt))    #
  p1_cc_sd = sd(p1_cc_Nt, na.rm=TRUE)
  
  bias1_cc_mean = mean(bias1_cc_Nt, na.rm=TRUE)
  bias1_cc_sd = sd(bias1_cc_Nt, na.rm=TRUE)
  
  p0_single_imp_mean=mean(p0_single_imp_Nt, na.rm=TRUE)
  p0_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),sum(is.na(p0_single_imp_Nt)), NA)
  p0_single_imp_sd=sd(p0_single_imp_Nt, na.rm=TRUE)
  
  bias0_single_imp_mean=mean(bias0_single_imp_Nt, na.rm=TRUE)
  bias0_single_imp_sd=sd(bias0_single_imp_Nt, na.rm=TRUE)
  
  p1_single_imp_mean=mean(p1_single_imp_Nt, na.rm=TRUE)
  p1_single_imp_na=ifelse(missing_approach %in% c("single_impute_current", "single_impute"),sum(is.na(p1_single_imp_Nt)), NA)    #
  p1_single_imp_sd=sd(p1_single_imp_Nt, na.rm=TRUE)
  
  bias1_single_imp_mean=mean(bias1_single_imp_Nt, na.rm=TRUE)
  bias1_single_imp_sd=sd(bias1_single_imp_Nt, na.rm=TRUE)
  
  p0_MI_mean=mean(p0_MI_Nt, na.rm=TRUE)
  p0_MI_na=ifelse(missing_approach=="MI", sum(is.na(p0_MI_Nt)), NA)                 #  
  p0_MI_sd=sd(p0_MI_Nt, na.rm=TRUE)
  
  bias0_MI_mean=mean(bias0_MI_Nt, na.rm=TRUE)
  bias0_MI_sd=sd(bias0_MI_Nt, na.rm=TRUE)
  
  p1_MI_mean=mean(p1_MI_Nt, na.rm=TRUE)
  p1_MI_na=ifelse(missing_approach=="MI", sum(is.na(p1_MI_Nt)), NA)                 #  
  p1_MI_sd=sd(p1_MI_Nt, na.rm=TRUE)
  
  bias1_MI_mean=mean(bias1_MI_Nt, na.rm=TRUE)
  bias1_MI_sd=sd(bias1_MI_Nt, na.rm=TRUE)
  
  
  
  results_Nt <- data.frame(data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
                                      alpha_setting=as.character(setting_i$alpha),
                                      missing_approach=missing_approach, MI_emp=MI_emp,controlled=controlled,
                                      pstar_mean=pstar_mean[Nt], pstar_sd=pstar_sd[Nt],
                                      pstar_obs_mean = pstar_obs_mean[Nt], pstar_obs_sd = pstar_obs_sd[Nt],
                                      ONS_mean=ONS_mean[Nt], ONS_sd=ONS_sd[Nt], 
                                      p0_cc_mean=p0_cc_mean, p0_cc_sd=p0_cc_sd, 
                                      p0_cc_na = p0_cc_na,
                                      p1_cc_mean=p1_cc_mean,p1_cc_sd=p1_cc_sd,
                                      p1_cc_na = p1_cc_na,
                                      bias0_cc_mean=bias0_cc_mean, bias0_cc_sd=bias0_cc_sd,
                                      bias1_cc_mean=bias1_cc_mean, bias1_cc_sd=bias1_cc_sd,
                                      p0_single_imp_mean=p0_single_imp_mean, p0_single_imp_sd =p0_single_imp_sd,
                                      p0_single_imp_na = p0_single_imp_na, 
                                      p1_single_imp_mean=p1_single_imp_mean, p1_single_imp_sd=p1_single_imp_sd,
                                      p1_single_imp_na = p1_single_imp_na,
                                      bias0_single_imp_mean=bias0_single_imp_mean, bias0_single_imp_sd=bias0_single_imp_sd,
                                      bias1_single_imp_mean=bias1_single_imp_mean, bias1_single_imp_sd=bias1_single_imp_sd,
                                      p0_MI_mean=p0_MI_mean, p0_MI_sd=p0_MI_sd,
                                      p0_MI_na = p0_MI_na,
                                      p1_MI_mean=p1_MI_mean,p1_MI_sd=p1_MI_sd,
                                      p1_MI_na = p1_MI_na,
                                      bias0_MI_mean=bias0_MI_mean, bias0_MI_sd=bias0_MI_sd,
                                      bias1_MI_mean=bias1_MI_mean, bias1_MI_sd=bias1_MI_sd))
  
  results_pstar <- data.frame(data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
                                         alpha_setting=as.character(setting_i$alpha),
                                         missing_approach=missing_approach, MI_emp=MI_emp, controlled=controlled, pstar=pstar[,300]))
  
  
  return(list(results_all, results_Nt, results_pstar))
  
}


