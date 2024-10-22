#---------------------------------------------------------
#This script contains functions to implement missing data strategies
#Last updated 22 October 2024
#-----------------------------------------------------------





handle_missing <- function(y_t, arm_t, p0_t, p1_t, approach, ...){
#Purpose: handle missing values either by complete cases, single imputation or Multiple imputation 
#Input: y_t: vector of binary outcomes, possibly with missing values 
#       arm_t: vector of binary treatment assignment 
#       p0_t: scalar for current MLE of success probability for arm 0 
#       p1_t: scalar for current MLE of success probability for arm 1
#       approach: "cc" for complete case, "mean_impute_current" for single imputation of current value, 
#                 "mean_impute" for single imputation of all missing values, "MI" for Multiple Imputation
#Output: y: new vector of binary outcomes, possibly with imputations 
#        p0: updated MLE of success probability for arm 0 
#        p1: updated MLE of success probability for arm 1
#        s_obs0: number of successes for arm 0
#        s_obs1: number of successes for arm 1  
#        f_obs0: number of failures for arm 0 
#        f_obs1: number of failures for arm 1
  
  
    #if no missing values, proceed with complete case analysis. 
    if(sum(is.na(y_t))==0 ){
      approach="cc"
    } 

  
    #single imputation
    if(approach=="single_impute_current" | approach=="single_impute"){
      
      if(is.na(p0_t)){
        p0_t <- 0.5
      }
      
      if(is.na(p1_t)){
        p1_t <- 0.5
      }
        
      #perform single imputation  
      impute_results <- impute_step(y_t, arm_t, p0_t, p1_t)
      
      #update numbers of successes and failures (this includes imputations)
      s_obs0 = impute_results$s_obs0
      f_obs0 = impute_results$f_obs0
      s_obs1 = impute_results$s_obs1
      f_obs1 = impute_results$f_obs1
      
      p0_t = impute_results$p0_t
      p1_t = impute_results$p1_t
      
        
      #if single_impute_current is used: include imputation in y vector.
      #Otherwise, retain the NAs in the y vector as they will be imputed again later. 
        if(approach=="single_impute_current"){
          y_t <- impute_results$y}
      
    }
  
   
    
    if(approach=="MI"){
      
      if(is.na(p0_t)){
        p0_t <- 0.5
      }
      
      if(is.na(p1_t)){
        p1_t <- 0.5
      }
      
      #repeat imputation step 10 times 
      
      impute_results <- replicate(10, impute_step(y_t, arm_t, p0_t, p1_t)) 
      
      #update numbers of successes and failures (this includes imputations)
      s_obs0 = unlist(impute_results[1,])
      f_obs0 = unlist(impute_results[2,])
      s_obs1 = unlist(impute_results[3,])
      f_obs1 = unlist(impute_results[4,])
      
     
      p0_t=mean(unlist(impute_results[5,])) 
      p1_t=mean(unlist(impute_results[6,])) 
          
    }
  
    
      
      if(approach=="impute_zero"){
        
        y_t[is.na(y_t)] <- 0
        
        s_obs0 = sum(y_t[arm_t==0]==1)
        f_obs0 = sum(y_t[arm_t==0]==0)
        s_obs1 = sum(y_t[arm_t==1]==1)
        f_obs1 = sum(y_t[arm_t==1]==0)
        
      }
  
  
    if(approach=="impute_one"){
      
      y_t[is.na(y_t)] <- 1
      
      s_obs0 = sum(y_t[arm_t==0]==1)
      f_obs0 = sum(y_t[arm_t==0]==0)
      s_obs1 = sum(y_t[arm_t==1]==1)
      f_obs1 = sum(y_t[arm_t==1]==0)
      
    }
  

      if(approach=="cc"){
        
      
       #if all observations so far are NA 
        if(sum(is.na(y_t))==length(y_t)){
          
          s_obs0 = f_obs0 = s_obs1 = f_obs1 = 0
          
          
        }else{
          
           
            s_obs0 = sum(y_t[arm_t==0]==1, na.rm=T)
            f_obs0 = sum(y_t[arm_t==0]==0, na.rm=T)
            s_obs1 = sum(y_t[arm_t==1]==1, na.rm=T)
            f_obs1 = sum(y_t[arm_t==1]==0, na.rm=T)
            
            p0_t = s_obs0/(s_obs0+f_obs0)
            p1_t = s_obs1/(s_obs1+f_obs1)
            
        }
            
    }
  
  
return(list(y=y_t, p0=p0_t, p1=p1_t,                      
            s_obs0=s_obs0, s_obs1=s_obs1, 
            f_obs0=f_obs0, f_obs1=f_obs1))

}



impute_step <- function(new_y, arm_t, p0_t, p1_t){
  #impute-step needed if performing multiple imputation (not explored in paper)
  
  num_na_arm1 <- length(new_y[is.na(new_y)& arm_t==1])   #number of missing values arm 1
  new_y[is.na(new_y)& arm_t==1] <- rbinom(num_na_arm1, 1, p1_t) #impute with mean 
  
  num_na_arm0 <- length(new_y[is.na(new_y)& arm_t==0])  #number of missing values arm0
  new_y[is.na(new_y)& arm_t==0] <- rbinom(num_na_arm0, 1, p0_t)  #impute with mean 
  
  
  #update numbers of successes and failures (this includes imputations)
  s_obs0 = sum(new_y[arm_t==0]==1)
  f_obs0 = sum(new_y[arm_t==0]==0)
  s_obs1 = sum(new_y[arm_t==1]==1)
  f_obs1 = sum(new_y[arm_t==1]==0)
  
  p0_new=s_obs0/(s_obs0+f_obs0)
  
  p1_new=s_obs1/(s_obs1+f_obs1)
  
  
  return(list(s_obs0=s_obs0,
              f_obs0=f_obs0,
              s_obs1=s_obs1, 
              f_obs1=f_obs1, 
              p0_t=p0_new, 
              p1_t=p1_new, 
              y = new_y))
}


