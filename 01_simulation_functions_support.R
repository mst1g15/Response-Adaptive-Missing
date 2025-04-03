#--------------------------------------------------------------------------------
#this script contains functions to support the main simulation function 
#--------------------------------------------------------------------------------

#load common functions
prob_success <- function(beta, x, t){ 
  #purpose: compute probability of success 
  #inputs: beta: data generating mechanism
  #        x: a value for the covariate 
  #        t: value of treatment arm 
  #output: computed probability 
  
  prob <- 1/(1 + exp(-(beta[1]+ beta[2]*x+beta[3]*t+beta[4]*x*t)))   
  
  return(prob)
}


prob_missing<- function(alpha, beta, x, t){
  #purpose: compute probability of missing outcome 
  #inputs: alpha: missing data mechanism
  #        beta: data generating mechanism
  #        x: a value for the covariate 
  #        t: value of treatment arm 
  #output: computed probability 
  
  y <- ifelse(prob_success(beta, x, t) >= 0.5, 1, 0)
  prob <- 1/(1+exp(-(alpha[1]+alpha[2]*x +alpha[3]*t+ alpha[4]*y)))
  return(prob)
  
}

prob_missing_linear <- function(alpha, beta, x, t){
  #purpose: compute probability of successful outcome given covariate and treatment 
  #inputs:  alpha: missing data mechanism
  #         beta: data generating mechanism
  #         x: a value for the covariate 
  #         t: value of treatment arm  
  #output: computed probability 
  
  prob <- 1/(1+exp(-((alpha[1]+alpha[4]*beta[1])+(alpha[2]+alpha[4]*beta[2])*x +
                       (alpha[3]+alpha[4]*beta[3])*t  + alpha[4]*beta[4]*t*x)))
  
  return(prob)
  
}



setting_summary <- function(alpha=c(0, 0, 0, 0), beta){ 
  #create table showing success probabilities for all combinations of t and x 
  
  arm0 <- c(prob_success(beta, 0, 0), prob_success(beta, 1, 0)) 
  arm1 <- c(prob_success(beta, 0, 1), prob_success(beta, 1, 1))
  
  if(is.character(alpha)==TRUE){
    p_arm0=p_arm1=c(0, 0)
  }else{
    
    p_arm0 <- c(prob_missing(alpha,beta, 0, 0), prob_missing(alpha,beta, 1, 0)) 
    p_arm1 <- c(prob_missing(alpha,beta, 0, 1), prob_missing(alpha,beta, 1, 1)) 
    
  }
  
  result <- rbind(c(arm0, mean(arm0)), 
                  c(arm1, mean(arm1)), 
                  c(p_arm0, mean(p_arm0)), 
                  c(p_arm1, mean(p_arm1)))
  colnames(result) <- c("conditional on X=x0", "conditional on X=x1", "marginal")
  rownames(result) <- c("p(success for arm0)",
                        "p(success for arm1)",
                        "p(missing value for arm0)",
                        "p(missing value for arm1)")
  
  return(result)
  
}




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


compute_prob <- function(x, k, beta){
  
  des_vec <- c(1, x, k, x*k) #one row of design matrix 
  f <- beta %*% des_vec  #linear predictor
  prob= 1/(1+exp(-f))
  
  return(prob)
}

comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY=FALSE)
}



analysis_step <- function(y, arm_t, obs_vec, p0, p1, method, missing_approach, allocate_all){
  #Purpose: compute p0 and p1 at the end of the trial with complete case
  # and other approaches (single imputation, MI) depending on missingness approach taken 
  #Input: y - total responses, obs_vec: indicator for observed/missing, 
  #       p0 - final value of p0, p1 - final value of p1, missing approach - approach 
  #       taken for treatment assignment, i.e. cc, mean_impute 
  #Output: estimates for p0 and p1 under cc, single imputation and MI 
  #        (NA where approach not applicable)
  
  
  #results for complete case analysis-------------------------------------------
  comp_y <- y*obs_vec
  
  s_obs0 = sum(comp_y[arm_t==0]==1, na.rm=TRUE)
  f_obs0 = sum(comp_y[arm_t==0]==0, na.rm=TRUE)
  s_obs1 = sum(comp_y[arm_t==1]==1, na.rm=TRUE)
  f_obs1 = sum(comp_y[arm_t==1]==0, na.rm=TRUE)
  
  cc_p0=s_obs0/(s_obs0+f_obs0)
  cc_p1=s_obs1/(s_obs1+f_obs1)
  
  #IPw estimators using complete cases 
  if(method_name %in% c("RTS", "Neyman")){
    HT_est0 <- (1/sum(obs_vec, na.rm=T)) * sum(as.numeric(comp_y==1 & arm_t==0)/allocate_all, na.rm=TRUE)  
    HT_est1 <- (1/sum(obs_vec, na.rm=T)) * sum(as.numeric(comp_y==1 & arm_t==1)/(1-allocate_all) ,na.rm=TRUE)
    
    IPW_est0_cc <- sum(obs_vec, na.rm=T) * HT_est0/ sum (as.numeric(arm_t==0&obs_vec==1)/allocate_all, na.rm=TRUE)
    IPW_est1_cc <- sum(obs_vec, na.rm=T) * HT_est1/ sum (as.numeric(arm_t==1&obs_vec==1)/(1-allocate_all), na.rm=TRUE)
  }else{
    IPW_est0_cc=IPW_est1_cc=NA   #not applicable for other methods 
  }
  
  #results including single imputations-----------------------------------------
  
  if(missing_approach %in% c("single_impute_current", "single_impute", "impute_zero", "impute_one")){
    
    single_imp_p0=p0
    single_imp_p1=p1
    
    #IPw estimators using imputations 
    if(method_name %in% c("RTS", "neyman")){
      HT_est0 <- (1/length(arm_t)) * sum(as.numeric(y==1 & arm_t==0)/allocate_all, na.rm=TRUE)  
      HT_est1 <- (1/length(arm_t)) * sum(as.numeric(y==1 & arm_t==1)/(1-allocate_all), na.rm=TRUE)
      
      IPW_est0_imp <- length(arm_t) * HT_est0/ sum (as.numeric(arm_t==0)/allocate_all)
      IPW_est1_imp <- length(arm_t) * HT_est1/ sum (as.numeric(arm_t==1)/(1-allocate_all))
    }else{
    IPW_est0_imp=IPW_est1_imp=NA
    }
  }else{
    single_imp_p0=single_imp_p1=NA
    IPW_est0_imp=IPW_est1_imp=NA
  }
  
  
  return(list(cc_p0=cc_p0, cc_p1=cc_p1, 
              single_imp_p0=single_imp_p0, single_imp_p1=single_imp_p1,
              IPW_est0_cc=IPW_est0_cc, IPW_est1_cc=IPW_est1_cc,
              IPW_est0_imp=IPW_est0_imp, IPW_est1_imp=IPW_est1_imp))
  
}
