#--------------------------------------------------------------------------------
#this script contains functions to for treatment allocation to 
#support the main simulation function 
#--------------------------------------------------------------------------------

create_permuted_block <- function(block_size) {
  #purpose: create a permuted block 
  #input: block size
  #output: vector for block 
  block <- sample(c(rep(0, block_size / 2), rep(1, block_size / 2))) # Equal numbers of 0s and 1s
  return(block)
}

calc_index <- function(method_name, sk0, skt, fk0, fkt, K, num_t=NULL){
  #purpose: calculate index for deterministic and semi-randomised designs
  #input: method_name: 
  #       sk0: prior number of successes
  #       skt: number of successes
  #       fk0: prior number of failures
  #       fkt: number of failures
  #       K: number of arms
  #       num_t: update to denominator for semi-randomised designs
  #output: calculated index 
  
  if(method_name=="CB"){
    
    index=(sk0 + skt)/(sk0 + fk0 + skt + fkt)
    
  }else if (method_name=="RBI"){
    
    index=(sk0 + skt)/(sk0 + fk0 + skt + fkt)+ rexp(1, 1/K) * K / (sk0 + fk0 + skt + fkt)
    
  }else if (method_name=="RBI_nt"){
    
    index=(sk0 + skt)/(sk0 + fk0 + skt + fkt)+ rexp(1, 1/K) * K / (sk0 + fk0 + num_t)
    
  }else if (method_name=="GI"){
    
    index=GIindex[sk0+skt, fk0+fkt]  
    
    
  }else if (method_name=="RGI"){
    
    index=GIindex[sk0+skt, fk0+fkt] + rexp(1, 1/K) * K / (sk0 + fk0 + skt + fkt) 

  }else if (method_name=="RGI_nt"){
    
    index=GIindex[sk0+skt, fk0+fkt] + rexp(1, 1/K) * K / (sk0 + fk0 + num_t) 

  }else if (method_name=="UCB"){
    
    index = (sk0 + skt) / (sk0 + fk0 + skt + fkt) + sqrt(2 * log(max(1,(skt+fkt)))/(sk0 + fk0 + skt + fkt))

  }
  
  return(index)
}


