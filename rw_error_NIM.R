rw_error_model <- function() {

 for (i in 1:N_trees){   
    for (n in 1:N_cores[i]){
      for (t in 1:T){
        logy[idx[i] + n-1,t] ~ dnorm(mu[i,t], sd = sigma);
      }
    }  
  }

    sigma ~ dunif(0, 1000)
    mu ~ dunif(-10,10)
}


