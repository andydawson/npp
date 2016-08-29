ring_model_avg <- function() {
  
  beta0     ~ dnorm(0, .00001)
  sig_d_obs ~ dunif(0, 1000)
  sig_x_obs ~ dunif(0, 1000)
  
  sig_x     ~ dunif(0, 1000)
  beta_sd   ~ dunif(0, 1000)
  beta_t_sd ~ dunif(0, 1000)
  
  # dbh likelihood
  for(j in 1:N_dbh) {
    logDobs[j] ~ dnorm(log(D[meas2d[j]]), sd = sig_d_obs)
  }
  
  # increment likelihood for trees with rings
  for(i in 1:N_inc) {
    logXobs[i] ~ dnorm(log(X[meas2x[i]]), sd = sig_x_obs)
  }
  
  # process evolution
  for (i in 1:N_trees){
    D0[i] ~ dunif(-200, 300)
    D[first_ti[i]] <- D0[i] + X[first_ti[i]]/10
    for (t in (first_ti[i]+1):cs_last_ti[i]){
      D[t] <- D[t-1] + X[t]/10
    }
  }
  
  for(i in 1:N_X) {
    X[i] ~ dlnorm(beta[x2tree[i]] + beta_t[x2year[i]], sd = sig_x)
  }
  
  for(i in 1:N_trees) {
    beta[i] ~ dnorm(beta0, sd = beta_sd)
  }
  
  for(j in 1:N_years) {
    beta_t[j] ~ dnorm(0, sd = beta_t_sd)
  }
  
}

