ring_model_t_date_sapl_test <- function() {
  
  beta0     ~ dnorm(0, sd = 1.0/0.00001)

  sig_d_obs ~ dunif(1e-6, 1000)
  sig_x_obs ~ dunif(1e-6, 1000)
  sig_x     ~ dunif(1e-6, 1000)
  
  beta_sd   ~ dunif(0, 1000)
  beta_t_sd ~ dunif(0, 1000)

  b0 ~ dnorm(0, 1)
  b1 ~ dunif(5, 25)

  
  # dbh likelihood
  for(j in 1:N_dbh) {
    logDobs[j] ~ dt(log(D[meas2d[j]-1] + expit(b0 + b1*(dbh_day_id[j]-0.5)) * X[meas2d[j]]/10.0),  1/(sig_d_obs * sig_d_obs), 3)
  }
  
  # ring-width likelihood
  for (i in 1:N_inc) {
    logXobs[i] ~ dnorm(log(X[meas2x[i]]), sd=sig_x_obs)
  }

  # constraint on trees not yet appearing in a census
  for(i in 1:N_saplings) {
    sapling_ind[i] <- D[sap2x[i]] < max_size[i]
    sapling_constraint[i] ~ dconstraint(sapling_ind[i])
  }

  # # constraint on trees not yet appearing in a census
  # for(i in 1:N_saplings) {
  #   not_cens[i] ~ dinterval(D[sap2x[i]], max_size[i])
  # }

  # process evolution
  for (i in 1:N_trees){
    D0[i] ~ dunif(-200, 300)
    D[first_ti[i]] <- D0[i] + 2.0 * X[first_ti[i]] / 10.0
    for (t in (first_ti[i]+1):cs_last_ti[i]){
      D[t] <- D[t-1] + 2.0 * X[t] / 10.0
    }
  }

 for (i in 1:N_X){
    X[i] ~ dlnorm(beta[x2tree[i]] + beta_t[x2year[i]], sd = sig_x)
  }
  
  for(i in 1:N_trees) {
    beta[i] ~ dnorm(beta0, sd = beta_sd)
  }
  
  for(j in 1:N_years) {
    beta_t[j] ~ dnorm(0, sd = beta_t_sd)
  }
}