ring_model_t_date_sapl_size_constraint <- function() {
  
  tau2 ~ dunif(-1, 1)       # 2-core correlation
  tau3 ~ dunif(-1.0/2.0, 1) # 3-core correlation; keep it pos def
  tau4 ~ dunif(-1.0/3.0, 1) # 4-core correlation; keep it pos def
  
  beta0      ~ dnorm(0, sd = 1.0/0.00001)
  beta_slope ~ dunif(0, 0.3)	

  sig_d_obs ~ dunif(1e-6, 1000)
  sig_x_obs ~ dunif(1e-6, 0.6)#1000)
  sig_x     ~ dunif(1e-6, 1000)
  
  beta_sd   ~ dunif(0, 1000)
  beta_t_sd ~ dunif(0, 1000)

  b0 ~ dnorm(0, 1)
  b1 ~ dunif(5, 25)

  # build cov matrices	
  for (i in 1:2) {
    cov2[i,i] <- sig_x_obs * sig_x_obs
  }
  cov2[1,2] <- tau2 * sig_x_obs * sig_x_obs
  cov2[2,1] <- cov2[1,2]
  
  for (i in 1:3) {
    cov3[i,i] <- sig_x_obs * sig_x_obs
  }
  cov3[1,2] <- tau3 * sig_x_obs * sig_x_obs
  cov3[1,3] <- tau3 * sig_x_obs * sig_x_obs
  cov3[2,3] <- tau3 * sig_x_obs * sig_x_obs
  
  cov3[2,1] <- cov3[1,2]
  cov3[3,1] <- cov3[1,3]
  cov3[3,2] <- cov3[2,3]
  
  for (i in 1:4) {
    cov4[i,i] <- sig_x_obs * sig_x_obs
  }
  cov4[1,2] <- tau4 * sig_x_obs * sig_x_obs
  cov4[1,3] <- tau4 * sig_x_obs * sig_x_obs
  cov4[1,4] <- tau4 * sig_x_obs * sig_x_obs
  cov4[2,3] <- tau4 * sig_x_obs * sig_x_obs
  cov4[2,4] <- tau4 * sig_x_obs * sig_x_obs
  cov4[3,4] <- tau4 * sig_x_obs * sig_x_obs
  
  cov4[2,1] <- cov4[1,2]
  cov4[3,1] <- cov4[1,3]
  cov4[4,1] <- cov4[1,4]
  cov4[3,2] <- cov4[2,3]
  cov4[4,2] <- cov4[2,4]
  cov4[4,3] <- cov4[3,4]


  # constraint tau2
  tau2_constraint ~ dconstraint(1/2 * sig_x_obs * sig_x_obs * (1 - tau2) == 0.278)
  tau3_constraint ~ dconstraint(2/3 * sig_x_obs * sig_x_obs * (1 - tau3) == 0.535)

  # ring-width likelihood
  
  # one core trees
  for (i in 1:n1cores) {
    logXobs[i1core2m[i]] ~ dnorm(log(X[meas2x[i1core2m[i]]]), sd=sig_x_obs)
  }
  
  # two core trees
  for (i in 1:n2cores) {
    mn_expr1[1:2,i] <- log(X[meas2x[i2core2m[i]]])*ones[1:2]
    logXobs[i2core2m[i]:(i2core2m[i]+1)] ~ dmnorm(mn_expr1[1:2,i], cov=cov2[1:2,1:2])
  }
  
  # three core trees
  for (i in 1:n3cores) {
    mn_expr2[1:3, i] <- log(X[meas2x[i3core2m[i]]])*ones[1:3]
    logXobs[i3core2m[i]:(i3core2m[i]+2)] ~ dmnorm(mn_expr2[1:3, i], cov=cov3[1:3,1:3])
  }
  
  # four core trees (only 2 in lyford)
  for (i in 1:n4cores) {
    mn_expr3[1:4, i] <- log(X[meas2x[i4core2m[i]]])*ones[1:4]
    logXobs[i4core2m[i]:(i4core2m[i]+3)] ~ dmnorm(mn_expr3[1:4,i], cov=cov4[1:4,1:4])
  }
  # # two core trees
  # for (i in 1:n2cores) {
  #   mn_expr1[1:2] <- log(X[meas2x[i2core2m[i]]])*ones[1:2]
  #   logXobs[i2core2m[i]:(i2core2m[i]+1)] ~ dmnorm(mn_expr1[1:2], cov=cov2[1:2,1:2])
  # }
  # 
  # # three core trees
  # for (i in 1:n3cores) {
  #   mn_expr2[1:3] <- log(X[meas2x[i3core2m[i]]])*ones[1:3]
  #   logXobs[i3core2m[i]:(i3core2m[i]+2)] ~ dmnorm(mn_expr2[1:3], cov=cov3[1:3,1:3])
  # }
  # 
  # # four core trees (only 2 in lyford)
  # for (i in 1:n4cores) {
  #   mn_expr3[1:4] <- log(X[meas2x[i4core2m[i]]])*ones[1:4]
  #   logXobs[i4core2m[i]:(i4core2m[i]+3)] ~ dmnorm(mn_expr3[1:4], cov=cov4[1:4,1:4])
  # }


  # # increment likelihood for trees with rings
  # for(i in 1:N_inc) {
  #     logXobs[i] ~ dnorm(log(X[meas2x[i]]), sd = sig_x_obs)
  # }
  
  # process evolution
  for (i in 1:N_trees){
    D0[i] ~ dunif(-200, 300)
    D[first_ti[i]] <- D0[i] + 2.0 * X[first_ti[i]] / 10.0
    for (t in (first_ti[i]+1):cs_last_ti[i]){
      D[t] <- D[t-1] + 2.0 * X[t] / 10.0
    }
  }
  
  # dbh likelihood
  for(j in 1:N_dbh) {
    logDobs[j] ~ dt(log(D[meas2d[j]-1] + expit(b0 + b1*(dbh_day_id[j]-0.5)) * X[meas2d[j]]/10.0),  1/(sig_d_obs * sig_d_obs), 3)
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

  # for (i in 1:N_X){
  #   X[i] ~ dlnorm(beta[x2tree[i]] + beta_t[x2year[i]], sd = sig_x)
  # }

  for(i in 1:N_trees){
    mnX[first_ti[i]] <- beta[i] + beta_slope * (D0[i] - open_dbh) * (D0[i] > open_dbh) + beta_t[1]
   for (j in (first_ti[i]+1):cs_last_ti[i]){
     mnX[j] <- beta[i] + beta_slope * (D[j-1] - open_dbh) * (D[j-1] > open_dbh) + beta_t[x2year[j]]
   }
  }

 for (i in 1:N_X){
    X[i] ~ dlnorm(mnX[i], sd = sig_x)
  }

  for(i in 1:N_trees) {
    beta[i] ~ dnorm(beta0, sd = beta_sd)
  }
  
  for(j in 1:N_years) {
    beta_t[j] ~ dnorm(0, sd = beta_t_sd)
  }
}