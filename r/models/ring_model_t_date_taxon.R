ring_model_t_date_taxon <- function() {
  
  tau2 ~ dunif(-1, 1)       # 2-core correlation
  tau3 ~ dunif(-1.0/2.0, 1) # 3-core correlation; keep it pos def
  tau4 ~ dunif(-1.0/3.0, 1) # 4-core correlation; keep it pos def
  
  beta0     ~ dnorm(0, sd = 1.0/0.00001)
  beta_slope ~ dunif(0, 0.3)

  sig_d_obs ~ dunif(1e-6, 1000)
  sig_x_obs ~ dunif(1e-6, 1000)
  sig_x     ~ dunif(1e-6, 1000)
  
  beta_sd   ~ dunif(0, 1000)
  beta_t_sd ~ dunif(0, 1000)
  beta_spp_sd ~ dunif(0, 1000)

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
  
  # ring-width likelihood
  
  # one core trees
  for (i in 1:n1cores) {
    logXobs[i1core2m[i]] ~ dnorm(log(X[meas2x[i1core2m[i]]]), sd=sig_x_obs)
  }
  
  # two core trees
  for (i in 1:n2cores) {
    logXobs[i2core2m[i]:(i2core2m[i]+1)] ~ dmnorm(log(X[meas2x[i2core2m[i]]])*ones[1:2], cov=cov2[1:2,1:2])
  }
  
  # three core trees
  for (i in 1:n3cores) {
    logXobs[i3core2m[i]:(i3core2m[i]+2)] ~ dmnorm(log(X[meas2x[i3core2m[i]]])*ones[1:3], cov=cov3[1:3,1:3])
  }
  
  # four core trees (only 2 in lyford)
  for (i in 1:n4cores) {
    logXobs[i4core2m[i]:(i4core2m[i]+3)] ~ dmnorm(log(X[meas2x[i4core2m[i]]])*ones[1:4], cov=cov4[1:4,1:4])
  }
  
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
    beta[i] ~ dnorm(beta_spp[taxon[i]], sd = beta_sd)
  }

  for (i in 1:N_taxa){
    beta_spp[i] ~ dnorm(beta0, sd = beta_spp_sd)
  }
  
  for(j in 1:N_years) {
    beta_t[j] ~ dnorm(0, sd = beta_t_sd)
  }
}

# ring_model_full <- function() {
  
#   # dbh likelihood
#   for(j in 1:N_dbh) {
#     logDobs[j] ~ dnorm(log(D[dbh_tree_id[j], dbh_year_id[j]]), sd = sig_d_obs)
#   }
  
#   # build cov matrices	
#   tau2 ~ dunif(-1, 1)
#   for (i in 1:2) {
#     cov2[i,i] <- sig_x_obs * sig_x_obs
#   }
#   cov2[1,2] <- tau2 * sig_x_obs * sig_x_obs
#   cov2[2,1] <- cov2[1,2]
  
#   tau3 ~ dunif(-1/2, 1) # keep it pos def
#   for (i in 1:4) {
#     cov3[i,i] <- sig_x_obs * sig_x_obs
#   }
  
#   cov3[1,2] <- tau3 * sig_x_obs * sig_x_obs
#   cov3[1,3] <- tau3 * sig_x_obs * sig_x_obs
#   cov3[2,3] <- tau3 * sig_x_obs * sig_x_obs
  
#   cov3[2,1] <- cov3[1,2]
#   cov3[3,1] <- cov3[1,3]
#   cov3[3,2] <- cov3[2,3]
  
#   tau4 ~ dunif(-1/3, 1) # keep it pos def
#   for (i in 1:4) {
#     cov4[i,i] <- sig_x_obs * sig_x_obs
#   }
  
#   cov4[1,2] <- tau4 * sig_x_obs * sig_x_obs
#   cov4[1,3] <- tau4 * sig_x_obs * sig_x_obs
#   cov4[1,4] <- tau4 * sig_x_obs * sig_x_obs
#   cov4[2,3] <- tau4 * sig_x_obs * sig_x_obs
#   cov4[2,4] <- tau4 * sig_x_obs * sig_x_obs
#   cov4[3,4] <- tau4 * sig_x_obs * sig_x_obs
  
#   cov4[2,1] <- cov4[1,2]
#   cov4[3,1] <- cov4[1,3]
#   cov4[4,1] <- cov4[1,4]
#   cov4[3,2] <- cov4[2,3]
#   cov4[4,2] <- cov4[2,4]
#   cov4[4,3] <- cov4[3,4]
  
#   # ring-width likelihood
  
#   # one core trees
#   for (i in 1:n1cores) {
#     logXobs[i1core2m[i]] ~ dnorm(log(X[m2tree[i1core2m[i]], m2ti[i1core2m[i]]]), sd=sig_x_obs)
#   }
  
#   # two core trees
#   for (i in 1:n2cores) {
#     logXobs[i2core2m[i]:(i2core2m[i]+1)] ~ dmnorm(log(X[m2tree[i2core2m[i]], m2ti[i2core2m[i]]])*ones[1:2], cov=cov2[1:2,1:2])
#   }
  
#   # three core trees
#   for (i in 1:n3cores) {
#     logXobs[i3core2m[i]:(i3core2m[i]+2)] ~ dmnorm(log(X[m2tree[i3core2m[i]], m2ti[i3core2m[i]]])*ones[1:3], cov=cov3[1:3,1:3])
#   }
  
#   # four core trees (only 2 in lyford)
#   for (i in 1:n4cores) {
#     logXobs[i4core2m[i]:(i4core2m[i]+3)] ~ dmnorm(log(X[m2tree[i4core2m[i]], m2ti[i4core2m[i]]])*ones[1:4], cov=cov4[1:4,1:4])
#   }
  
#   # process evolution
#   for(i in 1:N_trees) {
#     D0[i] ~ dunif(1e-6, 300)
#     D[i, 1] <- D0[i] + X[i, 1]/10
    
#     for(j in 1:last_ti[i]) {
#       X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
#     }
#     for(j in (last_ti[i]+1):N_years) {
#       # BUGS (and NIMBLE by default) should ignore if
#       # indexing is decreasing
#       X[i, j] <- 0
#     }
    
#     for(j in 2:N_years) {
#       D[i, j] <- D[i, j-1] + X[i, j]/10
#     }
#     beta[i] ~ dnorm(beta0, sd = beta_sd)
    
#   }
  
#   for(j in 1:N_years) {
#     beta_t[j] ~ dnorm(0, sd = beta_t_sd)
#   }
  
#   beta0     ~ dnorm(0, 0.00001)
#   sig_d_obs ~ dunif(0, 1000)
#   sig_x_obs ~ dunif(0, 1000)
  
#   sig_x     ~ dunif(0, 1000)
#   beta_sd   ~ dunif(0, 1000)
#   beta_t_sd ~ dunif(0, 1000)
# }

