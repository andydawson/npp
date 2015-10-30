ring_model_full <- function() {

    # dbh likelihood
    for(j in 1:N_dbh) {
        logDobs[j] ~ dnorm(log(D[dbh_tree_id[j], dbh_year_id[j]]), sd = sig_d_obs)
    }

    # build cov matrices	
    tau2 ~ dunif(-1, 1)
    for (i in 1:2) {
      cov2[i,i] <- sigma_x_obs * sigma_x_obs
    }
    cov2[1,2] <- tau2 * sigma_x_obs * sigma_x_obs
    cov2[2,1] <- cov2[1,2]

    tau3 ~ dunif(-1, 1)
    for (i in 1:4) {	
      cov3[i,i] <- sigma_x_obs * sigma_x_obs
    }
    
    cov3[1,2] <- tau3 * sigma_x_obs * sigma_x_obs
    cov3[1,3] <- tau3 * sigma_x_obs * sigma_x_obs
    cov3[1,4] <- tau3 * sigma_x_obs * sigma_x_obs
    cov3[2,3] <- tau3 * sigma_x_obs * sigma_x_obs
    cov3[2,4] <- tau3 * sigma_x_obs * sigma_x_obs
    cov3[3,4] <- tau3 * sigma_x_obs * sigma_x_obs

    cov3[2,1] <- cov3[1,2]
    cov3[3,1] <- cov3[1,3]
    cov3[4,1] <- cov3[1,4]
    cov3[3,2] <- cov3[2,3]
    cov3[4,2] <- cov3[2,4]
    cov3[4,3] <- cov3[3,4]

    # ring-width likelihood
    
    # one core trees
    for (i in 1:n1cores) {
      logXobs[i1core2m[i]] ~ dnorm(log(X[m2tree[i1core2m[i]], m2ti[i1core2m[i]]]), sigma_x_obs)
    }

    # two core trees
    for (i in 1:n2cores) {
      logXobs[i2core2m[i]:(i2core2m[i]+1)] ~ dmnorm(log(X[m2tree[i2core2m[i]], m2ti[i2core2m[i]]])*ones[1:2], cov2[1:2,1:2])
    }

    # three core trees
    for (i in 1:n3cores) {
      logXobs[i3core2m[i]:(i3core2m[i]+2)] ~ dmnorm(log(X[m2tree[i3core2m[i]], m2ti[i3core2m[i]]]), cov3[1:3,1:3])
    }

    # four core trees (only 2 in lyford)
    for (i in 1:n4cores) {
      logXobs[i4core2m[i]:(i4core2m[i]+2)] ~ dmnorm(log(X[m2tree[i4core2m[i]], m2ti[i4core2m[i]]]), cov3[1:4,1:4])
    }

    # process evolution
    for(i in 1:N_trees) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10
        
        for(j in 1:last_ti[i]) {
            X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
        }
        for(j in (last_ti[i]+1):T) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:T) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta0, sd = beta_sd)
         
    }

    for(j in 1:T) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    
    beta0     ~ dnorm(0, 0.00001)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x     ~ dunif(0, 1000)
    beta_sd   ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
}
	