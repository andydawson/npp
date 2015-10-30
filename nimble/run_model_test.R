library(nimble)

source('nimble/model_test.R')

load('data/dump/tree_full.rdata')

# full model
m <- nimbleModel(body(ring_model_full),
                 constants = list(N_inc = N_inc,
                                  N_dbh = N_dbh,
                                  N_trees = N_trees,
                                  T = T,
                                  dbh_tree_id = dbh_tree_id,
                                  dbh_year_id = dbh_year_id,
                                  m2tree = m2tree,
                                  m2ti = m2ti,
                                  m2nc = m2nc,
                                  n1cores = n1cores,
                                  n2cores=n2cores,
                                  n3cores=n3cores,
                                  n4cores=n4cores,
                                  i1core2m = i1core2m,
                                  i2core2m = i2core2m,
                                  i3core2m = i3core2m,
                                  i4core2m = i4core2m,
                                  ones = ones,
                                  meas2d = meas2d,
                                  meas2x = meas2x, 
                                  last_ti = last_ti
                 ), check=FALSE)


m$setData(list(logDobs = logDobs, logXobs = logXobs))

m$setInits(list(X = matrix(1, N_trees, T), 
                D0 = rep(5, N_trees),
                D = matrix(.1, N_trees, T), 
                beta0 = .1, 
                beta   = rep(.1, N_trees),
                beta_t = rep(.1, T), 
                sig_x_obs = .1,
                sig_d_obs = .1,
                beta_sd   = .1, 
                beta_t_sd = .1, 
                sig_x     = .1,
                tau2      = .1,
                tau3      = .1))

# spec <- configureMCMC(m, thin = 50)
# #spec$samplerSpecs[[206]]$control$scale=.01
# spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))
# 
# t4 <- proc.time()
# Rmcmc <- buildMCMC(spec)
# t5 <- proc.time()
# cm <- compileNimble(m)
# t6 <- proc.time()
# Cmcmc <- compileNimble(Rmcmc, project = m)
# t8 <- proc.time()
# Cmcmc$run(10000) #25000
# t9 <- proc.time()
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# # save(out, file = '~/Documents/projects/npp/nimble/output/model_test_output.Rdata')
