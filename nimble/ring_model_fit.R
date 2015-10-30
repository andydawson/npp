### code to fit the statistical model using NIMBLE

# note - to create an MCMC for a second model within the same R session you may need to give all the related model and MCMC objects new names as overwriting old names may cause problems

source('config')
source('nimble/rw_full_model.R')

load('data/dump/rw_data.rdata')

#tsplot <- function(x) plot(seq_along(x), x, type = 'l')

library(nimble)

# full model
m <- nimbleModel(body(ring_model_full),
                 constants = list(N_m = N_m,
		                  N_dbh = N_dbh,
		                  N_trees = N_trees,
				  T = T,
				  N_ncores = N_ncores,
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
				  ones     = ones
				  ))

# nT = nT, nDBH = nDBH,
#                                   nWidths = nWidths, dbh_tree_id = dbh_tree_id, 
#                                   dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
#                                   incr_year_id = incr_year_id, last_time = last_time ))

m$setData(list(logDobs = logDobs, logXobs = logXobs))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
                     D = matrix(.1, N_, T), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .1))


m <- nimbleModel(body(ringModel_plot),
                 constants = list(N_trees = N_trees, 
                                  T = T, 
                                  N_tot = N_tot,
                                  N_cores = N_cores))


# initial model with plot effects
m <- nimbleModel(body(ringModel_plot),
                 constants = list(N_trees = N_trees, 
                                  T = T, 
                                  N_tot = N_tot,
                                  N_cores = N_cores))

m$setData(list(logy = logy))
m$setInits(list(mu    = matrix(.1, N_trees, T),
                sigma = .1))

spec <- configureMCMC(m, thin = 50)
#spec$samplerSpecs[[206]]$control$scale=.01
#spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd', 'beta_k', 'beta_k_sd'))

t1 <- proc.time()
Rmcmc <- buildMCMC(spec)
t2 <- proc.time()

cm <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)

# some timing from NIMBLE 0.3-1
# 505 sec for Rmodel (8.5 mins)
# 300 for configureMCMC
# 140 for buildMCMC
# 270 to compile model
# 325 to compile MCMC

Cmcmc$run(10000)#25000)

out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

save(out, file = 'output/output_ring_PE.Rda')
