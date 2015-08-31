### code to fit the statistical model using NIMBLE

# note - to create an MCMC for a second model within the same R session you may need to give all the related model and MCMC objects new names as overwriting old names may cause problems

source("config")

source("model.R")

#setwd(dataDir)
#setwd('lyford')
load('data/lyford/data.Rda')

tsplot <- function(x) plot(seq_along(x), x, type = 'l')

library(nimble)

# initial model with plot effects
m <- nimbleModel(body(ringModel_plot),
                 constants = list(n = n, 
                                  nT = nT, 
                                  nDBH = nDBH,
                                  nWidths = nWidths,
				  nSites = nSites, 
                                  dbh_tree_id = dbh_tree_id, 
                                  dbh_year_id = dbh_year_id, 
                                  incr_tree_id = incr_tree_id,
                                  incr_year_id = incr_year_id,
				  tree_site_id = tree_site_id, 
                                  last_time    = last_time ))


m$setData(list(logDobs = logDobs, logXobs = logXobs))
m$setInits(list(X  = matrix(1, n, nT), 
                D0 = rep(5, n),
                D  = matrix(.1, n, nT), 
                beta0  = .1, 
                beta   = rep(.1, n),
                beta_t = rep(.1, nT), 
                sig_x_obs = .1,
                sig_d_obs = .1,
                beta_sd   = .1, 
                beta_t_sd = .1, 
                sig_x     = .1))

spec <- configureMCMC(m, thin = 50)
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))

Rmcmc <- buildMCMC(spec)
cm <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)

# some timing from NIMBLE 0.3-1
# 505 sec for Rmodel (8.5 mins)
# 300 for configureMCMC
# 140 for buildMCMC
# 270 to compile model
# 325 to compile MCMC

Cmcmc$run(25000)

out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

save(out, file = 'output/output_ring_plot.Rda')
