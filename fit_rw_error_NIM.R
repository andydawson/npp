source('rw_error_NIM.R')
load('data/dump/rw_data.rdata')

idx = array(NA, N_trees)
idx[1] <- 1;
for (i in 2:N_trees){
  idx[i] <- idx[i-1] + N_cores[i-1];
}

library(nimble)

m <- nimbleModel(body(rw_error_model),
                 constants = list(N_trees = N_trees, 
                                  T = T, 
                                  N_tot = N_tot,
                                  N_cores = N_cores,
                                  idx = idx))

m$setData(list(logy = logy))
m$setInits(list(mu    = matrix(.1, N_trees, T),
                sigma = .1))

spec <- configureMCMC(m, thin = 50)
spec$addMonitors(c('mu', 'sigma'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(25000)

out <- as.matrix(Cmcmc$mvSamples)
save(out, file = 'output/output_rw_error_NIM.Rda')
