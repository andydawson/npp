library(nimble)

machine = 'pine'
dvers = 'v5'
mvers = 'v6'

# load('data/dump/tree_full.rdata')
load(paste0('tree_full_20_no_census_', dvers, '.rdata'))
# load(paste0('data/dump/tree_full_20_', dvers, '.rdata'))


# check initial values
D_init = rep(0, N_X)
for (i in 1:N_trees){
  # for (j in 1:last_ti[i]){
  D_init[first_ti[i]:cs_last_ti[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti[i])
  # }
}

# D_init[sap2x] < max_size
# all(D_init[sap2x] < max_size)

fname_model = 'ring_model_t_date_sapl_size_pdbh_nc'

x2idx = match(x2tree, unique(x2tree))

source(paste0('r/models/', fname_model, '.R'))
# full model
m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  #N_dbh = N_dbh,
                                  N_pdbh = N_pdbh,
                                  N_trees = N_trees,
                                  N_years = N_years,
                                  N_X = N_X,
                                  #N_taxa = N_taxa,
                                  #dbh_tree_id = dbh_tree_id,
                                  #dbh_year_id = dbh_year_id,
                                  #dbh_day_id = dbh_day_id,
                                  pdbh_day_id = pdbh_day_id,
                                  open_dbh = open_dbh,
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
                                  #x2tree = x2idx, 
                                  x2year = x2year,
                                  ones = ones,
                                  #meas2d = meas2d,
                                  meas2x = meas2x, 
                                  pdbh2d = pdbh2d,
                                  last_ti = last_ti,
                                  first_ti = first_ti,
                                  cs_last_ti = cs_last_ti,
                                  #taxon = taxon,
                                  #sap2x = sap2x,
                                  max_size = max_size
                                  #N_saplings = N_saplings,
                                  #tree_site_id = tree_site_id
                                  
                 ), check=FALSE)


#debug(m$setData)

# m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s

m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s

# m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s



inits = list(X = rep(0.5, N_X), 
             D0 = rep(0, N_trees),
             D = rep(.1, N_X), 
             beta0 = .1, 
             beta   = rep(.1, N_trees),
             beta_t = rep(.1, N_years), 
             sig_x_obs = .5,
             sig_d_obs = .1,
             beta_sd   = .1, 
             beta_t_sd = .1, 
             sig_x     = .3,
             tau2      = -.1,
             tau3      = -.1,
             tau4      = -.1)
inits = c(inits, b0 = 0.5, b1 = 10)
#inits = c(inits, list(beta_slope=0.1))

m$setInits(inits)

spec <- configureMCMC(m, thin = 5) 
spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'tau3', 'tau4', 'sig_d_obs', 'tmp'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(10000) #25000
out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

save(out, file = paste0('output/', fname_model, '_nimble_', machine, '_', mvers, '.Rdata'))
