library(nimble)

machine = 'pine'
vers = 'v2'

# load('data/dump/tree_full.rdata')
load(paste0('data/dump/tree_full_20_', vers, '.rdata'))

# check initial values
D_init = rep(0, N_X)
for (i in 1:N_trees){
  # for (j in 1:last_ti[i]){
  D_init[first_ti[i]:cs_last_ti[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti[i])
  # }
}

D_init[sap2x] < max_size
all(D_init[sap2x] < max_size)

fname_model = 'ring_model_t_date_sapl_size_constraint'

source(paste0('r/models/', fname_model, '.R'))
# full model
m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  N_dbh = N_dbh,
                                  N_trees = N_trees,
                                  N_years = N_years,
                                  N_X = N_X,
                                  N_taxa = N_taxa,
                                  dbh_tree_id = dbh_tree_id,
                                  dbh_year_id = dbh_year_id,
                                  dbh_day_id = dbh_day_id,
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
                                  x2tree = x2tree, 
                                  x2year = x2year,
                                  ones = ones,
                                  meas2d = meas2d,
                                  meas2x = meas2x, 
                                  last_ti = last_ti,
                                  first_ti = first_ti,
                                  cs_last_ti = cs_last_ti,
                                  taxon = taxon,
                                  sap2x = sap2x,
                                  max_size = max_size,
                                  N_saplings = N_saplings,
                                  tree_site_id = tree_site_id
                                  
                 ), check=FALSE)


#debug(m$setData)

# m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s

m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings)))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s

# m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s



inits = list(X = rep(0.5, N_X), 
             D0 = rep(0, N_trees),
             # D = rep(.1, N_X), 
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


if (fname_model == 'ring_model_t_date_taxon') {
  inits = c(inits, list(b0 = 0.5, b1 = 10, beta_spp_sd = 0.1, beta_slope=0.1, beta_spp=rep(0.1, N_taxa)))
}

if ((fname_model == 'ring_model_t_date_sapl_size')|(fname_model == 'ring_model_t_date_sapl_size_constraint')) {
  inits = c(inits, list(beta_slope=0.1))
}


m$setInits(inits)

spec <- configureMCMC(m, thin = 5) 
spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'tau3', 'tau4', 'sig_d_obs'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(10000) #25000
out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

save(out, file = paste0('output/', fname_model, '_nimble_', machine, '_', vers, '.Rdata'))



# #################################################################################################################################
# ## tau2 and tau3 set to zero
# #################################################################################################################################
# 
# load('data/dump/tree_full.rdata')
# source('nimble/model_test.R')
# 
# tau2 = 0
# tau3 = 0
# 
# # full model
# m <- nimbleModel(body(ring_model_full),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   last_ti = last_ti,
#                                   tau2 = tau2,
#                                   tau3 = tau3
#                  ), check=FALSE)
# 
# #debug(m$setData)
# m$setData(list(logDobs = logDobs, logXobs = logXobs))
# 
# 
# m$setInits(list(X = matrix(1, N_trees, N_years), 
#                 D0 = rep(5, N_trees),
#                 D = matrix(.1, N_trees, N_years), 
#                 beta0 = .1, 
#                 beta   = rep(.1, N_trees),
#                 beta_t = rep(.1, N_years), 
#                 sig_x_obs = .1,
#                 sig_d_obs = .1,
#                 beta_sd   = .1, 
#                 beta_t_sd = .1, 
#                 sig_x     = .1))#,
# #                tau2      = .1,
# #                tau3      = .1))
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(1000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# # 
# save(out, file = '~/Documents/projects/npp/output/rw_full_nimble_tau0.Rdata')
# 
# col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
# 
# colnames(out)[which(col_names == 'X')]
# 
# X = matrix(colMeans(out[,which(col_names == 'X')]), nrow=N_trees, ncol=N_years, byrow=FALSE)
# 
# pdf('figures/rw_full_model_nimble_tau0.pdf', width=12, height=10)
# for (i in 1:N_trees){
#   
#   print(paste0('Tree ', i))
#   
#   tree = trees[i]
#   
#   #   tree_idx = which(x2tree == tree)
#   #   tree_years = x2year[tree_idx]
#   
#   # estimated increment
#   mu     = X[tree,]
#   
#   # raw data
#   dat_idx   = which(m2tree == tree)
#   core_nums = unique(m2orient[dat_idx])
#   rws = exp(logXobs[dat_idx])
#   
#   # mu = post[,1,start:(start+N_cores[i])]
#   plot(c(0,0), xlim=c(min(years),max(years)), ylim=c(min(c(rws,mu)), max(c(rws,mu))), type='n', main=paste0('Tree ', tree))
#   
#   lines(years, mu, col='blue')#, lty=2)
#   
#   #   if (nrow(dat) > 1) {
#   #     dat_mean = colMeans(dat)
#   #     lines(seq(1,T), dat_mean, col='red')
#   #   }
#   
#   for (core in core_nums){
#     idx = which((m2tree == tree) & (m2orient == core))
#     yrs = m2t[idx] 
#     lines(yrs, exp(logXobs[idx]), col='black', lty=2)
#   }
#   
# }
# dev.off()
# #################################################################################################################################
# ## average ring widths
# #################################################################################################################################
# 
# load('data/dump/tree_full.rdata')
# 
# cs_last_ti = cumsum(last_ti)
# first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)
# 
# fname_model = 'ring_model_avg'
# source('nimble/model_test.R')
# # full model
# m <- nimbleModel(body(ring_model_avg),
#                  constants = list(N_inc = N_inc_a,
#                                   N_dbh = N_dbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   m2tree = m2tree_a,
#                                   m2ti = m2ti_a,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti
#                                   
#                  ), check=FALSE)
# 
# #debug(m$setData)
# m$setData(list(logDobs = logDobs, logXobs = logXobs_a))
# 
# 
# # m$setInits(list(X = matrix(1, N_trees, N_years), 
# #                 D0 = rep(5, N_trees),
# #                 D = matrix(.1, N_trees, N_years), 
# #                 beta0 = .1, 
# #                 beta   = rep(.1, N_trees),
# #                 beta_t = rep(.1, N_years), 
# #                 sig_x_obs = .1,
# #                 sig_d_obs = .1,
# #                 beta_sd   = .1, 
# #                 beta_t_sd = .1, 
# #                 sig_x     = .1,
# #                 tau2      = .1,
# #                 tau3      = .1,
# #                 tau4      = .1))
# m$setInits(list(X = rep(1, N_X), 
#                 D0 = rep(5, N_trees),
#                 D = rep(.1, N_X), 
#                 beta0 = .1, 
#                 beta   = rep(.1, N_trees),
#                 beta_t = rep(.1, N_years), 
#                 sig_x_obs = .1,
#                 sig_d_obs = .1,
#                 beta_sd   = .1, 
#                 beta_t_sd = .1, 
#                 sig_x     = .1))
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# # 
# save(out, file = paste0('~/Documents/projects/npp/output/', fname_model, '_nimble.Rdata'))
# 
# col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
# colnames(out)[which(col_names == 'X')]
# post = out
# 
# pdf(paste0('figures/', fname_model, '_nimble.pdf'), width=12, height=10)
# for (i in 1:N_trees){
#   
#   print(paste0('Tree ', i))
#   
#   tree = trees[i]
#   tree_idx = which(x2tree == tree)
#   
#   tree_site = tree_site_id[i]
#   
#   # estimated increment
#   X_mu     = colMeans(post[,which(col_names == 'X')[tree_idx]])
#   
#   tree_years = x2year[tree_idx]
#   
#   # raw data
#   dat_idx   = which(m2tree == tree)
#   core_nums = unique(m2orient[dat_idx])
#   rws = exp(logXobs[dat_idx])
#   
#   # mu = post[,1,start:(start+N_cores[i])]
#   par(mfrow=c(2,1))
#   plot(c(0,0), xlim=c(min(tree_years),max(tree_years)), ylim=c(min(c(rws,X_mu)), max(c(rws,X_mu))), xlab='Year', ylab='ring-width',
#        type='n', main=paste0('Tree ', tree, '; site ', tree_site))
#   
#   lines(tree_years, X_mu, col='blue')#, lty=2)
#   
#   #   if (nrow(dat) > 1) {
#   #     dat_mean = colMeans(dat)
#   #     lines(seq(1,T), dat_mean, col='red')
#   #   }
#   
#   for (core in core_nums){
#     idx = which((m2tree == tree) & (m2orient == core))
#     yrs = m2ti[idx] 
#     lines(yrs, exp(logXobs[idx]), col='black', lty=2)
#   }
#   
#   # now plot D!
#   
#   D_mu = colMeans(post[,which(col_names=="D")[tree_idx]])
#   
#   idx_dbh = which(dbh_tree_id == tree)
#   yrs = dbh_year_id[idx_dbh]
#   
#   D_dat = exp(logDobs[idx_dbh])
#   
#   plot(c(0,0), xlim=c(min(tree_years),max(tree_years)), ylim=c(min(c(D_dat,D_mu)), max(c(D_dat,D_mu))), 
#        xlab='Year', ylab='dbh', type='n', main=paste0('Tree ', tree))
#   
#   lines(tree_years, D_mu, col='blue')#, lty=2)
#   points(yrs, D_dat, pch=19, col='black')  
#   
#   
#   
# }
# dev.off()
# 
# sig_x_obs = out[,which(col_names == 'sig_x_obs')]
# mean(sig_x_obs)
# 
# tau2 = out[,which(col_names == 'tau2')]
# tau3 = out[,which(col_names == 'tau3')]
# tau4 = out[,which(col_names == 'tau4')]
# 
# mean(tau2)
# mean(tau3)
# mean(tau4)
# 
# 
# #################################################################################################################################
# ## original base model from paciorek
# #################################################################################################################################
# 
# load('nimble/test/data.Rda')
# source('model.R')
# 
# # full model
# m <- nimbleModel(body(ringModel),
#                  constants = list(n = n, nT = nT, nDBH = nDBH,
#                                   nWidths = nWidths, dbh_tree_id = dbh_tree_id, 
#                                   dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
#                                   incr_year_id = incr_year_id, last_time = last_time ))
# 
# #debug(m$setData)
# m$setData(list(logDobs = logDobs, logXobs = logXobs))
# m$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
#                 D = matrix(.1, n, nT), beta0 = .1, beta = rep(.1, n),
#                 beta_t = rep(.1, nT), sig_x_obs = .1,
#                 sig_d_obs = .1,
#                 beta_sd = .1, beta_t_sd = .1, sig_x = .1))
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(1000) #25000
# out_avg <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# # 
# save(out_avg, file = '~/Documents/projects/npp/output/rw_full_nimble_avg.Rdata')
# 
# source(file = '~/Documents/projects/npp/output/rw_full_nimble_avg.Rdata')
# 
# col_names = sapply(strsplit(colnames(out_avg), '\\['), function(x) x[[1]])
# 
# colnames(out_avg)[which(col_names == 'X')]
# 
# X = matrix(colMeans(out_avg[,which(col_names == 'X')]), nrow=N_trees, ncol=N_years, byrow=FALSE)
# 
# pdf('figures/rw_full_model_nimble_avg.pdf', width=12, height=10)
# for (i in 1:N_trees){
#   
#   print(paste0('Tree ', i))
#   
#   tree = trees[i]
#   
#   #   tree_idx = which(x2tree == tree)
#   #   tree_years = x2year[tree_idx]
#   
#   # estimated increment
#   mu     = X[tree,]
#   
#   # raw data
#   dat_idx   = which(m2tree_a == tree)
#   # core_nums = unique(m2orient[dat_idx])
#   rws = exp(logXobs_a[dat_idx])
#   
#   # mu = post[,1,start:(start+N_cores[i])]
#   plot(c(0,0), xlim=c(min(years),max(years)), ylim=c(min(c(rws,mu)), max(c(rws,mu))), type='n', main=paste0('Tree ', tree))
#   
#   lines(years, mu, col='blue')#, lty=2)
#   
#   #   if (nrow(dat) > 1) {
#   #     dat_mean = colMeans(dat)
#   #     lines(seq(1,T), dat_mean, col='red')
#   #   }
#   
# #   for (core in core_nums){
# #     idx = which((m2tree_a == tree) & (m2orient == core))
# #     yrs = m2t[idx] 
# #     lines(yrs, exp(logXobs[idx]), col='black', lty=2)
# #   }
#   # for (core in core_nums){
#     idx = which((m2tree_a == tree))
#     yrs = m2t_a[idx] 
#     lines(yrs, exp(logXobs_a[idx]), col='black', lty=2)
#   # }
#   
# }
# dev.off()
# 
# 
# 
# #######################################################################################################################################
# # full model
# source('nimble/model_test.R')
# m <- nimbleModel(body(ring_model_full2),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_X = N_X,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   last_ti = last_ti
#                  ), check=FALSE)
# 
# #debug(m$setData)
# m$setData(list(logDobs = logDobs, logXobs = logXobs))
# 
# 
# m$setInits(list(X = matrix(1, N_trees, N_years), 
#                 D0 = rep(5, N_trees),
#                 D = matrix(.1, N_trees, N_years), 
#                 beta0 = .1, 
#                 beta   = rep(.1, N_trees),
#                 beta_t = rep(.1, N_years), 
#                 sig_x_obs = .1,
#                 sig_d_obs = .1,
#                 beta_sd   = .1, 
#                 beta_t_sd = .1, 
#                 sig_x     = .1,
#                 tau2      = .1,
#                 tau3      = .1,
#                 tau4      = .1))
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'tau3', 'tau4', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(2000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# # 
# save(out, file = '~/Documents/projects/npp/output/ring_model_nimble.Rdata')
