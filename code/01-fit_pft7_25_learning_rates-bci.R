library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

# adjust Epochs for 1 patch and 25 patches + 7 epochs and 35 epochs
Nepochs_vec = 10000

batchsize_vec = 10L

lrs = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001)

paths = c(
  "data/BCI/noSplits/pft-period7-25patches"
)

#paths = paths[!grepl("genus", paths)]

cl = parallel::makeCluster(6)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
})

res = parallel::parLapply(cl, 1:length(lrs), function(i) {

  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

  p = paths
  Nepochs = Nepochs_vec
  batchsize = batchsize_vec
  learning_rate = lrs[i]
  Npatches = ifelse(grepl("25patches", p), 25L, 1L)

  env_dt = fread(paste0(p, "/env_dt.csv"))
  obs_dt = fread(paste0(p, "/obs_dt.csv"))
  cohorts_dt =  fread(paste0(p, "/initial_cohorts1985.csv"))

  Nspecies = max(obs_dt$species)
  Nenv = ncol(env_dt) - 2

  ### get init parameters ###
  env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))

    if(grepl("genus",p)){
      env_obs$species_fac <- as.factor(env_obs$species)
      growth_init_model = glmmTMB(log(growth+1)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = gaussian())
      init_growth_env = coef(growth_init_model)[[1]][[1]]
      init_growth_env[is.na(init_growth_env)] = 0
    }else if(grepl("pft", p)){
      growth_init_model = lm(log(growth+1)~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs)
      init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
      init_growth_env[is.na(init_growth_env)] = 0
    }
    ### mort ####
    if(grepl("genus",p)){
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      env_obs$species_fac <- as.factor(env_obs$species)
      mort_init_model = glmmTMB(Fmort~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = beta_family())
      init_mort_env = coef(mort_init_model)[[1]][[1]]
      init_mort_env[is.na(init_mort_env)] = 0
    }else if(grepl("pft", p)){
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      mort_init_model = glmmTMB(Fmort~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
      init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_mort_env[is.na(init_mort_env)] = 0
    }
    ### reg ####
    if(grepl("genus",p)){
      reg_init_model = glmmTMB(ceiling(reg)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = "nbinom1")
      init_reg_env = coef(reg_init_model)[[1]][[1]]
      init_reg_env[is.na(init_reg_env)] = 0
    }else if(grepl("pft", p)){
      reg_init_model = glmmTMB(reg~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
      init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_reg_env[is.na(init_reg_env)] = 0
    }

    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
      growth_process = createProcess(~., initEnv = list(init_growth_env |> as.matrix()), func = FINN::growth, optimizeSpecies = TRUE),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env |> as.matrix()), func = FINN::regeneration, optimizeSpecies = TRUE),
      mortality_process = createProcess(~., initEnv = list(init_mort_env  |> as.matrix()), func = FINN::mortality, optimizeSpecies = TRUE)
    )

  cat("\ncreate cohort array")
  if(grepl("1patch",p)) {
    cohorts_dt_in <- cohorts_dt[,.(siteID, patchID = 1, cohortID, species, dbh = round(dbh,4), trees)]
    cohorts_dt_in_p = lapply(1:10, function(p) {
      patch = copy(cohorts_dt_in)
      patch$patchID = p
      return(patch)
    })
    cohort1 <- FINN::CohortMat(obs_df = data.table::rbindlist(cohorts_dt_in_p), sp = Nspecies)
    Npatches = 10L
  } else {
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)
  }


  m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = learning_rate, checkpoints = 5L,
         optimizer = torch::optim_ignite_adam, device = "gpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = FALSE,
         loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
  )

  torch::torch_save(m1, path = paste0("results/01_learning_rate/",basename(p), "_",learning_rate,".pt"))

  rm(m1)
  gc()
  torch::cuda_empty_cache()

})
