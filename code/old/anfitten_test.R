library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

# adjust Epochs for 1 patch and 25 patches + 7 epochs and 35 epochs
Nepochs_vec = c(6000, 3000, 6000, 3000, 10000, 6000, 10000, 6000)

# stop("define batchsizes for each model")
batchsize_vec = c(10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L)

paths = c(
  "data/BCI/noSplits/genus-period35-1patch",
  "data/BCI/noSplits/genus-period7-1patch",
  "data/BCI/noSplits/pft-period35-1patch",
  "data/BCI/noSplits/pft-period7-1patch",
  "data/BCI/noSplits/genus-period35-25patches",
  "data/BCI/noSplits/genus-period7-25patches",
  "data/BCI/noSplits/pft-period35-25patches",
  "data/BCI/noSplits/pft-period7-25patches"
)
i = 7
# cl = parallel::makeCluster(4L)
# parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
# parallel::clusterEvalQ(cl, {
#   library(data.table)
#   library(FINN)
#   library(torch)
#   library(glmmTMB)
# })
#
# res = parallel::parLapply(cl, 1:4, function(i) {
#
#   Sys.setenv(CUDA_VISIBLE_DEVICES=i-1)

  p = paths[i]
  Nepochs = Nepochs_vec[i]
  batchsize = batchsize_vec[i]

  env_dt = fread(paste0(p, "/env_dt.csv"))
  obs_dt = fread(paste0(p, "/obs_dt.csv"))
  cohorts_dt =  fread(paste0(p, "/initial_cohorts1985.csv"))

  Nspecies = max(obs_dt$species)
  Nenv = ncol(env_dt) - 2

  ### get init parameters ###
  env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))

  if(i < 3) {

    cat("estimating initial growth parameters")
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
      growth_process = createProcess(~., initEnv = list(init_growth_env), func = FINN::growth, optimizeSpecies = TRUE),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env), func = FINN::regeneration, optimizeSpecies = TRUE),
      mortality_process = createProcess(~., initEnv = list(init_mort_env), func = FINN::mortality, optimizeSpecies = TRUE)
    )

  } else {
    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
      growth_process = createProcess(~., func = FINN::growth, optimizeSpecies = TRUE),
      regeneration_process = createProcess(~., func = FINN::regeneration, optimizeSpecies = TRUE),
      mortality_process = createProcess(~.,  func = FINN::mortality, optimizeSpecies = TRUE)
    )
  }

  cat("\ncreate cohort array")
  cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)

  m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = 25L, lr = 0.01, checkpoints = 5L,
         optimizer = torch::optim_adam, device = "cpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = T,
         loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
         # start_time = 1L
  )

  # torch::torch_save(m1, path = paste0("results/01_full/",strsplit(p, "/")[[1]][2], "_full.pt"))
  # torch::torch_save(m1, path = paste0("results/01_full/01_pft-period7-25patches_full.pt"))
  torch::torch_save(m1, path = paste0("results/01_full/01_pft-period35-25patches_full.pt"))

  rm(m1)
  gc()
  torch::cuda_empty_cache()

})
