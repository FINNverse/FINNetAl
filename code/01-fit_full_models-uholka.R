library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

# adjust Epochs for 1 patch and 25 patches + 7 epochs and 35 epochs
Nepochs_vec = c(1000L, 2000L, 1000L, 2000L)

batchsize_vec = c(10L, 90L, 10L, 90L)
paths = c(
  "data/Uholka/noSplits/species-period3-9patches/",
  "data/Uholka/noSplits/species-period3-1patch/",
  "data/Uholka/noSplits/species-period15-9patches/",
  "data/Uholka/noSplits/species-period15-1patch/"
)

paths = paths[!grepl("genus", paths)]

cl = parallel::makeCluster(4L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
})

res = parallel::parLapply(cl, 1:length(paths), function(i) {

  myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
  dist = cbind(nodes,0:3)
  dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
  Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

  p = paths[i]
  Nepochs = Nepochs_vec[i]
  batchsize = batchsize_vec[i]
  Npatches = ifelse(grepl("9patches", p), 9L, 1L)

  env_dt = fread(paste0(p, "/env_dt.csv"))
  obs_dt = fread(paste0(p, "/obs_dt.csv"))
  cohorts_dt =  fread(paste0(p, "/initial_cohorts2000.csv"))

  Nspecies = max(obs_dt$species)
  # env_dt = env_dt[,-c("Psd","Psum")]
  Nenv = ncol(env_dt) - 2

  ### get init parameters ###
  env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))
  summary(env_obs)

    growth_init_model = lm(log(growth+1)~(Tmax+Tmin+Psum+Psd+awc):as.factor(species)+as.factor(species)+0, data = env_obs)
    init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
    init_growth_env[is.na(init_growth_env)] = 0
    init_growth_env[abs(init_growth_env) > 5 & init_growth_env < 0] = -5
    init_growth_env[abs(init_growth_env) > 5 & init_growth_env > 0] = 5

    ### mort ####
    env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
    mort_init_model = glmmTMB(Fmort~(Tmax+Tmin+Psum+Psd+awc):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
    init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
    init_mort_env[is.na(init_mort_env)] = 0
    init_mort_env[abs(init_mort_env) > 5 & init_mort_env < 0] = -5
    init_mort_env[abs(init_mort_env) > 5 & init_mort_env > 0] = 5

    ### reg ####
    reg_init_model = glmmTMB(reg~(Tmax+Tmin+Psum+Psd++awc):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
    init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
    init_reg_env[is.na(init_reg_env)] = 0
    init_reg_env[abs(init_reg_env) > 5 & init_reg_env < 0] = -5
    init_reg_env[abs(init_reg_env) > 5 & init_reg_env > 0] = 5

    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
      growth_process = createProcess(~., initEnv = list(init_growth_env |> as.matrix()), func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = T),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env |> as.matrix()), func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = T),
      mortality_process = createProcess(~., initEnv = list(init_mort_env  |> as.matrix()), func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = T)
    )

  cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)
  m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = 5L,
         optimizer = torch::optim_adam, device = "gpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = F,
         loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
  )

  if(!dir.exists(("results/Uholka/01_full/"))) dir.create("results/Uholka/01_full/", recursive = T)
  torch::torch_save(m1, path = paste0("results/Uholka/01_full/",basename(p), "_full.pt"))

  rm(m1)
  gc()
  torch::cuda_empty_cache()
})
