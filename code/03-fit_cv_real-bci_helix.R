library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

# logfile <- "log_03-fit_cv_real-bci.txt"
# if(file.exists(logfile)) file.remove(logfile)

Nepochs = 2
overwrite = F
lossvars_comb = "ba.trees.dbh.growth.mort.reg"
all_lossvars = c("ba", "trees", "dbh", "growth", "mort", "reg")

T_folds <- paste0("T",c(0, 1:2))
S_folds <- paste0("S", c(0, 1:5))
fold_names = expand.grid(list(S_folds,T_folds))
fold_names = paste0(fold_names$Var1, "_", fold_names$Var2)
cv_variants <- paste0(rep(fold_names,each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

directories <- list.files("data/BCI/CVsplits-realdata", full.names = T, pattern = "genus-")
directories = rev(directories)

cat("\nscript started")
all_variants = list()
for(i_dir in directories){
  for(i_cv in cv_variants){
  # all_variants = c(all_variants, list(list(i_dir = i_dir, cv_variants = cv_variants)))
  all_variants = c(all_variants, list(list(i_dir = i_dir, cv_variants = i_cv)))
  }
}

cl = parallel::makeCluster(32L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
})
i_var = all_variants[[1]]
# remove all objects in the environment but "overwrite" and "i_var"

# rm(list = setdiff(ls(), c("overwrite", "all_lossvars", "Nepochs", "i_var")))
parallel::clusterExport(cl, varlist = list("overwrite", "all_lossvars", "Nepochs"), envir = environment())
.null = parLapply(cl, all_variants, function(i_var){
    i_dir = i_var$i_dir
    i_cv = i_var$cv_variants
    cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
    cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
    response = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]
    i_name = basename(i_dir)
    out_dir = paste0("results/","02_realdata/",i_name,"_",i_cv,".pt")
    if(!file.exists(out_dir) | overwrite){

      # cat(paste("Starting", i_dir, i_cv, "at", Sys.time()), "\n", file = logfile, append = TRUE)

      myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
      dist = cbind(nodes,0:3)
      dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))

      Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

      cat("\nread data:", i_name)
      cat("\nCV variant:", i_cv)
      env_dt = fread(paste0(i_dir,"/env_dt_",cv_S,"_",cv_T,"_train.csv"))
      env_dt = env_dt[,.(siteID, year, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp)]

      obs_dt = fread(paste0(i_dir,"/obs_dt_",cv_S,"_",cv_T,"_train.csv"))
      # only pick columns siteID species  year   dbh    ba trees growth  mort   reg r_mean_ha
      obs_dt = obs_dt[,.(siteID, species, year, dbh, ba, trees, growth, mort, reg)]

      cohorts_dt = fread(paste0(i_dir,"/initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))
      # only pick columns siteID   dbh species patchID census trees cohortID
      cohorts_dt = cohorts_dt[,.(siteID, dbh, species, patchID, census, trees, cohortID)]

      Nspecies = max(obs_dt$species)
      Nenv = ncol(env_dt) - 2

      env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))
      ## get init parameters ####
      ### growth ####
      if(grepl("genus",i_name)){
        env_obs$species_fac <- as.factor(env_obs$species)
        growth_init_model = glmmTMB(log(growth+1)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = gaussian())
        init_growth_env = as.matrix(coef(growth_init_model)[[1]][[1]])
        init_growth_env[is.na(init_growth_env)] = 0
      }else if(grepl("pft", i_name)){
        growth_init_model = lm(log(growth+1)~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs)
        init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
        init_growth_env[is.na(init_growth_env)] = 0
      }
      ### mort ####
      if(grepl("genus",i_name)){
        env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
        env_obs$species_fac <- as.factor(env_obs$species)
        mort_init_model = glmmTMB(Fmort~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = beta_family())
        init_mort_env = as.matrix(coef(mort_init_model)[[1]][[1]])
        init_mort_env[is.na(init_mort_env)] = 0
      }else if(grepl("pft", i_name)){
        env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
        mort_init_model = glmmTMB(Fmort~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
        init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
        init_mort_env[is.na(init_mort_env)] = 0
      }
      ### reg ####
      if(grepl("genus",i_name)){
        reg_init_model = glmmTMB(ceiling(reg)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = "nbinom1")
        init_reg_env = as.matrix(coef(reg_init_model)[[1]][[1]])
        init_reg_env[is.na(init_reg_env)] = 0
      }else if(grepl("pft", i_name)){
        reg_init_model = glmmTMB(reg~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
        init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
        init_reg_env[is.na(init_reg_env)] = 0
      }

      max_val = 5
      init_growth_env[init_growth_env > max_val] = max_val
      init_growth_env[init_growth_env < -max_val] = -max_val
      init_mort_env[init_mort_env > max_val] = max_val
      init_mort_env[init_mort_env < -max_val] = -max_val
      init_reg_env[init_reg_env > max_val] = max_val
      init_reg_env[init_reg_env < -max_val] = -max_val

      obsNA = unlist(strsplit(response, ".", fixed = TRUE))
      for(i_lossvar in all_lossvars[!(all_lossvars %in% obsNA)]){
        obs_dt[[i_lossvar]] = NA_real_
      }

      ## create model ####
      m1 = finn(
        N_species = Nspecies,
        competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE, optimizeEnv = TRUE),
        growth_process = createProcess(~., initEnv = list(as.matrix(init_growth_env)), func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
        regeneration_process = createProcess(~., initEnv = list(init_reg_env), func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
        mortality_process = createProcess(~., initEnv = list(init_mort_env), func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE)
      )

      ## fit model ####
      cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)

      Nsites = length(unique(obs_dt$siteID))
      batchsize = ceiling((Nsites)/2) # TODO change back
      Npatches = uniqueN(cohorts_dt$patchID)

      m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = 5L,
             optimizer = torch::optim_ignite_adam, device = "gpu", record_gradients = FALSE, weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = FALSE,
             loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
      )

      if(!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir), recursive = T)
      torch::torch_save(m1, out_dir)

      #cat(paste("Finished", i_dir, i_cv, "at", Sys.time()), "\n", file = logfile, append = TRUE)
      rm(m1)
      gc()
      torch::cuda_empty_cache()
    }
})
# }
# }
parallel::stopCluster(cl)
cat("\nscript finished")
