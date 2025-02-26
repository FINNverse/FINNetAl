
library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

Nepochs = 2000
# timestamp_dt <- fread("input-timestamp.csv")
# timestamp_dt <- NULL
# timestamp <- timestamp_dt$timestamp

# cl = makeCluster(42)
# parallel::clusterEvalQ(cl, {library(data.table); library(FINN); library(torch); library(glmmTMB)})
# nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
# parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

lossvars_rates = c("growth", "mort", "reg")
lossvars_stand = c("ba", "trees")
all_lossvars = c(lossvars_stand, lossvars_rates)
# get all combinations of lossvars
lossvars_comb = lapply(1:length(lossvars_rates), function(x) unlist(combn(lossvars_rates, x, simplify = T)))
lossvars_comb = lapply(lossvars_comb, function(x) apply(x,2, function(y) paste0(y, collapse=".")))
lossvars_comb <- lapply(lossvars_comb, function(x) paste0(c("ba.","ba.trees."),rep(x,each = 2)))

T_folds <- paste0("T",c(0, 1:2))
S_folds <- paste0("S", c(0, 1:5))
fold_names = expand.grid(list(S_folds,T_folds))
fold_names = paste0(fold_names$Var1, "_", fold_names$Var2)
cv_variants <- paste0(rep(fold_names,each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

i_cv = cv_variants[1]

cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
response = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]
folders <- list.files("simulated_data/sim-split/CVfolds", pattern = "period")
i_folder = folders[1]

# only genus for now
folders = rev(folders[grepl("pft", folders)])


cl = parallel::makeCluster(20L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
})

for(i_folder in folders){
  parallel::clusterExport(cl, varlist = list("i_folder"))

  .null = parLapply(cl, cv_variants, function(i_cv){
    cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
    cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
    response = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]


    myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    dist = cbind(nodes,0:3)
    dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))

    Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

    cat("\nread data:", i_folder)
    cat("\nCV variant:", i_cv)
    env_dt = fread(paste0("simulated_data/sim-split/CVfolds/",i_folder,"/env_dt_",cv_S,"_",cv_T,"_train.csv"))[,-(9:11)]
    obs_dt = fread(paste0("simulated_data/sim-split/CVfolds/",i_folder,"/obs_dt_",cv_S,"_",cv_T,"_train.csv"))[,-(11:13)]
    cohorts_dt = fread(paste0("simulated_data/sim-split/CVfolds/",i_folder,"/initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))[,-(8:10)]

    Nspecies = max(obs_dt$species)
    Nenv = ncol(env_dt) - 2

    env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))
    ## get init parameters ####
    ### growth ####
    if(grepl("genus",i_folder)){
      env_obs$species_fac <- as.factor(env_obs$species)
      growth_init_model = glmmTMB(log(growth+1)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = gaussian())
      init_growth_env = coef(growth_init_model)[[1]][[1]]
      init_growth_env[is.na(init_growth_env)] = 0
    }else if(grepl("pft", i_folder)){
      growth_init_model = lm(log(growth+1)~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs)
      init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
      init_growth_env[is.na(init_growth_env)] = 0
    }
    ### mort ####
    if(grepl("genus",i_folder)){
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      env_obs$species_fac <- as.factor(env_obs$species)
      mort_init_model = glmmTMB(Fmort~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = beta_family())
      init_mort_env = coef(mort_init_model)[[1]][[1]]
      init_mort_env[is.na(init_mort_env)] = 0
    }else if(grepl("pft", i_folder)){
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      mort_init_model = glmmTMB(Fmort~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
      init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_mort_env[is.na(init_mort_env)] = 0
    }
    ### reg ####
    if(grepl("genus",i_folder)){
      reg_init_model = glmmTMB(ceiling(reg)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = "nbinom1")
      init_reg_env = coef(reg_init_model)[[1]][[1]]
      init_reg_env[is.na(init_reg_env)] = 0
    }else if(grepl("pft", i_folder)){
      reg_init_model = glmmTMB(reg~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
      init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_reg_env[is.na(init_reg_env)] = 0
    }


    obsNA = unlist(strsplit(response, ".", fixed = TRUE))
    for(i_lossvar in all_lossvars[!(all_lossvars %in% obsNA)]){
      obs_dt[[i_lossvar]] = NA_real_
    }

    ## create model ####
    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      growth_process = createProcess(~., initEnv = list(init_growth_env), func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env), func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      mortality_process = createProcess(~., initEnv = list(init_mort_env), func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE)
    )

    ## fit model ####
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt[,.(siteID, patchID, cohortID, species, dbh = round(dbh_cm,4), trees)], sp = Nspecies)

    Nsites = length(unique(obs_dt$siteID))
    batchsize = ceiling((Nsites)) # TODO change back
    Npatches = uniqueN(cohorts_dt$patchID)

    m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = 5L,
           optimizer = torch::optim_ignite_adam, device = "gpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = FALSE,
           loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
    )

    out_dir = paste0("results/","02_simulated/",i_folder,"_",i_cv,".pt")
    torch::torch_save(m1, out_dir)

    rm(m1)
    gc()
    torch::cuda_empty_cache()
  })
}
