
library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

Nepochs = 8000
# timestamp_dt <- fread("input-timestamp.csv")
# timestamp_dt <- NULL
# timestamp <- timestamp_dt$timestamp

all_lossvars = c("growth", "mort", "reg", "trees", "dbh", "ba")
lossvars_rates <- c("growth", "mort", "reg")
# All combinations of 1, 2, or 3 elements from lossvars_rates:
combo_rates <- unlist(lapply(1:length(lossvars_rates), function(n) apply(combn(lossvars_rates, n), 2, paste, collapse=".")))
lossvars_comb <- c(paste("ba.trees.dbh", combo_rates, sep="."),
                   "ba.growth.mort.reg", "ba.trees.dbh.growth.mort.reg",
                   lossvars_rates,
                   "ba",
                   "ba.trees.dbh"
)
lossvars_comb <- lossvars_comb[grepl("ba", lossvars_comb) | grepl("trees.dbh", lossvars_comb)]
lossvars_comb

T_folds <- paste0("T",c(0, 1:2))
S_folds <- paste0("S", c(0, 1:5))
fold_names = expand.grid(list(S_folds,T_folds))
fold_names = paste0(fold_names$Var1, "_", fold_names$Var2)
cv_variants <- paste0(rep(fold_names,each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

directories <- list.files("data/Uholka/CVsplits-simdata", full.names = T)
directories <- rev(directories)

cl = parallel::makeCluster(8L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
})

i_cv= cv_variants[1]
i_dir = directories[1]
for(i_dir in directories){
  parallel::clusterExport(cl, varlist = list("i_folder"))
  .null = parLapply(cl, cv_variants, function(i_cv){
    cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
    cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
    response = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]
    i_name = basename(i_dir)


    myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    dist = cbind(nodes,0:3)
    dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))

    Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

    cat("\nread data:", i_name)
    cat("\nCV variant:", i_cv)
    env_dt = fread(paste0(i_dir,"/env_dt_",cv_S,"_",cv_T,"_train.csv"))
    env_dt = env_dt[,.(siteID, year, Tmax, Tmin, Psum, Psd, awc)]

    obs_dt = fread(paste0(i_dir,"/obs_dt_",cv_S,"_",cv_T,"_train.csv"))
    # only pick columns siteID species  year   dbh    ba trees growth  mort   reg r_mean_ha
    obs_dt = obs_dt[,.(siteID, species, year, dbh, ba, trees, growth, mort, reg, r_mean_ha)]

    cohorts_dt = fread(paste0(i_dir,"/initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))
    # only pick columns siteID   dbh species patchID census trees cohortID
    cohorts_dt = cohorts_dt[,.(siteID, dbh, species, patchID, census, trees, cohortID)]

    Nspecies = max(obs_dt$species)
    Nenv = ncol(env_dt) - 2

    # Prec+SR_kW_m2+RH_prc+T_max+T_min+swp
    # Tmax+Tmin+Psum+Psd+awc
    env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))
    ## get init parameters ####
    ### growth ####
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
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)

    Nsites = length(unique(obs_dt$siteID))
    batchsize = ceiling((Nsites)) # TODO change back
    Npatches = uniqueN(cohorts_dt$patchID)

    m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = 5L,
           optimizer = torch::optim_ignite_adam, device = "cpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = FALSE,
           loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
    )

    out_dir = paste0("results/","02_simulated/",i_name,"_",i_cv,".pt")
    if(!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir), recursive = T)
    torch::torch_save(m1, out_dir)

    rm(m1)
    gc()
    torch::cuda_empty_cache()
  })
}
