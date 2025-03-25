#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## get index of array ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No batch index provided. Please provide a batch index as an argument.")
}
# batch_index <- 1
batch_index <- as.integer(args[1])
# Each batch contains 4 indices; for example, batch 1 -> 1:4, batch 2 -> 5:8, etc.
jobs_per_process = 6L
subset_indices <- (((batch_index - 1) * jobs_per_process) + 1):(batch_index * jobs_per_process)
cat("Running batch index:", batch_index, "\n")
cat("Processing tasks with indices:", paste(subset_indices, collapse = ", "), "\n")

library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)

Nepochs = 8000
overwrite = F
simmodel_path = "results/01_full/pft-period7-25patches_full.pt"
reps = 1:100

cat("\nscript started")
all_jobs = list()
for(i_rep in reps){
  for(i_simmodel_path in simmodel_path){
    all_jobs = c(all_jobs, list(list(rep = i_rep, simmodel_path = i_simmodel_path)))
  }
}

if(min(subset_indices) > length(all_jobs)) {
  stop("Batch index exceeds available task indices (1:length(all_jobs)).")
}
subset_indices <- subset_indices[subset_indices <= length(all_jobs)]

# Process only the subset for this batch
.selected_variants <- all_jobs[subset_indices]

cl = parallel::makeCluster(12L)
nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
  torch::torch_set_num_threads(8L)
})


# i_dir = directories[3]
# i_cv = cv_variants[15]
cat("\nscript started")
# for(i_dir in directories){
# for(i_cv in cv_variants){
i_var = .selected_variants[[2]]
parallel::clusterExport(cl, varlist = c(ls(envir = .GlobalEnv)), envir = environment())
.null = parLapply(cl, .selected_variants, function(i_var){
  i_rep = i_var$rep
  i_simmodel_path = i_var$simmodel_path
  i_name = gsub(".pt","",basename(i_simmodel_path))
  i_folder = gsub("_full","",i_name)
  out_dir = paste0("results/03_full100reps/",i_name,"_rep",i_rep,".pt")
  if(!file.exists(out_dir) | overwrite){

    myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    dist = cbind(nodes,0:3)
    dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))

    Sys.setenv(CUDA_VISIBLE_DEVICES=dev)

    obs_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/obs_dt.csv"))
    env_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/env_dt.csv"))
    cohorts_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/initial_cohorts1985.csv"))
    m <- torch::torch_load(i_simmodel_path)

    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = uniqueN(obs_dt$species))
    pred = m$simulate(env = env_dt, init_cohort = cohort1, patches = Npatches)
    pred_dt = pred$wide$site
    pred_dt[,reg := reg/0.1,]

    Nspecies = max(obs_dt$species)
    Nenv = ncol(env_dt) - 2

    env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))
    ## get init parameters ####
    if(grepl("genus",i_name)){
      ### growth ####
      env_obs$species_fac <- as.factor(env_obs$species)
      growth_init_model = glmmTMB(log(growth+1)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = gaussian())
      init_growth_env <- as.matrix(coef(growth_init_model)[[1]][[1]])
      init_growth_env[is.na(init_growth_env)] = 0
      ### mort ####
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      env_obs$species_fac <- as.factor(env_obs$species)
      mort_init_model = glmmTMB(Fmort~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = beta_family())
      init_mort_env = as.matrix(coef(mort_init_model)[[1]][[1]])
      init_mort_env[is.na(init_mort_env)] = 0
      ### reg ####
      reg_init_model = glmmTMB(ceiling(reg)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = "nbinom1")
      init_reg_env = as.matrix(coef(reg_init_model)[[1]][[1]])

      ### fix missing species ####
      init_reg_env[is.na(init_reg_env)] = 0
      missing_species_growth <- setdiff(1:Nspecies, as.integer(rownames(init_growth_env)))
      init_growth_env <- rbind(
        init_growth_env,
        matrix(0, length(missing_species_growth), ncol(init_growth_env), dimnames = list(missing_species_growth, colnames(init_growth_env)))
      )
      init_growth_env <- init_growth_env[order(as.numeric(rownames(init_growth_env))), ]

      missing_species_mort <- setdiff(1:Nspecies, as.integer(rownames(init_mort_env)))
      init_mort_env <- rbind(
        init_mort_env,
        matrix(0, length(missing_species_mort), ncol(init_mort_env), dimnames = list(missing_species_mort, colnames(init_mort_env)))
      )
      init_mort_env <- init_mort_env[order(as.numeric(rownames(init_mort_env))), ]

      missing_species_reg <- setdiff(1:Nspecies, as.integer(rownames(init_reg_env)))
      init_reg_env <- rbind(
        init_reg_env,
        matrix(0, length(missing_species_reg), ncol(init_reg_env), dimnames = list(missing_species_reg, colnames(init_reg_env)))
      )
      init_reg_env <- init_reg_env[order(as.numeric(rownames(init_reg_env))), ]
    }else if(grepl("pft", i_name)){
      ### growth ####
      growth_init_model = lm(log(growth+1)~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs)
      init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
      init_growth_env[is.na(init_growth_env)] = 0
      ### mort ####
      env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
      mort_init_model = glmmTMB(Fmort~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
      init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_mort_env[is.na(init_mort_env)] = 0
      ### reg ####
      reg_init_model = glmmTMB(reg~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
      init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
      init_reg_env[is.na(init_reg_env)] = 0
    }

    max_vals = 1
    init_growth_env[init_growth_env > max_vals] = max_vals
    init_growth_env[init_growth_env < -max_vals] = -max_vals
    init_mort_env[init_mort_env > max_vals] = max_vals
    init_mort_env[init_mort_env < -max_vals] = -max_vals
    init_reg_env[init_reg_env > max_vals] = max_vals
    init_reg_env[init_reg_env < -max_vals] = -max_vals

    ## create model ####
    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      growth_process = createProcess(~., initEnv = list(init_growth_env), func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env), func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      mortality_process = createProcess(~., initEnv = list(init_mort_env), func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE)
    )

    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)

    Nsites = length(unique(obs_dt$siteID))
    batchsize = ceiling(((Nsites)/2)) # TODO change back
    Npatches = uniqueN(cohorts_dt$patchID)

    m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = Inf,
           optimizer = torch::optim_ignite_adam, device = "gpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 10.0, 1, 1), plot_progress = FALSE,
           loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
    )

    if(!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir), recursive = T)
    m1$history = NULL
    torch::torch_save(m1, out_dir)

    rm(m1)
    gc()
    # torch::cuda_empty_cache()
  }else{
    cat("\nModel already exists:", out_dir)
  }
})
# }
# }
parallel::stopCluster(cl)
cat("\nscript finished")
