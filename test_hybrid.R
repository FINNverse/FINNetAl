#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## get index of array ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No batch index provided. Please provide a batch index as an argument.")
}
# batch_index <- 10
batch_index <- as.integer(args[1])
# Each batch contains 4 indices; for example, batch 1 -> 1:4, batch 2 -> 5:8, etc.
jobs_per_process = 2L
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
lossvars_comb = "ba.trees.dbh.growth.mort.reg"
all_lossvars = c("ba", "trees", "dbh", "growth", "mort", "reg")

# T_folds <- paste0("T",c(0, 1:2))
# S_folds <- paste0("S", c(0, 1:5))
T_folds <- paste0("T",c(0, 2))
S_folds <- paste0("S", c(0, 3))
transformer_folds = c("TF0","TF1")

fold_names = expand.grid(list(S_folds,T_folds,transformer_folds))
fold_names = paste0(fold_names$Var1, "_", fold_names$Var2, "_", fold_names$Var3)
cv_variants <- paste0(rep(fold_names,each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

directories <- list.files("data/BCI/CVsplits-realdata", full.names = T)

# directories <- directories[grepl("1patch", directories) & grepl("genus", directories)]
directories = rev(directories)

# for(i_dir in directories){
#   for(i_cv in cv_variants){
#     cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
#     cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
#     response = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]
#     i_name = basename(i_dir)
#     out_dir = paste0("results/","02_realdata/",i_name,"_",i_cv,".pt")
#     if(!file.exists(out_dir)) stop()
#   }
# }

cat("\nscript started")
all_variants = list()
for(i_dir in directories){
  for(i_cv in cv_variants){
    all_variants = c(all_variants, list(list(i_dir = i_dir, cv_variants = i_cv)))
  }
}

if(min(subset_indices) > length(all_variants)) {
  stop("Batch index exceeds available task indices (1:length(all_variants)).")
}
subset_indices <- subset_indices[subset_indices < length(all_variants)]

# Process only the subset for this batch
.selected_variants <- all_variants[subset_indices]

cl = parallel::makeCluster(jobs_per_process)
# nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(FINN)
  library(torch)
  library(glmmTMB)
  torch::torch_set_num_threads(8)
})


# i_dir = directories[3]
# i_cv = cv_variants[15]
cat("\nscript started")
# for(i_dir in directories){
# for(i_cv in cv_variants){
# i_var = .selected_variants[[2]]
i_var = all_variants[[2]]
parallel::clusterExport(cl, varlist = c(ls(envir = .GlobalEnv)), envir = environment())
.null = parLapply(cl, .selected_variants, function(i_var){
  i_dir = i_var$i_dir
  i_cv = i_var$cv_variants
  cv_S = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
  cv_T = tstrsplit(i_cv, "_", fixed = TRUE)[[2]][1]
  tf = tstrsplit(i_cv, "_", fixed = TRUE)[[3]][1]
  if(tf == "TF0") transformer = FALSE
  if(tf == "TF1") transformer = TRUE
  response = tstrsplit(i_cv, "_", fixed = TRUE)[[4]][1]
  obsNA = unlist(strsplit(response, ".", fixed = TRUE))
  
  i_name = basename(i_dir)
  out_dir = paste0("results/","02_realdata_hybrid",tf,"/",i_name,"_",cv_S,"_",cv_T,".pt")
  stopifnot(length(overwrite) > 0)
  stopifnot(length(all_lossvars) > 0)
  stopifnot(sum(all_lossvars %in% obsNA) > 0)
  
  
  
  if(!file.exists(out_dir) | overwrite){
    
    # cat(paste("Starting", i_dir, i_cv, "at", Sys.time()), "\n", file = logfile, append = TRUE)
    
    # myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    # dist = cbind(nodes,0:3)
    # dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
    
    # Sys.setenv(CUDA_VISIBLE_DEVICES=dev)
    
    cat("\nread data:", i_name)
    cat("\nCV variant:", i_cv)
    env_dt = fread(paste0(i_dir,"/env_dt_",cv_S,"_",cv_T,"_train.csv"))
    env_dt = env_dt[,.(siteID, year, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp )]
    
    obs_dt = fread(paste0(i_dir,"/obs_dt_",cv_S,"_",cv_T,"_train.csv"))
    # only pick columns siteID species  year   dbh    ba trees growth  mort   reg r_mean_ha
    obs_dt = obs_dt[,.(siteID, species, year, dbh, ba, trees, growth, mort, reg)]
    
    cohorts_dt = fread(paste0(i_dir,"/initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))
    # only pick columns siteID   dbh species patchID census trees cohortID
    cohorts_dt = cohorts_dt[,.(siteID, dbh, species, patchID, census, trees, cohortID)]
    
    obs_dt[dbh == 0, mort := NA_real_]
    obs_dt[dbh == 0, growth := NA_real_]
    obs_dt[dbh == 0, dbh := NA_real_]
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
    
    for(i_lossvar in all_lossvars[!(all_lossvars %in% obsNA)]){
      obs_dt[[i_lossvar]] = NA_real_
    }
    
    hybrid_growth = createHybrid(~., transformer = TRUE, emb_dim = 6L,dropout = 0.2, encoder_layers = 1L, dim_feedforward = 100L)
    hybrid_mort = createHybrid(~., transformer = TRUE, emb_dim = 6L,dropout = 0.2, encoder_layers = 1L, dim_feedforward = 100L)
    
    ## create model ####
    m1 = finn(
      N_species = Nspecies,
      competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      growth_process = hybrid_growth,
      #growth_process = createProcess(~., initEnv = list(init_growth_env), func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      regeneration_process = createProcess(~., initEnv = list(init_reg_env), func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      #mortality_process = createProcess(~., initEnv = list(init_mort_env), func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
      mortality_process = hybrid_mort
    )
    
    gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
      g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
      return(g)
    }
    m1$growth_func = m1$.__enclos_env__$private$set_environment(gh)
    
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)
    
    Nsites = length(unique(obs_dt$siteID))
    batchsize = ceiling(((Nsites))*1) # TODO change back
    Npatches = uniqueN(cohorts_dt$patchID)
    
    source("code/99_cohort500fix.R")
    m1$fit(data = obs_dt, batchsize = batchsize, env = env_dt, init_cohort = cohort1,  epochs = Nepochs, patches = Npatches, lr = 0.01, checkpoints = Inf,
           optimizer = torch::optim_ignite_adam, device = "gpu", record_gradients = FALSE,weights = c(0.1, 10, 1.0, 1, 1, 1), plot_progress = F,
           loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
    )
    
    if(!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir), recursive = T)
    torch::torch_save(m1, out_dir)
    
    
    gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
      self$nn_growth$train()
      g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
      #print(g$var()$item())
      return(g)
    }
    m1$growth_func = m1$.__enclos_env__$private$set_environment(gh)
    
    pp = m1$simulate(env_dt, init_cohort = cohort1, patches = Npatches)
    plot(pp$wide$site[, .(growth = mean(growth)), by = .(siteID, species)][order(siteID, species)][siteID %in% 1:15]$growth,
         obs_dt[, .(growth = mean(growth)), by = .(siteID, species)][order(siteID, species)][siteID %in% 1:15]$growth, method = "spearman")
    
    
    #cat(paste("Finished", i_dir, i_cv, "at", Sys.time()), "\n", file = logfile, append = TRUE)
    rm(m1)
    gc()
    torch::cuda_empty_cache()
  }else{
    cat("\nModel already exists:", out_dir)
  }
})
# }
# }
parallel::stopCluster(cl)
cat("\nscript finished")


gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
  #self$nn_growth$train()
  g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
  #print(g$var()$item())
  return(g)
}
m1$growth_func = m1$.__enclos_env__$private$set_environment(gh)
mh = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F) {
  #self$nn_mortality$train()
  m = self$nn_mortality(dbh = dbh, growth = self$g, trees = trees, light = light, species = species, env = pred)$sigmoid()
  return(m)
}
m1$mortality_func = m1$.__enclos_env__$private$set_environment(mh)



pred = m1$simulate(env = env_dt, patches = patches, init_cohort = cohort1)
pred_dt <- pred$wide$site
pred_dt$pred <- "pred"

# Correct the regeneration response
pred_dt[, reg := reg / 0.1]
obs_dt$pred <- "obs"

# Combine observed & predicted
comp_dt1 <- rbind(pred_dt, obs_dt, fill = TRUE)
# Add "all" species row
comp_dt1 <- rbind(
  comp_dt1,
  comp_dt1[
    ,
    .(
      species = "all",
      dbh     = mean(dbh, na.rm = TRUE),
      ba      = sum(ba),
      trees   = sum(trees),
      reg     = mean(reg, na.rm = TRUE),
      mort    = mean(mort, na.rm = TRUE),
      growth  = mean(growth, na.rm = TRUE)
    ),
    by = .(siteID, year, pred)
  ],
  fill = TRUE
)

# Melt for ggplot
melt_dt <- melt(comp_dt1,
                id.vars = c("siteID", "year", "species", "pred"),
                variable.name = "variable")
p_dat <- melt_dt[, .(value = mean(value)), by = .(species, year, pred, variable)]

library(ggplot2)
# Plot with ggplot
p1 <- ggplot(p_dat[variable != "r_mean_ha"],
             aes(x = year, y = value, color = factor(species))) +
  geom_point(aes(shape = pred, size = pred)) +
  geom_line(aes(linetype = pred)) +
  scale_size_manual(values = c(pred = 2, obs = 3)) +
  facet_wrap(~variable, scales = "free_y") +
  ggtitle(i_name)
print(p1)

df <- merge.data.table(obs_dt, pred_dt,
                       by = c("siteID", "year", "species"),
                       suffixes = c(".obs", ".pred"))

comp_allspecies_dt <- df[
  ,
  .(
    ba.obs     = sum(ba.obs) / uniqueN(siteID),
    ba.pred    = sum(ba.pred) / uniqueN(siteID),
    trees.obs  = sum(trees.obs) / uniqueN(siteID),
    trees.pred = sum(trees.pred) / uniqueN(siteID),
    dbh.obs    = mean(dbh.obs, na.rm = TRUE),
    dbh.pred   = mean(dbh.pred, na.rm = TRUE),
    reg.obs    = mean(reg.obs, na.rm = TRUE),
    reg.pred   = mean(r_mean_ha, na.rm = TRUE),
    mort.obs   = mean(mort.obs, na.rm = TRUE),
    mort.pred  = mean(mort.pred, na.rm = TRUE),
    growth.obs = mean(growth.obs, na.rm = TRUE),
    growth.pred= mean(growth.pred, na.rm = TRUE)
  ),
  by = .(siteID, year, species)
]
comp_allspecies_dt[,.(ba = cor(ba.obs, ba.pred, use = "complete.obs", method = "spearman"),
                      dbh = cor(dbh.obs, dbh.pred, use = "complete.obs", method = "spearman"),
                      trees = cor(trees.obs, trees.pred, use = "complete.obs", method = "spearman"),
                      growth = cor(growth.obs, growth.pred, use = "complete.obs", method = "spearman"),
                      mort = cor(mort.obs, mort.pred, use = "complete.obs", method = "spearman"),
                      reg = cor(reg.obs, reg.pred, use = "complete.obs", method = "spearman"))] %>% colMeans()


hist(obs_dt$growth)
