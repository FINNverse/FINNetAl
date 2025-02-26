library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
# library(parallel)
# setwd("~/working-directory/FINN-BCI/cluster/BCI-calibration/20250123_171017")
setwd("~/working-directory/FINN-BCI/cluster/20250219_101440/pft-period35")
# setwd("~/working-directory/FINN-BCI/cluster/BCI-stand-data-pft/20250218_100130")
Nepochs = 5000
timestamp_dt <- fread("input-timestamp.csv")
# timestamp_dt <- NULL
timestamp <- timestamp_dt$timestamp

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
cv_variants <- paste0(rep(c("S","T","ST"),each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

# cv_variants <- c("S_all","T_all","ST_all","S_norates", "T_norates","ST_norates")
i_cv = "F_ba.trees.growth.mort.reg"
# i_cv = "S_ba.trees.growth.reg"
# .null =
#   parallel::parSapplyLB(cl, cv_variants, function(i_cv) {
# for(i_cv in cv_variants){

# who am I
# myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
# dist = cbind(nodes,0:7)
# dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
#
# Sys.setenv(CUDA_VISIBLE_DEVICES=dev)
all_trees = fread(paste0("processed_data/calibration-data/BCI-stand-data-pft/all_trees.csv"))

all_trees[,.(gdbh_cm)]

cat("\nCV variant:", i_cv)
cat("\nread data")
cv_data = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
env_dt = fread(paste0("processed_data/calibration-data/BCI/pft-period35/blocked/env_dt_train",cv_data,".csv"))
round(cor(env_dt),2)
obs_dt = fread(paste0("processed_data/calibration-data/BCI/pft-period35/blocked/obs_dt_train",cv_data,".csv"))
# obs_dt$growth = obs_dt$growth*5
# obs_dt[year != 1, period_length := 5]
# obs_dt$AL <- as.numeric(obs_dt$AL)
# obs_dt[,trees_before := data.table::shift(trees), by = .(siteID, species)]
# obs_dt[,mort := mort*5*trees_before,]
cohorts_dt = fread(paste0("processed_data/calibration-data/BCI/pft-period35/blocked/initial_cohorts_train",cv_data,".csv"))

i_lossvars_vec = unlist(tstrsplit(tstrsplit(i_cv, "_", fixed = TRUE)[[2]], ".", fixed = TRUE))
# for(i_lossvar in all_lossvars[!(all_lossvars %in% i_lossvars_vec)]){
#   obs_dt[[i_lossvar]] = NA_real_
# }

Nspecies = max(obs_dt$species)
Nenv = ncol(env_dt) - 2

### get init parameters ###
env_obs = merge(obs_dt, env_dt, by = c('year', "siteID"))

cat("estimating initial growth parameters")
if(grepl("growth",i_cv)){
  # growth
  growth_init_model = lm(log(growth+1)~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs)
  # growth_init_model = glmmTMB(log(growth+1)~1+Prec+SR_kW_m2+RH_prc+T_max+T_min+swp+(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp||species_fac), data = env_obs, family = gaussian())
  # write model to file
  saveRDS(growth_init_model, file = paste0("BCI-linear_growth_init_model_",i_cv,".RDS"))
  # init_growth_env = matrix(summary(growth_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
  init_growth_env = matrix(coefficients(growth_init_model), Nspecies, Nenv+1)
  init_growth_env[is.na(init_growth_env)] = 0
}else{
  init_growth_env = NULL
}

cat("estimating initial mort parameters")
if(grepl("mort",i_cv)){
  # mort
  env_obs$Fmort = scales::rescale(env_obs$mort, c(0.0001, 1-0.0001))
  env_obs$species_fac <- as.factor(env_obs$species)
  mort_init_model = glmmTMB(Fmort~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = beta_family())
  # write model to file
  saveRDS(mort_init_model, file = paste0("BCI-linear_mort_init_model_",i_cv,".RDS"))
  init_mort_env = matrix(summary(mort_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
  init_mort_env[is.na(init_mort_env)] = 0
}else{
  init_mort_env = NULL
}

cat("estimating initial reg parameters")
# if(grepl("reg",i_cv)){
if(T){
  reg_init_model = glmmTMB(reg~(Prec+SR_kW_m2+RH_prc+T_max+T_min+swp):as.factor(species)+as.factor(species)+0, data = env_obs, family = "nbinom1")
  # write model to file
  saveRDS(reg_init_model, file = paste0("BCI-linear_reg_init_model_",i_cv,".RDS"))
  init_reg_env = matrix(summary(reg_init_model)$coefficients$cond[,1], Nspecies, Nenv+1)
  init_reg_env[is.na(init_reg_env)] = 0
}else{
  init_reg_env = NULL
}

cat("\ncreate cohort array")
cohorts_dt_in <- data.table()
for(i in 1:1){
  cohorts_dt_in <- rbind(
    cohorts_dt_in,
    cohorts_dt[,.(siteID, patchID = i, cohortID, species, dbh = round(dbh_cm,4), trees)]
  )
}
cohort1 <- FINN::CohortMat(obs_df = cohorts_dt_in, sp = Nspecies)

cat("\nrun model for 1 epoch")
## Simulate
mortality2 = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F) {
  environment = torch::torch_sigmoid(pred+self$g*parMort[,3][species]+light*parMort[,1][species] +  parMort[,2][species]*(dbh / ( 100)))
  predM = environment
  self$predM = predM
  return(predM)
}

m1 = finn(
  N_species = Nspecies,
  competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
  growth_process = createProcess(~., initEnv = init_growth_env, func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(
    ~., initEnv = init_reg_env, func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process = createProcess(
    ~., initEnv = init_mort_env, func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
m1$register_buffer("par_mortality_lower", torch::torch_tensor(c(-4, -4, -4)*2))
m1$register_buffer("par_mortality_upper", torch::torch_tensor(c(4, 4, 4)*2))
m1$par_mortality_unconstrained = torch::nn_parameter(cbind(rep(0.0, Nspecies), rep(0.0, Nspecies), rep(0.0, Nspecies)) %>% torch_tensor())

cohorts_dt_in <- cohorts_dt[,.(siteID, patchID = 1, cohortID, species, dbh = round(dbh_cm,4), trees)]
cohorts_dt_in_p = lapply(1:1, function(p) {
  patch = copy(cohorts_dt_in)
  patch$patchID = p
  return(patch)
})

cohort1 <- FINN::CohortMat(obs_df = data.table::rbindlist(cohorts_dt_in_p), sp = Nspecies)
#tmp = copy(obs_dt)
#tmp[year==2]$growth = NA_real_
m1$fit(data = obs_dt, batchsize = 200L, env = env_dt, init_cohort = cohort1,  epochs = 1000L, patches = 1L, lr = 0.01, checkpoints = 5L,
       optimizer = torch::optim_adam, device = "cpu", record_gradients = FALSE,
       # weights = c(0.1, 1, 1.0, 5, 10, 1),
       weights = c(1, 1, 1, 1, 1, 1),
       loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
)
torch::torch_save(m1,"BCI-linear-model-test.torch")
# m1 = saveRDS(m1,"BCI-linear_final_model.RDS")
# m1 = readRDS("BCI-linear_final_model.RDS")
#
# m2 = finn(
#   N_species = Nspecies,
#   competition_process = createProcess(~0, func = FINN::competition, initSpecies = m1$parameters_r),
#   growth_process = createProcess(~., initEnv = init_growth_env, func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
#   regeneration_process = createProcess(
#     ~., initEnv = init_reg_env, func = FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE, sample_regeneration = T),
#   mortality_process = createProcess(
#     ~., initEnv = init_mort_env, func = mortality2, optimizeSpecies = TRUE, optimizeEnv = TRUE,
#     initSpecies = cbind(
#       rep(0.5, Nspecies), # mortality rate
#       rep(0.2, Nspecies), # dbh threshold
#       rep(0, Nspecies) # g
#     ),
#     lower = c(-4, -4, -4), upper = c(4, 4, 4) )
# )

table(round(m1$pred$wide$site$reg), m1$pred$wide$site$year)
# true_reg = function(species, parReg, pred, light, debug = F) {
#   reg = y[,,,6]$mean(2) # site, Arten
#   return(reg$unsqueeze(2L))
# # }
#
# # save model for back transformation of parameter list
# saveRDS(m1, file = paste0("BCI-linear_parameters_1epoch-model",i_cv,".RDS"))
# rm(m1)
# gc()
# torch::cuda_empty_cache()
#
# saveRDS(m1, file = paste0("BCI-linear_final_model_",timestamp,i_cv,".RDS"))
#
# cat("\nscript finished")
# # })

pred = m1$simulate(env_dt, init_cohort = cohort1, patches = 100)

library(ggplot2)

# pred = m1$pred

# plot(obs_dt[,.(mean_reg = mean(reg, na.rm = TRUE)), by = .(species, siteID)][order(species)]$mean_reg,
# data.table(pred$wide$site)[,.(mean_reg = mean(reg, na.rm = TRUE)*10), by = .(species, siteID)][order(species)]$mean_reg)


ss = pred$wide$site
# obs_dt$year = as.integer(as.factor(obs_dt$year))
# obs_dt$siteID = as.integer(as.factor(obs_dt$siteID))
# obs_dt$species = as.integer(as.factor(obs_dt$species))
ss$year = as.integer(as.factor(ss$year))
ss$siteID = as.integer(as.factor(ss$siteID))
ss$species = as.integer(as.factor(ss$species))
df = merge.data.table(obs_dt, ss, by = c("siteID", "year", "species"), suffixes = c(".obs", ".pred"))
comp_allspecies_dt <- df[,.(
  ba.obs = sum(ba.obs, na.rm = T)/uniqueN(siteID),
  ba.pred = sum(ba.pred, na.rm = T)/uniqueN(siteID),
  trees.obs = sum(trees.obs, na.rm = T)/uniqueN(siteID),
  trees.pred = sum(trees.pred, na.rm = T)/uniqueN(siteID),
  dbh.obs = mean(dbh.obs, na.rm = TRUE),
  dbh.pred = mean(dbh.pred, na.rm = TRUE),
  reg.obs = mean(reg.obs, na.rm = TRUE),
  reg.pred = mean(reg.pred, na.rm = TRUE),
  mort.obs = mean(mort.obs, na.rm = TRUE),
  mort.pred = mean(mort.pred, na.rm = TRUE),
  growth.obs = mean(growth.obs, na.rm = T),
  growth.pred = mean(growth.pred, na.rm = T)
), by = .(siteID,species,year)]

p_dat = melt(comp_allspecies_dt[,.(
  ba.obs = sum(ba.obs, na.rm = T)/uniqueN(siteID),
  trees.obs = sum(trees.obs, na.rm = T)/uniqueN(siteID),
  dbh.obs = mean(dbh.obs, na.rm = TRUE),
  reg.obs = mean(reg.obs, na.rm = TRUE),
  mort.obs = mean(mort.obs, na.rm = TRUE),
  growth.obs = mean(growth.obs, na.rm = T),
  ba.pred = sum(ba.pred, na.rm = T)/uniqueN(siteID),
  trees.pred = sum(trees.pred, na.rm = T)/uniqueN(siteID),
  dbh.pred = mean(dbh.pred, na.rm = TRUE),
  reg.pred = mean(reg.pred, na.rm = TRUE),
  mort.pred = mean(mort.pred, na.rm = TRUE),
  growth.pred = mean(growth.pred, na.rm = T)
), by = .(species,year)], id.vars = c("year","species"))
p_dat[, pred_obs := tstrsplit(variable, ".", fixed = TRUE)[[2]]]
p_dat[, variable := tstrsplit(variable, ".", fixed = TRUE)[[1]]]
ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
  geom_point(aes(shape = factor(species))) +
  geom_line(aes(linetype = pred_obs, ))+
  facet_wrap(~variable, scales = "free_y", ncol = 2)


par(mfrow = c(3, 2))
# plot with correlation in title
plot(
  growth.obs~growth.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$growth.obs, comp_allspecies_dt$growth.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)
plot(reg.obs~reg.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$reg.obs, comp_allspecies_dt$reg.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)
plot(mort.obs~mort.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$mort.obs, comp_allspecies_dt$mort.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)
plot(ba.obs~ba.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$ba.obs, comp_allspecies_dt$ba.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)
plot(trees.obs~trees.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$trees.obs, comp_allspecies_dt$trees.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)
plot(dbh.obs~dbh.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0("Correlation: ", round(cor(comp_allspecies_dt$dbh.obs, comp_allspecies_dt$dbh.pred, use = "complete.obs", method = "spearman"), 2)))
abline(0, 1)

comp_allspecies_dt[,.(
  growth_cor = cor(growth.obs,growth.pred, use = "complete.obs", method = "spearman"),
  mort_cor = cor(mort.obs,mort.pred, use = "complete.obs", method = "spearman"),
  reg_cor = cor(reg.obs,reg.pred, use = "complete.obs", method = "spearman"),
  trees_cor = cor(trees.obs,trees.pred, use = "complete.obs", method = "spearman"),
  dbh_cor = cor(dbh.obs,dbh.pred, use = "complete.obs", method = "spearman"),
  ba_cor = cor(ba.obs,ba.pred, use = "complete.obs", method = "spearman")
),by=species]

