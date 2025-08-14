library(data.table)
library(FINN)
library(torch)
library(glmmTMB)
library(parallel)
library(iml)

# load data
i_dir = "data/BCI/CVsplits-realdata/pft-period7-25patches/"
cv_S = "S0"
cv_T = "T0"
env_dt = fread(paste0(i_dir,"/env_dt_",cv_S,"_",cv_T,"_train.csv"))
env_dt = env_dt[,.(siteID, year, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp )]
obs_dt = fread(paste0(i_dir,"/obs_dt_",cv_S,"_",cv_T,"_train.csv"))
obs_dt = obs_dt[,.(siteID, species, year, dbh, ba, trees, growth, mort, reg)]
cohorts_dt = fread(paste0(i_dir,"/initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))
cohorts_dt = cohorts_dt[,.(siteID, dbh, species, patchID, census, trees, cohortID)]
obs_dt[dbh == 0, mort := NA_real_]
obs_dt[dbh == 0, growth := NA_real_]
obs_dt[dbh == 0, dbh := NA_real_]
Nspecies = max(obs_dt$species)
Nenv = ncol(env_dt) - 2
cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)

# load models
model_hybrid = torch::torch_load("results/02_realdata_hybrid_mortTF0/pft-period7-25patches_S0_T0.pt")
m_process = torch::torch_load("results/02_realdata/pft-period7-25patches_S0_T0_ba.trees.dbh.growth.mort.reg.pt")

raw_debug = model_hybrid$simulate(env = env_dt, init_cohort = cohort1, patches = 25L, debug = F, device = "cpu")
raw_debug_process = m_process$simulate(env = env_dt, init_cohort = cohort1, patches = 25L, debug = F, device = "cpu")


################ ALE Hybrid und process #################

model_hybrid = torch::torch_load("results/02_realdata_hybridTF0/pft-period7-25patches_S0_T0.pt")
model_hybrid$raw_g = NULL
model_hybrid$environment = NULL
model_hybrid$eval()
gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
  
  self$nn_growth$train()
  g = (self$nn_growth(dbh = dbh, light = light, trees = trees, species = species, env = pred) - exp(1))$exp()
  #print(g$var()$item())
  self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4), 
                                                                   light$unsqueeze(4), 
                                                                   trees$unsqueeze(4),
                                                                   species$unsqueeze(4)$float(),
                                                                   g$unsqueeze(4)
  ), dim = 4)) ))
  
  self$environment = c(self$environment, list(as.matrix(pred)))
  return(g)
}
unlockBinding(sym = "growth_func", model_hybrid$.__enclos_env__$self)
model_hybrid$growth_func = model_hybrid$.__enclos_env__$private$set_environment(gh)
cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)
model_hybrid$eval()
raw_debug = model_hybrid$simulate(env = env_dt, init_cohort = cohort1, patches = 25L, debug = F, device = "cpu")
site = 1
patch = 1
time = 1
df = data.frame()

for(time in 1:6) {
  for(site in 1:20) {
    for(patch in 1:25) {
      tmp = model_hybrid$raw_g[[time]][site, patch,,]
      tmp = cbind(tmp, matrix(c(1, as.matrix(env_dt[siteID==site & year == time])[1,-(1:2)]), nrow = nrow(tmp), ncol = 7, byrow = T))
      tmp = data.frame(tmp)
      colnames(tmp) = c("dbh", "light", "trees", "species", "growth", "intercept", colnames(env_dt)[3:8])
      tmp$time = time
      tmp = tmp[tmp$trees > 0.5,]
      df = rbind(df, tmp)
    }
  }
}
df = as.data.table(df)
df = df[df$dbh > 0.5,]

df_results_ale = data.frame()
for(i in 1:5) {
  for(J in 1:100) {
    model_hybrid$train()
    model_hybrid$nn_growth$train()
    SPECIES = i
    df_sp1 = df[df$species == SPECIES,] %>% dplyr::select(-species)
    indices = sample.int(nrow(df_sp1), 20000)
    predict_function = function(model, newdata) {
      
      # inputs   dbh     light trees species Intercept Prec   SR_kW_m2   RH_prc    T_max    T_min         swp 
      sites = nrow(newdata)
      dbh = torch_tensor(array(newdata$dbh, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_hybrid$device)
      trees = torch_tensor(array(newdata$trees, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_hybrid$device)
      light = torch_tensor(array(newdata$light, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_hybrid$device)
      species = torch_tensor(array(SPECIES, dim = c(sites, 1, 1)), dtype = torch_long(), device = model_hybrid$device)
      pred = torch_tensor(as.matrix(newdata |> dplyr::select(intercept, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp)), device = model_hybrid$device)
      
      g = (model_hybrid$nn_growth(dbh = dbh, light = light,trees=trees, species = species, env = pred) - exp(1))$exp()
      return((g$squeeze() %>% as.matrix())[,1])
    }
    
    
    set.seed(1)
    expl = DALEX::explain(model_hybrid, data= (df_sp1 %>% dplyr::select(-time, -growth))[indices,],y = (df_sp1 %>% dplyr::select(-time))[indices,]$growth, predict_function =predict_function )
    ale <- DALEX::model_profile(expl,
                                #variables = "light",
                                method = "ale",
                                center = FALSE, grid_points = 10) 
    # model_parts(expl) |> plot()
    #plot(ale)
    df_tmp = ale$agr_profiles 
    df_tmp$species = SPECIES
    df_tmp$sample = J
    df_results_ale = rbind(df_results_ale, df_tmp)
  }
}


colnames(df_results_ale) = c("features", "label", "x", "y", "ids", "species", "sample")

#### Process
model_process = torch::torch_load("results/02_realdata/pft-period7-25patches_S0_T0_ba.trees.dbh.growth.mort.reg.pt")
model_process$eval()
gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
  self$nn_growth$train()
  g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
  #print(g$var()$item())
  self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4), 
                                                                   light$unsqueeze(4), 
                                                                   trees$unsqueeze(4),
                                                                   species$unsqueeze(4)$float(),
                                                                   g$unsqueeze(4)), dim = 4)) ))
  self$environment = c(self$environment, list(as.matrix(pred)))
  return(g)
}

gh = function (dbh, species, parGrowth, pred, light, light_steepness = 10, 
               debug = F, trees = NULL) 
{
  shade = ((1/(1 + torch::torch_exp(-light_steepness * (light - 
                                                          parGrowth[, 1][species]))) - 1/(1 + torch::torch_exp(light_steepness * 
                                                                                                                 parGrowth[, 1][species])))/(1/(1 + torch::torch_exp(-light_steepness * 
                                                                                                                                                                       (1 - parGrowth[, 1][species]))) - 1/(1 + torch::torch_exp(light_steepness * 
                                                                                                                                                                                                                                   parGrowth[, 1][species]))))
  environment = torch::torch_exp(pred)
  growth = shade * environment * (torch::torch_exp(-parGrowth[, 
                                                              2][species] * dbh))
  if (debug == TRUE) 
    out = list(shade = shade, light = light, environment = environment, 
               growth = growth)
  else out = growth
  self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4), 
                                                                   light$unsqueeze(4), 
                                                                   trees$unsqueeze(4),
                                                                   species$unsqueeze(4)$float(),
                                                                   out$unsqueeze(4)), dim = 4)) ))
  return(out)
}

unlockBinding(sym = "growth_func", model_process$.__enclos_env__$self)
model_process$growth_func = model_process$.__enclos_env__$private$set_environment(gh)
cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = Nspecies)
raw_debug = model_process$simulate(env = env_dt, init_cohort = cohort1, patches = 25L, debug = F, device = "cpu")
site = 1
patch = 1
time = 1
df = data.frame()
for(time in 1:6) {
  for(site in 1:20) {
    for(patch in 1:25) {
      tmp = model_process$raw_g[[time]][site, patch,,]
      tmp = cbind(tmp, matrix(c(1, as.matrix(env_dt[siteID==site & year == time])[1,-(1:2)]), nrow = nrow(tmp), ncol = 7, byrow = T))
      tmp = data.frame(tmp)
      colnames(tmp) = c("dbh", "light", "trees", "species", "growth", "intercept", colnames(env_dt)[3:8])
      tmp$time = time
      tmp = tmp[tmp$trees > 0.5,]
      df = rbind(df, tmp)
    }
  }
}
model_process$raw_g[[1]]

df = as.data.table(df)
df = df[df$dbh > 0.5,]


df_results_ale_process = data.frame()
library(ingredients)
for(i in 1:5) {
  
  SPECIES = i
  df_sp1 = df[df$species == SPECIES,] %>% dplyr::select(-species)
  predict_function = function(model, newdata) {
    # inputs   dbh     light trees species Intercept Prec   SR_kW_m2   RH_prc    T_max    T_min         swp 
    sites = nrow(newdata)
    dbh = torch_tensor(array(newdata$dbh, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_process$device)
    trees = torch_tensor(array(newdata$trees, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_process$device)
    light = torch_tensor(array(newdata$light, dim = c(sites, 1, 1)), dtype = torch_float32(), device = model_process$device)
    species = torch_tensor(array(SPECIES, dim = c(sites, 1, 1)), dtype = torch_long(), device = model_process$device)
    pred = torch_tensor(as.matrix(newdata |> dplyr::select(intercept, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp)), device = model_process$device)
    
    pred = model_process$nn_growth(pred)
    pred = FINN::index_species(pred, species)
    g = model_process$growth_func(dbh = dbh, trees = trees, light = light, species = species, pred = pred, parGrowth = model_process$par_growth)
    
    #g = (model_process$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
    return(g$squeeze() %>% as.matrix())
  }
  
  
  
  expl = DALEX::explain(model_process, data= (df_sp1 %>% dplyr::select(-time, -growth))[indices,],y = (df_sp1 %>% dplyr::select(-time))[indices,]$growth, predict_function =predict_function )
  ale <- DALEX::model_profile(expl,
                              #variables = "light",
                              method = "ale",
                              center = FALSE, grid_points = 500) 
  #plot(ale)
  df_tmp = ale$agr_profiles 
  df_tmp$species = SPECIES
  df_results_ale_process = rbind(df_results_ale_process, df_tmp)
}

colnames(df_results_ale_process) = c("features", "label", "x", "y", "ids", "species")


save(df_results_ale, df_results_ale_process, file = "results/03_xAI/results_xAI.RData")
