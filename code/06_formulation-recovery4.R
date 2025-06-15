# This script visualizes how FINN can recover complex formulations of processes
library(FINN)
library(data.table)
library(ggplot2)
library(torch)


## define default parameter settings and simulation setup for FINN ####

Ntimesteps = 150  # number of timesteps
Nsites = 10 # number of sites
Npatches = 20
patch_size = 0.1
Nsp = 1 # number of species
sp_id = 1

FINN.seed(1234)
shadeSP = 0.7

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  1, # intercept regulating the overall effect size
  0
),Nsp, 2))

# # growth parameters
# parGrowth = matrix(c(
#   shadeSP, # see above
#   0.06 # the second growth parameter modulates the size dependent growth
# ),Nsp, 2) 

# parGrowthEnv = list(matrix(c(
#   -0.5, # intercept regulating the overall effect size
#   1 # the second parameter modulates the effect of the environmental variable
# ),Nsp, 2))

# scale0 <- function(x){
#   if(length(x) == 1){
#     return(x)
#   } else {
#     return((scale(x)))
#   }
# }

# mortality parameters
parMort = matrix(c(
  0.5, # see above
  -1, # the second growth parameter modulates the size dependent mortality
  0.3 # the third mort parameter modulates the growth dependent mortality
),Nsp, 3)
parMortEnv = list(matrix(c(
  -4, # intercept regulating the overall effect size
  0# the second parameter modulates the effect of the environmental variable
), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
parComp = matrix(c(
  c(0.6), # parHeight
  0.5 # compStr
),Nsp, 2)

# we first generate a data.table with all combinations of site and timestep.
env_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)

dist_dt <- env_dt

# for this very simple model we will have a constant environment for all sites and timesteps
# env_dt$env1 = rep(0, Ntimesteps*Nsites)
env_dt[,env1 := as.numeric(scale(siteID)),]
# env_dt[,env1 := 0,]
env_dt[,env1 := env1+rnorm(.N,0,0.3),]
env_dt[,env1 := -env1,]
# hist(env_dt$env1)
disturbance_frequency = 0.0
disturbance_intensity = rbinom(Ntimesteps*Nsites,1,0.2)*runif(Ntimesteps*Nsites, 0.5, 1)
dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity

# Environmental effect ####

## models ####
###  M1 (true model) ####

parGrowth_m1env = cbind(
  c(0.1),# growthLight
  c(0.07)# growthSize
)


parGrowthEnv_m1env = list(matrix(c(
  0,#-0.5, # intercept regulating the overall effect size
  1 # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# parGrowthEnv_m1env = list(matrix(c(
#   -2.5, # intercept regulating the overall effect size
#   1 # the second parameter modulates the effect of the environmental variable
# ),Nsp, 2))
# 
# parGrowthEnv_m1env = list(matrix(c(
#   0.1, # intercept regulating the overall effect size
#   1 # the second parameter modulates the effect of the environmental variable
# ),Nsp, 2))


growth_m1 <- function (dbh, species, parGrowth, pred, light, light_steepness = 10, 
                       debug = F, trees = NULL) {
  shade = ((1/(1 + torch::torch_exp(-light_steepness * (light - 
                                                          parGrowth[, 1][species]))) - 1/(1 + torch::torch_exp(light_steepness * 
                                                                                                                 parGrowth[, 1][species])))/(1/(1 + torch::torch_exp(-light_steepness * 
                                                                                                                                                                       (1 - parGrowth[, 1][species]))) - 1/(1 + torch::torch_exp(light_steepness * 
                                                                                                                                                                                                                                   parGrowth[, 1][species]))))
  environment = torch::torch_exp(-pred^2)
  growth = shade * environment * (torch::torch_exp(-parGrowth[, 
                                                              2][species] * dbh))
  if (debug == TRUE) 
    out = list(shade = shade, light = light, environment = environment, 
               growth = growth)
  else out = growth
  return(out)
}
growth_m1env = createProcess(~1+env1, initEnv = parGrowthEnv_m1env,initSpecies = parGrowth_m1env, func = growth_m1)
reg_m1env = createProcess(~1, initEnv = list(matrix(parRegEnv[[1]][,1])),initSpecies = parReg, func = FINN::regeneration)
mort_m1env = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality)

m1env = finn(N_species = Nsp, 
             competition_process = createProcess(~0, func = FINN::competition),
             growth_process = growth_m1env,
             regeneration_process = reg_m1env,
             mortality_process = mort_m1env
)

###  M2 (wrong process model) ####

growth_m2env = createProcess(~1+env1, initEnv = parGrowthEnv_m1env,initSpecies = parGrowth_m1env, func = FINN::growth)

m2env = finn(N_species = Nsp, 
             competition_process = createProcess(~0, func = FINN::competition),
             growth_process = growth_m2env,
             regeneration_process = reg_m1env,
             mortality_process = mort_m1env
)

###  M3 (hybrid model) ####

growth_m3env = createHybrid(~., transformer = FALSE, dropout = 0.2)

m3env = finn(
  N_species = Nsp,
  competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = F, optimizeEnv = F),
  growth_process = growth_m3env,
  regeneration_process = reg_m1env,
  mortality_process = mort_m1env
)

gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
  g = (self$nn_growth(dbh = dbh, light = light, species = species, env = pred) - exp(1))$exp()
  return(g)
}

m3env$growth_func = m3env$.__enclos_env__$private$set_environment(gh)
source("code/99_cohort500fix2.R")


##  simulate from M1 ####

preds_env = 
  m1env$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = Npatches, debug = T)

### plot timeseries ####
p_dat <- preds_env$long$site[, .(value = mean(value)), by = .(year, species, variable, siteID)]
p_dat[, variable2 := factor(
  variable,
  levels = c("dbh", "ba", "trees", "AL", "growth", "mort", "r_mean_ha"),
  labels =  c("avg. DBH [cm]", "Basal Area [m²/ha]", "Trees [N/ha]",
              "Available Light [%]", "Growth [%/100]", "Mortality [%/100]",
              "Reg. Mean [N/ha]")
),]
p_dat <- p_dat[!is.na(variable2)]

# p_dat[variable %in% c("ba", "trees", "reg"), value := value/patch_size,]
p_dat[variable %in% c("ba", "trees"), value := value/patch_size,]
p <- ggplot(p_dat, aes(x = year, y = value, color = factor(siteID, ordered = T))) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Value") +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left") +  # Remove label_parsed
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    strip.background = element_blank(),
    legend.position = "bottom"
  )
# scale_color_gradient(low = "blue", mid = "orange", "red")
p

### plot simulated response pattern ####
p_dat2 <- preds_env$wide$site
p_dat2 <- merge(env_dt, p_dat2, by = c("year", "siteID"))
dbhcm_2_bam2 <- function(x) {
  x <- x/100
  return(x^2 * pi)
}
plot(p_dat2$env1, dbhcm_2_bam2(p_dat2$growth*p_dat2$dbh))

##  calibrate models ####

obs_dt <- preds_env$wide$site
cohorts_dt <- data.table(preds_env$wide$cohort)
# cohorts_dt <- predictions[["patches_100"]]$wide$cohort

obs_dt <- obs_dt[year >= 125 & year <= 150]
env_dt_calib <- env_dt[year >= 125 & year <= 150]

cohorts_dt <- cohorts_dt[dbh != 0]

init_cohort <- FINN::CohortMat$new(cohorts_dt[year == 124 & dbh >= 1], sp = 1)

obs_dt[dbh == 0, mort := NA_real_]
obs_dt[dbh == 0, growth := NA_real_]
obs_dt[dbh == 0, dbh := NA_real_]
obs_dt[,reg := reg/0.1,]
obs_dt[,year := as.integer(as.factor(year)),]
env_dt_calib[,year := as.integer(as.factor(year)),]

obs_dt_p <- merge(obs_dt, env_dt, by = c('year', "siteID"))
plot(obs_dt_p$env1, obs_dt_p$growth*obs_dt_p$dbh, xlab = "env1", ylab = "growth")
ggplot(obs_dt_p, aes(x = env1, y = growth*dbh)) +
  geom_point() +
  geom_smooth()+
  labs(x = "env1", y = "growth") +
  theme_classic()

# -1 oder 0 + 

m3env$fit(data = obs_dt, batchsize = 10, env = env_dt_calib, init_cohort = init_cohort,  epochs = 1000, patches = Npatches, lr = 0.001,
          optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
          weights = c(0.1, 10, 1.0, 10.0, 1, 1),
          loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
)

m2env$fit(data = obs_dt, batchsize = 10, env = env_dt_calib, init_cohort = init_cohort,  epochs = 1000, patches = Npatches, lr = 0.001,
          optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
          weights = c(0.1, 10, 1.0, 10.0, 1, 1),
          loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
)



## plot simulations from calibrated models ####

pred_m1env <- m1env$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 20, debug = F)
pred_m2env <- m2env$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 20, debug = F)
pred_m3env <- m3env$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 20, debug = F)


pred_m1env_dt <- pred_m1env$wide$site[dbh >= 1]
pred_m2env_dt <- pred_m2env$wide$site[dbh >= 1]
pred_m3env_dt <- pred_m3env$wide$site[dbh >= 1]
# cohorts_dt <- predictions[["patches_100"]]$wide$cohort
pred_m1env_dt <- merge(pred_m1env_dt, env_dt, by = c('year', "siteID"))
pred_m2env_dt <- merge(pred_m2env_dt, env_dt, by = c('year', "siteID"))
pred_m3env_dt <- merge(pred_m3env_dt, env_dt, by = c('year', "siteID"))

par(mfrow = c(2,2))
plot(pred_m1env_dt$env1, pred_m1env_dt$growth*pred_m1env_dt$dbh)
plot(pred_m2env_dt$env1, pred_m2env_dt$growth*pred_m2env_dt$dbh)
plot(pred_m3env_dt$env1, pred_m3env_dt$growth*pred_m3env_dt$dbh)


m3env$fit(data = obs_dt, batchsize = 10, env = env_dt_calib, init_cohort = init_cohort,  epochs = 8000, patches = Npatches, lr = 0.001,
          optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
          weights = c(0.1, 10, 1.0, 10.0, 1, 1),
          loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
)

m2env$fit(data = obs_dt, batchsize = 10, env = env_dt_calib, init_cohort = init_cohort,  epochs = 8000, patches = Npatches, lr = 0.001,
          optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
          weights = c(0.1, 10, 1.0, 10.0, 1, 1),
          loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
)



# Recover complex environmental effect ####
## Recover complex environmental effect ####

