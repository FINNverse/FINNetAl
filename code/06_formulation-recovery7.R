# This script visualizes how FINN can recover complex formulations of processes
library(FINN)
library(data.table)
library(ggplot2)
library(torch)

#' Draw a GAM smooth (±SE) without plotting points
#'
#' @param x,y        numeric vectors (explanatory, response) of equal length.
#' @param n          grid size for predictions (controls curve resolution).
#' @param se.mult    SE multiplier (≈2 for 95 % band).
#' @param add        logical; overlay on current plot (`TRUE`) or start a new
#'                   blank frame (`FALSE`, default).
#' @param formula    GAM formula (default = y ~ s(x, bs = "cs")).
#' @param family     glm family; if NULL and all(y > 0) defaults to
#'                   `Gamma(link = "log")`, otherwise `gaussian()`.
#' @param gam_args   named list of extra arguments for `mgcv::gam()` (optional).
#' @param col,lwd    colour / line-width for the central smooth.
#' @param se.col,se.lty  colour / lty for the ±SE curves.
#' @param ...        further arguments passed **only** to `plot()` when
#'                   `add = FALSE` (axis labels, limits, …).
#'
#' @return invisibly a data.frame of x, fit, upr, lwr.
#'
#' @examples
#' set.seed(1)
#' x <- runif(120, 0, 10)
#' y <- exp(0.5 * sin(x)) + rlnorm(120, 0, 0.1)
#'
#' ## (1) default smooth chosen for positive data
#' lines_gam(x, y)
#'
#' ## (2) a *smoother* curve with REML + higher gamma
#' plot(x, y, pch = 19, col = "#00000044")
#' lines_gam(
#'   x, y,
#'   add       = TRUE,
#'   formula   = y ~ s(x, k = 30),     # more basis functions
#'   gam_args  = list(method = "REML", gamma = 1.6),
#'   col       = "tomato"
#' )
lines_gam <- function(x, y,
                      n        = 100,
                      se.mult  = 2,
                      add      = FALSE,
                      formula  = NULL,
                      family   = NULL,
                      gam_args = list(),
                      col      = "#6c9fca",
                      lwd      = 3,
                      se.col   = col,
                      se.lty   = 2,
                      se = F,
                      ...) {
  
  ## --------- sanity checks & setup -----------------------------------------
  stopifnot(length(x) == length(y))
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required; install it with install.packages('mgcv').")
  
  dat <- data.frame(x = x, y = y)
  
  ## choose family if none supplied
  if (is.null(family))
    family <- if (all(y > 0)) Gamma(link = "log") else gaussian()
  
  ## default formula if none supplied
  if (is.null(formula))
    formula <- y ~ s(x, bs = "cs")
  
  ## --------- fit model (extra args allowed) ---------------------------------
  gam_call <- c(
    list(formula = formula, data = dat, family = family),
    gam_args
  )
  mod <- do.call(mgcv::gam, gam_call)
  
  ## --------- predictions on link scale -------------------------------------
  xs  <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n)
  pr  <- predict(mod, newdata = data.frame(x = xs), se.fit = TRUE)
  
  fit.lin <- pr$fit
  se.lin  <- pr$se.fit
  
  ## back-transform to response scale
  linkinv <- family$linkinv
  fit <- linkinv(fit.lin)
  upr <- linkinv(fit.lin + se.mult * se.lin)
  lwr <- linkinv(fit.lin - se.mult * se.lin)
  
  ## --------- draw -----------------------------------------------------------
  if (!add) {
    plot(x, y, type = "n", ...)       # blank plotting region
  } else if (is.null(dev.list())) {
    stop("No active graphics device – call plot() first or set add = FALSE.")
  }
  
  lines(xs, fit, col = col,      lwd = lwd)
  if(se == T){
    lines(xs, upr, col = se.col,   lty = se.lty)
    lines(xs, lwr, col = se.col,   lty = se.lty)
  }
  
  invisible(data.frame(x = xs, fit = fit, upr = upr, lwr = lwr))
}
## define default parameter settings and simulation setup for FINN ####

Ntimesteps = 20  # number of timesteps
Ntimesteps_nodebug = 500  # number of timesteps
Nsites = 20 # number of sites
Npatches = 10
patch_size = 0.1
Nsp = 1 # number of species
sp_id = 1

FINN.seed(1235)
shadeSP = 0.9

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  2, # intercept regulating the overall effect size
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
  0.3, # see above
  -3, # the second growth parameter modulates the size dependent mortality
  0.5 # the third mort parameter modulates the growth dependent mortality
),Nsp, 3)
parMortEnv = list(matrix(c(
  -3, # intercept regulating the overall effect size
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
      year = 1:Ntimesteps_nodebug
    )
  )
)

dist_dt <- env_dt

# for this very simple model we will have a constant environment for all sites and timesteps
# env_dt$env1 = rep(0, Ntimesteps*Nsites)
env_dt[,env1 := as.numeric(scale(siteID))*1,]
# hist(env_dt$env1)
# env_dt[,env1 := 0,]
env_dt[,env1 := env1+rnorm(.N,0,0.2),]
# hist(env_dt$env1)
env_dt[,env1 := -env1,]
# hist(env_dt$env1)
disturbance_frequency = 0.0
disturbance_intensity = rbinom(Ntimesteps_nodebug*Nsites,1,0.2)*runif(Ntimesteps_nodebug*Nsites, 0.5, 1)
dist_dt$intensity = rbinom(Ntimesteps_nodebug*Nsites, 1, disturbance_frequency)*disturbance_intensity

# Environmental effect ####

## models ####
###  M1 (true model) ####

parGrowth_m1env = cbind(
  c(0.1),# growthLight
  c(0.03)# growthSize
)


parGrowthEnv_m1env = list(matrix(c(
  0,#-0.5, # intercept regulating the overall effect size
  1.3 # the second parameter modulates the effect of the environmental variable
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

competition2 = function (dbh, species, trees, parComp, h = NULL, patch_size_ha, 
                         ba = NULL, cohortHeights = NULL) 
{
  parHeight = parComp[, 1]
  parCompStr = parComp[, 2]
  if (is.null(ba)) 
    ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha) * 
    parCompStr[species] * 0.1
  if (is.null(cohortHeights)) 
    cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
  if (is.null(h)) {
    h = cohortHeights
    ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - 
                                               h$permute(c(1, 2, 4, 3)) + 0.01)/0.01)$sigmoid_()))$sum(-2)
  }
  else {
    ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 
                                                0.1)/0.01)$sigmoid_()))$sum(-2)
  }
  light = 1 - ba_height
  light = torch_clamp(light, min = 0)
  return(light)
}

comp_m1env = createProcess(~0, func = FINN::competition, optimizeSpecies = FALSE, optimizeEnv = FALSE)
growth_m1env = createProcess(~1+env1, initEnv = parGrowthEnv_m1env,initSpecies = parGrowth_m1env, func = growth_m1)
reg_m1env = createProcess(~1, initEnv = list(matrix(parRegEnv[[1]][,1])),initSpecies = parReg, func = FINN::regeneration, sample_regeneration = F, optimizeSpecies = FALSE, optimizeEnv = FALSE)
mort_m1env = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality, optimizeSpecies = FALSE, optimizeEnv = FALSE)

if(!dir.exists("results/06_formrecov")){
  dir.create("results/06_formrecov", recursive = T)
}

if(file.exists("results/06_formrecov/C_m1env.pt")){
  m1env = torch::torch_load("results/06_formrecov/C_m1env.pt")
  }else{
    m1env = finn(N_species = Nsp, 
                 competition_process = comp_m1env,
                 growth_process = growth_m1env,
                 regeneration_process = reg_m1env,
                 mortality_process = mort_m1env
    )
  torch::torch_save(m1env,"results/06_formrecov/C_m1env.pt")
  }

# -1 oder 0 + 


###  M2 (wrong process model) ####

growth_m2env = createProcess(~1+env1, initEnv = parGrowthEnv_m1env,initSpecies = parGrowth_m1env, func = FINN::growth, optimizeSpecies = T, optimizeEnv = T)

m2env = finn(N_species = Nsp, 
             competition_process = comp_m1env,
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

env_dt_nodebug <- env_dt
env_dt <- env_dt[year <= Ntimesteps]
dist_dt_nodebug <- dist_dt
dist_dt <- dist_dt[year <= Ntimesteps]
preds_env = 
  m1env$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = Npatches, debug = T)

preds_env_nodebug = 
  m1env$simulate(init_cohort = NULL, env = env_dt_nodebug, disturbance = dist_dt_nodebug, device = "cpu", patches = Npatches, debug = F)

### plot timeseries ####
p_dat <- preds_env_nodebug$long$site[, .(value = mean(value)), by = .(year, species, variable, siteID)]
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

p_dat <- merge(p_dat, env_dt[, .(avg_env = mean(env1)), by = .(siteID)], by = "siteID")

p <- ggplot(p_dat, aes(x = year, y = value, group = siteID, color = avg_env)) +
  geom_line(alpha = 0.5) +
  labs(x = "Year", y = "Value") +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left") +  # Remove label_parsed
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    strip.background = element_blank(),panel.background = element_rect(fill = "grey70"),
    legend.position = "bottom"
  )+
  scale_color_gradient2(name = "environment", low = "darkblue", mid = "white", high = "darkred")

pdf("figures/C_timeseries_m1env.pdf", width = 7, height = 7)
p
dev.off()

### plot simulated response pattern ####
p_dat2 <- preds_env$wide$site
p_dat2 <- merge(env_dt, p_dat2, by = c("year", "siteID"))
dbhcm_2_bam2 <- function(x) {
  x <- x/100
  return(x^2 * pi)
}
plot(p_dat2$env1, dbhcm_2_bam2(p_dat2$growth*p_dat2$dbh))
plot(p_dat2$env1, p_dat2$growth*p_dat2$dbh)

ggplot(p_dat2, aes(x = env1, y = growth)) +
  geom_point() +
  geom_smooth()+
  facet_grid(~cut(year, 5))+
  labs(x = "env1", y = "growth") +
  theme_classic()

##  calibrate models ####

obs_dt <- preds_env$wide$site
cohorts_dt <- data.table(preds_env$wide$cohort)
# cohorts_dt <- predictions[["patches_100"]]$wide$cohort

obs_dt <- obs_dt[year >= 10 & year <= 20]
env_dt_calib <- env_dt[year >= 10 & year <= 20]

cohorts_dt <- cohorts_dt[dbh != 0]
plot(V1~year,cohorts_dt[patchID == 1,uniqueN(cohortID), by = .(year, siteID, patchID)])

cohorts_dt[year == 1 & patchID == 1]
cohorts_dt[,uniqueN(cohortID), by = .(year)]
cohorts_dt <- cohorts_dt[year == 9 & dbh >= 1]
summary(cohorts_dt$dbh)
cohorts_dt2 <- cohorts_dt[,.(
  trees = sum(round(trees))
  # ), by = .(siteID, patchID, species, dbh = cut(dbh, seq(1,max(cohorts_dt$dbh)+.1,0.1)))]
), by = .(siteID, patchID, species, dbh = cut(dbh, seq(1,max(cohorts_dt$dbh, na.rm = T)+.1,0.1), labels = seq(1,max(cohorts_dt$dbh),0.1), include.lowest=T))]
cohorts_dt2[,dbh := as.numeric(as.character(dbh)),]
cohorts_dt2[,cohortID := 1:.N,]
init_cohort <- FINN::CohortMat$new(cohorts_dt2, sp = 1)

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


if(file.exists("results/06_formrecov/C_m3env1000.pt")){
  m3env = torch::torch_load("results/06_formrecov/C_m3env1000.pt")
}else{
  m3env$fit(data = obs_dt, batchsize = 20, env = env_dt_calib, init_cohort = init_cohort,  epochs = 500, patches = Npatches, lr = 0.001,
            optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
            weights = c(0.1, 10, 1.0, 10.0, 1, 1),
            loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
  )
  torch::torch_save(m3env,"results/06_formrecov/C_m3env1000.pt")
}

if(file.exists("results/06_formrecov/C_m2env1000.pt")){
  m2env = torch::torch_load("results/06_formrecov/C_m2env1000.pt")
  }else{
    m2env$fit(data = obs_dt, batchsize = 5, env = env_dt_calib, init_cohort = init_cohort,  epochs = 500, patches = Npatches, lr = 0.001,
              optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
              weights = c(0.1, 10, 1.0, 10.0, 1, 1),
              loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
    )
    torch::torch_save(m2env,"results/06_formrecov/C_m2env1000.pt")
  }
## plot simulations from calibrated models ####

pred_m1env <- m1env$simulate(init_cohort = init_cohort, env = env_dt_calib, disturbance = dist_dt, device = "cpu", patches = Npatches, debug = T)
pred_m2env <- m2env$simulate(init_cohort = init_cohort, env = env_dt_calib, disturbance = dist_dt, device = "cpu", patches = Npatches, debug = T)
pred_m3env <- m3env$simulate(init_cohort = init_cohort, env = env_dt_calib, disturbance = dist_dt, device = "cpu", patches = Npatches, debug = T)

# pred_m1env_dt <- pred_m1env$wide$site[dbh >= 1]
# pred_m2env_dt <- pred_m2env$wide$site[dbh >= 1]
# pred_m3env_dt <- pred_m3env$wide$site[dbh >= 1]
# pred_m1env_dt <- merge(pred_m1env_dt, env_dt, by = c('year', "siteID"))
pred_m1env_dt <- pred_m1env$wide$cohort[dbh >= 1 & !is.na(g)][,growth := g*0.5,]
pred_m2env_dt <- pred_m2env$wide$cohort[dbh >= 1 & !is.na(g)][,growth := g*0.5,]
pred_m3env_dt <- pred_m3env$wide$cohort[dbh >= 1 & !is.na(g)][,growth := g*0.5,]
pred_m1env_dt <- merge(pred_m1env_dt, env_dt, by = c('year', "siteID"))
pred_m2env_dt <- merge(pred_m2env_dt, env_dt, by = c('year', "siteID"))
pred_m3env_dt <- merge(pred_m3env_dt, env_dt, by = c('year', "siteID"))

pdf("figures/Cform-recovery-abcd.pdf", width = 6, height = 6)
par(mfrow = c(1,1), mar = c(4,4,1,1), pty="s", las = 1)
formula = y ~ s(x,bs="cr",  k= 5)
# plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment")
# plot with rotaed y lab and as perfect square
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment")
# lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
# legend("topleft", legend = c("M1 (true form)"), col = c("black"), lty = 1, bty = "n", lwd = 3)
# legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
# dev.off()

# pdf("figures/Cform-recovery-b.pdf", width = 6, height = 6)
# par(mfrow = c(1,1))
formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
# legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
legend("topleft", legend = c("M1 (true form)"), col = c("black"), lty = 1, bty = "n", lwd = 3)

formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = NA, ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
# legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
legend("topleft", legend = c("M1 (true form)"), col = c("black"), lty = 1, bty = "n", lwd = 3)
formula = y ~ s(x,bs="cr",  k= 5)

plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = alpha("black",0.01), ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
# legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
legend("topleft", legend = c("M1 (true form)"), col = c("black"), lty = 1, bty = "n", lwd = 3)
# dev.off()

# pdf("figures/Cform-recovery-c.pdf", width = 6, height = 6)
# par(mfrow = c(1,1))
formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)"), col = c("black", "#fe944e"), lty = 1, bty = "n", lwd = 3)
# dev.off()
formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = NA, ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)"), col = c("black", "#fe944e"), lty = 1, bty = "n", lwd = 3)

formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = alpha("black",0.01), ylab = "dbh growth", xlab = "environment")
points(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = alpha("#fe944e",0.01))
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
# lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)"), col = c("black", "#fe944e"), lty = 1, bty = "n", lwd = 3)
# dev.off()

# pdf("figures/Cform-recovery-d.pdf", width = 6, height = 6)
# par(mfrow = c(1,1))
formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)

formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = NA, ylab = "dbh growth", xlab = "environment")
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)

formula = y ~ s(x,bs="cr",  k= 5)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = alpha("black",0.01), ylab = "dbh growth", xlab = "environment")
points(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = alpha("#fe944e",0.01))
points(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = alpha("#6c9fca",0.01))
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
dev.off()

pdf("figures/Cform-recovery_all.pdf", width = 10, height = 10)
par(mfrow = c(2,2), mar = c(4,4,2,1), pty="s", las = 1)
formula = y ~ s(x, k = 5)
ylim = c(0,8.7)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment", ylim = ylim)
mtext("a", side = 3, line = 0.5, adj = -0.3, font = 2)  # 'a' at top-left
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)"), col = c("black"), lty = 1, bty = "n", lwd = 3)
plot(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment", ylim = ylim)
mtext("b", side = 3, line = 0.1, adj = -0.3, font = 2)  # 'a' at top-left
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M2 (wrong form)"), col = c("#fe944e"), lty = 1, bty = "n", lwd = 3)
plot(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "grey70", ylab = "dbh growth", xlab = "environment", ylim = ylim)
mtext("c", side = 3, line = 0.1, adj = -0.3, font = 2)  # 'a' at top-left
lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M3 (hybrid model)"), col = c("#6c9fca"), lty = 1, bty = "n", lwd = 3)
plot(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = NA, ylab = "dbh growth", xlab = "environment", ylim = ylim)
mtext("d", side = 3, line = 0.1, adj = -0.3, font = 2)  # 'a' at top-left
lines_gam(pred_m1env_dt$env1, (pred_m1env_dt$growth*pred_m1env_dt$dbh), col = "black", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m2env_dt$env1, (pred_m2env_dt$growth*pred_m2env_dt$dbh), col = "#fe944e", add = T, formula = formula, ylim = c(0,3))
lines_gam(pred_m3env_dt$env1, (pred_m3env_dt$growth*pred_m3env_dt$dbh), col = "#6c9fca", add = T, formula = formula, ylim = c(0,3))
legend("topleft", legend = c("M1 (true form)", "M2 (wrong form)", "M3 (hybrid model)"), col = c("black", "#fe944e", "#6c9fca"), lty = 1, bty = "n", lwd = 3)
dev.off()

# m3env$fit(data = obs_dt, batchsize = 5, env = env_dt_calib, init_cohort = init_cohort,  epochs = 8000, patches = Npatches, lr = 0.001,
#           optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
#           weights = c(0.1, 10, 1.0, 10.0, 1, 1),
#           loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
# )
# torch::torch_save(m3env,"results/06_formrecov/C_m3env1000.pt")
# 
# m2env$fit(data = obs_dt, batchsize = 10, env = env_dt_calib, init_cohort = init_cohort,  epochs = 8000, patches = Npatches, lr = 0.001,
#           optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
#           weights = c(0.1, 10, 1.0, 10.0, 1, 1),
#           loss= c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom")
# )
# torch::torch_save(m2env,"results/06_formrecov/C_m2env1000.pt")



# Recover complex environmental effect ####
## Recover complex environmental effect ####
