library(FINN)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(data.table)

# Function to extract the legend as a grob
get_legend <- function(ggplot_obj) {
  # Build the gtable
  plot_gtable <- ggplot_gtable(ggplot_build(ggplot_obj))
  
  # Find the index of the guide box (the legend)
  guide_index <- which(sapply(plot_gtable$grobs, function(x) x$name) == "guide-box")
  
  # Return the legend grob
  if (length(guide_index) > 0) {
    return(plot_gtable$grobs[[guide_index]])
  } else {
    warning("No legend found.")
    return(NULL)
  }
}

label_grob_f = function(x){
  textGrob(x, x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
}

Ntimesteps = 100  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
Nsp = 1 # number of species
sp_id = 1

FINN.seed(1234)
shadeSP = 0.2

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  1, # intercept regulating the overall effect size
  0
),Nsp, 2))

# growth parameters
parGrowth = matrix(c(
  shadeSP, # see above
  0.06 # the second growth parameter modulates the size dependent growth
),Nsp, 2) 

parGrowthEnv = list(matrix(c(
  1, # intercept regulating the overall effect size
  0 # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

scale0 <- function(x){
  if(length(x) == 1){
    return(x)
  } else {
    return((scale(x)))
  }
}

# mortality parameters
parMort = matrix(c(
  0.2, # see above
  -0.7, # the second growth parameter modulates the size dependent mortality
  0 # the third mort parameter modulates the growth dependent mortality
),Nsp, 3)
parMortEnv = list(matrix(c(
  -2.5, # intercept regulating the overall effect size
  -2.5 # the second parameter modulates the effect of the environmental variable
), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
parComp = matrix(c(
  c(0.6), # parHeight
  0.5 # compStr
),Nsp, 2)

# Create a wide-format data.table with one row per species
pars_dt <- data.table(
  speciesID   = 1:Nsp,
  reg         = parReg,
  growth1      = parGrowth[, 1],
  growth2      = parGrowth[, 2],
  mort1        = parMort[, 1],
  mort2        = parMort[, 2],
  mort3        = parMort[, 3],
  compHeight  = parComp[, 1],
  compStrength= parComp[, 2],
  regEnv1     = sapply(parRegEnv, function(x) x[1, 1]),
  regEnv2     = sapply(parRegEnv, function(x) x[1, 2]),
  growthEnv1  = sapply(parGrowthEnv, function(x) x[1, 1]),
  growthEnv2  = sapply(parGrowthEnv, function(x) x[1, 2]),
  mortEnv1    = sapply(parMortEnv, function(x) x[1, 1]),
  mortEnv2    = sapply(parMortEnv, function(x) x[1, 2])
)
pars_dt

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
env_dt$env1 = rep(0, Ntimesteps)
disturbance_frequency = 0.0
disturbance_intensity = rbinom(Ntimesteps*Nsites,1,0.2)*runif(Ntimesteps*Nsites, 0.5, 1)
dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity

predictions <- list()
simulationModel = finn(N_species = Nsp, 
                       competition_process = createProcess(~0, func = FINN::competition),
                       growth_process = createProcess(~1+env1, initEnv = parGrowthEnv,initSpecies = parGrowth, func = FINN::growth),
                       regeneration_process = createProcess(~1+env1, initEnv = parRegEnv,initSpecies = parReg, func = FINN::regeneration),
                       mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality),
)

predictions[["patches_1"]] = 
  simulationModel$simulate(init_cohort = NULL, env = env_dt, disturbance= dist_dt, device = "cpu", patches = 1)

simulationModel = finn(N_species = Nsp, 
                       competition_process = createProcess(~0, func = FINN::competition),
                       growth_process = createProcess(~1+env1, initEnv = parGrowthEnv,initSpecies = parGrowth, func = FINN::growth),
                       regeneration_process = createProcess(~1+env1, initEnv = parRegEnv,initSpecies = parReg, func = FINN::regeneration),
                       mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality)
)

predictions[["patches_100"]] = 
  simulationModel$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 20, debug = T)

p_list <- list()
for(i in c("patches_1", "patches_100")){
  p_dat <- predictions[[i]]$long$site[, .(value = mean(value)), by = .(year, species, variable)]
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
  p <- ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
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
    ) +
    scale_color_discrete(name = fifelse(grepl("pred_", i),"PFT","Random Species"))
  legend_p <- get_legend(p)
  p_list[[i]] <- p+theme(legend.position = "none")
  # }
}

grid.arrange(
  grobs = c(
    list(label_grob_f("a"), label_grob_f("b")),
    p_list[1:2], list(legend_p)),
  ncol = 2,
  layout_matrix = rbind(c(1, 2), c(3, 4), c(5,5)),
  heights = c(0.05, 0.95,0.1)
)


obs_dt <- predictions[["patches_100"]]$wide$site
cohorts_dt <- predictions[["patches_100"]]$wide$cohort

obs_dt <- obs_dt[year >= 80]
env_dt_calib <- env_dt[year >= 80]
init_cohort <- FINN::CohortMat$new(cohorts_dt[year == 79,], sp = 1)


simulationModel2 = finn(N_species = Nsp, 
                        competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = F, optimizeEnv = F),
                        growth_process = createProcess(~1+env1, initEnv = NULL,initSpecies = NULL, func = FINN::growth),
                        regeneration_process = createProcess(~1+env1, initEnv = parRegEnv,initSpecies = parReg, func = FINN::regeneration, optimizeSpecies = F, optimizeEnv = F),
                        mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality, optimizeSpecies = F, optimizeEnv = F),
)

simulationModel2$fit(data = obs_dt, batchsize = 1, env = env_dt_calib, init_cohort = init_cohort,  epochs = 500, patches = 100, lr = 0.01,
                     optimizer = torch::optim_adam, device = "cpu", plot_progress = T,
                     loss= c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "mse", regeneration = "nbinom")
)

