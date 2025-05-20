library(FINN)
library(data.table)
library(ggplot2)

Ntimesteps = 500  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 5 # number of species

FINN.seed(1234)
# we draw the same shade parameters for each process for simplicity
# shade parameters correspond to the fraction of light a species needs to succesfully grow, regenerate, or survive.
shadeSP = c(0.1,0.2,0.5,0.5,0.7)

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  c(1,2,3,3,5), # intercept regulating the overall effect size
  # runif(Nsp, 1, 4), # intercept regulating the overall effect size
  runif(Nsp, -2, 2) # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# growth parameters
parGrowth = matrix(c(
  shadeSP, # see above
  c(0.04,0.05,0.05,0.06,0.1) # the second growth parameter modulates the size dependent growth
  # runif(Nsp, 0, 0.1) # the second growth parameter modulates the size dependent growth
),Nsp, 2) 

parGrowthEnv = list(matrix(c(
  c(0.2,0.3,0.5,1,1)*0.5, # intercept regulating the overall effect size
  # runif(Nsp, -2, -0.5), # intercept regulating the overall effect size
  runif(Nsp, -2, -0.5) # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# mortality parameters
parMort = matrix(c(
  as.numeric(scale(shadeSP)), # see above
  as.numeric(scale(parGrowth[,2])), # the second growth parameter modulates the size dependent mortality
  # runif(Nsp, -1, -0.1), # the second mort parameter modulates the size dependent mortality
  rep(0,Nsp) # the third mort parameter modulates the growth dependent mortality
  # runif(Nsp, -1, -.1) # the third mort parameter modulates the growth dependent mortality
),Nsp, 3)
parMortEnv = list(matrix(c(
  runif(Nsp, -3, -2), # intercept regulating the overall effect size
  runif(Nsp, -3, -2) # the second parameter modulates the effect of the environmental variable
), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
# parHeight = runif(Nsp, 0.3, 0.7)
# growth parameters
parComp = matrix(c(
  # runif(Nsp, 0.3, 0.7), # parHeight
  c(0.5,0.5,0.4,0.7,0.6), # parHeight
  # runif(Nsp, 0.2, 0.2) # Competition strength
  c(0.3,0.2,0.2,0.2,0.1) # parHeight
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

# we can also specify the intensity of disturbances for each timestep

# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_frequency = 0.05 

# the disturbance intensity at each timestep is the fraction of patches that is disturbed at that time step
# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_intensity = rbinom(Ntimesteps*Nsites,1,0.2)*runif(Ntimesteps*Nsites, 0.5, 1)

# this will result in 0 to 20 % of the patches being disturbed at each timestep with a change of 1% that a disturbance occurs at the timestep
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
                       mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality),
)

predictions[["patches_100"]] = 
  simulationModel$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 100)

i = "patches_100"
for(i in c("patches_1", "patches_100")){
  # p_dat <- predictions[[i]]$long$site[variable == "ba", .(value = mean(value)), by = .(year, species, variable)]
  p_dat <- predictions[[i]]$long$site[, .(value = mean(value)), by = .(year, species, variable)]
  p_dat[, variable2 := factor(
    variable,
    levels = c("dbh", "ba", "trees", "AL", "growth", "mort", "reg", "r_mean_ha"),
    labels =  c("avg. DBH [cm]", "Basal Area [m²/ha]", "Trees [N/ha]", 
                "Available Light [%]", "Growth [cm/yr]", "Mortality [%]", 
                "Reg. Count [N/ha]", "Reg. Mean [N/ha]")
  ),]
  # p_dat[variable %in% c("ba", "trees", "reg"), value := value/patch_size,]
  p_dat[variable %in% c("ba", "trees"), value := value/patch_size,]
  p <- ggplot(p_dat[year <= 100], aes(x = year, y = value, color = factor(species))) +
    geom_line() +
    theme_minimal() +
    labs(x = "Year", y = "Value") +
    coord_cartesian(ylim = c(0, NA)) +
    facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left") +  # Remove label_parsed
    theme(
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90)
    ) +
    guides(color = guide_legend(title = "Species", override.aes = list(linewidth = 5), ncol = 2, title.position = "top")) +
    scale_color_discrete(name = "Species") +
    ggtitle(paste0("patches = ",gsub("patches_", "", i)))
  
  print(p)
}

