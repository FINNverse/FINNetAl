#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## test implication of biased rates ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(FINN)
library(data.table)
library(ggplot2)

Ntimesteps = 50  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 5 # number of species
FINN.seed(1234)
# we draw the same shade parameters for each process for simplicity
# shade parameters correspond to the fraction of light a species needs to succesfully grow, regenerate, or survive.
shadeSP = runif(Nsp, 0.01, 0.99)

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  runif(Nsp, 3, 4), # intercept regulating the overall effect size
  runif(Nsp, -2, 2) # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# growth parameters
parGrowth = matrix(c(
  shadeSP, # see above
  runif(Nsp, 0, 0.1) # the second growth parameter modulates the size dependent growth
),Nsp, 2)

parGrowthEnv = list(matrix(c(
  runif(Nsp, -1, 2), # intercept regulating the overall effect size
  runif(Nsp, -1, 2) # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# mortality parameters
parMort = matrix(c(
  shadeSP, # see above
  runif(Nsp, 1, 4) # the second growth parameter modulates the size dependent mortality
),Nsp, 2)
parMortEnv = list(matrix(c(
  runif(Nsp, 0, 0.5), # intercept regulating the overall effect size
  runif(Nsp, 0, .5) # the second parameter modulates the effect of the environmental variable
), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
# parHeight = runif(Nsp, 0.3, 0.7)
# growth parameters
parComp = matrix(c(
  runif(Nsp, 0.3, 0.7), # parHeight
  runif(Nsp, 0.2, 0.2) # Competition strength
),Nsp, 2)

## Environment and disturbances {#sec-env-dist}

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
disturbance_frequency = 0.1

# the disturbance intensity at each timestep is the fraction of patches that is disturbed at that time step
# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_intensity = rbinom(Ntimesteps*Nsites,1,0.2)*runif(Ntimesteps*Nsites, 0.5, 1)

# this will result in 0 to 20 % of the patches being disturbed at each timestep with a change of 1% that a disturbance occurs at the timestep
dist_dt$intensity = dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity



mort2 = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F) {
  mort = torch::torch_zeros_like(dbh) + 0.01
  return(mort)
}


## Simulate
simulationModel = finn(N_species = Nsp,
                       competition_process = createProcess(~0, func = FINN::competition),
                       growth_process = createProcess(~1+env1, initEnv = parGrowthEnv,initSpecies = parGrowth, func = FINN::growth),
                       regeneration_process = createProcess(~1+env1, initEnv = parRegEnv,initSpecies = parReg, func = FINN::regeneration),
                       # mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = FINN::mortality),
                       mortality_process = createProcess(~1+env1, initEnv = parMortEnv,initSpecies = parMort, func = mort2),
)

predictions =
  simulationModel$simulate(init_cohort = NULL, env = env_dt, device = "cpu", patches = 50, debug = T)

hist(predictions$wide$cohort$trees)
hist(predictions$wide$cohort$m)
summary(predictions$wide$site$mort)
hist(predictions$wide$site$mort) # 0.01
median(predictions$wide$site$mort) # 0.01

site_dt_full <- copy(predictions$wide$site)
cohort_dt_full <- copy(predictions$wide$cohort)

# check g.finn vs g.derived
cohort_dt_full[,":="(
  dbh_before = data.table::shift(dbh, 1, type = "lag"),
  year_before = data.table::shift(year, 1, type = "lag")
), by = .(patchID, cohortID)]

cohort_dt_full[is.na(m), reg_year := year,]

# 1) For each cohort, order by year
# 2) If m is NA, keep its reg_year in a temporary column; otherwise NA
# 3) Fill-forward (locf) that temporary column
# 4) Wherever m is not NA, overwrite reg_year with the filled-forward value
# 5) Drop the temporary column
cohort_dt_full[
  order(patchID, cohortID, year),
  `:=`(tmp_reg_year = ifelse(is.na(m), reg_year, NA_real_)),
  by = cohortID][
  , tmp_reg_year := nafill(tmp_reg_year, type="locf"),
  by = cohortID][
  !is.na(m), reg_year := tmp_reg_year][
  , tmp_reg_year := NULL]
cohort_dt_full[,cohortID := as.numeric(as.factor(paste0(cohortID,reg_year))),]

cohort_dt_full[,":="(
  dbh_before = data.table::shift(dbh, 1, type = "lag")
), by = cohortID]
cohort_dt_full[trees > 0 & dbh != 0, g_rel := (dbh-dbh_before)/dbh_before,]
plot(cohort_dt_full$g_rel, cohort_dt_full$g)


cohort_dt_thin <- copy(predictions$wide$cohort)[year %in% seq(1, Ntimesteps, 5)]
cohort_dt_thin[,":="(
  dbh_before = data.table::shift(dbh, 1, type = "lag"),
  year_before = data.table::shift(year, 1, type = "lag"),
  trees_before = data.table::shift(trees, 1, type = "lag")
), by = cohortID]

cohort_dt_thin[trees > 0 & dbh != 0, g_rel := ((dbh-dbh_before)/dbh_before)/5,]

ggplot()+
  geom_line(data = cohort_dt_thin, aes(x = year, y = g_rel, color = "thin", group = cohortID))+
  geom_line(data = cohort_dt_full, aes(x = year, y = g_rel, color = "full", group = cohortID))+
  facet_wrap(~species)
growth_dt_cohorts_thin = cohort_dt_thin[,.(
  g_rel = sum(g_rel*trees,na.rm = T)/sum(trees, na.rm = T)
), by = .(species, year)]
growth_dt_cohorts_full = cohort_dt_full[,.(
  g_rel = sum(g_rel*trees,na.rm = T)/sum(trees, na.rm = T)
), by = .(species, year)]

ggplot()+
  geom_line(data = growth_dt_cohorts_thin, aes(x = year, y = g_rel, color = "thin"))+
  geom_line(data = growth_dt_cohorts_full, aes(x = year, y = g_rel, color = "full"))+
  geom_line(data = site_dt_full, aes(x = year, y = growth, color = "finn"))+
  facet_wrap(~species)


# mortality
cohort_dt_full[is.na(year_before), reg := trees,]
mort_cohort_dt_full <- cohort_dt_full[,.(
  trees_alive = sum(trees, na.rm = T),
  trees_survived = sum(trees, na.rm = T)-sum(reg, na.rm = T),
  trees_before = sum(trees_before, na.rm = T),
  reg = sum(reg, na.rm = T)
), by = .(species, year, patchID)]
mort_cohort_dt_full[, trees_died := trees_before - trees_survived,]
mort_cohort_dt_full[trees_before > 0, ":="(
  m2 = mean(trees_died/trees_before,na.rm = T)
), by = .(species, year)]
mort_cohort_dt_full[, ":="(
  reg = mean(reg, na.rm = T),
  trees_before = mean(trees_before, na.rm = T)
), by = .(species, year)]
mort_cohort_dt_full[m2 < 0, m2 := 0]

cohort_dt_thin[is.na(year_before), reg := trees,]
mort_cohort_dt_thin <- cohort_dt_thin[,.(
  trees_alive = sum(trees, na.rm = T),
  trees_survived = sum(trees, na.rm = T)-sum(reg, na.rm = T),
  trees_before = sum(trees_before, na.rm = T),
  reg = sum(reg, na.rm = T)
), by = .(species, year, patchID)]
mort_cohort_dt_thin[trees_before > 0, trees_died := trees_before - trees_survived,]
mort_cohort_dt_thin[, ":="(
  m2 = mean(trees_died/trees_before, na.rm = T)/5
  ), by = .(species, year)]
mort_cohort_dt_thin[, ":="(
  reg = mean(reg, na.rm = T)/5,
  trees_before = mean(trees_before, na.rm = T)
  ), by = .(species, year)]
mort_cohort_dt_thin[m2 < 0, m2 := 0]


ggplot()+
  geom_line(data = mort_cohort_dt_thin, aes(x = year, y = m2, color = "thin"), alpha = 0.3)+
  geom_line(data = mort_cohort_dt_full, aes(x = year, y = m2, color = "full"), alpha = 0.3)+
  geom_line(data = site_dt_full, aes(x = year, y = mort, color = "finn"), alpha = 0.7)+
  geom_line(data = site_dt_full, aes(x = year, y = trees_dead/trees_before, color = "finn2"), alpha = 0.3)+
  facet_wrap(~species, scales = "free")+
  coord_cartesian(ylim = c(0,0.02))
  # geom_hline(yintercept = 0.01, linetype = "dashed", color = "red")

ggplot()+
  geom_line(data = mort_cohort_dt_thin, aes(x = year, y = reg, color = "thin"),alpha = 0.4)+
  geom_line(data = mort_cohort_dt_full, aes(x = year, y = reg, color = "full"),alpha = 0.4)+
  geom_line(data = site_dt_full, aes(x = year, y = reg, color = "finn"),alpha = 0.4)+
  facet_wrap(~species)

mort_cohort_dt_compfinn_thin <- merge(mort_cohort_dt_thin, site_dt_full, by = c("species", "year"), suffixes = c(".thin", ".finn"))
mort_cohort_dt_compfinn_full <- merge(mort_cohort_dt_full, site_dt_full, by = c("species", "year"), suffixes = c(".thin", ".finn"))
growth_cohort_dt_compfinn_thin <- merge(growth_dt_cohorts_thin, site_dt_full, by = c("species", "year"), suffixes = c(".thin", ".finn"))
growth_cohort_dt_compfinn_full <- merge(growth_dt_cohorts_full, site_dt_full, by = c("species", "year"), suffixes = c(".thin", ".finn"))

par(mfrow=c(3,2), mar=c(4,4,2,0.5))
plot(mort_cohort_dt_compfinn_thin$m2, mort_cohort_dt_compfinn_thin$mort, xlab = "thin", ylab = "finn", main = "mortality")
abline(0,1, col = "red")
plot(mort_cohort_dt_compfinn_full[m2 < 1]$m2, mort_cohort_dt_compfinn_full[m2 < 1]$mort, xlab = "full", ylab = "finn", main = "mortality")
abline(0,1, col = "red")

plot(mort_cohort_dt_compfinn_thin$reg.thin, mort_cohort_dt_compfinn_thin$reg.finn, xlab = "thin", ylab = "finn", main = "regeneration")
abline(0,1, col = "red")
plot(mort_cohort_dt_compfinn_full$reg.thin, mort_cohort_dt_compfinn_full$reg.finn , xlab = "full", ylab = "finn", main = "regeneration")
abline(0,1, col = "red")

plot(growth_cohort_dt_compfinn_thin$g_rel, growth_cohort_dt_compfinn_thin$growth, xlab = "thin", ylab = "finn", main = "growth")
abline(0,1, col = "red")
plot(growth_cohort_dt_compfinn_full$g_rel, growth_cohort_dt_compfinn_full$growth, xlab = "full", ylab = "finn", main = "growth")
abline(0,1, col = "red")




ggplot()+
  geom_line(data = cohort_dt_thin[!is.na(m),.(m2 = sum(m2*trees, na.rm=T)/sum(trees,na.rm=T)), .(species, year)], aes(x = year, y = g_rel, color = "thin"))+
  geom_line(data = cohort_dt_full[!is.na(m),.(m2 = sum(m2*trees, na.rm=T)/sum(trees,na.rm=T)), .(species, year)], aes(x = year, y = g_rel, color = "full"))+
  geom_line(data = site_dt_full, aes(x = year, y = mort, color = "finn"))+
  facet_wrap(~species)


ggplot(cohort_dt_full[cohortID %in% cohort_dt_full[g_rel < 0]$cohortID], aes(x = year, y = dbh))+
  geom_point()+
  facet_wrap(~cohortID, ncol=1, scales = "free")


summary(cohort_dt_full$g_rel)
summary(cohort_dt_full$g)


year_seq = seq(min(site_dt_full$year, na.rm = T), max(site_dt_full$year), 5)
site_dt_incomplete <- site_dt_full[year %in% year_seq, -c("growth", "mort", "reg")]
cohort_dt_incomplete <- cohort_dt_full[year %in% year_seq, -c("m", "g", "trees_before")]

# check growth
cohort_dt_full[,":="(
  dbh_before = data.table::shift(dbh, 1, type = "lag")
), by = cohortID]
#  with relative growth
cohort_dt_full[, g_rel := (dbh-dbh_before)/dbh_before,]

growth_comp_dt <- merge(cohort_dt_full, site_dt_full[,.(species, year, growth)], by = c("species", "year"))

plot(growth_comp_dt$g_rel, growth_comp_dt$growth)


cohort_dt_full[, g_abs := (dbh-dbh_before),]
growth_dt_incomplete = cohort_dt_full[,.(
  g_rel = sum(g_rel*trees,na.rm = T)/sum(trees, na.rm = T),
  g_abs = (sum(g_abs*trees,na.rm = T)/sum(trees, na.rm = T)/5)
  ), by = .(species, year)]
growth_dt_full = site_dt_full[,.(g = mean(growth,na.rm = T)), by = .(species, year)]
ggplot() +
  geom_line(data = cohort_dt_incomplete, aes(x = year, y = g_rel, color = "incomplete rel"))+
  geom_line(data = cohort_dt_incomplete, aes(x = year, y = g_abs, color = "incomplete abs"))+
  geom_line(data = growth_dt_full, aes(x = year, y = g, color = "full rel"))+
  facet_wrap(~species)+
  ylab("")
# with absolute growth
cohort_dt_incomplete[, g := (dbh-dbh_before),]
growth_dt_incomplete = cohort_dt_incomplete[,.(g = sum(g*trees,na.rm = T)/sum(trees, na.rm = T)), by = .(species, year)]
growth_dt_full = site_dt_full[,.(g = mean(growth,na.rm = T)), by = .(species, year)]
ggplot() +
  geom_line(data = growth_dt_cohorts, aes(x = year, y = g, color = "incomplete"))+
  geom_line(data = growth_dt_full, aes(x = year, y = g, color = "full"))+
  facet_wrap(~species)+
  ylab("absolute growth")





merge(
  growth_dt_cohorts,
  site_dt_full[,.(species, year, growth)], by = c("species", "year"))


trees_incomplete =
cohort_dt_incomplete[,.(
  trees = sum(trees),
), by = .(species, year)][order(year)]

trees_incomplete[,trees_before, by = species]

obs_dt_incomplete







