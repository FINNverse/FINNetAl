#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# plot simulated data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(data.table)
library(ggplot2)
library(FINN)
seed = 123

models <- list.files("results/02_realdata/", full.names = T, pattern = "pft-period7-25patches_S0_T0")
# models <- list.files("results/01_full/", pattern = "01", full.names = T)
m_list = lapply(models, function(x) torch::torch_load(x))
m_names <- sapply(strsplit(models, "/"), tail, 1)
names(m_list) <- sapply(strsplit(m_names, "_"), function(x) x[1])
# pdf("figures/01-results.pdf", width = 10, height = 7)
# for(i in 1:length(m_list)){
FINN.seed(seed)

i=1
i_name = names(m_list)[i]
m = m_list[[i]]
if(grepl("1patch$", i_name)){
  Npatches = 1L
} else if(grepl("25patches$", i_name)){
  Npatches = 25L
}

# read full raw data
obs_dt = fread(paste0("data/BCI/noSplits/", i_name,"/obs_dt.csv"))
env_dt = fread(paste0("data/BCI/noSplits/", i_name,"/env_dt.csv"))
cohorts_dt = fread(paste0("data/BCI/noSplits/", i_name,"/initial_cohorts1985.csv"))

# predict for model
cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = uniqueN(obs_dt$species))

env_dt_long <- data.table()
start_t = 0
for(k in 1:20){
  env_dt_temp = copy(env_dt)
  # randomly shuffle years in env_dt
  env_dt_temp$year = env_dt_temp$year+start_t
  env_dt_temp$year = as.integer(as.character(factor(env_dt_temp$year, levels = sample(unique(env_dt_temp$year)))))
  env_dt_long = rbind(env_dt_long, env_dt_temp)
  start_t = max(env_dt_long$year)
  max(env_dt_long$year) == uniqueN(env_dt_long$year)
}

# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_frequency = 1

Ntimesteps = uniqueN(env_dt_long$year)
Ntimesteps = max(env_dt_long$year)
Nsites = uniqueN(env_dt_long$siteID)
# the disturbance intensity at each timestep is the fraction of patches that is disturbed at that time step
# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_intensity = runif(Ntimesteps*Nsites, 0.0043, 0.016)
length(disturbance_intensity)
dist_dt <- env_dt_long[,.(siteID, year)]
# this will result in 0 to 20 % of the patches being disturbed at each timestep with a change of 1% that a disturbance occurs at the timestep
dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity
m2 = finn(
  N_species = m$N_species,
  competition_process = createProcess(~0, func = FINN::competition, initSpecies = m$par_competition_r),
  growth_process = createProcess(~., initEnv = m$parameters_r$nn_growth.0.weight, func = FINN::growth, initSpecies = m$par_growth_r),
  regeneration_process = createProcess(~., initEnv = m$parameters_r$nn_regeneration.0.weight, func = FINN::regeneration, dispersion_parameter = m$par_theta_recruits, initSpecies = m$par_regeneration_r),
  mortality_process = createProcess(~., initEnv = m$parameters_r$nn_mortality.0.weight, func = FINN::mortality, m$par_mortality_r)
)
m3 = finn(
  N_species = m$N_species,
  competition_process = m$process_competition,
  growth_process = m$process_growth,
  regeneration_process = m$process_regeneration,
  mortality_process = m$process_mortality
)
# pred = m2$simulate(env = env_dt_long, disturbance = dist_dt, init_cohort = NULL, patches = Npatches)
pred = m3$simulate(env = env_dt_long, disturbance = dist_dt, init_cohort = NULL, patches = Npatches)
# pred = m$simulate(env = env_dt, init_cohort = cohort1, patches = Npatches)
pred_dt = pred$wide$site
pred_dt$pred = "pred"
pred_dt[,reg := reg/0.1,] # correct regeneration response
obs_dt$pred = "obs"

# create comparison data
comp_dt1 = rbind(pred_dt, obs_dt, fill = T)
comp_dt1 <- rbind(comp_dt1, comp_dt1[, .(
  species = "all",
  dbh = mean(dbh, na.rm = T),
  ba = sum(ba),
  trees = sum(trees),
  reg = mean(reg, na.rm = T),
  mort = mean(mort, na.rm = T),
  growth = mean(growth, na.rm = T)
) , by = .(siteID, year, pred, period_length)], fill = T)

# pick 10 most abundant species for plotting
top10_species = obs_dt[,.(ba = sum(ba)),by=species][order(ba, decreasing = T)]$species[1:(min(c(10, uniqueN(obs_dt$species))))]
# comp_dt1 = comp_dt1[species %in% c(top10_species, "all")]
# comp_dt1 = comp_dt1[species %in% c(top10_species, "all")]
melt_dt = melt(comp_dt1, id.vars = c("siteID", "year", "species", "pred", "period_length"), variable.name = "variable")
p_dat = melt_dt[,.(value = mean(value)), by = .(species,year, pred,variable)]
p1 = ggplot()+
  geom_point(data = p_dat[variable != "r_mean_ha"],
             aes(x = year, y = value, color = factor(species), shape = pred), alpha = 0.3)+
  geom_line(data = p_dat[variable != "r_mean_ha"],
            aes(x = year, y = value, color = factor(species), shape = pred, linetype = pred),alpha = 0.3)+
  geom_hline(data = p_dat[variable != "r_mean_ha",.(
    yintercept = mean(value[year > 0 | pred == "obs"], na.rm = T)
  ), by = .(species, pred, variable)], aes(yintercept = yintercept, color = factor(species), linetype = pred),linewidth= 1, alpha = 0.7)+
  scale_size_manual(values = c(pred = 2, obs = 3))+
  facet_wrap(~variable, scales = "free_y")+
  ggtitle(i_name)
print(p1)

comp_dt1_means = comp_dt1[year > 0 | pred == "obs",.(
  ba = mean(ba, na.rm = T),
  trees = mean(trees, na.rm = T),
  dbh = mean(dbh, na.rm = T),
  reg = mean(reg, na.rm = T),
  mort = mean(mort, na.rm = T),
  growth = mean(growth, na.rm = T)
), by = .(pred, species)]
compd_dt1_means_dcast <- dcast(comp_dt1_means, species~pred, value.var = c("ba", "trees", "dbh", "reg", "mort", "growth"), sep = ".")
par(mfrow = c(2, 3), pty = "s")
for(j in c("ba", "trees", "dbh", "reg", "mort", "growth")){
  plot(
    compd_dt1_means_dcast[[paste0(j,".obs")]]~compd_dt1_means_dcast[[paste0(j,".pred")]], col = as.factor(compd_dt1_means_dcast[["species"]]), main = paste0(j,"\nCorrelation: ", round(cor(compd_dt1_means_dcast[[paste0(j,".obs")]], compd_dt1_means_dcast[[paste0(j,".pred")]], use = "complete.obs", method = "spearman"), 2)),
    xlim = range(c(compd_dt1_means_dcast[[paste0(j,".obs")]],compd_dt1_means_dcast[[paste0(j,".pred")]]), na.rm = T), ylim = range(c(compd_dt1_means_dcast[[paste0(j,".obs")]],compd_dt1_means_dcast[[paste0(j,".obs")]]), na.rm = T),
    asp = 1, cex = 2, pch = 8, xlab = "prediction", ylab = "observation")
  # add legend for colors outside of plot
  abline(0, 1)
}
legend("topright", legend = levels(as.factor(compd_dt1_means_dcast[["species"]])), col = 1:length(levels(as.factor(compd_dt1_means_dcast[["species"]]))) , pch = 8)
plot(
  compd_dt1_means_dcast[[paste0(j,".obs")]]~compd_dt1_means_dcast[[paste0(j,".pred")]], col = as.factor(compd_dt1_means_dcast[["species"]]), main = paste0(i_name,"\nCorrelation: ", round(cor(compd_dt1_means_dcast[[paste0(j,".obs")]], compd_dt1_means_dcast[[paste0(j,".pred")]], use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(compd_dt1_means_dcast[[paste0(j,".obs")]],compd_dt1_means_dcast[[paste0(j,".pred")]]), na.rm = T), ylim = range(c(compd_dt1_means_dcast[[paste0(j,".obs")]],compd_dt1_means_dcast[[paste0(j,".obs")]]), na.rm = T),
  asp = 1)
abline(0, 1)

all_years = unique(pred_dt$year)
all_obs_years = unique(obs_dt$year)
period_length = all_obs_years[2]-all_obs_years[1]
intervals = rep(unique(obs_dt$year), each = period_length)
growth_mort_pred_dt =
  pred_dt[,.(
    growth = mean(growth),
    mort = mean(mort)
  ), by = .(
    siteID,
    year = as.integer(as.character(factor(year, levels = all_years, labels = intervals))),
    species)]

df1 = merge.data.table(obs_dt[,-c("growth", "mort")], pred_dt[,-c("growth", "mort")], by = c("siteID", "year", "species"), suffixes = c(".obs", ".pred"))
df2 = merge.data.table(obs_dt[,.(siteID, year, species, growth, mort)], growth_mort_pred_dt[,.(siteID, year, species, growth, mort)], by = c("siteID", "year", "species"), suffixes = c(".obs", ".pred"))
df = merge.data.table(df1, df2, by = c("siteID", "year", "species"))
comp_allspecies_dt <- df[,.(
  ba.obs = sum(ba.obs)/uniqueN(siteID),
  ba.pred = sum(ba.pred)/uniqueN(siteID),
  trees.obs = sum(trees.obs)/uniqueN(siteID),
  trees.pred = sum(trees.pred)/uniqueN(siteID),
  dbh.obs = mean(dbh.obs, na.rm = TRUE),
  dbh.pred = mean(dbh.pred, na.rm = TRUE),
  reg.obs = mean(reg.obs, na.rm = TRUE),
  reg.pred = mean(r_mean_ha, na.rm = TRUE),
  mort.obs = mean(mort.obs, na.rm = TRUE),
  mort.pred = mean(mort.pred, na.rm = TRUE),
  growth.obs = mean(growth.obs, na.rm = T),
  growth.pred = mean(growth.pred, na.rm = T)
), by = .(siteID,year, species)]

par(mfrow = c(2, 3), pty = "s")
# plot with correlation in title
plot(
  growth.obs~growth.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$growth.obs, comp_allspecies_dt$growth.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$growth.obs,comp_allspecies_dt$growth.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$growth.obs,comp_allspecies_dt$growth.pred), na.rm = T),
  asp = 1)
abline(0, 1)
plot(
  reg.obs~reg.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$reg.obs, comp_allspecies_dt$reg.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$reg.obs,comp_allspecies_dt$reg.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$reg.obs,comp_allspecies_dt$reg.pred), na.rm = T),
  asp = 1)
abline(0, 1)
plot(
  mort.obs~mort.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$mort.obs, comp_allspecies_dt$mort.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$mort.obs,comp_allspecies_dt$mort.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$mort.obs,comp_allspecies_dt$mort.pred), na.rm = T),
  asp = 1)
abline(0, 1)
plot(
  ba.obs~ba.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$ba.obs, comp_allspecies_dt$ba.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$ba.obs,comp_allspecies_dt$ba.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$ba.obs,comp_allspecies_dt$ba.pred), na.rm = T),
  asp = 1)
abline(0, 1)
plot(
  trees.obs~trees.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$trees.obs, comp_allspecies_dt$trees.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$trees.obs,comp_allspecies_dt$trees.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$trees.obs,comp_allspecies_dt$trees.pred), na.rm = T),
  asp = 1)
abline(0, 1)
plot(
  dbh.obs~dbh.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$dbh.obs, comp_allspecies_dt$dbh.pred, use = "complete.obs", method = "spearman"), 2)),
  xlim = range(c(comp_allspecies_dt$dbh.obs,comp_allspecies_dt$dbh.pred), na.rm = T), ylim = range(c(comp_allspecies_dt$dbh.obs,comp_allspecies_dt$dbh.pred), na.rm = T),
  asp = 1)
abline(0, 1)
# same plots with ranges
par(mfrow = c(2, 3))
# plot with correlation in title
plot(
  growth.obs~growth.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$growth.obs, comp_allspecies_dt$growth.pred, use = "complete.obs", method = "spearman"), 2))
  ,xlim = range(comp_allspecies_dt$growth.obs, na.rm = T), ylim = range(comp_allspecies_dt$growth.obs, na.rm = T))
abline(0, 1)
plot(
  reg.obs~reg.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$reg.obs, comp_allspecies_dt$reg.pred, use = "complete.obs", method = "spearman"), 2))
  ,xlim = range(comp_allspecies_dt$reg.obs, na.rm = T), ylim = range(comp_allspecies_dt$reg.obs, na.rm = T))
abline(0, 1)
plot(
  mort.obs~mort.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$mort.obs, comp_allspecies_dt$mort.pred, use = "complete.obs", method = "spearman"), 2))
  , xlim = range(comp_allspecies_dt$mort.obs, na.rm = T), ylim = range(comp_allspecies_dt$mort.obs, na.rm = T))
abline(0, 1)
plot(ba.obs~ba.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$ba.obs, comp_allspecies_dt$ba.pred, use = "complete.obs", method = "spearman"), 2))
     , xlim = range(comp_allspecies_dt$ba.obs, na.rm = T), ylim = range(comp_allspecies_dt$ba.obs, na.rm = T))
abline(0, 1)
plot(trees.obs~trees.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$trees.obs, comp_allspecies_dt$trees.pred, use = "complete.obs", method = "spearman"), 2)),
     xlim = range(comp_allspecies_dt$trees.obs, na.rm = T), ylim = range(comp_allspecies_dt$trees.obs, na.rm = T))
abline(0, 1)
plot(dbh.obs~dbh.pred, data = comp_allspecies_dt, col = comp_allspecies_dt$species, main = paste0(i_name,"\nCorrelation: ", round(cor(comp_allspecies_dt$dbh.obs, comp_allspecies_dt$dbh.pred, use = "complete.obs", method = "spearman"), 2)),
     xlim = range(comp_allspecies_dt$dbh.obs, na.rm = T), ylim = range(comp_allspecies_dt$dbh.obs, na.rm = T))
abline(0, 1)
}
dev.off()


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# make splits ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
i = 1
# "pft-period7-25patches_S1_T0"
# i_name = "pft-period7-25patches"
m_list = c(m_list, rep("realdata",8))
names(m_list)[m_list == "realdata"] = c(
  "pft-period7-1patch-realdata",
  "pft-period7-25patches-realdata",
  "pft-period35-1patch-realdata",
  "pft-period35-25patches-realdata",
  "genus-period7-1patch-realdata",
  "genus-period7-25patches-realdata",
  "genus-period35-1patch-realdata",
  "genus-period35-25patches-realdata"
)
# genus-period7-25patches
i=7
for(i in 1:length(m_list)){
  FINN.seed(seed)
  i_name = names(m_list)[i]
  i_folder = gsub("-realdata","", i_name)
  m = m_list[[i]]

  Npatches = ifelse(grepl("1patch", i_name), 1L, 25L)
  Nepochs = ifelse(grepl("period7", i_name), 7L, 35L)

  period_length = Nepochs/7

  # read full raw data
  obs_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/obs_dt.csv"))
  env_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/env_dt.csv"))
  cohorts_dt = fread(paste0("data/BCI/noSplits/", i_folder,"/initial_cohorts1985.csv"))

  # simulate or assign realdata
  if(!grepl("realdata", i_name)){
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = uniqueN(obs_dt$species))
    pred = m$simulate(env = env_dt, init_cohort = cohort1, patches = Npatches)
    pred_dt = pred$wide$site
    pred_dt[,reg := reg/0.1,]

    out_dir0 = paste0("data/BCI/CVsplits-simdata/")
  }else if(grepl("realdata", i_name)){
    pred_dt = copy(obs_dt)

    out_dir0 = paste0("data/BCI/CVsplits-realdata/")
  }

  # create temporal folds
  years <- c(1990, 1995, 2000, 2005, 2010, 2015)
  if(grepl("period7", i_name)){
    all_years = seq(1985,2015, 5)
  }else{
    all_years = 1985:2015
  }
  obs_dt_idx = seq(min(obs_dt$year), max(obs_dt$year), by = max(obs_dt$period_length, na.rm = T))
  env_dt_idx = seq(min(env_dt$year), max(env_dt$year), by = 1)
  temporal_folds_list <- list(
    list(test = c(1,2), train = c(3,6)),
    list(test = c(5,6), train = c(1,4))
  )

  temporal_folds_dt = data.table()
  k = 1
  holdout = "train"
  for(k in 1:length(temporal_folds_list)){
    for(holdout in c("test", "train")){
      init_year = years[temporal_folds_list[[k]][[holdout]][1]]-5
      env_start = which(init_year == all_years)
      if(grepl("period7", i_name)) env_start = env_start+1
      temporal_folds_dt = rbind(
        temporal_folds_dt,
        data.table(
          splitID = k,
          start = years[temporal_folds_list[[k]][[holdout]][1]], end = years[temporal_folds_list[[k]][[holdout]][2]],
          obs_start = obs_dt_idx[temporal_folds_list[[k]][[holdout]][1]], obs_end = obs_dt_idx[temporal_folds_list[[k]][[holdout]][2]],
          env_start = env_start, env_end = env_dt_idx[temporal_folds_list[[k]][[holdout]][2]*period_length],
          init_year = init_year, holdout = holdout
        )
      )
    }
  }

  # read spatial holdouts
  spatial_folds_dt = fread(paste0("data/BCI/noSplits/",i_folder,"/spatial_folds_dt.csv"))[,.(siteID, fold = spatial_fold)]

  if(!dir.exists(paste0(out_dir0, i_folder, "/")))
    dir.create(paste0(out_dir0, i_folder), recursive = T)

  fwrite(temporal_folds_dt,paste0("data/BCI/noSplits/",i_folder,"/temporal_folds_dt_applied.csv"))
  fwrite(spatial_folds_dt,paste0("data/BCI/noSplits/",i_folder,"/spatial_folds_dt_applied.csv"))

  temporal_splits = unique(temporal_folds_dt$splitID)
  spatial_splits = sort(unique(spatial_folds_dt$fold))
  for(i_T in c(0,temporal_splits)){
    for(i_S in c(0, spatial_splits)){
      # define years of initial cohorts for train and test
      if(i_T == 0){
        train_cohort_init_year = 1985
        test_cohort_init_year = 1985
        train_years_obs = seq(min(temporal_folds_dt$obs_start), max(temporal_folds_dt$obs_end), period_length)
        test_years_obs = seq(min(temporal_folds_dt$obs_start), max(temporal_folds_dt$obs_end), period_length)
        train_years_env = seq(min(temporal_folds_dt$env_start), max(temporal_folds_dt$env_end), 1)
        test_years_env = seq(min(temporal_folds_dt$env_start), max(temporal_folds_dt$env_end), 1)
      }else{
        Tfolds_dt = temporal_folds_dt[splitID == i_T]
        test_years = seq(Tfolds_dt[holdout == "test"]$start, Tfolds_dt[holdout == "test"]$end, period_length)
        train_years = seq(Tfolds_dt[holdout == "train"]$start, Tfolds_dt[holdout == "train"]$end, period_length)
        train_years_obs = seq(Tfolds_dt[holdout == "train"]$obs_start, Tfolds_dt[holdout == "train"]$obs_end, period_length)
        test_years_obs = seq(Tfolds_dt[holdout == "test"]$obs_start, Tfolds_dt[holdout == "test"]$obs_end, period_length)
        train_years_env = seq(Tfolds_dt[holdout == "train"]$env_start, Tfolds_dt[holdout == "train"]$env_end, 1L)
        test_years_env = seq(Tfolds_dt[holdout == "test"]$env_start, Tfolds_dt[holdout == "test"]$env_end, 1)
        train_cohort_init_year = Tfolds_dt[holdout == "train"]$init_year
        test_cohort_init_year = Tfolds_dt[holdout == "test"]$init_year
      }

      # read correct initial cohorts for train and test
      initial_cohort_train = fread(paste0("data/BCI/noSplits/",i_folder,"/initial_cohorts",train_cohort_init_year,".csv"))
      initial_cohort_test = fread(paste0("data/BCI/noSplits/",i_folder,"/initial_cohorts",test_cohort_init_year,".csv"))

      # get siteIDs of spatial folds for train and test
      if(i_S == 0){
        train_sites = unique(obs_dt$siteID)
        test_sites = unique(obs_dt$siteID)
      }else{
        train_sites = spatial_folds_dt[fold != i_S, siteID]
        test_sites = spatial_folds_dt[fold == i_S, siteID]
      }
      Sfolds_dt = rbind(
        data.table(siteID_full = train_sites, splitID = i_S, holdout = "train"),
        data.table(siteID_full = test_sites, splitID = i_S, holdout = "test")
      )
      # add new siteID to ensure siteID starting from 1
      Sfolds_dt[, siteID_holdout := 1:.N, by = holdout]

      # get sites for training and testing from fully predicted data set
      obs_train_dt = pred_dt[siteID %in% train_sites & year %in% train_years_obs]
      obs_test_dt = pred_dt[siteID %in% test_sites & year %in% test_years_obs]
      env_train_dt = env_dt[siteID %in% train_sites & year %in% train_years_env]
      env_test_dt = env_dt[siteID %in% test_sites & year %in% test_years_env]
      init_cohort_train = initial_cohort_train[siteID %in% train_sites]
      init_cohort_test = initial_cohort_test[siteID %in% test_sites]

      # replace original site id with new siteIDs that start from 1
      obs_train_dt = merge(obs_train_dt, Sfolds_dt[holdout == "train"], by.x = "siteID", by.y = "siteID_full")
      obs_train_dt[,siteID := siteID_holdout,]
      obs_test_dt = merge(obs_test_dt, Sfolds_dt[holdout == "test"], by.x = "siteID", by.y = "siteID_full")
      obs_test_dt[,siteID := siteID_holdout,]
      env_train_dt = merge(env_train_dt, Sfolds_dt[holdout == "train"], by.x = "siteID", by.y = "siteID_full")
      env_train_dt[,siteID := siteID_holdout,]
      env_test_dt = merge(env_test_dt, Sfolds_dt[holdout == "test"], by.x = "siteID", by.y = "siteID_full")
      env_test_dt[,siteID := siteID_holdout,]
      init_cohort_test = merge(init_cohort_test, Sfolds_dt[holdout == "test"], by.x = "siteID", by.y = "siteID_full")
      init_cohort_test[,siteID := siteID_holdout,]
      init_cohort_train = merge(init_cohort_train, Sfolds_dt[holdout == "train"], by.x = "siteID", by.y = "siteID_full")
      init_cohort_train[,siteID := siteID_holdout,]

      obs_train_dt[,year := year-min(year)+period_length]
      # obs_train_dt$year = as.integer(as.factor(obs_train_dt$year))
      obs_test_dt[,year := year-min(year)+period_length]
      # obs_test_dt$year = as.integer(as.factor(obs_test_dt$year))
      # env_train_dt[,year := year-min(year)+min_year]
      env_train_dt$year = as.integer(as.factor(env_train_dt$year))
      env_test_dt$year = as.integer(as.factor(env_test_dt$year))

      # save final files
      cv_label = paste0("S", i_S, "_T", i_T)
      fwrite(obs_train_dt, paste0(out_dir0, i_folder, "/obs_dt_", cv_label, "_train.csv"))
      fwrite(obs_test_dt, paste0(out_dir0, i_folder, "/obs_dt_", cv_label, "_test.csv"))
      fwrite(env_train_dt, paste0(out_dir0, i_folder, "/env_dt_", cv_label, "_train.csv"))
      fwrite(env_test_dt, paste0(out_dir0, i_folder, "/env_dt_", cv_label, "_test.csv"))
      fwrite(init_cohort_test, paste0(out_dir0, i_folder, "/initial_cohorts_", cv_label, "_test.csv"))
      fwrite(init_cohort_train, paste0(out_dir0, i_folder, "/initial_cohorts_", cv_label, "_train.csv"))
    }
  }
}


