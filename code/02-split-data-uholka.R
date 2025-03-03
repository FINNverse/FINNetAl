#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# plot simulated data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(data.table)
library(ggplot2)
library(FINN)
seed = 123

models <- list.files("results/01_full/", full.names = T)
# models <- list.files("results/01_full/", pattern = "01", full.names = T)
m_list = lapply(models, function(x) torch::torch_load(x))
m_names <- sapply(strsplit(models, "/"), tail, 1)
names(m_list) <- sapply(strsplit(m_names, "_"), function(x) x[1])
i=1
source("code/bci-plot-functions.R")
# --- Example usage ---
# > models <- list.files("results/01_full/", full.names = TRUE)
out = plot_simulated_data(models, seed = 123, pdf_path = "figures/01-results.pdf")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# make splits ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
i = 1
# "pft-period3-9patches_S1_T0"
# i_name = "pft-period3-9patches"
m_list = c(m_list, rep("realdata",8))
names(m_list)[m_list == "realdata"] = c(
  "pft-period3-1patch-realdata",
  "pft-period3-9patches-realdata",
  "pft-period15-1patch-realdata",
  "pft-period15-9patches-realdata",
  "genus-period3-1patch-realdata",
  "genus-period3-9patches-realdata",
  "genus-period15-1patch-realdata",
  "genus-period15-9patches-realdata"
)
# genus-period3-9patches
i=8
for(i in 1:length(m_list)){
  FINN.seed(seed)
  i_name = names(m_list)[i]
  i_folder = gsub("-realdata","", i_name)
  m = m_list[[i]]

  Npatches = ifelse(grepl("1patch", i_name), 1L, 9L)
  Nepochs = ifelse(grepl("period3", i_name), 3L, 15L)

  period_length = Nepochs/3

  # read full raw data
  obs_dt = fread(paste0("data/Uholka/noSplits/", i_folder,"/obs_dt.csv"))
  env_dt = fread(paste0("data/Uholka/noSplits/", i_folder,"/env_dt.csv"))
  cohorts_dt = fread(paste0("data/Uholka/noSplits/", i_folder,"/initial_cohorts1985.csv"))

  # simulate or assign realdata
  if(!grepl("realdata", i_name)){
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = uniqueN(obs_dt$species))
    pred = m$simulate(env = env_dt, init_cohort = cohort1, patches = Npatches)
    pred_dt = pred$wide$site
    pred_dt[,reg := reg/0.1,]

    out_dir0 = paste0("data/Uholka/CVsplits-simdata/")
  }else if(grepl("realdata", i_name)){
    pred_dt = copy(obs_dt)

    out_dir0 = paste0("data/Uholka/CVsplits-realdata/")
  }

  # create temporal folds
  years <- c(1990, 1995, 2000, 2005, 2010, 2015)
  if(grepl("period3", i_name)){
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
      # if(grepl("period3", i_name)) env_start = env_start+1
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
  spatial_folds_dt = fread(paste0("data/Uholka/noSplits/",i_folder,"/spatial_folds_dt.csv"))[,.(siteID, fold = spatial_fold)]

  if(!dir.exists(paste0(out_dir0, i_folder, "/")))
    dir.create(paste0(out_dir0, i_folder), recursive = T)

  fwrite(temporal_folds_dt,paste0("data/Uholka/noSplits/",i_folder,"/temporal_folds_dt_applied.csv"))
  fwrite(spatial_folds_dt,paste0("data/Uholka/noSplits/",i_folder,"/spatial_folds_dt_applied.csv"))

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
      initial_cohort_train = fread(paste0("data/Uholka/noSplits/",i_folder,"/initial_cohorts",train_cohort_init_year,".csv"))
      initial_cohort_test = fread(paste0("data/Uholka/noSplits/",i_folder,"/initial_cohorts",test_cohort_init_year,".csv"))

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


