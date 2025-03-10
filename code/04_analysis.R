library(data.table)
library(ggplot2)
library(ggh4x)
library(FINN)
library(torch)

result_folders = c("02_realdata", "02_simulated")
# get current git hash
git_hash <- system("git rev-parse --short HEAD", intern = TRUE)

if(!file.exists(paste0("results/pred_list_",git_hash,".rds"))){
  pred_list <- list()
  i_result_folder = result_folders[1]
  for(i_result_folder in result_folders){
    cat("\nstarting ", i_result_folder,"\n")
    if(i_result_folder == "02_realdata"){
      split_location = "CVsplits-realdata/"
      simreal = "real"
    }else if(i_result_folder == "02_simulated"){
      split_location = "CVsplits-simdata/"
      simreal = "sim"
    }

    files_list = list.files(path = paste0("results/", i_result_folder,"/"), pattern = "*.pt", full.names = TRUE)

    # files_list = files_list[grepl("T0", files_list)]

    # Read all files
    m_list <- lapply(files_list, torch::torch_load)
    names(m_list) <- gsub(".pt","", basename(files_list))

    # i=4
    i = 1
    for(i in 1:length(m_list)){
      m = m_list[[i]]
      name = names(m_list)[i]
      folder = strsplit(name,"_")[[1]][1]
      cat("\r", i, "of", length(m_list), "start model:", name, "                                          ")
      if(!grepl("species", folder)) dataset = "BCI"
      if(grepl("species", folder)) dataset = "Uholka"
      files_dir = paste0("data/",dataset,"/",split_location,folder,"/")
      cv_S = tstrsplit(name, "_", fixed = TRUE)[[2]][1]
      cv_T = tstrsplit(name, "_", fixed = TRUE)[[3]][1]
      cv = paste0(cv_S, cv_T)
      response = tstrsplit(name, "_", fixed = TRUE)[[4]][1]
      obs_dt_train <- fread(paste0(files_dir, "obs_dt_",cv_S,"_",cv_T,"_train.csv"))
      obs_dt_test <- fread(paste0(files_dir, "obs_dt_",cv_S,"_",cv_T,"_test.csv"))
      env_dt_train <- fread(paste0(files_dir, "env_dt_",cv_S,"_",cv_T,"_train.csv"))
      env_dt_train <- env_dt_train[,-c("splitID", "holdout", "siteID_holdout")]
      env_dt_test <- fread(paste0(files_dir, "env_dt_",cv_S,"_",cv_T,"_test.csv"))
      env_dt_test <- env_dt_test[,-c("splitID", "holdout", "siteID_holdout")]
      cohorts_dt_train <- fread(paste0(files_dir, "initial_cohorts_",cv_S,"_",cv_T,"_train.csv"))
      cohorts_dt_test <- fread(paste0(files_dir, "initial_cohorts_",cv_S,"_",cv_T,"_test.csv"))

      Nspecies = max(obs_dt_train$species)
      Npatches = max(cohorts_dt_train$patchID)

      # cohorts_dt_train <- cohorts_dt_train[,.(siteID, patchID, cohortID, species, dbh = round(dbh_cm,4), trees)]
      cohorts_train = FINN::CohortMat(cohorts_dt_train, sp = Nspecies)
      # cohorts_dt_test <- cohorts_dt_test[,.(siteID, patchID, cohortID, species, dbh = round(dbh_cm,4), trees)]
      cohorts_test = FINN::CohortMat(cohorts_dt_test, sp = Nspecies)

      pred_train = m$simulate(env = env_dt_train, init_cohort = cohorts_train, patches = Npatches)
      pred_test = m$simulate(env = env_dt_test, init_cohort =  cohorts_test, patches = Npatches)

      pred_list[[dataset]][[simreal]][[name]] <- list(
        pred = list(
          train = pred_train,
          test = pred_test
        ),
        obs = list(
          train = obs_dt_train,
          test = obs_dt_test
          )
        )
    }
  }
  saveRDS(pred_list, paste0("results/pred_list_",git_hash,".rds"))
}else{
  pred_list <- readRDS(paste0("results/pred_list_",git_hash,".rds"))
}

# pred_l = pred_list[[1]]
# i=1
all_dt_list <- list()
for(dataset in unique(names(pred_list))){
  pred_dataset = pred_list[[dataset]]
  for(simreal in names(pred_dataset)){
    cat("\nstart", dataset, simreal, "\n")
    pred_dataset_simreal = pred_dataset[[simreal]]
    for(i in 1:length(pred_dataset_simreal)){
      pred_l = pred_dataset_simreal[[i]]
      name = names(pred_dataset_simreal)[i]

      if (i %% 100 == 0) {
        cat("\r", i, "of", length(pred_dataset_simreal), "combining", name, "                                          ")
      }

      pred_train_temp = pred_l$pred$train$long$site
      pred_train_temp <- pred_train_temp[variable != "reg"]
      pred_train_temp[variable == "r_mean_ha", variable := "reg",]
      pre_test_temp = pred_l$pred$test$long$site
      pre_test_temp <- pre_test_temp[variable != "reg"]
      pre_test_temp[variable == "r_mean_ha", variable := "reg",]

      dt <-
        rbindlist(
          list(
            pred_train = data.table(pred_train_temp, pred_obs = "pred", test_train = "train", simreal = simreal, dataset = dataset),
            pred_test = data.table(pre_test_temp, pred_obs = "pred", test_train = "test", simreal = simreal, dataset = dataset),
            obs_dt_train =
              data.table(
                melt(pred_l$obs$train[,.(siteID, year, species, ba, dbh, trees, growth, mort, reg)], id.vars = c("siteID", "species", "year")),
                pred_obs = "obs", test_train = "train", simreal = simreal, dataset = dataset),
            obs_dt_test =
              data.table(
                melt(pred_l$obs$test[,.(siteID, year, species, ba, dbh, trees, growth, mort, reg)], id.vars = c("siteID", "species", "year")),
                pred_obs = "obs", test_train = "test", simreal = simreal, dataset = dataset)
          ),
          use.names=TRUE
        )
      parts <- strsplit(name, "_")[[1]]
      dt[, scale := parts[1]]
      dt[, cv := paste0(parts[2], parts[3])]
      dt[, response := parts[4]]
      # all_dt <- rbind(all_dt, dt)
      all_dt_list[[length(all_dt_list) + 1]] = dt
    }
  }
}
all_dt <- rbindlist(all_dt_list)
all_dt[pred_obs == "pred" & is.na(value) & variable == "reg", value := 0,]


# r2 = function(pred, obs, na.rm = T) {
#   SS_res <- sum((obs - pred)^2, na.rm = na.rm)  # Sum of squared residuals
#   SS_tot <- sum((obs - mean(obs))^2, na.rm = na.rm) # Total sum of squares
#   R_squared <- 1 - (SS_res / SS_tot)
#   R_squared[R_squared<=0] = 0
#   return(R_squared)
# }

# dcast by pred_obs
pred_dt1 <- dcast(all_dt, siteID+year+species+test_train+simreal+dataset+scale+cv+response~pred_obs+variable, value.var = "value", sep = ".")
pred_dt1[,pred.trees_before :=  data.table::shift(pred.trees), by = .(siteID,species,test_train,simreal,dataset,scale,cv,response)]
pred_dt1[,obs.trees_before := data.table::shift(obs.trees), by = .(siteID,species,test_train,simreal,dataset,scale,cv,response)]

pred_dt1[(pred.trees_before == 0) & pred.trees == 0, ":="(pred.mort = NA, pred.growth = NA),]
pred_dt1[(obs.trees_before == 0) & obs.trees == 0, ":="(obs.mort = NA, obs.growth = NA),]
pred_dt1[obs.trees == 0, obs.dbh := NA]
pred_dt1[pred.trees == 0, pred.dbh := NA]

pred_dt2 <- melt(pred_dt1[,-c("obs.trees_before","pred.trees_before")], id.vars = c("siteID","year","species","test_train","simreal","dataset","scale","cv","response"))
pred_dt2[,pred_obs := tstrsplit(variable, ".", fixed = T)[[1]],]
pred_dt2[,variable := tstrsplit(variable, ".", fixed = T)[[2]],]
pred_dt <- dcast(pred_dt2, siteID+year+species+test_train+simreal+dataset+scale+cv+response+variable~pred_obs, value.var = "value", sep = ".")

pred_dt <- pred_dt[!is.na(obs),]
pred_dt <- pred_dt[cv != "S0T0"]

# pred_dt[grepl("T0", cv), spatial_holdout := T,]
pred_dt[!grepl("S0", cv), spatial_holdout := T,]
pred_dt[grepl("S0", cv), spatial_holdout := F,]
pred_dt[!grepl("T0", cv), temporal_holdout := T,]
pred_dt[grepl("T0", cv), temporal_holdout := F,]

cors_dt1_all =
  pred_dt[!is.na(obs) & !is.na(pred),.(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    spearmans_r = cor(pred, obs, method = "spearman"),
    # pearson_r = cor(pred, obs, method = "pearson"),
    obs_center = sum(range(obs, na.rm = T))/2,
    pred_center = sum(range(pred, na.rm = T)/2),
    N = .N
  ), by = .(variable, cv, spatial_holdout, temporal_holdout, test_train, scale, response, simreal, dataset)]

cors_dt1_all[is.infinite(pearson_r)]
cors_dt1_all[is.na(pearson_r)]
cors_dt1_all[is.na(spearmans_r)]

cors_dt1_species =
  pred_dt[!is.na(obs) & !is.na(pred),.(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    spearmans_r = cor(pred, obs, method = "spearman"),
    obs_center = sum(range(obs, na.rm = T))/2,
    pred_center = sum(range(pred, na.rm = T)/2),
    N = .N
  ), by = .(variable, cv, spatial_holdout, temporal_holdout, test_train, scale, response, simreal, dataset, species)]

cors_dt <-
  rbind(
    cors_dt1_all[,.(
    rmse = mean(rmse, na.rm = T),
    spearmans_r = mean(spearmans_r, na.rm = T),
    obs_center = mean(obs_center, na.rm = T),
    pred_center = mean(pred_center, na.rm = T),
    N = sum(N),
    species_cor = "overall"
    ), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset)],
    cors_dt1_species[,.(
    rmse = mean(rmse, na.rm = T),
    spearmans_r = mean(spearmans_r, na.rm = T),
    obs_center = mean(obs_center, na.rm = T),
    pred_center = mean(pred_center, na.rm = T),
    N = sum(N),
    species_cor = "species"
    ), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset)]
  )

cors_dt_onlyspecies <- cors_dt1_species[,.(
  rmse = mean(rmse, na.rm = T),
  spearmans_r = mean(spearmans_r, na.rm = T),
  obs_center = mean(obs_center, na.rm = T),
  pred_center = mean(pred_center, na.rm = T),
  N = sum(N)
), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset, species)]

# count number of "." in response
cors_dt[, N_dots := stringr::str_count(response, "\\.")]

cors_dt[, ylabels := paste0(response, " (",simreal,")"),]
cors_dt[, N_dots := stringr::str_count(ylabels, "\\.")]
cors_dt[simreal == "real", N_dots := 50,]
cors_dt[, ylabels := forcats::fct_reorder(ylabels, N_dots),]

cors_dt[,variable := factor(variable, levels = c("ba","trees", "dbh", "growth", "mort", "reg"), labels = c("ba","trees", "dbh", "growth", "mort", "reg")),]

for(i_species in c("overall", "species")){
  for(i_dataset in c("BCI", "Uholka")){
    # p_dat_s = cors_dt[spatial_holdout == T & temporal_holdout == F & dataset == i_dataset]
    p_dat_s = cors_dt[spatial_holdout == T & temporal_holdout == F & dataset == i_dataset & species_cor == i_species]
    p_dat_s$holdout = "spatial"
    # p_dat_s[scale == "species-period3-1patch" & variable == "dbh" & simreal == "real" & ylabels == "ba.trees.dbh.growth.mort.reg (real)" & holdout == "spatial" & test_train == "test"]
    # p_dat_t = cors_dt[temporal_holdout == T & spatial_holdout == F & dataset == i_dataset]
    p_dat_t = cors_dt[temporal_holdout == T & spatial_holdout == F & dataset == i_dataset & species_cor == i_species]
    p_dat_t$holdout = "temporal"
    p_dat = rbind(p_dat_s, p_dat_t)
    p_all = ggplot(
      p_dat,
      aes(y = forcats::fct_reorder(ylabels, N_dots), x = (as.factor(variable)))
      )+
      geom_tile(aes(fill = spearmans_r))+
      facet_grid(gsub("-","\n",scale)~paste0(holdout,"\n",test_train), scales = "fixed")+
      geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
      # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
      scale_fill_gradientn(colors = c("blue", "white", "red"), limits = c(0,1))+
      ggtitle(paste(unique(p_dat_s$dataset)))+
      theme_classic()+
      theme(legend.position = "none")+
      # rotate facet labels
      theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
        )+
      # rotate x labels
      xlab("predicted variable")+
      ylab("response variables")
    p_all
    ggsave(paste0("figures/fits_",i_dataset,"_",i_species,".png"), p_all, width = 12, height = 10)
  }
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## compare parameters ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
par_list <- list()
i_result_folder = result_folders[1]
for(i_result_folder in result_folders){
  cat("\nstarting ", i_result_folder,"\n")
  if(i_result_folder == "02_realdata"){
    split_location = "CVsplits-realdata/"
    simreal = "real"
  }else if(i_result_folder == "02_simulated"){
    split_location = "CVsplits-simdata/"
    simreal = "sim"
  }
  files_list = list.files(path = paste0("results/", i_result_folder,"/"), pattern = "*.pt", full.names = TRUE)
  # Read all files
  m_list <- lapply(files_list, torch::torch_load)
  names(m_list) <- gsub(".pt","", basename(files_list))
  i = 1
  for(i in 1:length(m_list)){
    m = m_list[[i]]
    name = names(m_list)[i]
    cat("\r", i, "of", length(m_list), "start model:", name, "                                          ")
    folder = strsplit(name,"_")[[1]][1]
    # get full true model
    true_model = torch_load(paste0("results/01_full/", folder, "_full.pt"))
    par_list[[paste0(name,"_",simreal)]] <- list(
      true = true_model$parameters_r,
      fit = m$parameters_r
    )
  }
}
j = 1
pars_dt <- data.table()
cat("\n")
for(j in 1:length(par_list)){
  cat("\r", j, "of", length(par_list))
  j_model = names(par_list)[j]
  pars_fit = par_list[[j_model]]$fit
  pars_true = par_list[[j_model]]$true
  i = 1
  for(i in 1:length(pars_true)){
    i_par = names(pars_true)[i]
    dt_true = data.table(pars_true[[i_par]])
    dt_fit = data.table(pars_fit[[i_par]])
    colnames(dt_true) = paste0(i_par, 1:ncol(dt_true))
    colnames(dt_fit) = paste0(i_par, 1:ncol(dt_fit))
    if(nrow(dt_true) > 1) dt_true$species = 1:nrow(dt_true) else dt_true$species = 0
    if(nrow(dt_fit) > 1) dt_fit$species = 1:nrow(dt_fit) else dt_fit$species = 0
    melt_dt_true = melt(dt_true, id.vars = "species")
    melt_dt_fit = melt(dt_fit, id.vars = "species")
    pars_dt = rbind(
      pars_dt,
      rbind(
        data.table(melt_dt_true, par = i_par, model = j_model, truefit = "true"),
        data.table(melt_dt_fit, par = i_par, model = j_model, truefit = "fit")
        )
      )
  }
}

pars_dt2 <- dcast(pars_dt, species+par+model+variable~truefit, value.var = "value")

pars_dt2[, scale := tstrsplit(model,"_",fixed = T)[[1]],]
pars_dt2[, cv := paste0(tstrsplit(model,"_",fixed = T)[[2]],tstrsplit(model,"_",fixed = T)[[3]]),]
pars_dt2[, response := tstrsplit(model,"_",fixed = T)[[4]],]
pars_dt2[, simreal := tstrsplit(model,"_",fixed = T)[[5]],]
pars_dt2[grepl("species", model), dataset := "Uholka",]
pars_dt2[!grepl("species", model), dataset := "BCI",]


pars_cor_dt <-
pars_dt2[,.(
  var = var(true - fit, na.rm = T),
  bias = mean(true - fit, na.rm = T)
), by = .(par, scale, cv, response, dataset, simreal)]

pars_cor_dt[grepl("S0", cv), S_test_train := "train",]
pars_cor_dt[!grepl("S0", cv), S_test_train := "test",]
pars_cor_dt[grepl("T0", cv), T_test_train := "train",]
pars_cor_dt[!grepl("T0", cv), T_test_train := "test",]

pars_cor_dt2 <- pars_cor_dt[,.(
  var = mean(var,na.rm = T),
  bias = mean(bias, na.rm = T)
), by = .(par, scale,S_test_train, T_test_train, response, dataset, simreal)]

pars_cor_dt2[, N_dots := stringr::str_count(response, "\\.")]

for(i_dataset in c("BCI", "Uholka")){
  p_dat_s = pars_cor_dt2[T_test_train == "train" & dataset == i_dataset & simreal == "sim"]
  p_dat_t = pars_cor_dt2[S_test_train == "train" & dataset == i_dataset & simreal == "sim"]
  p_s_bias <- ggplot(p_dat_s, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
    geom_tile(aes(fill = bias))+
    # geom_point(aes(size = bias, color = bias))+
    facet_grid(scale~S_test_train+simreal, scales = "fixed")+
    geom_text(aes(label = round(bias,2)), size = 3, color = "black")+
    # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    ggtitle(paste0(i_dataset, " spatial holdout"))+
    theme_classic()+
    # theme(legend.position = "none")+
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # change to multicolor scale that is easier to read
    scale_fill_gradient2()+
    xlab("parameter")+
    ylab("response variables")

  p_t_bias <- ggplot(p_dat_t, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
    geom_tile(aes(fill = bias))+
    # geom_point(aes(size = bias, color = bias))+
    facet_grid(scale~T_test_train+simreal, scales = "fixed")+
    geom_text(aes(label = round(bias,2)), size = 3, color = "black")+
    # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    ggtitle(paste0(i_dataset, " temporal holdout"))+
    theme_classic()+
    # theme(legend.position = "none")+
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # change to multicolor scale that is easier to read
    scale_fill_gradient2()+
    xlab("parameter")+
    ylab("response variables")

  p_all_bias = gridExtra::grid.arrange(p_s_bias, p_t_bias, ncol = 2)

  p_s_var <- ggplot(p_dat_s, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
    geom_tile(aes(fill = var))+
    # geom_point(aes(size = var, color = var))+
    facet_grid(scale~S_test_train+simreal, scales = "fixed")+
    geom_text(aes(label = round(var,2)), size = 3, color = "black")+
    # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    ggtitle(paste0(i_dataset, " spatial holdout"))+
    theme_classic()+
    # theme(legend.position = "none")+
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # change to multicolor scale that is easier to read
    scale_fill_viridis_c(option = "C")+
    xlab("parameter")+
    ylab("response variables")

  p_t_var <- ggplot(p_dat_t, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
    geom_tile(aes(fill = var))+
    # geom_point(aes(size = var, color = var))+
    facet_grid(scale~T_test_train+simreal, scales = "fixed")+
    geom_text(aes(label = round(var,2)), size = 3, color = "black")+
    # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    ggtitle(paste0(i_dataset, " temporal holdout"))+
    theme_classic()+
    # theme(legend.position = "none")+
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # change to multicolor scale that is easier to read
    scale_fill_viridis_c(option = "C")+
    xlab("parameter")+
    ylab("response variables")

  p_all_var = gridExtra::grid.arrange(p_s_var, p_t_var, ncol = 2)

  ggsave(paste0("figures/parameters_var_",i_dataset,".png"), p_all_var, width = 30, height = 10)
  ggsave(paste0("figures/parameters_bias_",i_dataset,".png"), p_all_bias, width = 30, height = 10)
}
