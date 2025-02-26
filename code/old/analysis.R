library(data.table)
library(FINN)
library(torch)
library(forcats)
library(dplyr)
library(ggplot2)

folder = "20250123_171017"
# folder = "20250121_182554"
figure_folder = paste0("figures/",folder,"/")
if(!dir.exists(figure_folder)) dir.create(figure_folder)
lossvars_rates = c("growth", "mort", "reg")
lossvars_stand = c("ba", "trees")
all_lossvars = c(lossvars_stand, lossvars_rates)
# get all combinations of lossvars
lossvars_comb = lapply(1:length(lossvars_rates), function(x) unlist(combn(lossvars_rates, x, simplify = T)))
lossvars_comb = lapply(lossvars_comb, function(x) apply(x,2, function(y) paste0(y, collapse=".")))
lossvars_comb <- lapply(lossvars_comb, function(x) paste0(c("ba.","ba.trees."),rep(x,each = 2)))
cv_variants <- paste0(rep(c("S","T","ST"),each = length(unlist(lossvars_comb))),"_",unlist(lossvars_comb))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## read fitted model ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
datapath = paste0("/Users/yannekkaber/working-directory/FINN-BCI/processed_data/cluster-output/",folder,"/processed_data/calibration-data/BCI-stand-data/blocked/")
parpath = paste0("/Users/yannekkaber/working-directory/FINN-BCI/processed_data/cluster-output/",folder,"/")
files = list.files(paste0(parpath))
parfiles = files[grepl("parameters", files) & !grepl("1epoch-model", files) & !grepl("final_model", files) & !grepl("_test-parameters", files)]
mfiles = files[grepl("final_model", files)]
mtransfiles = files[grepl("1epoch-model", files)]

datalist = list()
# cv_variants = c("S","T","ST")
# loss_vars = c("all","norates")
# hybrid_variants = c("BCI-linear", "BCI-DNN", "BCI-GrowthDNN")
hybrid_variants = c("BCI-linear")
i_cv = cv_variants[1]
i_hybrid = hybrid_variants[1]
# cv_variants = "S"
# cv_variants = cv_variants[unlist(lapply(cv_variants, function(x) any(grepl(paste0(x, "\\.RDS$"), parfiles))))]
for(i_cv in cv_variants){
  # read calibration data
  i_datname = tstrsplit(i_cv, "_", fixed = TRUE)[[1]][1]
  datalist[[i_cv]][["train"]][["obs_dt"]] = fread(paste0(datapath,"obs_dt_train",i_datname,".csv"))
  datalist[[i_cv]][["train"]][["env_dt"]] = fread(paste0(datapath,"env_dt_train",i_datname,".csv"))
  datalist[[i_cv]][["train"]][["cohorts_dt"]] = fread(paste0(datapath,"initial_cohorts_train",i_datname,".csv"))

  # read test data
  datalist[[i_cv]][["test"]][["obs_dt"]] = fread(paste0(datapath,"obs_dt_test",i_datname,".csv"))
  datalist[[i_cv]][["test"]][["env_dt"]] = fread(paste0(datapath,"env_dt_test",i_datname,".csv"))
  datalist[[i_cv]][["test"]][["cohorts_dt"]] = fread(paste0(datapath,"initial_cohorts_test",i_datname,".csv"))

  # read calibration output
  for(i_hybrid in hybrid_variants){
    i_parspath = tail(sort(paste0(parpath,parfiles[grepl(paste0("[0-9]+", i_cv, "\\.RDS$"), parfiles) & grepl(i_hybrid, parfiles)])),1)
    datalist[[i_cv]][[i_hybrid]][["pars"]] = readRDS(i_parspath)
    i_mpath = tail(sort(paste0(parpath,mfiles[grepl(paste0("[0-9]+", i_cv, "\\.RDS$"), mfiles) & grepl(i_hybrid, mfiles)])),1)
    datalist[[i_cv]][[i_hybrid]][["m"]] = readRDS(i_mpath)
    i_mpath = tail(sort(paste0(parpath,mtransfiles[grepl(paste0("model", i_cv, "\\.RDS$"), mtransfiles) & grepl(i_hybrid, mtransfiles)])),1)
    datalist[[i_cv]][[i_hybrid]][["m"]] = readRDS(i_mpath)
    # i_mtrans = tail(sort(paste0(parpath,mtransfiles[grepl(paste0("model", i_cv,"_",i_lossvars, "\\.RDS$"), mtransfiles) & grepl(i_hybrid, mtransfiles)])),1)
    # datalist[[i_cv]][[i_hybrid]][["mtrans"]] = readRDS(paste0(parpath,mtransfiles[grepl(paste0("model"', i_cv,"_",i_lossvars, "\\.RDS$"), mtransfiles) & grepl(i_hybrid, mtransfiles)]))
  }
}

i_cv = cv_variants[1]
i_hybrid = hybrid_variants[1]

parlist_out = list()
out_dt = data.table()
i_train = "train"
for(i_train in c("train", "test")){
  for(i_cv in cv_variants){
    for(i_hybrid in hybrid_variants){
      obs_dt = datalist[[i_cv]][["train"]]$obs_dt
      env_dt = datalist[[i_cv]][["train"]]$env_dt
      cohorts_dt = datalist[[i_cv]][["train"]]$cohorts_dt

      obs_dt_test = datalist[[i_cv]][[i_train]]$obs_dt
      env_dt_test = datalist[[i_cv]][[i_train]]$env_dt
      cohorts_dt_test = datalist[[i_cv]][[i_train]]$cohorts_dt

      pars = datalist[[i_cv]][[i_hybrid]]$pars
      m = datalist[[i_cv]][[i_hybrid]]$m

      # all_pars[[1]]
      # create plots for alll elements in pars
      all_valid_epochs = c()
      for(j in 1:length(pars[[1]])){
        checkNA = sapply(1:length(pars), function(i) any(is.na(pars[[i]][[j]])))
        if(all(!checkNA)) valid_epochs = length(pars) else valid_epochs = min(which(checkNA))-100
        all_valid_epochs = c(all_valid_epochs, valid_epochs)
      }
      valid_epoch = all_valid_epochs[which.min(all_valid_epochs)]
      pdf(paste0(figure_folder,"traceplots_", i_cv,"_", i_hybrid, "_", i_train, ".pdf"), width = 8, height = 6)
      par(mfrow=c(2,2))
      for(j in names(pars[[1]])){
        for(k in 1:dim(pars[[1]][[j]])[2]){
          sapply(1:valid_epochs, function(i) pars[[i]][[j]][,k]) %>% t() %>% matplot(type = "l", main = paste0(j, " ", k))
        }
      }
      dev.off()



      # all_pars <- list()
      # i=1
      # for(j in 1:length(pars)){
      fitted_pars = pars[[valid_epoch]]
      par_names = names(fitted_pars)
      parlist <- list()
      # i = par_names[5]
      for(i in par_names){
        if(is.null(m$model$speciesPars_ranges[[i]])){
          parlist[[i]] = as.matrix(fitted_pars[[i]])
        }else{
          parlist[[i]] = as.matrix(m$model$getPars(fitted_pars[[i]], m$model$speciesPars_ranges[[i]]))
        }
      }

      parlist_out[[i_cv]][[i_hybrid]][[i_train]] = parlist
        # all_pars <- c(all_pars,parlist)
      # }
      # get all parametrs for each timestep in length(pars)

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## read and format parameters ####
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

      # pars <- readRDS("data/calibration-output/BCI_parameters_01ha_9_9_wo_T_v2.RDS")
      # stand_dt = fread("cluster/BCI-calibration/20250120_233316/processed_data/calibration-data/BCI-stand-data/stand_dt.csv")

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## continue simulation with fitted parameters ####
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      Nspecies = nrow(parlist$parComp)
      # cohorts_dt_in <- cohorts_dt[,.(siteID, patchID = 1, cohortID, species, dbh = round(dbh_cm,4), trees)]
      cohorts_dt_in <- cohorts_dt_test[,.(siteID, patchID = 1, cohortID, species, dbh = round(dbh_cm,4), trees)]
      cohort1 <- FINN::CohortMat$new(obs_df = cohorts_dt_in, sp = Nspecies)

      if(i_hybrid == "BCI-GrowthDNN"){
        predictions =
          simulateForest(env = env_dt_test,
                         sp = Nspecies,
                         init = cohort1,
                         patch_size = 0.1,
                         patches = 1, speciesPars_ranges = m$model$speciesPars_ranges,
                         competitionProcess = createProcess(~., func = competition,initSpecies = parlist$parComp),
                         growthProcess = m$growthProcess,
                         mortalityProcess = createProcess(~., func = mortality, initEnv = parlist$nnMort.0.weight,initSpecies = parlist$parMort),
                         regenerationProcess = createProcess(~., func = regeneration, initEnv = parlist$nnReg.0.weight,initSpecies = as.vector(parlist$parReg)),
                         device = "cpu")
      }else{
        predictions =
          simulateForest(env = env_dt_test,
                         sp = Nspecies,
                         init = cohort1,
                         patch_size = 0.1,
                         patches = 1, speciesPars_ranges = m$model$speciesPars_ranges,
                         competitionProcess = createProcess(~., func = competition,initSpecies = parlist$parComp),
                         growthProcess = createProcess(~., func = growth, initEnv = parlist$nnGrowth.0.weight,initSpecies = parlist$parGrowth),
                         mortalityProcess = createProcess(~., func = mortality, initEnv = parlist$nnMort.0.weight,initSpecies = parlist$parMort),
                         regenerationProcess = createProcess(~., func = regeneration, initEnv = parlist$nnReg.0.weight,initSpecies = as.vector(parlist$parReg)),
                         device = "cpu")
      }
      pred_dt_test = predictions$wide$site

      obs_dt_test$year = as.integer(factor(obs_dt_test$year))
      comp_dt = merge(pred_dt_test, obs_dt_test, by = c("siteID", "year", "species"), suffixes = c(".pred", ".obs"))

      comp_allspecies_dt <- comp_dt[,.(
        ba.obs = sum(ba.obs)/uniqueN(siteID),
        ba.pred = sum(ba.pred)/uniqueN(siteID),
        trees.obs = sum(trees.obs)/uniqueN(siteID),
        trees.pred = sum(trees.pred)/uniqueN(siteID),
        dbh.obs = mean(dbh.obs, na.rm = TRUE),
        dbh.pred = mean(dbh.pred, na.rm = TRUE),
        reg.obs = mean(reg.obs),
        reg.pred = mean(r_mean_ha),
        mort.obs = mean(mort.obs, na.rm = TRUE),
        mort.pred = mean(mort.pred, na.rm = TRUE),
        growth.obs = mean(growth.obs, na.rm = T),
        growth.pred = mean(growth.pred, na.rm = T)
      ), by = .(year, species)]
      plot(comp_allspecies_dt$growth.obs, comp_allspecies_dt$growth.pred)
      # comp_allspecies_dt[,":="(
      #   ba.diff = ba.obs-ba.pred,
      #   trees.diff = trees.obs-trees.pred,
      #   dbh.diff = dbh.obs-dbh.pred,
      #   reg.diff = reg.obs-reg.pred,
      #   growth.diff = growth.obs-growth.pred,
      #   mort.diff = mort.obs-mort.pred
      # ),]

      out_dt = rbind(
        out_dt,
        data.table(
          cv = i_cv,
          hybrid = i_hybrid,
          train = i_train,
          epoch = valid_epoch,
          comp_allspecies_dt
        )
      )
      cat("finished:", i_cv, i_hybrid, i_train, "\n")
    }
  }
}

plot(growth.obs~growth.pred,data = out_dt)
out_dt_long = melt.data.table(out_dt, id.vars = c("cv", "hybrid", "train", "epoch", "year", "species"))
# out_dt_long = melt.data.table(out_dt, id.vars = c("cv", "hybrid", "train", "epoch", "siteID", "year", "species"))
# out_dt_long = melt.data.table(out_dt_all, id.vars = c("cv", "hybrid", "train", "epoch", "siteID", "year", "species"))
out_dt_long[, c("variable", "sim_obs") := tstrsplit(variable, "\\.", fixed = FALSE)]
#dcast to have value.obs and value.pred
out_dt_long = dcast(out_dt_long, cv + hybrid + train + epoch + year + species + variable ~ sim_obs, value.var = "value")

plot(obs~pred, out_dt_long[variable == "growth",])

r2 = function(pred, obs) {
  SS_res <- sum((obs - pred)^2)  # Sum of squared residuals
  SS_tot <- sum((obs - mean(obs))^2) # Total sum of squares
  R_squared <- 1 - (SS_res / SS_tot)
  R_squared[R_squared<=0] = 0
  return(R_squared)
}
# out_dt_long = copy(out_dt_long_orig)
# out_dt_long_orig = copy(out_dt_long)
# out_dt_long = out_dt_long[,.(
#   obs = mean(obs, na.rm = T),
#   pred = mean(pred, na.rm = T),
# ), by = .(cv, hybrid, train, epoch, year, species, variable)]

# plot(obs~pred, out_dt_long[variable == "growth",])

out_dt_long[!is.na(obs),":="(
  rmse = sqrt(mean((obs-pred)^2, na.rm = T)),
  spearmanr = cor(obs, pred, method = "spearman"),
  r2 = r2(pred, obs)
), by = .(cv, hybrid, train, epoch, variable)]

out_dt_long[, c("cv2", "lossvars") := tstrsplit(cv, "_", fixed = TRUE)]

p = ggplot(out_dt_long, aes(y = (as.factor(lossvars)), x = (as.factor(variable))))+
  geom_tile(aes(fill = spearmanr))+
  facet_grid(cv2~train, scales = "free")+
  # add epochs as text to each tile
  geom_text(aes(label = round(spearmanr,2)), size = 3, color = "black")+
  # pick heat color scale with earthy colors
  # scale_fill_viridis_c(option = "D", direction = 1,limits = c(-1,1))
  # choose appropriate scale for values -1 to 1
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))
  # set color range from 0 to 1
p
ggsave(paste0(figure_folder,"spearmanr_heatmap.pdf"), p, width = 8, height = 8)
# tstrsplit cv

# ggplot(out_dt_long, aes(x = obs, y = pred, color = factor(species)))+
#   geom_point()+
#   facet_wrap(cv~variable, scales = "free")+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_minimal()
#
# # plot spearmanr
# ggplot(out_dt_long, aes(x = cv, y = r2, fill = lossvars))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   facet_grid(variable~train, scales = "fixed")+
#   theme_minimal()
#
# ggplot(out_dt_long, aes(x = cv, y = rmse, fill = lossvars))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   facet_grid(variable~train, scales = "free")+
#   theme_minimal()
#
# ggplot(out_dt_long[train == "test"], aes(x = cv, y = spearmanr, fill = lossvars))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   facet_wrap(~variable, scales = "free")+
#   theme_minimal()
#
# ggplot(out_dt_long[train == "test"], aes(x = cv, y = rmse))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   facet_wrap(~variable, scales = "free")+
#   theme_minimal()

pdf(paste0(figure_folder,"correlation_plots.pdf"), width = 8, height = 6)
par(mfrow=c(2,3))
i_cv = "S_ba.trees.growth.mort.reg"
i_hybrid = "BCI-linear"
for(i_cv in cv_variants){
  for(i_hybrid in hybrid_variants){
    for(i_train in c("train", "test")){
      # for(i_lossvars in c("all","norates")){
        for(i in c("reg", "mort", "growth", "ba", "trees", "dbh")) {
          # plot
          obs = out_dt_long[cv == i_cv & hybrid == i_hybrid & train == i_train & variable == i]$obs
          pred = out_dt_long[cv == i_cv & hybrid == i_hybrid & train == i_train & variable == i]$pred
          species = out_dt_long[cv == i_cv & hybrid == i_hybrid & train == i_train & variable == i]$species
          noNAs = !is.na(obs) & !is.na(pred)
          obs = obs[noNAs]
          pred = pred[noNAs]
          plot(
            obs, pred, col = species,
            ylim = c(0,max(c(obs, pred), na.rm = T)), xlim = c(0,max(c(obs, pred), na.rm = T)),
            main = paste(i,"\n", i_cv,"\n", i_hybrid, i_train), xlab = "observation", ylab = "prediction")
          # add spearman correlation to the plot as text make sure to place it in the plot on the top left
          spearmanr = round(cor(obs,pred, method = "spearman"), 2)
          rmse = round(sqrt(mean((obs-pred)^2, na.rm = T)), 2)
          text(0.2*max(obs, pred, na.rm = T), max(obs, pred, na.rm = T)-0.2*max(obs, pred, na.rm = T), paste0("R=", spearmanr,"\nrmse=", rmse), cex = 1.5, font = 2, col = "red")
        }
      }
    # }
  }
}
dev.off()





