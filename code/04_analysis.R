library(data.table)
library(ggplot2)
library(ggh4x)
library(FINN)
library(torch)
source("code/plot-functions.R")

result_folders = c("02_realdata", "02_simulated", "02_realdata_hybrid")

all_models <- c(
    list.files("results/02_realdata/", full.names = T, recursive = T),
    list.files("results/02_realdata_hybrid/", full.names = T, recursive = T)
)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## create pdf for S0_T0 variants####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

noCV_models = all_models[grepl("S0_T0", all_models)]
pdf("figures/02-results_noCV.pdf", width = 15, height = 10)
all_cors_dt <- data.table()
all_cors_dt_plot <- data.table()
i = noCV_models[1]
for(i in noCV_models){
  cat("\nstarting", i,"\n", which(i == noCV_models), "of", length(noCV_models))
  pred_dt = build_model_dt(i)
  all_dt = pred_dt[[1]]
  name = names(pred_dt)[1]
  period_length = unique(all_dt[!is.na(period_length)]$period_length)
  all_dt <- all_dt[,-c("period_length")]
  # aggregate species by pft
  if(grepl("genus",name)){
    # species_dt <- fread("data/BCI/data-cleaning/genus/species_assigned.csv")
    # all_dt2 <- merge(all_dt, unique(species_dt[,.(species = speciesID, pft = PFT_2axes)]), by = "species", all.x = T, allow.cartesian = T)
    # all_dt2[,species := pft,]
    # all_dt2 <- all_dt2[,-c("pft")]
    # all_dt2[,.(
    #   ba.pred = sum(ba.pred, na.rm = T),
    #   ba.obs = sum(ba.obs, na.rm = T),
    #   trees.pred = sum(trees.pred, na.rm = T),
    #   trees.obs = sum(trees.obs, na.rm = T),
    #   dbh.pred = sum(dbh.pred*trees.pred, na.rm = T)/sum(trees.pred, na.rm = T),
    #   dbh.obs = sum(dbh.obs*trees.obs, na.rm = T)/sum(trees.obs, na.rm = T),
    #   growth.pred = sum(growth.pred*trees.pred, na.rm = T)/sum(trees.pred, na.rm = T),
    #   growth.obs = sum(growth.obs*trees.obs, na.rm = T)/sum(trees.obs, na.rm = T),
    #   mort.pred = sum(mort.pred*trees.pred, na.rm = T)/sum(trees.pred, na.rm = T),
    #   mort.obs = sum(mort.obs*trees.obs, na.rm = T)/sum(trees.obs, na.rm = T),
    #   reg.pred = sum(reg.pred*trees.pred, na.rm = T)/sum(trees.pred, na.rm = T),
    #   reg.obs = sum(reg.obs*trees.obs, na.rm = T)/sum(trees.obs, na.rm = T)
    # ), by = .(siteID,year,species,test_train)]
    top10_ba_species <- all_dt[,.(ba = sum(ba.obs,na.rm = T)), by = .(species)][order(-ba)][1:10]$species
    all_dt2 <- all_dt
  }else{
    top10_ba_species = unique(all_dt$species)
    all_dt2 <- all_dt
  }
  all_dt2 <- melt(all_dt2, id.vars = c("siteID","year","species","test_train"))
  all_dt2[grepl("pred",variable), pred_obs := "pred",]
  all_dt2[grepl("obs",variable), pred_obs := "obs",]
  all_dt2[, variable := gsub(".obs|.pred","",variable),]

  p_dat <- all_dt2[species %in% top10_ba_species,.(value = mean(value,na.rm = T)), by = .(year,species,test_train, pred_obs, variable)]
  p1 = ggplot(p_dat[test_train == "train"], aes(x = year, y = value))+
    geom_point(aes(shape = pred_obs, color = factor(species), size = pred_obs), alpha = 0.5)+
    geom_line(data = p_dat[test_train == "train" & !is.na(value)],  aes(linetype = pred_obs, color = factor(species)))+
    geom_line(aes(linetype = pred_obs, color = factor(species)))+
    facet_wrap(~variable,scales = "free_y")+
    ggtitle(name)+
    scale_size_manual(values = c(3,1))
  print(p1)

  cor_agg <- all_dt2[test_train == "train"]
  if(period_length > 1){
    cor_agg[,periodID := as.integer(as.factor(year))]
    cor_agg[is.na(value)]
    years = unique(cor_agg$year)
    periods = rep(seq(1, uniqueN(cor_agg$year)/period_length, by = 1), each = period_length)
    periods_year = rep(seq(period_length, max(cor_agg$year), by = period_length), each = period_length)
    cor_agg[,periods_year := as.integer(as.character(factor(year, levels = years, labels = periods_year)))]
    cor_agg[,periodID := as.integer(factor(year, levels = years, labels = periods))]
    cor_agg2 <-
      rbindlist(
        list(
        cor_agg[variable %in% c("ba", "dbh", "trees") & periods_year == year,.(value = mean(value,na.rm = T)), by= .(siteID,periodID,species,test_train, pred_obs, variable)],
        cor_agg[variable %in% c("reg"), .(value = sum(value, na.rm = T)), .(siteID,periodID,species,test_train, pred_obs,variable)],
        cor_agg[variable %in% c("mort","growth"), .(value = mean(value, na.rm = T)), .(siteID,periodID,species,test_train, pred_obs, variable)]
        )
      )
  }else{
    cor_agg2 <- cor_agg[, .(value = mean(value,na.rm = T)),.(siteID,periodID = year,species,test_train, pred_obs, variable)]
  }
  cor_dt_in <- dcast(cor_agg2,
    siteID+periodID+species+test_train+variable~pred_obs, value.var = "value")


  cor_dt_in_plot <- dcast(cor_agg2[,.(
    value = mean(value,na.rm = T)
    ), by = .(periodID,species,test_train,variable,pred_obs)],
    periodID+species+test_train+variable~pred_obs, value.var = "value")
  cor_dt <-
    cor_dt_in[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_species <- cor_dt_in[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable,species, test_train)]
  cor_dt_plot <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_plot_species <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable,species, test_train)]


  p2a = ggplot(cor_dt_in[species %in% top10_ba_species], aes(x = pred, y = obs, color = factor(species)))+
    geom_point(alpha = 0.1)+
    geom_smooth(method = "lm")+
    facet_wrap(~variable, scales = "free", ncol = 2)+
    ggtitle(paste0("per site: ",name))+
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")+
    theme_minimal()+
    theme(legend.position = "bottom")
  # print(p2a)
  cor_dt_species <- rbind(
    data.table(cor_dt, species = "overall"),
    cor_dt_species, fill = T
  )
  p2b = ggplot(cor_dt_species[species %in% top10_ba_species], aes(y = species, x = variable))+
    geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
    facet_wrap(~test_train)+
    theme_classic()+
    ggtitle(paste0("per site and species: ",name))
  # print(p2b)
  # arrange p2a and p2b
  gridExtra::grid.arrange(p2a, p2b, ncol = 2)
  p3a = ggplot(cor_dt_in_plot[species %in% top10_ba_species], aes(x = pred, y = obs, color = factor(species)))+
    geom_point(alpha = 0.8)+
    geom_smooth(method = "lm")+
    facet_wrap(~variable, scales = "free")+
    ggtitle(paste0("whole plot: ",name))+
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")+
    theme_minimal()+
    theme(legend.position = "bottom")
  # print(p3)
  cor_dt_plot_species <- rbind(
    data.table(cor_dt_plot, species = "overall"),
    cor_dt_plot_species, fill = T
  )
  p3b = ggplot(cor_dt_plot_species[species %in% top10_ba_species], aes(y = species, x = variable))+
    geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
    geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
    facet_wrap(~test_train)+
    theme_classic()+
    ggtitle(paste0("whole plot and species: ",name))
  # print(p3b)
  # arrange p3a and p3b
  gridExtra::grid.arrange(p3a, p3b, ncol = 2)
  ## per species

  cor_dt[,":="(name = basename(name), data = basename(dirname(name)))]
  cor_dt_plot[,":="(name = basename(name), data = basename(dirname(name)))]
  all_cors_dt <- rbind(all_cors_dt, cor_dt)
  all_cors_dt_plot <- rbind(all_cors_dt_plot, cor_dt_plot)
}
p4 = ggplot(all_cors_dt, aes(y = paste(data,name), x = variable))+
  geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
  facet_wrap(~test_train)+
  theme_classic()+
  ggtitle("per site")
print(p4)
p5 = ggplot(all_cors_dt_plot, aes(y = paste(data,name), x = variable))+
  geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
  facet_wrap(~test_train)+
  theme_classic()+
  ggtitle("whole plot")
print(p5)
dev.off()

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## correlations for all combinations ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
gc()
all_models <- c(
  list.files("results/02_simulated/", full.names = T, recursive = T),
  list.files("results/02_realdata/", full.names = T, recursive = T),
  list.files("results/02_realdata_hybrid/", full.names = T, recursive = T)
)

all_cors_dt <- data.table()
all_cors_dt_plot <- data.table()
all_cors_dt_species <- data.table()
all_cors_dt_plot_species <- data.table()
for(i in all_models){
  cat("\n", which(i == all_models), "of", length(all_models),"\nstarting", i, "\n")
  pred_dt = build_model_dt(i)
  all_dt = pred_dt[[1]]
  name = names(pred_dt)[1]
  if(grepl("period35",name) | grepl("period15",name)){
    period_length = 5
  }else{
    period_length = 1
  }
  all_dt <- all_dt[,-c("period_length")]

  all_dt <- melt(all_dt, id.vars = c("siteID","year","species","test_train"))
  all_dt[grepl("pred",variable), pred_obs := "pred",]
  all_dt[grepl("obs",variable), pred_obs := "obs",]
  all_dt[, variable := gsub(".obs|.pred","",variable),]

  cor_agg <- all_dt
  if(period_length > 1){
    cor_agg[,periodID := as.integer(as.factor(year))]
    cor_agg[is.na(value)]
    years = unique(cor_agg$year)
    periods = rep(seq(1, uniqueN(cor_agg$year)/period_length, by = 1), each = period_length)
    periods_year = rep(seq(period_length, max(cor_agg$year), by = period_length), each = period_length)
    cor_agg[,periods_year := as.integer(as.character(factor(year, levels = years, labels = periods_year)))]
    cor_agg[,periodID := as.integer(factor(year, levels = years, labels = periods))]
    cor_agg2 <-
      rbindlist(
        list(
          cor_agg[variable %in% c("ba", "dbh", "trees") & periods_year == year,.(value = mean(value,na.rm = T)), by= .(siteID,periodID,species,test_train, pred_obs, variable)],
          cor_agg[variable %in% c("reg"), .(value = sum(value, na.rm = T)), .(siteID,periodID,species,test_train, pred_obs,variable)],
          cor_agg[variable %in% c("mort","growth"), .(value = mean(value, na.rm = T)), .(siteID,periodID,species,test_train, pred_obs, variable)]
        )
      )
  }else{
    cor_agg2 <- cor_agg[, .(value = mean(value,na.rm = T)),.(siteID,periodID = year,species,test_train, pred_obs, variable)]
  }
  cor_dt_in <- dcast(cor_agg2,siteID+periodID+species+test_train+variable~pred_obs, value.var = "value")

  cor_dt_in_plot <- dcast(cor_agg2[,.(
    value = mean(value,na.rm = T)
    ),by = .(periodID,species,test_train,variable,pred_obs)],
    periodID+species+test_train+variable~pred_obs, value.var = "value")
  cor_dt <-
    cor_dt_in[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_species <- cor_dt_in[!is.na(obs) & !is.na(pred) & test_train == "train",.(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    spearmans_r = cor(pred, obs, method = "spearman"),
    obs_center = sum(range(obs, na.rm = T))/2,
    pred_center = sum(range(pred, na.rm = T)/2),
    N = .N
  ), by = .(variable,species, test_train)]
  cor_dt_plot <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_plot_species <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred) & test_train == "train",.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable,species, test_train)]

  cor_dt[,":="(name = basename(name), data = basename(dirname(name)))]
  cor_dt_plot[,":="(name = basename(name), data = basename(dirname(name)))]
  cor_dt_species[,":="(name = basename(name), data = basename(dirname(name)))]
  cor_dt_plot_species[,":="(name = basename(name), data = basename(dirname(name)))]
  all_cors_dt <- rbind(all_cors_dt, cor_dt)
  all_cors_dt_plot <- rbind(all_cors_dt_plot, cor_dt_plot)
  all_cors_dt_species <- rbind(all_cors_dt_species, cor_dt_species)
  all_cors_dt_plot_species <- rbind(all_cors_dt_plot_species, cor_dt_plot_species)
  rm(cor_dt, cor_dt_plot, cor_dt_species, cor_dt_plot_species, all_dt, pred_dt, cor_agg, cor_agg2)
  gc()
}
fwrite(all_cors_dt, "results/all_cors_dt.csv")
fwrite(all_cors_dt_plot, "results/all_cors_dt_plot.csv")
fwrite(all_cors_dt_species, "results/all_cors_dt_species.csv")
fwrite(all_cors_dt_plot_species, "results/all_cors_dt_plot_species.csv")

# p4 = ggplot(all_cors_dt, aes(y = paste(data,name), x = variable))+
#   geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
#   scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#   geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
#   facet_wrap(~test_train)+
#   theme_classic()+
#   ggtitle("per site")
# print(p4)
# p5 = ggplot(all_cors_dt_plot, aes(y = paste(data,name), x = variable))+
#   geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
#   scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#   geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
#   facet_wrap(~test_train)+
#   theme_classic()+
#   ggtitle("whole plot")
# print(p5)




#
# all_dt2 <- dcast(
#   all_dt,
#   siteID+year+species+test_train+variable~pred_obs, value.var = "value")
#
# all_dt_fullPlot <- all_dt2[,.(
#   obs = mean(obs, na.rm = T),
#   pred = mean(pred, na.rm = T)
# ), by = .(year, species, test_train, variable)]
#
# cors_dt1_all =
#   all_dt2[!is.na(obs) & !is.na(pred),.(
#     rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
#     spearmans_r = cor(pred, obs, method = "spearman"),
#     obs_center = sum(range(obs, na.rm = T))/2,
#     pred_center = sum(range(pred, na.rm = T)/2),
#     N = .N
#   ), by = .(variable, test_train)]
# cors_dt1_species =
#   all_dt2[!is.na(obs) & !is.na(pred),.(
#     rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
#     spearmans_r = cor(pred, obs, method = "spearman"),
#     obs_center = sum(range(obs, na.rm = T))/2,
#     pred_center = sum(range(pred, na.rm = T)/2),
#     N = .N
#   ), by = .(variable, test_train, species)]
#
# cors_dt1_comb <-
#   rbind(data.table(cors_dt1_all, species = "overall"), data.table(cors_dt1_species))
#
# ggplot(cors_dt1_comb[test_train == "test"], aes(y = variable, x = species))+
#   geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 3)+
#   scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#   geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
#   facet_wrap(~test_train)+
#   theme_classic()
#
# ggplot(cors_dt1_species, aes(x = spearmans_r^2, y = forcats::fct_reorder(factor(species), spearmans_r^2)))+
#   geom_bar(stat = "identity", aes(fill = test_train), position = position_dodge2())+
#   facet_wrap(~variable)+
#   theme_minimal()
#   # coord_cartesian(xlim = c(0,1))
#
#
#
# pred_dt <- pred_dt[!is.na(obs),]
# pred_dt <- pred_dt[cv != "S0T0"]
#
# pred_dt[!grepl("S0", cv), spatial_holdout := T,]
# pred_dt[grepl("S0", cv), spatial_holdout := F,]
# pred_dt[!grepl("T0", cv), temporal_holdout := T,]
# pred_dt[grepl("T0", cv), temporal_holdout := F,]
#
# cors_dt1_all =
#   pred_dt[!is.na(obs) & !is.na(pred),.(
#     rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
#     spearmans_r = cor(pred, obs, method = "spearman"),
#     # pearson_r = cor(pred, obs, method = "pearson"),
#     obs_center = sum(range(obs, na.rm = T))/2,
#     pred_center = sum(range(pred, na.rm = T)/2),
#     N = .N
#   ), by = .(variable, cv, spatial_holdout, temporal_holdout, test_train, scale, response, simreal, dataset)]
#
# cors_dt1_species =
#   pred_dt[!is.na(obs) & !is.na(pred),.(
#     rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
#     spearmans_r = cor(pred, obs, method = "spearman"),
#     obs_center = sum(range(obs, na.rm = T))/2,
#     pred_center = sum(range(pred, na.rm = T)/2),
#     N = .N
#   ), by = .(variable, cv, spatial_holdout, temporal_holdout, test_train, scale, response, simreal, dataset, species)]
#
# cors_dt <-
#   rbind(
#     cors_dt1_all[,.(
#     rmse = mean(rmse, na.rm = T),
#     spearmans_r = mean(spearmans_r, na.rm = T),
#     obs_center = mean(obs_center, na.rm = T),
#     pred_center = mean(pred_center, na.rm = T),
#     N = sum(N),
#     species_cor = "overall"
#     ), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset)],
#     cors_dt1_species[,.(
#     rmse = mean(rmse, na.rm = T),
#     spearmans_r = mean(spearmans_r, na.rm = T),
#     obs_center = mean(obs_center, na.rm = T),
#     pred_center = mean(pred_center, na.rm = T),
#     N = sum(N),
#     species_cor = "species"
#     ), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset)]
#   )
#
# cors_dt_onlyspecies <- cors_dt1_species[,.(
#   rmse = mean(rmse, na.rm = T),
#   spearmans_r = mean(spearmans_r, na.rm = T),
#   obs_center = mean(obs_center, na.rm = T),
#   pred_center = mean(pred_center, na.rm = T),
#   N = sum(N)
# ), by = .(variable, spatial_holdout, temporal_holdout, test_train, scale,response, simreal, dataset, species)]
#
# # count number of "." in response
# cors_dt[, N_dots := stringr::str_count(response, "\\.")]
#
# cors_dt[, ylabels := paste0(response, " (",simreal,")"),]
# cors_dt[, N_dots := stringr::str_count(ylabels, "\\.")]
# cors_dt[simreal == "real", N_dots := 50,]
# cors_dt[, ylabels := forcats::fct_reorder(ylabels, N_dots),]
#
# cors_dt[,variable := factor(variable, levels = c("ba","trees", "dbh", "growth", "mort", "reg"), labels = c("ba","trees", "dbh", "growth", "mort", "reg")),]
#
# for(i_species in c("overall", "species")){
#   # for(i_dataset in c("BCI", "Uholka")){
#   for(i_dataset in c("BCI")){
#     # p_dat_s = cors_dt[spatial_holdout == T & temporal_holdout == F & dataset == i_dataset]
#     p_dat_s = cors_dt[spatial_holdout == T & temporal_holdout == F & dataset == i_dataset & species_cor == i_species]
#     p_dat_s$holdout = "spatial"
#     # p_dat_s[scale == "species-period3-1patch" & variable == "dbh" & simreal == "real" & ylabels == "ba.trees.dbh.growth.mort.reg (real)" & holdout == "spatial" & test_train == "test"]
#     # p_dat_t = cors_dt[temporal_holdout == T & spatial_holdout == F & dataset == i_dataset]
#     p_dat_t = cors_dt[temporal_holdout == T & spatial_holdout == F & dataset == i_dataset & species_cor == i_species]
#     p_dat_t$holdout = "temporal"
#     p_dat = rbind(p_dat_s, p_dat_t)
#     p_all = ggplot(
#       p_dat,
#       aes(y = forcats::fct_reorder(ylabels, N_dots), x = (as.factor(variable)))
#       )+
#       geom_tile(aes(fill = spearmans_r))+
#       facet_grid(gsub("-","\n",scale)~paste0(holdout,"\n",test_train), scales = "fixed")+
#       geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
#       # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#       scale_fill_gradientn(colors = c("blue", "white", "red"), limits = c(0,1))+
#       ggtitle(paste(unique(p_dat_s$dataset)))+
#       theme_classic()+
#       theme(legend.position = "none")+
#       # rotate facet labels
#       theme(
#         strip.text.y = element_text(angle = 0),
#         axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
#         )+
#       # rotate x labels
#       xlab("predicted variable")+
#       ylab("response variables")
#     p_all
#     # ggsave(paste0("figures/fits_",i_dataset,"_",i_species,".png"), p_all, width = 12, height = 10)
#     ggsave(paste0("figures/testfits_",i_dataset,"_",i_species,".png"), p_all, width = 12, height = 10)
#   }
# }
#
# #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# ## compare parameters ####
# #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# par_list <- list()
# i_result_folder = result_folders[1]
# for(i_result_folder in result_folders){
#   cat("\nstarting ", i_result_folder,"\n")
#   if(i_result_folder == "02_realdata"){
#     split_location = "CVsplits-realdata/"
#     simreal = "real"
#   }else if(i_result_folder == "02_simulated"){
#     split_location = "CVsplits-simdata/"
#     simreal = "sim"
#   }
#   files_list = list.files(path = paste0("results/", i_result_folder,"/"), pattern = "*.pt", full.names = TRUE)
#   # Read all files
#   m_list <- lapply(files_list, torch::torch_load)
#   names(m_list) <- gsub(".pt","", basename(files_list))
#   i = 1
#   for(i in 1:length(m_list)){
#     m = m_list[[i]]
#     name = names(m_list)[i]
#     cat("\r", i, "of", length(m_list), "start model:", name, "                                          ")
#     folder = strsplit(name,"_")[[1]][1]
#     # get full true model
#     true_model = torch_load(paste0("results/01_full/", folder, "_full.pt"))
#     par_list[[paste0(name,"_",simreal)]] <- list(
#       true = true_model$parameters_r,
#       fit = m$parameters_r
#     )
#   }
# }
# j = 1
# pars_dt <- data.table()
# cat("\n")
# for(j in 1:length(par_list)){
#   cat("\r", j, "of", length(par_list))
#   j_model = names(par_list)[j]
#   pars_fit = par_list[[j_model]]$fit
#   pars_true = par_list[[j_model]]$true
#   i = 1
#   for(i in 1:length(pars_true)){
#     i_par = names(pars_true)[i]
#     dt_true = data.table(pars_true[[i_par]])
#     dt_fit = data.table(pars_fit[[i_par]])
#     colnames(dt_true) = paste0(i_par, 1:ncol(dt_true))
#     colnames(dt_fit) = paste0(i_par, 1:ncol(dt_fit))
#     if(nrow(dt_true) > 1) dt_true$species = 1:nrow(dt_true) else dt_true$species = 0
#     if(nrow(dt_fit) > 1) dt_fit$species = 1:nrow(dt_fit) else dt_fit$species = 0
#     melt_dt_true = melt(dt_true, id.vars = "species")
#     melt_dt_fit = melt(dt_fit, id.vars = "species")
#     pars_dt = rbind(
#       pars_dt,
#       rbind(
#         data.table(melt_dt_true, par = i_par, model = j_model, truefit = "true"),
#         data.table(melt_dt_fit, par = i_par, model = j_model, truefit = "fit")
#         )
#       )
#   }
# }
#
# pars_dt2 <- dcast(pars_dt, species+par+model+variable~truefit, value.var = "value")
#
# pars_dt2[, scale := tstrsplit(model,"_",fixed = T)[[1]],]
# pars_dt2[, cv := paste0(tstrsplit(model,"_",fixed = T)[[2]],tstrsplit(model,"_",fixed = T)[[3]]),]
# pars_dt2[, response := tstrsplit(model,"_",fixed = T)[[4]],]
# pars_dt2[, simreal := tstrsplit(model,"_",fixed = T)[[5]],]
# pars_dt2[grepl("species", model), dataset := "Uholka",]
# pars_dt2[!grepl("species", model), dataset := "BCI",]
#
#
# pars_cor_dt <-
# pars_dt2[,.(
#   var = var(true - fit, na.rm = T),
#   bias = mean(true - fit, na.rm = T)
# ), by = .(par, scale, cv, response, dataset, simreal)]
#
# pars_cor_dt[grepl("S0", cv), S_test_train := "train",]
# pars_cor_dt[!grepl("S0", cv), S_test_train := "test",]
# pars_cor_dt[grepl("T0", cv), T_test_train := "train",]
# pars_cor_dt[!grepl("T0", cv), T_test_train := "test",]
#
# pars_cor_dt2 <- pars_cor_dt[,.(
#   var = mean(var,na.rm = T),
#   bias = mean(bias, na.rm = T)
# ), by = .(par, scale,S_test_train, T_test_train, response, dataset, simreal)]
#
# pars_cor_dt2[, N_dots := stringr::str_count(response, "\\.")]
#
# for(i_dataset in c("BCI", "Uholka")){
#   p_dat_s = pars_cor_dt2[T_test_train == "train" & dataset == i_dataset & simreal == "sim"]
#   p_dat_t = pars_cor_dt2[S_test_train == "train" & dataset == i_dataset & simreal == "sim"]
#   p_s_bias <- ggplot(p_dat_s, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
#     geom_tile(aes(fill = bias))+
#     # geom_point(aes(size = bias, color = bias))+
#     facet_grid(scale~S_test_train+simreal, scales = "fixed")+
#     geom_text(aes(label = round(bias,2)), size = 3, color = "black")+
#     # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#     ggtitle(paste0(i_dataset, " spatial holdout"))+
#     theme_classic()+
#     # theme(legend.position = "none")+
#     # rotate x axis labels
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#     # change to multicolor scale that is easier to read
#     scale_fill_gradient2()+
#     xlab("parameter")+
#     ylab("response variables")
#
#   p_t_bias <- ggplot(p_dat_t, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
#     geom_tile(aes(fill = bias))+
#     # geom_point(aes(size = bias, color = bias))+
#     facet_grid(scale~T_test_train+simreal, scales = "fixed")+
#     geom_text(aes(label = round(bias,2)), size = 3, color = "black")+
#     # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#     ggtitle(paste0(i_dataset, " temporal holdout"))+
#     theme_classic()+
#     # theme(legend.position = "none")+
#     # rotate x axis labels
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#     # change to multicolor scale that is easier to read
#     scale_fill_gradient2()+
#     xlab("parameter")+
#     ylab("response variables")
#
#   p_all_bias = gridExtra::grid.arrange(p_s_bias, p_t_bias, ncol = 2)
#
#   p_s_var <- ggplot(p_dat_s, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
#     geom_tile(aes(fill = var))+
#     # geom_point(aes(size = var, color = var))+
#     facet_grid(scale~S_test_train+simreal, scales = "fixed")+
#     geom_text(aes(label = round(var,2)), size = 3, color = "black")+
#     # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#     ggtitle(paste0(i_dataset, " spatial holdout"))+
#     theme_classic()+
#     # theme(legend.position = "none")+
#     # rotate x axis labels
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#     # change to multicolor scale that is easier to read
#     scale_fill_viridis_c(option = "C")+
#     xlab("parameter")+
#     ylab("response variables")
#
#   p_t_var <- ggplot(p_dat_t, aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(par))))+
#     geom_tile(aes(fill = var))+
#     # geom_point(aes(size = var, color = var))+
#     facet_grid(scale~T_test_train+simreal, scales = "fixed")+
#     geom_text(aes(label = round(var,2)), size = 3, color = "black")+
#     # scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
#     ggtitle(paste0(i_dataset, " temporal holdout"))+
#     theme_classic()+
#     # theme(legend.position = "none")+
#     # rotate x axis labels
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#     # change to multicolor scale that is easier to read
#     scale_fill_viridis_c(option = "C")+
#     xlab("parameter")+
#     ylab("response variables")
#
#   p_all_var = gridExtra::grid.arrange(p_s_var, p_t_var, ncol = 2)
#
#   ggsave(paste0("figures/parameters_var_",i_dataset,".png"), p_all_var, width = 30, height = 10)
#   ggsave(paste0("figures/parameters_bias_",i_dataset,".png"), p_all_bias, width = 30, height = 10)
# }
