#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## correlations for all combinations ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(torch)
library(FINN)
library(data.table)
library(ggplot2)
source("code/plot-functions.R")

gc()
all_models <- c(
  list.files("results/02_realdata", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridTF0", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridTF1", full.names = T, recursive = T),
  list.files("results/02_realdata_hybrid_mortTF0", full.names = T, recursive = T),
  list.files("results/02_realdata_hybrid_mortTF1", full.names = T, recursive = T),
  # list.files("results/02_realdata_hybridTF1fixed", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmall", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmallDropout", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmallDropoutfixed", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMedium", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMediumDropout", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMediumDropoutfixed", full.names = T, recursive = T)
)

overwrite = F
if(overwrite){
  all_cors_dt <- data.table()
  all_cors_dt_plot <- data.table()
  all_cors_dt_species <- data.table()
  all_cors_dt_plot_species <- data.table()
}else{
  all_cors_dt <- fread("results/all_cors_dt.csv")
  all_cors_dt_plot <- fread("results/all_cors_dt_plot.csv")
  all_cors_dt_species <- fread("results/all_cors_dt_species.csv")
  all_cors_dt_plot_species <- fread("results/all_cors_dt_plot_species.csv")

  existing_cors <- paste0("results/",all_cors_dt$data, "/", all_cors_dt$name)
  all_models <- all_models[!all_models %in% existing_cors]

  # all_models <- c(paste0("results/02_realdata/",gsub(".pt","",unique(basename(all_models))),"_ba.trees.dbh.growth.mort.reg.pt"), all_models)
}
all_models <- all_models[grepl("period7-25patches", all_models)]

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
    cor_dt_in[!is.na(obs) & !is.na(pred),.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_species <- cor_dt_in[!is.na(obs) & !is.na(pred),.(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    spearmans_r = cor(pred, obs, method = "spearman"),
    obs_center = sum(range(obs, na.rm = T))/2,
    pred_center = sum(range(pred, na.rm = T)/2),
    N = .N
  ), by = .(variable,species, test_train)]
  cor_dt_plot <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred),.(
      rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
      spearmans_r = cor(pred, obs, method = "spearman"),
      obs_center = sum(range(obs, na.rm = T))/2,
      pred_center = sum(range(pred, na.rm = T)/2),
      N = .N
    ), by = .(variable, test_train)]
  cor_dt_plot_species <-
    cor_dt_in_plot[!is.na(obs) & !is.na(pred),.(
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

cors_list = list(
  "perSite" = all_cors_dt,
  "fullPlot" = all_cors_dt_plot,
  "perSiteSpecies" = all_cors_dt_species,
  "fullPlotSpecies" =all_cors_dt_plot_species
)
i = 1
for(i in 1:length(cors_list)){
  # set dataset
  cors_list[[i]][grepl("species-",name),dataset := "Uholka",]
  cors_list[[i]][grepl("genus-",name) | grepl("pft-",name),dataset := "BCI",]
  # set simreal
  cors_list[[i]][grepl("realdata",data),simreal := "real",]
  cors_list[[i]][grepl("simulated",data),simreal := "simulated",]
  # set hybrid
  
  cors_list[[i]][grepl("hybrid", data),hybrid := gsub("hybrid","",tstrsplit(gsub("hybrid_","hybrid",data),"_",fixed = T)[[3]]),]
  cors_list[[i]][!grepl("hybrid", data),hybrid := "nohybrid",]

  # set cv
  cors_list[[i]][, cv := paste0(tstrsplit(name,"_")[[2]],tstrsplit(name,"_")[[3]]),]
  cors_list[[i]][, Scv := tstrsplit(name,"_")[[2]],]
  cors_list[[i]][, Tcv := tstrsplit(name,"_")[[3]],]
  cors_list[[i]][, cv := gsub(".pt","",cv),]
  cors_list[[i]][, Scv := gsub(".pt","",Scv),]
  cors_list[[i]][, Tcv := gsub(".pt","",Tcv),]
  # set response
  cors_list[[i]][, response := tstrsplit(name,"_")[[4]],]
  # set scale
  cors_list[[i]][, scale := tstrsplit(name,"_")[[1]],]
}

fwrite(cors_list[["perSite"]], "results/cors_perSite.csv")
fwrite(cors_list[["fullPlot"]], "results/cors_fullPlot.csv")
fwrite(cors_list[["perSiteSpecies"]], "results/cors_perSiteSpecies.csv")
fwrite(cors_list[["fullPlotSpecies"]], "results/cors_fullPlotSpecies.csv")

i=1
pdf("figures/03-hybrid_cors.pdf", width = 16, height = 10)
for(i in 1:length(cors_list)){
  dt = cors_list[[i]]
  i_name = names(cors_list)[i]
  if(!("species" %in% names(dt))) dt[, species := 1,]
  pdat <- dcast(dt, species+variable+test_train+dataset+cv+scale+simreal~hybrid, value.var = "spearmans_r")
  # difference for "Medium"        "MediumDropout" "Small"         "SmallDropout"  "TF0"           "TF1"           "mortTF0"       "mortTF1"       "no"      
  pdat[,":="(
    Medium.diff = Medium - nohybrid,
    MediumDropout.diff = MediumDropout - nohybrid,
    Small.diff = Small - nohybrid,
    SmallDropout.diff = SmallDropout - nohybrid,
    TF0.diff = TF0 - nohybrid,
    TF1.diff = TF1 - nohybrid,
    mortTF0.diff = mortTF0 - nohybrid,
    mortTF1.diff = mortTF1 - nohybrid
  ),]

  pdat2 <- melt(pdat[,-c("Medium", "MediumDropout", "Small", "SmallDropout", "TF0", "TF1", "mortTF0", "mortTF1")],
  # pdat2 <- melt(pdat[,-c("growth_NoTF", "growth_TF", "mort_NoTF", "mort_TF", "growth_TF1fixed")],
                id.vars = c("nohybrid","species","variable","test_train","dataset","cv","scale","simreal"), variable.name = "hybrid", value.name = "r_diff")

  pdat2 <- pdat2[,.(
    r_diff = mean(r_diff, na.rm = T)
  ), by = .(nohybrid,variable,test_train,dataset,cv,scale,hybrid)]


  # ggplot(pdat2, aes(x = nohybrid, y = r_diff,color = variable))+
  #   geom_point()+
  #   geom_smooth(method = "lm")+
  #   # facet_grid(scale~hybrid+test_train)+
  #   facet_wrap(~variable)+
  #   theme_classic()+
  #   theme(
  #     strip.text.y = element_text(angle = 0),
  #     strip.background = element_blank()
  #   )

  pdat2 <- pdat2[,hybrid := factor(
    hybrid, levels = sort(as.character(unique(pdat2$hybrid))),
    ordered = T),]
  p=ggplot(pdat2, aes(x = cv, y = variable))+
    geom_tile(aes(fill = nohybrid))+
    facet_grid(scale~test_train+hybrid)+
    theme_classic()+
    theme(
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank()
    )+
    #fill color scale that is centered at 0
    scale_fill_gradient2(name = "spearmans_r",low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))+
    ggtitle(paste0(i_name, " ", "nohybrid cor"))
  print(p)
  p=ggplot(dt[hybrid != "nohybrid"], aes(x = cv, y = variable))+
    geom_tile(aes(fill = spearmans_r))+
    facet_grid(scale~test_train+paste0(hybrid_process,"_TF",as.integer(transformer=="yes")))+
    theme_classic()+
    theme(
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank()
    )+
    #fill color scale that is centered at 0
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))+
    ggtitle(paste0(i_name, " ", "hybrid corr"))
  print(p)
  p=ggplot(pdat2, aes(x = cv, y = variable))+
    geom_tile(aes(fill = r_diff))+
    facet_grid(scale~test_train+hybrid)+
    theme_classic()+
    theme(
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank()
    )+
    #fill color scale that is centered at 0
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))+
    ggtitle(paste0(i_name, " ", "diff"))
  print(p)
}
dev.off()

all_models <- c(
  list.files("results/02_realdata", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridTF0", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridTF1", full.names = T, recursive = T),
  list.files("results/02_realdata_hybrid_mortTF0", full.names = T, recursive = T),
  list.files("results/02_realdata_hybrid_mortTF1", full.names = T, recursive = T),
  # list.files("results/02_realdata_hybridTF1fixed", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmall", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmallDropout", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridSmallDropoutfixed", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMedium", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMediumDropout", full.names = T, recursive = T),
  list.files("results/02_realdata_hybridMediumDropoutfixed", full.names = T, recursive = T)
)

exmple_models <- all_models[grepl("pft-period7-25patches_S0_T0", all_models)]
# exmple_models2 = gsub("results/02_realdata_hybridTF1","results/02_realdata_hybridTF1fixed",exmple_models)
pred_nohybrid = build_model_dt(exmple_models[1])[[1]][test_train == "train"]
pred_growthTF0 = build_model_dt(exmple_models[2])[[1]][test_train == "train"]
pred_growthTF1 = build_model_dt(exmple_models[3])[[1]][test_train == "train"]
pred_growthTF1fixed = build_model_dt(gsub("results/02_realdata_hybridTF1","results/02_realdata_hybridTF1fixed",exmple_models[3]))[[1]][test_train == "train"]

build_model_dt(exmple_models[1])
models = lapply(exmple_models, function(x) {
  build_model_dt(x)[[1]]
})
names(models) = exmple_models
for(i in 1:length(models)) models[[i]]$hybrid = basename(dirname(names(models)[i]))

pdat = rbindlist(models)
pdat[dbh.obs == 0]
pdat <- pdat[dbh.pred >= 1 & dbh.obs >= 1]
pdat[hybrid == "02_realdata_hybridTF1fixed"]$growth.pred
# pdat <- rbindlist(models)
p1 = ggplot(pdat, aes(x = growth.obs, y = growth.pred))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~hybrid, ncol = 1)+
  ggtitle("growth")+
  coord_cartesian(ylim = c(0,1))
p1
p2 = ggplot(pdat, aes(x = dbh.pred, y = mort.pred))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~hybrid, ncol = 1)+
  ggtitle("mort")+
  coord_cartesian(ylim = c(0,1), xlim = c(0,50))
p2
p1 = ggplot(rbindlist(list(
  data.table(pred_nohybrid, hybrid = "nohybrid"),
  data.table(pred_growthTF0, hybrid = "growthTF0"),
  data.table(pred_growthTF1, hybrid = "growthTF1"),
  data.table(pred_growthTF1fixed, hybrid = "growthTF1fixed")
), fill = T), aes(x = dbh.pred, y = growth.pred))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~hybrid, ncol = 1)+
  ggtitle("growth")

p1
p2 = ggplot(pred_growthTF0, aes(x = dbh.pred, y = growth.pred))+
  geom_point()+
  geom_smooth()+
  ggtitle("growth no transformer")

p3 = ggplot(pred_growthTF1, aes(x = dbh.pred, y = growth.pred))+
  geom_point()+
  geom_smooth()+
  ggtitle("growth transformer")
gridExtra::grid.arrange(p1,p2,p3, ncol = 1)


plot(growth.pred~dbh.pred, pred_nohybrid)
# fit model
fm = lm(growth.pred~dbh.pred, pred_nohybrid)
summary(fm)
# plot fitted line
abline(fm, col = "red")

dt[simreal == "real" & test_train == "test",.N, by = .(Tcv, Scv,variable, hybrid, scale)][N > 1]
dt[simreal == "real" & test_train == "test" & Tcv == "T2" & Scv == "S5" & hybrid == "nohybrid" & variable == "trees" & scale == "pft-period7-25patches"]
ggplot(dt[simreal == "real" & test_train == "test"], aes(x = Tcv, y = Scv))+
  geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 1)+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  geom_text(aes(label = round(spearmans_r,2)), size = 2, color = "black")+
  facet_grid(scale~variable+hybrid)+
  #rotate y facets
  theme_classic()+
  theme(
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank()
  )

# ggplot(dt[simreal == "real" & test_train == "test"], aes(x = Tcv, y = Scv))+
ggplot(dt[simreal == "real" & test_train == "test" & hybrid == "nohybrid"], aes(x = paste0(variable), y = scale))+
  geom_tile(aes(fill = spearmans_r), color = "white", linewidth = 0.1)+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  geom_text(aes(label = round(spearmans_r,2)), size = 2, color = "black")+
  # facet_grid(scale~variable+hybrid)+
  facet_grid(Scv~dataset+Tcv)+
  #rotate y facets
  theme_classic()+
  theme(
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank()
  )

pdt1 <- dt[test_train == "test" & !grepl("0",cv) & hybrid == "nohybrid",.(spearmans_r = mean(spearmans_r)), by = .(simreal,variable, scale)]
ggplot(pdt1, aes(y = spearmans_r, x = variable, fill = simreal))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(~scale)+
  theme_classic()+
  theme(
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank()
  )
