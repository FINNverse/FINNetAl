library(data.table)
library(ggplot2)

out_dir0 = "data/BCIsecf"

if(!dir.exists(out_dir0)) {
  dir.create(out_dir0,recursive = T)
}else{
  unlink(out_dir0, recursive = T)
  dir.create(out_dir0,recursive = T)
}

source("code/data_preparation/bci-data-prep-functions.R")
dbh_cmTOba_m2 <- function(dbh) {
  dbh = dbh/100
  return(pi*dbh^2/4)
}

secf_dt <- data.table(readxl::read_xls("data/raw-data/BCNM_secondary_forest_1ha_plots.xls", sheet = 1))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# assign species ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
species_codes <- fread("data/BCI/data-cleaning/pft/species_assigned.csv")

secf_dt$sp = tolower(secf_dt$sp_code)

secf_dt <- merge(secf_dt, species_codes, by = "sp", all.x = TRUE)

secf_dt[, species := speciesID,]

secf_dt[site_age == 40, siteID := as.integer(as.factor(site)),]

secf_dt[,year := (site_age-40)/5+1,]

all_trees <- secf_dt[,.(
  siteID = siteID,
  site_name = site,
  site_age,
  year,
  species,
  dbh_cm = dbh/10,
  ba = dbh_cmTOba_m2(dbh/10)
  ),]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# calculate stand variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

stand_dt <- all_trees[,.(
  dbh_mean = mean(dbh_cm),
  trees = .N,
  ba = sum(dbh_cmTOba_m2(dbh_cm))
  ), by = .(siteID, site_name, site_age, year, species)]

stand_dt_all <- stand_dt[,.(
  dbh_mean = mean(dbh_mean),
  trees = sum(trees),
  ba = sum(ba),
  species = "all"
  ), by = .(siteID,site_name, site_age, year)]

pdf(paste0(out_dir0, "/pdat_secundary_forest.pdf"), width = 10, height = 5)
p_list <- list(
  "individual species" = stand_dt,
  "all species" = stand_dt_all
)
for(i in 1:length(p_list)){
  i_name = names(p_list)[i]
  i_data = p_list[[i]]
  alpha = fifelse(i_name == "all species", 1, 0.3)
  p_ba <- ggplot(i_data, aes(x = site_age, y = ba, color = factor(species)))+
    geom_jitter(width = 0.5, alpha = alpha)+
    theme_minimal()+
    # theme(legend.position = "none")+
    ggtitle(paste("Basal area [m2] of", i_name))
  print(p_ba)
  p_trees <- ggplot(i_data, aes(x = site_age, y = trees, color = factor(species)))+
    geom_jitter(width = 0.5, alpha = alpha)+
    theme_minimal()+
    # theme(legend.position = "none")+
    ggtitle(paste("Number of trees of", i_name))
  print(p_trees)
  p_dbh <- ggplot(i_data, aes(x = site_age, y = dbh_mean, color = factor(species)))+
    geom_jitter(width = 0.5, alpha = alpha)+
    theme_minimal()+
    # theme(legend.position = "none")+
    ggtitle(paste("Mean dbh [cm] of", i_name))
  print(p_dbh)
}
dev.off()

# backtransform site_age in data.table
# stand_dt[,site_age2 := 40+(year-1)*5,]
fwrite(stand_dt, paste0(out_dir0,"/stand_dt_secf.csv"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# initial cohorts ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## patch size 1 ha ####
# i_age = 40
# for(i_age in c(40)){
#   # split cohorts by dbh size class for efficiency
#   dbh_cm_class_size = 0.1
#   initial_trees <- all_trees[
#     site_age == i_age, .(
#     trees = .N
#   ), by = .(dbh_cm = cut(
#     dbh_cm,
#     breaks = seq(1, 351, dbh_cm_class_size),
#     labels = seq(1, 351-dbh_cm_class_size, dbh_cm_class_size),
#     include.lowest = T), species, siteID = siteID, site_name)
#   ]
#   initial_trees <- initial_trees[order(siteID)]
#   initial_trees[,cohortID := 1:.N,by = .(siteID)]
#   initial_trees[, dbh_cm := as.numeric(as.character(dbh_cm)),]
# 
#   # check if basic stand properties in initial_trees are similar to all_trees and stand_dt
#   comp_stand_dt_from_initial_trees <- initial_trees[, .(
#     ba = sum(dbh_cmTOba_m2(dbh_cm)*trees, na.rm = T),
#     trees = sum(trees, na.rm = T),
#     dbh_mean = sum(dbh_cm*trees, na.rm = T)/sum(trees, na.rm = T)
#   ), by = .(siteID, species)]
#   compdt1 <- merge(
#     comp_stand_dt_from_initial_trees, stand_dt[site_age == i_age],
#     by = c("siteID", "species"), suffixes = c(".init", ".stand_dt")
#   )
# 
#   compdt1[, ba.diff := abs(ba.init - ba.stand_dt),]
#   compdt1[, dbh.diff := abs(dbh_mean.init - dbh_mean.stand_dt),]
#   compdt1[, trees.diff := abs(trees.init - trees.stand_dt),]
# 
#   if(any(compdt1$ba.diff > 0.1)) {
#     plot(compdt1$ba.init, compdt1$ba.stand_dt)
#     stop("ba difference is larger than 0.01")
#   }
# 
#   if(any(compdt1$dbh.diff > dbh_cm_class_size+0.001)) {
#     plot(compdt1$dbh_mean.init, compdt1$dbh_mean.stand_dt)
#     stop(paste("dbh difference is larger than", dbh_cm_class_size,"+0.001"))
#   }
# 
#   if(any(compdt1$trees.diff > 1)) {
#     plot(compdt1$trees.init, compdt1$trees.stand_dt)
#     stop("trees difference is larger than 1")
#   }
# 
#   ### initial_cohorts.csv ####
#   fwrite(initial_trees, paste0(out_dir,"/initial_cohorts_1patches_1ha_age",i_age,".csv"))
# }

## patch size 0.1 ha -> 10 patches ####

### define 0.1 ha patches ####

# split data in 10 patches with similar ba
patch_size = 0.1
Npatches = 5
Nrep = 20
all_patches_dt <- data.table()
set.seed(123)
for(i_site in unique(all_trees[site_age == 40]$siteID)){
  for(i_rep in 1:Nrep){
    dbh_classes <- cut(all_trees[siteID == i_site,]$dbh_cm, breaks = seq(0, 350, 10))
    site_dt <- copy(all_trees[siteID == i_site,])
    site_dt$repID = i_rep
    # site_dt$siteID = i_site
    site_dt[,patchID := sample(1:Npatches, .N, replace = T),]
    all_patches_dt <- rbind(all_patches_dt, site_dt)
  }
}
# all_trees[site_age == 40,.(ba = sum(dbh_cmTOba_m2(dbh_cm))), by = .(siteID, site_age,species)]
# new=all_patches_dt[,.(ba = sum(dbh_cmTOba_m2(dbh_cm))), by = .(siteID, repID, patchID, species)][,.(ba = max(ba)), by= .(species, siteID)]
# old=all_patches_dt[,.(ba = sum(dbh_cmTOba_m2(dbh_cm))), by = .(siteID, repID, patchID, species)][,.(ba = max(ba)), by= .(species, siteID)]
# 
# cohorts_dt[species == 1,.(ba = sum(dbh_cmTOba_m2(dbh_cm)*trees)), by = .(siteID, repID, species)]
# stand_dt50ha <- fread(paste0(dir0, "/data-cleaning/",i_species,"/stand_dt.csv"))
# # save histogram as pdf
# pdf(paste0(out_dir,"/ba_histogram_0.1ha_patches.pdf"), width = 5, height = 6)
# par(mfrow = c(2,1), mar = c(4,3,1,1), mgp =  c(2,0.7,0))
# all_1ha = all_patches_dt[,.(ba = sum(ba)), by = .(siteID,repID,patchID)]
# all_50ha = stand_dt50ha[,.(ba = sum(ba)), by = .(siteID,census)]
# ranges = range(c(all_1ha$ba, all_50ha$ba))
# hist(all_1ha$ba, xlim = c(0,ranges[2]), axes = F, breaks = seq(0, max(all_50ha$ba)+0.1, 0.1), main = NULL, xlab = "ba [m2] of 0.1 ha patches in 1 ha secundary forest")
# title("a)", adj = 0)
# axis(
#   side   = 1,
#   at     = seq(0, max(all_50ha$ba), by = 1),   # where to place the ticks
#   labels = seq(0, max(all_50ha$ba), by = 1)    # matching labels
# )
# axis(2)
# hist(all_50ha$ba, xlim = c(0,ranges[2]), axes = F, breaks = seq(0, max(all_50ha$ba)+.1, .1), main = NULL, xlab = "ba [m2] of 0.1 ha patches in 50 ha BCI forest plot")
# axis(
#   side   = 1,
#   at     = seq(0, max(all_50ha$ba), by = 1),   # where to place the ticks
#   labels = seq(0, max(all_50ha$ba), by = 1)    # matching labels
# )
# axis(2)
# title("b)", adj = 0)
# dev.off()

i_age = 40
for(i_age in c(40)){
  # split cohorts by dbh size class for efficiency
  dbh_cm_class_size = 0.1
  initial_trees <- all_patches_dt[
    site_age == i_age, .(
      trees = .N
    ), by = .(dbh_cm = cut(
      dbh_cm,
      breaks = seq(1, 351, dbh_cm_class_size),
      labels = seq(1, 351-dbh_cm_class_size, dbh_cm_class_size),
      include.lowest = T), species, siteID = siteID, site_name, patchID, repID)
  ]
  initial_trees <- initial_trees[order(siteID)]
  initial_trees[,cohortID := 1:.N,by = .(siteID)]
  initial_trees[, dbh_cm := as.numeric(as.character(dbh_cm)),]

  # check if basic stand properties in initial_trees are similar to all_trees and stand_dt
  comp_stand_dt_from_initial_trees <- initial_trees[, .(
    ba = sum(dbh_cmTOba_m2(dbh_cm)*trees, na.rm = T),
    trees = sum(trees, na.rm = T),
    dbh_mean = sum(dbh_cm*trees, na.rm = T)/sum(trees, na.rm = T)
  ), by = .(siteID, species, repID)]
  compdt1 <- merge(
    comp_stand_dt_from_initial_trees, stand_dt[site_age == i_age],
    by = c("siteID", "species"), suffixes = c(".init", ".stand_dt")
  )

  compdt1[, ba.diff := abs(ba.init - ba.stand_dt),]
  compdt1[, dbh.diff := abs(dbh_mean.init - dbh_mean.stand_dt),]
  compdt1[, trees.diff := abs(trees.init - trees.stand_dt),]

  if(any(compdt1$ba.diff > 0.1)) {
    plot(compdt1$ba.init, compdt1$ba.stand_dt)
    stop("ba difference is larger than 0.01")
  }

  if(any(compdt1$dbh.diff > dbh_cm_class_size+0.001)) {
    plot(compdt1$dbh_mean.init, compdt1$dbh_mean.stand_dt)
    stop(paste("dbh difference is larger than", dbh_cm_class_size,"+0.001"))
  }

  if(any(compdt1$trees.diff > 1)) {
    plot(compdt1$trees.init, compdt1$trees.stand_dt)
    stop("trees difference is larger than 1")
  }
  
  # hist(initial_trees[, .(ba = sum(dbh_cmTOba_m2(dbh_cm))), by = .(siteID, patchID, repID)]$ba)
  
  fwrite(initial_trees[repID == 1], paste0(out_dir0,"/initial_cohorts_10patches_0.1ha_age",i_age,".csv"))
  
  initial_trees[,cohortID := 1:.N,by = .(siteID, repID)]
  initial_trees[,patchID := as.integer(factor(paste0(repID,"_",patchID))), by=siteID]
  ### initial_cohorts.csv ####
  fwrite(initial_trees, paste0(out_dir0,"/initial_cohorts_100patches_0.1ha_age",i_age,".csv"))
}

# # add environmental data
# env_dt <- fread(paste0(dir0, i_species,"-", i_period,"-1patch/env_dt.csv"))
# env_dt <- env_dt[siteID == 1,][,-"siteID"]
# env_dt[,swp := 0,]
# 
# # get all simulation years
# sim_years <- seq(min(stand_dt$year), max(stand_dt$year), by = 1)
# 
# set.seed(123)
# sampled_idx <- sample(1:nrow(env_dt), length(sim_years), replace = T)
# 
# env_dt_secf <- env_dt[sampled_idx,]
# env_dt_secf[,year := sim_years,]
# 
# 
# all_env_dt <- data.table()
# for(i in unique(stand_dt[!is.na(siteID)]$siteID)){
#   all_env_dt <- rbind(
#     all_env_dt,
#     data.table(
#       siteID = i,
#       env_dt_secf
#       )
#     )
# }
# fwrite(all_env_dt, paste0(out_dir,"/env_dt_secf.csv"))
# }