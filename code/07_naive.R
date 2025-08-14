library(data.table)
library(FINN)
library(torch)
library(cito)
# obs_dt <- fread("data/BCI/noSplits/pft-period7-25patches/obs_dt.csv")
# env_dt <- fread("data/BCI/noSplits/pft-period7-25patches/env_dt.csv")

preds_all <- data.table()
true_all <- data.table()
k=1
for(k in 0:5){
  obs_dt_train <- fread(paste0("data/BCI/CVsplits-realdata/pft-period7-25patches/obs_dt_S",k,"_T0_train.csv"))
  obs_dt_train <- obs_dt_train[,.(siteID, year, species, dbh, ba, trees, growth, mort, reg)]
  obs_dt_test <- fread(paste0("data/BCI/CVsplits-realdata/pft-period7-25patches/obs_dt_S",k,"_T0_test.csv"))
  obs_dt_test <- obs_dt_test[,.(siteID, year, species, dbh,ba, trees, growth, mort, reg)]
  
  env_dt_train <- fread(paste0("data/BCI/CVsplits-realdata/pft-period7-25patches/env_dt_S",k,"_T0_train.csv"))
  env_dt_train = env_dt_train[,.(siteID, year, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp )]
  env_dt_test <- fread(paste0("data/BCI/CVsplits-realdata/pft-period7-25patches/env_dt_S",k,"_T0_test.csv"))
  env_dt_test = env_dt_test[,.(siteID, year, Prec, SR_kW_m2, RH_prc, T_max, T_min, swp )]
  
  
  
  env_obs_train = merge(obs_dt_train, env_dt_train, by = c('year', "siteID"))
  train_df = data.frame()
  for(i in 1:(max(env_obs_train$year)-1 ))  {
    df_train = merge(env_obs_train[year==i+1], env_obs_train[year==i],  by = c('siteID', "species"), all = TRUE)
    df_train[,c("Prec.x", "SR_kW_m2.x",  "RH_prc.x",    "T_max.x", "T_min.x", "swp.x", "year.y" ):=NULL]
    train_df = rbind(train_df, df_train)
  }
  names(train_df)
  train_df = as.data.frame(train_df)
  centers_train = attr(scale(train_df[,-(1:9)]), "scaled:center")
  scales_train = attr(scale(train_df[,-(1:9)]), "scaled:scale")
  train_df[,-(1:9)] = scale(train_df[,-(1:9)])
  train_df$species = as.factor(train_df$species)
  train = train_df#[train_df$siteID %in% 1:15, ]
  
  env_obs_test = merge(obs_dt_test, env_dt_test, by = c('year', "siteID"))
  test_df = data.frame()
  for(i in 1:(max(env_obs_test$year)-1 ))  {
    df_test = merge(env_obs_test[year==i+1], env_obs_test[year==i],  by = c('siteID', "species"), all = TRUE)
    df_test[,c("Prec.x", "SR_kW_m2.x",  "RH_prc.x",    "T_max.x", "T_min.x", "swp.x", "year.y" ):=NULL]
    test_df = rbind(test_df, df_test)
  }
  names(test_df)
  test_df = as.data.frame(test_df)
  centers_test = attr(scale(test_df[,-(1:9)]), "scaled:center")
  scales_test = attr(scale(test_df[,-(1:9)]), "scaled:scale")
  test_df[,-(1:9)] = scale(test_df[,-(1:9)])
  test_df$species = as.factor(test_df$species)
  test = test_df#[test_df$siteID %in% 1:15, ]
  
  
  
  custom = function(pred, true) { 
    torch::nnf_mse_loss(pred[,1], true[,1])*0.1 + 
      torch::nnf_mse_loss(pred[,2], true[,2])*10.0+
      torch::nnf_mse_loss(pred[,3]$exp(), true[,3])+
      torch::nnf_mse_loss(pred[,4], true[,4])*10.0+
      torch::nnf_mse_loss(pred[,5], true[,5])*1.0+
      torch::nnf_mse_loss(pred[,6]$exp(), true[,6])
  }
  
  m = dnn(cbind(dbh.x, ba.x, trees.x, growth.x, mort.x, reg.x)~dbh.y + ba.y +trees.y + growth.y + mort.y + reg.y + Prec.y + SR_kW_m2.y+RH_prc.y+T_max.y+T_min.y+swp.y+e(species, dim = 2),
          data = train , loss = custom, lr = 0.001,optimizer = "adam", burnin = Inf, plot = F, epochs = 200)
  
  tmp = test
  predictions = data.frame()
  for(i in 2:6) {
    pred = predict(m, newdata = tmp[tmp$year.x == i,])
    pred[,3] = exp(pred[,3])
    pred[,6] = exp(pred[,6])
    pred[pred<0] = 0
    colnames(pred) = c("dbh.x", "ba.x", "trees.x", "growth.x", "mort.x", "reg.x")
    if(i < 6) tmp[tmp$year.x == (i+1),10:15] = (pred - matrix(centers_train[1:6], nrow = nrow(pred), ncol = 6L, byrow = TRUE) )/ matrix(scales_train[1:6], nrow = nrow(pred), ncol = 6L, byrow = TRUE) 
    pred = data.frame(pred)
    pred$siteID = tmp[tmp$year.x == i,]$siteID
    pred$species = tmp[tmp$year.x == i,]$species
    pred$year = i
    predictions = rbind(predictions, pred)
  }
  preds = as.data.table(predictions)[,.(dbh.x = mean(dbh.x), ba.x = mean(ba.x), trees.x = mean(trees.x), growth.x = mean(growth.x), mort.x = mean(mort.x), reg.x = mean(reg.x)), by = .(species, year)]
  true = as.data.table(test)[,.(dbh.x = mean(dbh.x), ba.x = mean(ba.x), trees.x = mean(trees.x), growth.x = mean(growth.x), mort.x = mean(mort.x), reg.x = mean(reg.x)), by = .(species, year.x)]

  preds_all <- rbind(preds_all, data.table(preds, k = k))
  true_all <- rbind(true_all, data.table(true, k = k))
  
}



k=1
for(K in 1:5){
  preds = preds_all[k == K][,-"k"]
  true = true_all[k==K][,-"k"]
  par(mfrow = c(2, 3))
  counter = 1
  mins = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  maxs = c(13, 5.0, 350, 0.35, 0.35, 350)
  for(var in colnames(preds)[-(1:2)]) {
    # var="dbh.x"
    # preds[species==1][[var]]
    matplot(sapply(1:5, function(i) preds[species==i][[var]]), type = "l", lty = 1, lwd = 2.0, main = var, ylim = c(mins[counter], maxs[counter]))
    matplot(sapply(1:5, function(i) true[species==i][[var]]), type = "l", lty = 2, lwd = 2.0, add = TRUE)
    counter = counter +1
  }
}

names(true_all)[names(true_all) == "year.x"] = "year"
comb_dt <- merge(preds_all, true_all, by = c("species", "year", "k"), suffixes = c(".pred", ".true"))
names(comb_dt) = gsub(".x.",".",names(comb_dt))

comb_dt = melt(comb_dt, id.vars = c("species", "year", "k"))
comb_dt[, var := tstrsplit(variable,"\\.")[[1]],]
comb_dt[, simpred := tstrsplit(variable,"\\.")[[2]],]

comb_dt_allspecies <- rbind(
  comb_dt[var %in% c("ba", "trees"), .(value = sum(value)), by = .(year,k,var,simpred)],
  comb_dt[var %in% c("dbh", "growth", "mort", "reg"), .(value = mean(value)), by = .(year,k,var,simpred)]
  )
comb_dt_allspecies$species = "all"

comb_dt <- rbind(comb_dt[,-"variable"], comb_dt_allspecies)

comb_dt <- dcast(comb_dt, species + year + k + var ~ simpred, value.var = "value")
cors_dt <- comb_dt[,.(r = cor(pred,true)), by = .(var, species, k)]

library(ggplot2)
ggplot(cors_dt[k != 0], aes(x = factor(species), y = r))+
  geom_hline(yintercept = 0.0)+
  geom_boxplot()+
  # geom_point()+
  facet_wrap(~var)


fwrite(cors_dt[k != 0], "results/cors_naive.csv")


# 
# init = test[1:5,]
# init[,10:15] = 0.0
# init
# predictions_eq = data.frame()
# for(i in 1:500) {
#   pred = predict(m, newdata = init)
#   pred[,3] = exp(pred[,3])
#   pred[,6] = exp(pred[,6])
#   pred = jitter(pred, factor = 2)
#   pred[pred<0] = 0
#   colnames(pred) = c("dbh.x", "ba.x", "trees.x", "growth.x", "mort.x", "reg.x")
#   init[,10:15] = (pred - matrix(centers[1:6], nrow = nrow(pred), ncol = 6L, byrow = TRUE) )/ matrix(scales[1:6], nrow = nrow(pred), ncol = 6L, byrow = TRUE)
#   
#   init[,16:21] = as.matrix(test[sample.int(nrow(test), 1),16:21], nrow = 5, ncol = 6, byrow = TRUE)
#   
#   pred = data.frame(pred)
#   pred$siteID = 1
#   pred$species = init$species
#   pred$year = i
#   predictions_eq = rbind(predictions_eq, pred)
# }
# predictions_eq = as.data.table(predictions_eq)[,.(dbh.x = mean(dbh.x), ba.x = mean(ba.x), trees.x = mean(trees.x), growth.x = mean(growth.x), mort.x = mean(mort.x), reg.x = mean(reg.x)), by = .(species, year)]
# 
# par(mfrow = c(2, 3))
# counter = 1
# mins = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# maxs = c(13, 5.0, 350, 0.35, 0.35, 350)
# for(var in colnames(predictions_eq)[-(1:2)]) {
#   matplot(sapply(1:5, function(i) predictions_eq[species==i][[var]]), type = "l", lty = 1, lwd = 2.0, main = var, ylim = c(mins[counter], maxs[counter]))
#   counter = counter +1
# }

