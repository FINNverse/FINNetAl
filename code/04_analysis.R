library(data.table)
library(ggplot2)
library(ggh4x)
library(FINN)
library(torch)

files_list = list.files(path = "results/02_realdata/", pattern = "*.pt", full.names = TRUE)
# files_list = files_list[grepl("T0", files_list)]

# Read all files
m_list <- lapply(files_list, torch::torch_load)
names(m_list) <- gsub(".pt","", basename(files_list))

# i=4
pred_list <- list()
i = 1
for(i in 1:length(m_list)){
  m = m_list[[i]]
  name = names(m_list)[i]
  folder = strsplit(name,"_")[[1]][1]
  files_dir = paste0("data/BCI/CVsplits-realdata/",folder,"/")
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

  pred_list[[name]] <- list(
    pred = list(
      train = pred_train,
      test = pred_test
    ),
    obs = list(
      train = obs_dt_train,
      test = obs_dt_test
      )
    )
  cat("\n", i, "of", length(m_list), "models simulated")
}

# pred_l = pred_list[[1]]
# i=1
all_dt <- data.table()
for(i in 1:length(pred_list)){
  pred_l = pred_list[[i]]
  name = names(pred_list)[i]
  dt <-
    rbindlist(
      list(
        pred_train = data.table(pred_l$pred$train$long$site, pred_obs = "pred", test_train = "train"),
        pred_test = data.table(pred_l$pred$test$long$site, pred_obs = "pred", test_train = "test"),
        obs_dt_train =
          data.table(
            melt(pred_l$obs$test[,c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "mort", "reg")], id.vars = c("siteID", "species", "year")),
            pred_obs = "obs", test_train = "train"),
        obs_dt_test =
          data.table(
            melt(pred_l$obs$test[,c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "mort", "reg")], id.vars = c("siteID", "species", "year")),
            pred_obs = "obs", test_train = "test")
      ),
      use.names=TRUE
    )
  dt$scale = strsplit(name,"_")[[1]][1]
  dt$cv = paste0(strsplit(name,"_")[[1]][2], strsplit(name,"_")[[1]][3])
  dt$response = strsplit(name,"_")[[1]][4]
  all_dt <- rbind(all_dt, dt)
}
all_dt[variable == "reg" & pred_obs == "pred", value := value/0.1,]
all_dt <- all_dt[variable != "r_mean_ha"]

# r2 = function(pred, obs, na.rm = T) {
#   SS_res <- sum((obs - pred)^2, na.rm = na.rm)  # Sum of squared residuals
#   SS_tot <- sum((obs - mean(obs))^2, na.rm = na.rm) # Total sum of squares
#   R_squared <- 1 - (SS_res / SS_tot)
#   R_squared[R_squared<=0] = 0
#   return(R_squared)
# }
calculate_r2 <- function(observed, predicted) {
  # Remove NA values
  valid_indices <- complete.cases(observed, predicted)
  observed <- observed[valid_indices]
  predicted <- predicted[valid_indices]
  # Compute R-squared
  ss_total <- sum((observed - mean(observed))^2)
  ss_residual <- sum((observed - predicted)^2)
  r2 <- 1 - (ss_residual / ss_total)
  return(r2)
}

table(pred_dt$scale)
# dcast by pred_obs
pred_dt <- dcast(all_dt, siteID+year+species+variable+test_train+scale+cv+response~pred_obs, value.var = "value")

pred_dt[test_train == "test" & response == "ba.growth" & cv == "S1T0"]
pred_dt[test_train == "train" & response == "ba.growth" & cv == "S1T0"]

table(pred_dt$cv)
pred_dt[grepl("S0", cv), S_test_train := "train",]
pred_dt[!grepl("S0", cv), S_test_train := "test",]
pred_dt[grepl("T0", cv), T_test_train := "train",]
pred_dt[!grepl("T0", cv), T_test_train := "test",]
# pred_dt[grecv == "S0T0", full_test_train := "train",]
# pred_dt[cv != "S0T0", full_test_train := "test",]

cors_dt1 =
  pred_dt[,.(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    spearmans_r = cor(pred, obs, method = "spearman"),
    # r2 = r2(pred, obs),
    r2 = calculate_r2(pred, obs),
    obs_center = sum(range(obs, na.rm = T))/2,
    pred_center = sum(range(pred, na.rm = T)/2),
    N = .N
  ), by = .(variable, cv, S_test_train, T_test_train ,scale,response)]

cors_dt <- cors_dt1[,.(
  rmse = mean(rmse),
  spearmans_r = mean(spearmans_r),
  r2 = mean(r2),
  # r2_2 = mean(r2_2),
  obs_center = mean(obs_center),
  pred_center = mean(pred_center),
  N = sum(N)
  ), by = .(variable, S_test_train, T_test_train ,scale,response)]
# count number of "." in response
cors_dt[, N_dots := stringr::str_count(response, "\\.")]


p_s = ggplot(
  cors_dt[T_test_train == "train"],
  aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(variable)))
  )+
  geom_tile(aes(fill = r2))+
  # geom_tile(aes(fill = spearmans_r))+
  facet_grid(scale~S_test_train, scales = "free")+
  # geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
  geom_text(aes(label = round(r2,2)), size = 3, color = "black")+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  ggtitle("spatial holdout")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("predicted variable")+
  ylab("response variables")
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))
  # scale_fill_gradient(low = "blue", high = "red", )
p_t = ggplot(
  cors_dt[S_test_train == "train"],
  aes(y = forcats::fct_reorder(response, N_dots), x = (as.factor(variable)))
  )+
  geom_tile(aes(fill = r2))+
  # geom_tile(aes(fill = spearmans_r))+
  facet_grid(scale~T_test_train, scales = "free")+
  # geom_text(aes(label = round(spearmans_r,2)), size = 3, color = "black")+
  geom_text(aes(label = round(r2,2)), size = 3, color = "black")+
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1))+
  ggtitle("temporal holdout")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("predicted variable")+
  ylab("response variables")
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1))
  # scale_fill_gradient(low = "blue", high = "red", )
p_all = gridExtra::grid.arrange(p_s, p_t, ncol = 2)
ggsave("figures/02_simulated_r2.png", p_all, width = 15, height = 6)
#
# p2 <- ggplot(pred_dt1[test_train == "train" & response == "ba.growth"], aes(x = pred, y = obs, color = factor(species))) +
#   geom_point(alpha = 0.1) +
#   stat_smooth(method = "lm", se = FALSE) +
#   geom_abline(intercept = 0, slope = 1, color = "black") +
#   facet_grid2(cv+response ~ variable, scales = "free", independent = "all")+
#   geom_text(
#     data = cors_dt[test_train == "train" & response == "ba.growth"],
#     aes(label = paste0(round(r2,2)), x = pred_center, y = obs_center),
#     size = 5, color = "black", fontface = "bold" )+
#   theme(
#     panel.grid = element_blank(),     # Remove gridlines
#     panel.background = element_blank(), # Remove panel background
#     axis.text = element_blank(),      # Remove axis text
#     axis.ticks = element_blank(),     # Remove axis ticks
#     strip.background = element_blank(), # Remove strip background
#     strip.text.y = element_text(angle = 0) # Make facet labels horizontal
#   )
# p2
#
