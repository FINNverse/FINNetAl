library(torch)
library(FINN)

models = list.files("results/03_full100reps/", full.names = T, recursive = T)
m_true = torch::torch_load("results/01_full/pft-period7-25patches_full.pt")
m_list = sapply(models, torch::torch_load)
m_pars = lapply(m_list, function(x) x$parameters_r)

i = names(m_true$parameters_r)[1]
m_true$nn_growth$parameters

diffs_df = data.frame()
par_names = c("par_competition_r", "par_growth_r", "par_mortality_r", "par_regeneration_r")
for(i in par_names){
  for(c in 1:ncol(m_true[[i]])){
    for(m in m_list[1:length(m_list)]){
      diffs_df = rbind(
        diffs_df,
        data.frame(diff = m_true[[i]][,c] - m[[i]][,c], par = paste0(i,c))
        )
    }
  }
}

par_names_nn = c("nn_mortality.0.weight", "nn_growth.0.weight", "nn_regeneration.0.weight")
for(i in par_names_nn){
  for(c in 1:ncol(m_true$parameters_r[[i]])){
    for(m in m_list[1:length(m_list)]){
      diffs_df = rbind(
        diffs_df,
        data.frame(diff = m_true$parameters_r[[i]][,c] - m$parameters_r[[i]][,c], par = paste0(i,c))
      )
    }
  }
}

par(mfrow = c(1,1), mar = c(3,8,1,1))
boxplot(diff ~ par, data = diffs_df, horizontal = T, las = 1, cex.axis = 0.7)
abline(v = 0)



