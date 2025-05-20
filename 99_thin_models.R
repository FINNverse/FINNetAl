library(FINN)
library(torch)

all_models <- list.files("results", pattern = ".pt", recursive = T, full.names = T)
all_models <- all_models[!grepl("01_learning_rate", all_models) & !grepl("S0_T0", all_models) & !grepl("full", all_models)]

i = all_models[1]
for(i in all_models){
  m = torch::torch_load(i)
  if(!is.null(m$history)){
    history_idx = unique(ceiling(c(seq(1,length(m$history), length.out = ceiling(length(m$history)*0.1)), length(m$history))))
    m$history = m$history[history_idx]
  }
  if(!is.null(m$param_history)){
    param_idx = unique(ceiling(c(seq(1,length(m$param_history), length.out = ceiling(length(m$param_history)*0.1)), length(m$param_history))))
    m$param_history = m$param_history[param_idx]
  }
  torch::torch_save(m, i)
}


