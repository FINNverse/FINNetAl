library(torch)
if(!dir.exists("results/02_realdata_hybridTF1fixed")) dir.create("results/02_realdata_hybridTF1fixed")

transformer_models = list.files("results/02_realdata_hybridTF1")
sapply(transformer_models, function(p) {

  tmp_m = torch::torch_load(paste0(i_dir, p))
  tmp_m$param_history = NULL
  gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
    self$nn_growth$train()
    g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
    return(g)
  }
  tmp_m$growth_func = tmp_m$.__enclos_env__$private$set_environment(gh)
  torch::torch_save(tmp_m, paste0(i_dir,"fixed/", p))
})


# apply fix for multiple directories

dirs = c("results/02_realdata_hybridSmall","results/02_realdata_hybridSmallDropout","results/02_realdata_hybridMedium","results/02_realdata_hybridMediumDropout")
i_dir = dirs[2]
for(i_dir in dirs){
  transformer_models = list.files(paste0(i_dir))
  # p = transformer_models[1]
  sapply(transformer_models, function(p) {
    tmp_m = torch::torch_load(paste0(i_dir,"/", p))
    tmp_m$param_history = NULL
    gh = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {
      self$nn_growth$train()
      g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred) - exp(1))$exp()
      return(g)
    }
    tmp_m$growth_func = tmp_m$.__enclos_env__$private$set_environment(gh)
    out_dir = paste0(i_dir,"fixed/", p)
    if(!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir))
    torch::torch_save(tmp_m, out_dir)
  })
}





