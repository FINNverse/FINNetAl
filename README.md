

# Inferring processes within dynamic forest models using hybrid modeling

This repository contains the code to reproduce the results in
“**Inferring processes within dynamic forest models using hybrid
modeling**” by Maximilian Pichler & Yannek Käber

## Preprint

Preprint can be found on the [arXiv
server](https://arxiv.org/abs/2508.01228)

## Data

Results are based on simulated data and on empirical data. The original
BCI plot data is available from <https://doi.org/10.15146/5xcp-0d46> ,
meteorological data from <https://doi.org/10.15486/ngt/1771850>, water
potential map from <https://doi.org/10.6084/m9.figshare.c.4372898.v1>,
and species list from <https://doi.org/10.1126/science.aaz4797>. The
FINN model documentation is provided in the SI Appendix. FINN is
available as R package and can be installed from
<https://github.com/FINNverse/FINN>.

## Scripts

Scripts to reproduce the analysis can be found in the `code` folder. To
rerun everything, run in the following order:

Prepare data (from raw BCI data):

``` r
source("code/data_preparation/01-bci-data-preparation.R")
source("code/data_preparation/01-uholka-data-preparation.R")
source("code/data_preparation/01b-bci-secundary-forest-preparation.R")
source("code/data_preparation/02-spatial-holdout-folds-bci.R")
```

Run simulations and models:

- Run scripts from 01-07 in the code folder (GPU support is
  recommended!)

## Results

### Figure 2 - Recovery of functional form

<div id="fig-figure_2">

![](figures/fig-figure_2-1.png)

**Figure** 1: Hybrid-FINN recovers the functional form of a process and
improves predictive performance. To test whether Hybrid-FINN can recover
the functional form of a process, we simulated data from a Process-FINN
model (M1), in which the functional form of the growth process includes
an environmental unimodal niche curve (quadratic effect of the
environment). Then, we tested two models: M3, a Hybrid-FINN in which
growth was replaced by a deep neural network; and M2, a misspecified
Process-FINN model with an assumed linear relationship between
environment and growth. (a) compares the calibrated growth-environment
relationships of models M2 and M3 to the true model M1. (b) compares the
predictive performance of models M2 and M3 with the observed data,
measured by the Spearman correlation factor. Correlations are based on
equilibrium simulations and represent the relation for each variable and
site over 500 simulated timesteps.

</div>

### Figure 3 - Spearman correlation between predicted and observed

<div id="fig-figure_3">

![](figures/fig-figure_3-1.png)

**Figure** 2: Spearman correlation of observed and simulated variables
for five PFTs from the five-fold blocked spatial cross-validation of all
models calibrated/trained on the BCI data. ‘Process’ refers to
Process-FINN as a full process model with fixed functional forms,
‘Hybrid’ refers to Hybrid-FINN in which the growth process is replaced
by a DNN, ‘Naïve NN’ refers to a single DNN that was trained on the
response variables for each time step. The naïve NN received the derived
response variables of each previous time step as well as the current
environment as input. PFTs are according to Rüger et al. (2020)

</div>

### Figure 4 - Succession trajectories

<div id="fig-figure_4">

![](figures/fig-figure_4-1.png)

**Figure** 3: Simulated successional trajectories of stand
characteristics and demographic rates of five PFTs for Process-FINN (a)
and Hybrid-FINN (b), where the growth process is replaced by a DNN. Note
that rates are not annual and represent changes of five year intervals.
Both models were trained with the BCI forest data with 7 censuses.
Forest dynamics were simulated for 600 years. Disturbance regime of
Rüger et al. (2020) was used.

</div>

<!-- #### boxplots of cv corr -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_cv-cor-only -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 7 -->

<!-- #| fig-height: 5 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Calculated pearson correlation from the five-fold blocked spatial cross-validation. The first two rows show the observed and predicted values for basal area and tree density. The third and fourth row show the observed and predicted values for diameter at breast height (dbh), regeneration, growth, and mortality." -->

<!-- var_labels3 = c("Basal Area [m²/ha]", -->

<!--                 "Trees [N/ha]", -->

<!--                 "DBH [cm]", -->

<!--                 "Regeneration [N/ha]", -->

<!--                 "Growth [%/100]", -->

<!--                 "Mortality [%/100]") -->

<!-- all_cors[,variable2 := factor( -->

<!--     variable,  -->

<!--     levels = vars, -->

<!--     labels = var_labels3),] -->

<!-- p2_full =  -->

<!--   ggplot(all_cors, aes(x = type2, y = spearmans_r, color = species2))+ -->

<!--     geom_hline(yintercept = c(-1,0,1), linetype = "solid", color = "grey50", linewidth = 1)+ -->

<!--     geom_boxplot(position = position_dodge())+ -->

<!--     labs(x = "", y = "Spearmans Correlation")+ -->

<!--     theme_classic()+ -->

<!--     facet_wrap(~variable2, ncol = 3)+ -->

<!--     coord_cartesian(ylim = c(-1,1))+ -->

<!--     # fully remove stip -->

<!--     theme( -->

<!--           legend.position = "bottom")+ -->

<!--     #move y axis from left side to right side -->

<!--     scale_y_continuous(position = "left", breaks = c(-1,0,1))+ -->

<!--     scale_color_manual(name = "PFT", values = c(pft_cols, all = "black")) -->

<!-- p2_full -->

<!-- ``` -->

### Figure 5 - Explainable AI

<div id="fig-figure_5">

![](figures/fig-figure_5-1.png)

**Figure** 4: Simulated relative growth response to available light and
tree size (dbh) for 5 PFTs for the hybrid model (a) and the process
model (b). Visualizations are based on accumulated local effect plots.
Uncertainties for Hybrid-FINN (A) are based on Monte-Carlo Dropout. Swp
is the observed soil water potential.

</div>

<!-- ## Figure 6 (competition) -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_6 -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 10 -->

<!-- #| fig-height: 10 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "" -->

<!-- ba_dt =  -->

<!--   data.table( -->

<!--     expand.grid( -->

<!--       list( -->

<!--         trees = 1:1000, -->

<!--         dbh = seq(0, 100,1) -->

<!--       ) -->

<!--     ) -->

<!--   ) -->

<!-- ba_dt$rowID = 1:nrow(ba_dt) -->

<!-- ba_dt[, ba := BA_stand(dbh, trees,1), by = rowID] -->

<!-- p_ba = -->

<!--   ggplot(ba_dt, aes(x = dbh, y = trees, fill = ba)) + -->

<!--   geom_tile() + -->

<!--   scale_fill_gradientn( -->

<!--     colours = rev(viridis::viridis(100)), -->

<!--     limits = c(0, 100), -->

<!--     na.value = rev(viridis::viridis(100))[100], -->

<!--     breaks = seq(0, 100, 20), -->

<!--     labels = c(seq(0, 80, 20),">100"), -->

<!--     name = "Basal Area\n[m2/ha]" -->

<!--   )+ -->

<!--   xlab("dbh [cm]")+ -->

<!--   ylab("trees [ha]")+ -->

<!--   theme_classic() -->

<!-- height_dt = -->

<!--   data.table( -->

<!--     expand.grid( -->

<!--       list( -->

<!--         dbh = seq(0, 300,1), -->

<!--         parHeight = seq(0,1,0.2) -->

<!--       ) -->

<!--     ) -->

<!--   ) -->

<!-- height_dt$rowID = 1:nrow(height_dt) -->

<!-- height_dt[, height := height(dbh, parHeight), by = rowID] -->

<!-- p_height <- -->

<!--   ggplot(height_dt, aes(x = dbh, y = height, color = factor(parHeight, ordered = T)))+ -->

<!--   ylab("height(dbh,parHeight)")+ -->

<!--   xlab("dbh [cm]")+ -->

<!--   geom_line()+ -->

<!--   theme_classic()+ -->

<!--   scale_color_discrete(name = "parHeight") -->

<!-- ## Competition -->

<!-- parHeight_vec = c(0.7) -->

<!-- parCompStr = c(0.1,0.2,0.3) -->

<!-- parComp = cbind(parHeight_vec, parCompStr) -->

<!-- cohort_layers_dt <- data.table( -->

<!--   name = c("understory", "midstory", "overstory"), -->

<!--   trees = c(1000, 300, 50), -->

<!--   dbh = c(5, 30, 90), -->

<!--   dbh_shape = c(2, 4, 5) -->

<!-- ) -->

<!-- cohort_layers1_dt <- data.table() -->

<!-- for(i_c in 1:nrow(cohort_layers_dt)){ -->

<!--   dbh = cohort_layers_dt$dbh[i_c] -->

<!--   trees = cohort_layers_dt$trees[i_c] -->

<!--   dbh_shape = cohort_layers_dt$dbh_shape[i_c] -->

<!--   temp_cohort_df1 = data.table(FINN::rweibull_cohorts(trees = trees, dbh_shape = dbh_shape, dbh_scale = dbh, dbh_class_range = 0.1, species = 1, siteID = 1)) -->

<!--   temp_cohort_df1$cohortID = paste0(i_c,"_",temp_cohort_df1$cohortID) -->

<!--   temp_cohort_df1$trees_ha = trees -->

<!--   temp_cohort_df1$mean_dbh = dbh -->

<!--   temp_cohort_df1$dbh_shape = dbh_shape -->

<!--   temp_cohort_df1$name = cohort_layers_dt$name[i_c] -->

<!--   cohort_layers1_dt <- rbind(cohort_layers1_dt,temp_cohort_df1) -->

<!-- } -->

<!-- cohort_layers2_dt <- data.table() -->

<!-- counter = 0 -->

<!-- for(sp in 1:nrow(parComp)){ -->

<!--   counter = counter + 1 -->

<!--   temp_cohort_df2 = copy(cohort_layers1_dt) -->

<!--   temp_cohort_df2$siteID = counter -->

<!--   temp_cohort_df2$species = sp -->

<!--   temp_cohort_df2$parCompStr = parCompStr[sp] -->

<!--   cohort_layers2_dt <- rbind(cohort_layers2_dt, temp_cohort_df2) -->

<!-- } -->

<!-- cohort_layers2_dt[,cohortID := as.integer(as.factor(cohortID)), by = siteID] -->

<!-- cohort_layers2_dt$cohortID = as.integer(cohort_layers2_dt$cohortID) -->

<!-- cohort_layers2_dt = cohort_layers2_dt[order(cohortID)] -->

<!-- cohort = CohortMat$new(obs_df = cohort_layers2_dt) -->

<!-- basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = 1) -->

<!-- light = competition(dbh = cohort$dbh, species = cohort$species, -->

<!--                     trees = cohort$trees, parComp = torch::torch_tensor(parComp), -->

<!--                     h = NULL, patch_size_ha = 1) -->

<!-- for(i_site in unique(cohort_layers2_dt$siteID)){ -->

<!--   cohort_layers2_dt[siteID == i_site, basal_area2 := torch::as_array(basal_area)[i_site,1,1:nrow(cohort_layers2_dt[siteID == i_site])]] -->

<!--   cohort_layers2_dt[siteID == i_site, light2 := torch::as_array(light)[i_site,1,1:nrow(cohort_layers2_dt[siteID == i_site])]] -->

<!-- } -->

<!-- cohort_layers2_dt[,":="(ba = round(sum(basal_area2))), by = .(siteID, name)] -->

<!-- cohort_layers2_dt[,":="(dbh_density = dweibull(dbh,shape = dbh_shape, scale = mean_dbh)), by = .(name)] -->

<!-- cohort_layers2_dt[,name := factor(name, levels = c("understory", "midstory", "overstory"), ordered = T),] -->

<!-- label_dt <- cohort_layers2_dt[,.(mean_dbh = mean(dbh), ba = round(sum(basal_area2)), max_density = max(dbh_density)), by = .(siteID, name)] -->

<!-- p_cohorts = ggplot(cohort_layers2_dt)+ -->

<!--   geom_histogram(aes(x = dbh, fill = name, y = ..density..), alpha = 0.3)+ -->

<!--   geom_line(aes(x = dbh, y = dbh_density, color = name), linewidth = 1)+ -->

<!--   # facet_wrap(~siteID)+ -->

<!--   geom_text(data = label_dt, aes(x = mean_dbh, y = max_density+0.03, label = paste0(name,"\n"," ba = ", ba), color = name), size = 5)+ -->

<!--   guides( -->

<!--     color = guide_legend(title = "layer", override.aes = list(label = "")), -->

<!--     fill = guide_legend(title = "layer") -->

<!--     )+ -->

<!--   xlab("dbh [cm]") + -->

<!--   theme_bw()+ -->

<!--   coord_cartesian(ylim = c(0, max(label_dt$max_density)+0.05))+ -->

<!--   theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) -->

<!-- p_light = ggplot(cohort_layers2_dt)+ -->

<!--   geom_boxplot(aes(x = name, y = light2, color = name), linewidth = 1)+ -->

<!--   facet_wrap(~factor(siteID, levels = 1:3, labels = paste("parCompStr =", c(parComp[,2]))))+ -->

<!--   # add axis labels -->

<!--   xlab("layer") + -->

<!--   ylab("light [%]") + -->

<!--   theme_classic()+ -->

<!--   theme(legend.position = "none", axis.title.x = element_blank()) -->

<!-- gridExtra::grid.arrange( -->

<!--   grobs = list( -->

<!--     label_grob_f("a"), label_grob_f("b"), label_grob_f("c"), label_grob_f("d"), -->

<!--     p_ba, p_height, p_cohorts, p_light -->

<!--   ),  -->

<!--   layout_matrix = rbind( -->

<!--     c(1,2), -->

<!--     c(5,6), -->

<!--     c(3,NA), -->

<!--     c(7,7), -->

<!--     c(4,NA), -->

<!--     c(8,8) -->

<!--   ), -->

<!--   # include letters to each plot -->

<!--   widths = c(0.6,0.6),  -->

<!--   heights = c(0.05,0.95,0.05,0.95,0.05,0.95) -->

<!-- ) -->

<!-- ``` -->

<!-- ## Figure 7 (scaled sigmoids) -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_7 -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 11 -->

<!-- #| fig-height: 5 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "" -->

<!-- light = seq(0,1,0.01) -->

<!-- regLight = seq(0.1, 0.9, 0.2) -->

<!-- growthLight = regLight -->

<!-- regLight_f <- function(light,regLight, k = 10) { -->

<!--   regP = (1 / (1 + exp(-k * (light - regLight))) - 1 / (1 + exp(k * regLight))) / (1 - 1 / (1 + exp(k * (1 - regLight)))) -->

<!--   return(regP) -->

<!--   } -->

<!-- growthShade_f <- function(light, growthLight, k = 10) { -->

<!--   shade = ((1 / (1 + exp(-k * (light - growthLight))) - 1 / (1 + exp(k * growthLight))) / -->

<!--              (1 / (1 + exp(-k * (1 - growthLight))) - 1 / (1 + exp(k * growthLight)))) -->

<!--   return(shade) -->

<!--   } -->

<!-- reg_dt <-  -->

<!--   data.table( -->

<!--       process = "regeneration", -->

<!--       light = rep(light, length(regLight)), -->

<!--       regLight = factor(rep(regLight, each = length(light)), ordered = T), -->

<!--       response = as.vector(sapply(regLight, function(x) regLight_f(light, x), simplify = TRUE)) -->

<!--     ) -->

<!-- growth_dt <- data.table( -->

<!--       process = "growth", -->

<!--       light = rep(light, length(growthLight)), -->

<!--       growthLight = factor(rep(growthLight, each = length(light)), ordered = T), -->

<!--       response = as.vector(sapply(growthLight, function(x) growthShade_f(light, x), simplify = TRUE)) -->

<!--     ) -->

<!-- reg_p <- ggplot(reg_dt, aes(light, response, colour = regLight)) + -->

<!--   geom_line(size = 1) + -->

<!--   labs(x = "Light", -->

<!--        y = "Scaled light response of (regeneration)") + -->

<!--   theme_classic() -->

<!-- # reg_p -->

<!-- growth_p <- ggplot(growth_dt, aes(light, response, colour = growthLight)) + -->

<!--   geom_line(size = 1) + -->

<!--   labs(x = "Light", -->

<!--        y = "Scaled light response of (growth)") + -->

<!--   theme_classic() -->

<!-- # growth_p -->

<!-- gridExtra::grid.arrange(reg_p+ggtitle("a"), growth_p+ggtitle("b"), ncol = 2) -->

<!-- ``` -->

<!-- ## Appendix -->

<!-- ### emergent patterns -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_emerg_patterns -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 9 -->

<!-- #| fig-height: 9 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Annual growth and mortality rates of trees vs. dbh. Observed values (red line) and simulated values from the process model and the hybrid model (cyan line). The lines show fitted gams to the data and visualize how well the simulated patterns match the observed patterns. Note that the observed data only contains 0 and 1 values for mortality while the simulated values show the estimated mortality probability." -->

<!-- # load models -->

<!-- process_m <- torch::torch_load(fitted_models[1]) -->

<!-- hybrid_m <- torch::torch_load(fitted_models[2]) -->

<!-- if(file.exists("figures/emergent-patterns.csv")){ -->

<!--   p_dat <- fread("figures/emergent-patterns.csv") -->

<!-- }else{ -->

<!--   cohort_dt <- fread("data/BCI/noSplits/pft-period7-25patches/initial_cohorts1985.csv") -->

<!--   env_dt <- fread("data/BCI/noSplits/pft-period7-25patches/env_dt.csv") -->

<!--   cohort1 <- CohortMat$new(obs_df = cohort_dt) -->

<!--   plot(process_m$parameters_r$par_regeneration_unconstrained) -->

<!--   obs_trees <- fread("data/BCI/data-cleaning/pft/all_trees.csv") -->

<!--   proc_sim = process_m$simulate(env = env_dt, init_cohort = cohort1, patches = 25, patch_size = 0.1, debug = T) -->

<!--   hybr_sim = hybrid_m$simulate(env = env_dt, init_cohort = cohort1, patches = 25, patch_size = 0.1, debug = T) -->

<!--   proc_dat = data.table(proc_sim$wide$cohort, type = "process") -->

<!--   hybr_dat = data.table(hybr_sim$wide$cohort, type = "hybrid") -->

<!--   obs_dat = obs_trees[census > 1985 & status2 != "regeneration",.( -->

<!--     species, -->

<!--     g = relative_growth_5yr, -->

<!--     m = as.integer(status2 == "died"), -->

<!--     dbh = gdbh_cm_before, -->

<!--     dbh2 = gdbh_cm -->

<!--     )] -->

<!--   # obs_dat <- fread("data/BCI/noSplits/pft-period7-25patches/obs_dt.csv") -->

<!--   p_dat <- rbindlist(list( -->

<!--       data.table(proc_dat, dbh2 = proc_dat$dbh, obs = "simulated"), -->

<!--       data.table(hybr_dat, dbh2 = hybr_dat$dbh, obs = "simulated"), -->

<!--       data.table(obs_dat, obs = "observed", type = "hybrid"), -->

<!--       data.table(obs_dat, obs = "observed", type = "process") -->

<!--       ), fill = T -->

<!--       ) -->

<!--   fwrite(p_dat, "figures/emergent-patterns.csv") -->

<!-- } -->

<!-- p1 = ggplot()+ -->

<!--   geom_point(data = p_dat[dbh >= 1 & obs != "observed"], aes(dbh,(1+m)^(1/5)-1), alpha = 0.3)+ -->

<!--   # geom_rug(data = p_dat[dbh >= 1 & obs == "observed" & m == 1], aes(x = dbh, alpha = 0.001, color = "observed"))+ -->

<!--   geom_smooth( -->

<!--     data = p_dat[dbh >= 1], aes(dbh,(1+m)^(1/5)-1, color = obs), -->

<!--     method = "gam", method.args = list(family = "binomial") -->

<!--     )+ -->

<!--   facet_grid(factor(type, levels = c("observed", "hybrid", "process"), ordered = T)~paste0("PFT = ", species))+ -->

<!--   scale_x_log10()+ -->

<!--   theme_classic()+ -->

<!--   labs(x = "dbh [cm]", y = "Mortality [%/100]") -->

<!-- dbh_cmTOba_m2 <- function(dbh) { -->

<!--   dbh = dbh/100 -->

<!--   return(pi*dbh^2/4) -->

<!--   } -->

<!-- p2 = ggplot()+ -->

<!--   geom_point( -->

<!--     data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m)) & obs != "observed"], aes(dbh,(1+g)^(1/5)-1), alpha = 0.3)+ -->

<!--   # geom_rug(data = p_dat[dbh2 >= 1 & obs == "observed" & g > 0 & (m != 1|is.na(m))], aes(x = dbh, alpha = 0.001, color = "observed"))+ -->

<!--   geom_smooth( -->

<!--     data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m))], aes(dbh,(1+g)^(1/5)-1, color = obs), -->

<!--     method = "gam", method.args = list(family = "gaussian") -->

<!--     )+ -->

<!--   coord_cartesian(ylim = c(0,NA))+ -->

<!--   facet_grid(factor(type, levels = c("observed", "hybrid", "process"), ordered = T)~paste0("PFT = ", species))+ -->

<!--   scale_x_log10()+ -->

<!--   theme_classic()+ -->

<!--   labs(x = "dbh [cm]", y = "Growth [%/100]") -->

<!-- p3 = ggplot()+ -->

<!--   # geom_point( -->

<!--   #   data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m)) & obs != "observed"], -->

<!--   #   aes(dbh_cmTOba_m2(dbh),dbh_cmTOba_m2(dbh*((g)^(1/5)-1))), alpha = 0.3)+ -->

<!--   # geom_rug(data = p_dat[dbh2 >= 1 & obs == "observed" & g > 0 & (m != 1|is.na(m))], aes(x = dbh, alpha = 0.001, color = "observed"))+ -->

<!--   geom_smooth( -->

<!--     data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m))], -->

<!--     aes(dbh,dbh_cmTOba_m2(dbh*((g-1))), color = obs), -->

<!--     method = "gam", method.args = list(family = "gaussian"), se = F -->

<!--     )+ -->

<!--   coord_cartesian(ylim = c(0,NA))+ -->

<!--   facet_grid(factor(type, levels = c("observed", "hybrid", "process"), ordered = T)~paste0("PFT = ", species))+ -->

<!--   # scale_x_log10()+ -->

<!--   theme_classic()+ -->

<!--   labs(x = "dbh [cm]", y = "BAI [m2]") -->

<!-- max_dt <- p_dat[,.(max(dbh,na.rm = T)),by = .(species,obs)] -->

<!-- p4 = ggplot()+ -->

<!--   # geom_point( -->

<!--   #   data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m)) & obs != "observed"], -->

<!--   #   aes(dbh_cmTOba_m2(dbh),dbh_cmTOba_m2(dbh*((g)^(1/5)-1))), alpha = 0.3)+ -->

<!--   # geom_rug(data = p_dat[dbh2 >= 1 & obs == "observed" & g > 0 & (m != 1|is.na(m))], aes(x = dbh, alpha = 0.001, color = "observed"))+ -->

<!--   geom_smooth( -->

<!--     data = p_dat[dbh2 >= 1 & g > 0 & (m != 1|is.na(m))], -->

<!--     aes(dbh_cmTOba_m2(dbh),dbh_cmTOba_m2(dbh*((g-1))), color = obs), -->

<!--     method = "gam", method.args = list(family = "gaussian"), se = F -->

<!--     )+ -->

<!--   coord_cartesian(ylim = c(0,NA))+ -->

<!--   facet_grid(factor(type, levels = c("observed", "hybrid", "process"), ordered = T)~paste0("PFT = ", species))+ -->

<!--   # scale_x_log10()+ -->

<!--   theme_classic()+ -->

<!--   labs(x = "basal area [m2]", y = "BAI [m2]") -->

<!-- max_dt <- p_dat[,.(max(dbh,na.rm = T)),by = .(species,obs)] -->

<!-- leg_p <- get_legend(p1+theme(legend.title = element_blank())) -->

<!-- grid.arrange(grobs = list( -->

<!--   p1+ggtitle("a")+theme(legend.position = "none"), -->

<!--   p2+ggtitle("b")+theme(legend.position = "none"), -->

<!--   p3+ggtitle("b")+theme(legend.position = "none"), -->

<!--   leg_p -->

<!--   ), -->

<!--   layout_matrix = rbind( -->

<!--     c(1,4), -->

<!--     c(2,4), -->

<!--     c(3,4) -->

<!--   ), widths = c(0.9,0.1), heights = c(0.5,.5,0.5)) -->

<!-- ``` -->

<!-- ### response correlations -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_response-cors -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 10 -->

<!-- #| fig-height: 5 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Spearman correlations from five fold patial cross validation to asses how well a model can be trained with different combinations of response variables. To ensure a full coverage of each response variable, we simulated new data from the process model that was trained on the BCI forest inventory data. These simulations were used to calibrate models with different combinations of response variables. The correlations were derived as average from a spatially blocked five-fold cross validation." -->

<!-- cors_fullPlot = fread("results/cors_fullPlot.csv") -->

<!-- cors_fullPlotSpecies <- fread("results/cors_fullPlotSpecies.csv") -->

<!-- all_cors <- rbind( -->

<!--   data.table(cors_fullPlot, species = "all"), -->

<!--   cors_fullPlotSpecies -->

<!-- ) -->

<!-- p_cors <- cors_fullPlot[hybrid == "nohybrid"] -->

<!-- # p_cors <- cors_fullPlot[hybrid == "nohybrid"] -->

<!-- p_cors2 <- p_cors[cv != "T0S0",.( -->

<!--   r = mean(spearmans_r), -->

<!--   nresponses = length(tstrsplit(response, ".", fixed = TRUE)) -->

<!--   # nresponses = tstrsplit(response, ".", fixed = TRUE)[[1]] -->

<!--   ), by = .(variable,response=gsub(".pt","",response),scale, test_train,simreal)] -->

<!-- ggplot(p_cors2[grepl("pft", scale) & simreal == "simulated"],aes(x = variable, y = forcats::fct_reorder(response, -nresponses), fill = r)) + -->

<!--   geom_tile() + -->

<!--   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1)) + -->

<!--   labs(x = "Predicted Variable", y = "Loss Variables", fill = "Spearmans Correlation") + -->

<!--   theme_classic() + -->

<!--   facet_grid(test_train~scale, )+ -->

<!--   theme(axis.text.x = element_text(angle = 45, hjust = 1)) -->

<!-- ``` -->

<!-- ### hybrid model correlations -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_hybridvariant-cors -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 10 -->

<!-- #| fig-height: 7 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Performance of different hybrid modeling architecture and calibration setups. Each hybrid model was trained on the observed BCI forest inventory data. The spearman correlations represents the average from a spatially blocked five-fold cross validation." -->

<!-- cors_fullPlot = fread("results/cors_fullPlot.csv") -->

<!-- cors_fullPlotSpecies <- fread("results/cors_fullPlotSpecies.csv") -->

<!-- p_cors <- cors_fullPlot[response == "ba.trees.dbh.growth.mort.reg.pt" | response == ""] -->

<!-- p_cors2 <- p_cors[cv != "S0T0",.( -->

<!--   r = mean(spearmans_r, na.rm = T), -->

<!--   nresponses = length(tstrsplit(response, ".", fixed = TRUE)) -->

<!--   # nresponses = tstrsplit(response, ".", fixed = TRUE)[[1]] -->

<!--   ), by = .(variable, scale, test_train,simreal, hybrid)] -->

<!-- p_pft <- ggplot(p_cors2[grepl("pft", scale) & simreal == "real"],aes(x = variable, y = hybrid, fill = r)) + -->

<!--   geom_tile() + -->

<!--   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1)) + -->

<!--   labs(x = "Predicted Variable", y = "Loss Variables", fill = "Spearmans Correlation") + -->

<!--   theme_classic() + -->

<!--   facet_grid(test_train~scale)+ -->

<!--   theme(axis.text.x = element_text(angle = 45, hjust = 1)) -->

<!-- p_genus <- ggplot(p_cors2[grepl("genus", scale) & simreal == "real"],aes(x = variable, y = hybrid, fill = r)) + -->

<!--   geom_tile() + -->

<!--   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1)) + -->

<!--   labs(x = "Predicted Variable", y = "Loss Variables", fill = "Spearmans Correlation") + -->

<!--   theme_classic() + -->

<!--   facet_grid(test_train~scale)+ -->

<!--   theme(axis.text.x = element_text(angle = 45, hjust = 1)) -->

<!-- grid.arrange( -->

<!--   grobs = list( -->

<!--     p_pft+ggtitle("a"), p_genus+ggtitle("b") -->

<!--   ), -->

<!--   layout_matrix = rbind( -->

<!--     c(1), -->

<!--     c(2) -->

<!--   ), -->

<!--   heights = c(0.6,0.4), widths = c(1) -->

<!-- ) -->

<!-- ``` -->

<!-- ### fitted parameters -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_fitted_parameters -->

<!-- #| fig-format: png -->

<!-- #| fig-width: 10 -->

<!-- #| fig-height: 8 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Estimated parameters for process-FINN and hybrid-FINN calibrated on the BCI forest inventory data." -->

<!-- models_process = list.files("results/02_realdata/", full.names = T, recursive = T, pattern = "pft-period7-25patches_S0_T0") -->

<!-- models_hybrid = list.files("results/02_realdata_hybridTF0/", full.names = T, recursive = T, pattern = "pft-period7-25patches_S0_T0") -->

<!-- m_true_process = torch::torch_load(models_process) -->

<!-- m_true_hybrid = torch::torch_load(models_hybrid) -->

<!-- par_df = data.table() -->

<!-- par_names = c("par_competition_r", "par_growth_r", "par_mortality_r", "par_regeneration_r") -->

<!-- for(i in par_names){ -->

<!--   for(c in 1:ncol(m_true_process[[i]])){ -->

<!--     par_df = rbind( -->

<!--       par_df, -->

<!--       data.table( -->

<!--         true_process = m_true_process[[i]][,c], -->

<!--         true_hybrid = m_true_hybrid[[i]][,c], -->

<!--         par = paste0(i,c), -->

<!--         species = 1:length(m_true_process[[i]][,c]) -->

<!--         ) -->

<!--       ) -->

<!--   } -->

<!-- } -->

<!-- par_names_nn = c("nn_mortality.0.weight", "nn_growth.0.weight", "nn_regeneration.0.weight") -->

<!-- for(i in par_names_nn){ -->

<!--   for(c in 1:ncol(m_true_process$parameters_r[[i]])){ -->

<!--     par_df = rbind( -->

<!--       par_df, -->

<!--       data.table( -->

<!--         true_process = m_true_process$parameters_r[[i]][,c], -->

<!--         true_hybrid = m_true_hybrid$parameters_r[[i]][,c], -->

<!--         par = paste0(i,c), -->

<!--         species = 1:length(m_true_process$parameters_r[[i]][,c]) -->

<!--         ), fill = TRUE -->

<!--     ) -->

<!--   } -->

<!-- } -->

<!-- par_df[grepl("competition",par),process := "competition",] -->

<!-- par_df[grepl("competition_r1",par),name := "compHeight",] -->

<!-- par_df[grepl("competition_r2",par),name := "compStr",] -->

<!-- par_df[grepl("growth",par),process := "growth",] -->

<!-- par_df[grepl("growth_r1",par),name := "growthLight",] -->

<!-- par_df[grepl("growth_r2",par),name := "growthSize",] -->

<!-- par_df[grepl("mortality",par),process := "mortality",] -->

<!-- par_df[grepl("mortality_r1",par),name := "mortLight",] -->

<!-- par_df[grepl("mortality_r2",par),name := "mortSize",] -->

<!-- par_df[grepl("mortality_r3",par),name := "mortGrowth",] -->

<!-- par_df[grepl("regeneration",par),process := "regeneration",] -->

<!-- par_df[grepl("regeneration_r1",par),name := "regLight",] -->

<!-- par_df[grepl("par",par),type := "process parameter",] -->

<!-- par_df[grepl("nn_",par),type := "environment parameter",] -->

<!-- par_df[type == "environment parameter",name := gsub("weight","",tstrsplit(par, ".", fixed = T)[[3]]),] -->

<!-- # make "intercept" allways the first -->

<!-- par_df[,name := factor(name, levels = c("intercept", unique(par_df[name != "intercept"]$name))),] -->

<!-- p_dat <- melt(par_df, id.vars = c("name", "process", "species", "type"), -->

<!--                  measure.vars = c("true_process", "true_hybrid"), -->

<!--                  variable.name = "true", value.name = "value") -->

<!-- p_dat[,model := gsub("true_","",true),] -->

<!-- p_dat[,species2 := factor(species, levels = 1:5, labels = names(pft_cols)),] -->

<!-- p_proc_comp <- ggplot( -->

<!--   p_dat[type == "process parameter" & process == "competition"], -->

<!--   aes(x = name, y = value, color = factor(model), fill = factor(model)) -->

<!--   ) + -->

<!--   geom_hline(yintercept = 0)+ -->

<!--   facet_wrap(process~paste0(species2), ncol = 5, scales = "fixed") + -->

<!--   theme_classic()+ -->

<!--   geom_bar(width = 0.2, outlier.shape = NA, stat = "identity", position = position_dodge())+ -->

<!--   theme(axis.title.x = element_blank()) -->

<!-- p_proc_growth <- ggplot( -->

<!--   p_dat[type == "process parameter" & process == "growth" & model == "process"], -->

<!--   aes(x = name, y = value, -->

<!--       color = factor(model, levels = c("hybrid", "process")), -->

<!--       fill = factor(model, levels = c("hybrid", "process"))) -->

<!--   ) + -->

<!--     geom_hline(yintercept = 0)+ -->

<!--   scale_color_manual(values = c("process" = "#00BFC4"))+ -->

<!--   scale_fill_manual(values = c("process" = "#00BFC4"))+ -->

<!--   facet_wrap(process~paste0(species2), ncol = 5, scales = "fixed") + -->

<!--   theme_classic()+ -->

<!--   geom_bar(width = 0.2, outlier.shape = NA, stat = "identity", position = position_dodge())+ -->

<!--   theme(axis.title.x = element_blank()) -->

<!-- p_proc_mort <- ggplot( -->

<!--   p_dat[type == "process parameter" & process == "mortality"], -->

<!--   aes(x = name, y = value, color = factor(model), fill = factor(model)) -->

<!--   ) + -->

<!--     geom_hline(yintercept = 0)+ -->

<!--   facet_wrap(process~paste0(species2), ncol = 5, scales = "fixed") + -->

<!--   theme_classic()+ -->

<!--   geom_bar(width = 0.2, outlier.shape = NA, stat = "identity", position = position_dodge())+ -->

<!--   theme(axis.title.x = element_blank()) -->

<!-- p_proc_reg <- ggplot( -->

<!--   p_dat[type == "process parameter" & process == "regeneration"], -->

<!--   aes(x = name, y = value, color = factor(model), fill = factor(model)) -->

<!--   ) + -->

<!--     geom_hline(yintercept = 0)+ -->

<!--   facet_wrap(process~paste0(species2), ncol = 5, scales = "fixed") + -->

<!--   theme_classic()+ -->

<!--   geom_bar(width = 0.2, outlier.shape = NA, stat = "identity", position = position_dodge())+ -->

<!--   theme(axis.title.x = element_blank()) -->

<!-- p_proc_legend <- get_legend(p_proc_comp+theme(legend.position = "bottom")+guides(color = guide_legend(title = "Model"), fill = guide_legend(title = "Model"))) -->

<!-- # arrange in three columns -->

<!-- fac_levels = c("intercept","Prec","SR_kW_m2","RH_prc","T_max","T_min","swp") -->

<!-- p_env <- ggplot( -->

<!--   p_dat[type == "environment parameter"], -->

<!--        aes(x = factor(name, levels = 1:7, labels = fac_levels), y = value, color = factor(model), fill = factor(model)) -->

<!--   ) + -->

<!--   geom_hline(yintercept = 0)+ -->

<!--   facet_wrap(process~paste0(species2), ncol = 5, scales = "free_x") + -->

<!--   labs(x = "Parameter", y = "Difference")+ -->

<!--   theme_classic()+ -->

<!--   geom_bar(width = 0.7, outlier.shape = NA, stat = "identity", position = position_dodge2(width = 1.5, padding = 0.2, preserve = "single"))+ -->

<!--   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) -->

<!-- grid.arrange( -->

<!--   grobs = list( -->

<!--     p_proc_comp+theme(legend.position = "none", axis.title.y = element_blank()), -->

<!--     p_proc_growth+theme(legend.position = "none", axis.title.y = element_blank()), -->

<!--     p_proc_mort+theme(legend.position = "none", axis.title.y = element_blank()), -->

<!--     p_proc_reg+theme(legend.position = "none", axis.title.y = element_blank()), -->

<!--     p_proc_legend -->

<!--   ), -->

<!--   layout_matrix = rbind( -->

<!--     c(1), -->

<!--     c(2), -->

<!--     c(3), -->

<!--     c(4), -->

<!--     c(5) -->

<!--   ), -->

<!--   # include letters to each plot -->

<!--   heights = c(0.6,0.6,0.6,0.6,0.1), -->

<!-- ) -->

<!-- grid.arrange( -->

<!--   grobs = list( -->

<!--     p_env+theme(legend.position = "none", axis.title.y = element_blank()), -->

<!--     p_proc_legend -->

<!--   ), -->

<!--   layout_matrix = rbind( -->

<!--     c(1), -->

<!--     c(2) -->

<!--   ), -->

<!--   # include letters to each plot -->

<!--   heights = c(0.6*3.2,0.1), -->

<!-- ) -->

<!-- # gridExtra::grid.arrange( -->

<!-- #   grobs = list( -->

<!-- #     label_grob_f("a"), label_grob_f("b"), -->

<!-- #     p_proc, p_env -->

<!-- #   ), -->

<!-- #   layout_matrix = rbind( -->

<!-- #     c(1,2), -->

<!-- #     c(3,4) -->

<!-- #   ), -->

<!-- #   # include letters to each plot -->

<!-- #   widths = c(0.6,0.6), -->

<!-- #   heights = c(0.05,0.95) -->

<!-- # ) -->

<!-- ``` -->

<!-- ### Spatial holdout -->

<!-- ```{r} -->

<!-- #| label: fig-Fig_spatial_folds -->

<!-- #| fig-format: pdf -->

<!-- #| fig-width: 8 -->

<!-- #| fig-height: 3 -->

<!-- #| echo: false -->

<!-- #| silent: true -->

<!-- #| message: false -->

<!-- #| warning: false -->

<!-- #| fig-cap: "Spatial holdout folds used for the spatial cross validation." -->

<!-- spatial_folds_dt <- fread("data/BCI/noSplits/pft-period7-25patches/spatial_folds_dt.csv") -->

<!-- swp_dt <- fread("data/BCI/noSplits/pft-period7-25patches/swp_dt.csv") -->

<!-- spatial_folds_dt = merge(spatial_folds_dt, swp_dt, by = c("siteID")) -->

<!-- p_spatial2 = ggplot(spatial_folds_dt, aes(x = x_class, y = y_class, fill = factor(spatial_fold))) + -->

<!--   geom_tile()+ -->

<!--   theme_minimal() + -->

<!--   labs(title = "Spatial Folds", color = "Fold")+ -->

<!--   # facet_wrap(~factor(fold)) + -->

<!--   guides(fill = guide_legend(title = "fold"))+ -->

<!--   xlab("X")+ -->

<!--   ylab("Y")+ -->

<!--   coord_fixed() -->

<!-- p_spatial2 -->

<!-- ``` -->

## Session info

    R version 4.4.3 (2025-02-28)
    Platform: aarch64-apple-darwin20
    Running under: macOS Sequoia 15.5

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: Europe/Berlin
    tzcode source: internal

    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4      
     [5] purrr_1.1.0       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1     
     [9] tidyverse_2.0.0   gghalves_0.1.4    ggbeeswarm_0.7.2  gridExtra_2.3    
    [13] FINN_0.1.000      torch_0.15.1      data.table_1.17.8 ggplot2_3.5.2    

    loaded via a namespace (and not attached):
     [1] generics_0.1.3    lattice_0.22-6    stringi_1.8.4     hms_1.1.3        
     [5] digest_0.6.37     magrittr_2.0.3    timechange_0.3.0  evaluate_1.0.3   
     [9] fastmap_1.2.0     Matrix_1.7-2      jsonlite_1.9.1    processx_3.8.6   
    [13] safetensors_0.1.2 ps_1.9.0          mgcv_1.9-1        scales_1.3.0     
    [17] coro_1.1.0        cli_3.6.4         rlang_1.1.5       splines_4.4.3    
    [21] bit64_4.6.0-1     munsell_0.5.1     withr_3.0.2       yaml_2.3.10      
    [25] tools_4.4.3       tzdb_0.4.0        colorspace_2.1-1  vctrs_0.6.5      
    [29] R6_2.6.1          lifecycle_1.0.4   bit_4.5.0.1       vipor_0.4.7      
    [33] pkgconfig_2.0.3   beeswarm_0.4.0    callr_3.7.6       pillar_1.10.1    
    [37] gtable_0.3.6      glue_1.8.0        Rcpp_1.0.14       xfun_0.51        
    [41] tidyselect_1.2.1  rstudioapi_0.17.1 knitr_1.49        farver_2.1.2     
    [45] nlme_3.1-167      htmltools_0.5.8.1 labeling_0.4.3    rmarkdown_2.29   
    [49] compiler_4.4.3   
