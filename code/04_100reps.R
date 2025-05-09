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

log_mod <- function(x, base = exp(1)) {
  # Vorzeichenmatrix (0 bekommt hier schon das richtige Vorzeichen = 0)
  sgn <- sign(x)
  # Logarithmus von (|x| + 1); bei x = 0 ergibt das log(1) = 0
  log_part <- log(abs(x) + 1, base = base)
  # Kombination von Vorzeichen und Betrag
  result <- sgn * log_part
  # Nur zur Sicherheit: stellt sicher, dass exakt x == 0 auch exakt 0 ergibt
  result[x == 0] <- 0
  return(result)
}
log_mod(-10:10)
par(mfrow = c(1,1), mar = c(3,8,1,1))
boxplot(log_mod(diff) ~ par, data = diffs_df, horizontal = T, las = 1, cex.axis = 0.7)
abline(v = 0)


# create transformation for log modulu
library(scales) 
# -----------------------------------------------------------
# 1)  Transformation definieren
# -----------------------------------------------------------
log_mod_trans <- function(base = exp(1)) {
  trans_new(
    name      = paste0("log_mod-", format(base)),
    transform = function(x) sign(x) * log(abs(x) + 1, base = base),
    inverse   = function(y) sign(y) * (base^abs(y) - 1),
    breaks    = pretty_breaks(),      # Tick‑Positionen auf Originalskala
    domain    = c(-Inf, Inf)
  )
}

# -----------------------------------------------------------
# 2)  Beispielplot mit log_mod‐Skala, Original‑Labels
# -----------------------------------------------------------
ggplot(diffs_df, aes(y = par, x = diff)) +
  geom_boxplot() +
  labs(y = "Parameter", x = "Difference") +
  scale_x_continuous(
    trans  = log_mod_trans(),           # ← hier anwenden
    breaks = pretty_breaks(n = 7)       # (optional) schönere Tick‑Abstände
  ) +
  theme_minimal()

library(ggplot2)
library(ggbeeswarm)
library(gghalves)
library(data.table)
names(diffs_df)
diffs_dt <- data.table(diffs_df)
unique(diffs_df$par)
diffs_dt[grepl("competition",par),process := "competition",]
diffs_dt[grepl("competition_r1",par),name := "compHeight",]
diffs_dt[grepl("competition_r2",par),name := "compStr",]
diffs_dt[grepl("growth",par),process := "growth",]
# diffs_dt[grepl("nn_growth.0.weight1",par),process := "growth",]
# diffs_dt[grepl("nn_growth.0.weight1",par),name := "intercept",]
diffs_dt[grepl("growth_r1",par),name := "growthLight",]
diffs_dt[grepl("growth_r2",par),name := "growthSize",]
diffs_dt[grepl("mortality",par),process := "mortality",]
# diffs_dt[grepl("nn_mortality.0.weight1",par),process := "mortality",]
# diffs_dt[grepl("nn_mortality.0.weight1",par),name := "intercept",]
diffs_dt[grepl("mortality_r1",par),name := "mortLight",]
diffs_dt[grepl("mortality_r2",par),name := "mortSize",]
diffs_dt[grepl("mortality_r3",par),name := "mortGrowth",]
diffs_dt[grepl("regeneration",par),process := "regeneration",]
diffs_dt[grepl("regeneration_r1",par),name := "regLight",]
# diffs_dt[grepl("nn_regeneration.0.weight1",par),process := "regeneration",]
# diffs_dt[grepl("nn_regeneration.0.weight1",par),name := "intercept",]
# diffs_dt[grepl("par",par) | grepl("weight1", par),type := "process parameter",]
# diffs_dt[grepl("nn_",par) & !grepl("weight1", par),type := "environment parameter",]
diffs_dt[grepl("par",par),type := "process parameter",]
diffs_dt[grepl("nn_",par),type := "environment parameter",]

diffs_dt[type == "environment parameter",name := gsub("weight","",tstrsplit(par, ".", fixed = T)[[3]]),]

# make "intercept" allways the first
diffs_dt[,name := factor(name, levels = c("intercept", unique(diffs_dt[name != "intercept"]$name))),]

ggplot(diffs_dt[type == "process parameter"], aes(x = name, y = diff)) +
  # see::geom_violinhalf(scale = "width", flip = T)+
  geom_hline(yintercept = 0, color = "red") +
  geom_half_violin(scale = "width", fill = "blue")+
  geom_half_boxplot(side = "r", outlier.shape = NA)+
  facet_wrap(~process, ncol = 2, scales = "free") +
  labs(x = "Parameter", y = "Difference")+
  ggthemes::theme_base()

fac_levels = c("intercept","Prec","SR_kW_m2","RH_prc","T_max","T_min","swp")
ggplot(diffs_dt[type == "environment parameter"], 
       aes(x = factor(name, levels = 1:7, labels = fac_levels), y = diff)
       ) +
  # see::geom_violinhalf(scale = "width", flip = T)+
  geom_hline(yintercept = 0, color = "red") +
  geom_half_violin(scale = "width", fill = "blue")+
  geom_half_boxplot(side = "r", outlier.shape = NA)+
  facet_wrap(~process, ncol = 1, scales = "free") +
  labs(x = "Parameter", y = "Difference")+
  ggthemes::theme_base()


