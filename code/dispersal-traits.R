library(data.table)

traits <- fread("traits_yannek.csv")

genus_dt <- fread("data/BCI/data-cleaning/genus/species_assigned.csv")
pft_dt <- fread("data/BCI/data-cleaning/pft/species_assigned.csv")

traits <- merge(traits, genus_dt[,.(taxon = paste(Genus,Species), speciesID,PFT_2axes)], by = c("taxon"))

traits <- traits[,.(trait_name, trait_value, taxon, speciesID, PFT_2axes)]

library(torch)

model = torch::torch_load("results/02_realdata/genus-period7-1patch_S0_T0_ba.trees.dbh.growth.mort.reg.pt")

reg_par1 = model$par_regeneration_r
reg_par2 = model$parameters_r$nn_regeneration.0.weight

interc_dt <- data.table(
  interc = reg_par2[,1],
  speciesID = 1:nrow(reg_par2)
)
# define mean function that works with characters where either the most frequent value or the average is taken
mean_char = function(x){
  if(is.character(x)){
    if(length(unique(x)) == 1){
      return(unique(x))
    }else{
      return(base::mean(as.numeric(x)))
    }
  }else{
    return(base::mean(x))
  }
}
traits2 = traits[,.(trait_value = as.character(mean_char(trait_value))),by=.(speciesID,trait_name)]
traits2[,trait_value := as.numeric(trait_value)]
traits2[,NAs := all(!is.na(trait_value)), by=.(trait_name)]
traits2 = dcast(traits2[NAs == T], speciesID ~ trait_name, value.var = "trait_value")
dat = merge(traits2, interc_dt[,.(interc, speciesID)], by = "speciesID")
library(glmmTMB)


lm(formula = interc ~., data = data.frame(dat[,-1]))

library(gllvm)
library(sjSDM)
# sjSDM::install_sjSDM()
library(mvabund)
data(antTraits)
y <- as.matrix(antTraits$abund)
X <- as.matrix(antTraits$env) |> scale()
TR <- antTraits$traits
fitF1 <- gllvm(y = y, X = X, TR = TR, family = "negative.binomial")
coefplot(fitF1)

# Approximation via sjSDM
m = sjSDM(Y = y, env = X, family = "nbinom", se= TRUE)
summ = summary(m)
se = summ$se # standard errors
effs = coef(m)[[1]] |> t()

dim(se)
dim(effs)

library(metafor)
# Canopy is the second third row
effs[3,]
m_canopy_pilosity = rma(yi=effs[3,], sei = se[3,], mods = ~0+TR$Pilosity)
summary(m_canopy_pilosity)
# TR$Pilosity3 x Canopy is significant (as found by gllvm)

# bare ground against Polymorphism
m_bare_polymorphism = rma(yi=effs[2,], sei = se[2,], mods = ~0+TR$Polymorphism)
summary(m_bare_polymorphism)

