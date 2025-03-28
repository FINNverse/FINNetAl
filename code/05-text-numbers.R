library(data.table)
bci_genus_dt <- fread("data/BCI/data-cleaning/genus/species_assigned.csv")
bci_pft_dt <- fread("data/BCI/data-cleaning/pft/species_assigned.csv")
uholka_species_dt <- fread("data/Uholka/noSplits/species-period3-9patches/species_assigned.csv")

dt_text = data.table(
  "N of species BCI" = uniqueN(bci_genus_dt$sp),
  "N of genus BCI" = uniqueN(bci_genus_dt$speciesID),
  "N of pft from BCI" = uniqueN(bci_pft_dt$speciesID),
  "N of species Uholka" = uniqueN(uholka_species_dt$species)
)
dt_text

dt_text = data.frame(
  "N of species BCI" = uniqueN(genus_dt$sp),
  "N of genus BCI" = uniqueN(genus_dt$speciesID),
  "N of pft from BCI" = uniqueN(pft_dt$speciesID)
  )
t(dt_text)
