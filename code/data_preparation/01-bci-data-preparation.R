library(data.table)
library(ggplot2)
library(raster)

out_dir0 = "data/BCI"

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
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# raw data cleaning ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## tree records ####

# define vector with census years
census_years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

# define variables for calculating areas
total_area_m2 = 1000*500
total_area_ha = total_area_m2/10000
area_of_square_m2 = total_area_m2/500
area_of_square_ha = area_of_square_m2/10000
# BCI data was downloaded from https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
# load all tree data from censuses 1985-2015
dt_list <- list()
for (i in 2:8) {
  # load rdata file
  load(paste0("data/raw-data/doi_10_15146_5xcp-0d46__v20190607/bci.tree/bci.tree",i,".rdata"))
  dt <- data.table(get(paste0("bci.tree",i)))
  dt$census = census_years[i]
  dt_list[[i]] <- dt
}
all_trees <- rbindlist(dt_list, fill=TRUE)

### define species ####
# note that species are redifined based on the number of observations further below
# load species lookup table
load(paste0("data/raw-data/doi_10_15146_5xcp-0d46__v20190607/bci.spptable.rdata"))
spptable <- data.table(bci.spptable)

# read xlsx with valid species from ruger et al 2020
valid_species_dt <- data.table(readxl::read_xlsx("data/raw-data/ruger2020/aaz4797_ruger_data_s1.xlsx", sheet = 2))
valid_species <- tolower(unique(valid_species_dt$sp))

# View(valid_species_dt[order(SampleSize)])
all_trees <- merge(all_trees, spptable[,.(Genus,Species, sp)], by = "sp")
all_trees <- all_trees[sp %in% valid_species]
all_trees <- merge(all_trees, valid_species_dt[,.(sp = tolower(sp), PFT_2axes)], by = "sp")
all_trees <- all_trees[order(census)]
all_trees[,year := census,]

all_trees_clean <- clean_tree_level(all_trees, out_dir = paste0(out_dir0, "/data-cleaning"))

x_length = 40
y_length = 25
grid_list = create_grid(all_trees_clean, x_length = x_length, y_length = y_length)
grid_dt_1patch = grid_list$grid_dt_1patch
grid_dt_25patches = grid_list$grid_dt_25patches
all_trees_grid = grid_list$all_trees_grid

clean_site_dir = paste0(out_dir0, "/data-cleaning/site")
if(!dir.exists(clean_site_dir)) dir.create(clean_site_dir, recursive = T)
fwrite(grid_dt_1patch, paste0(clean_site_dir, "/grid_dt_1patch.csv"))
fwrite(grid_dt_25patches, paste0(clean_site_dir, "/grid_dt_25patches.csv"))
fwrite(all_trees_clean,paste0(clean_site_dir, "/all_trees_clean.csv"))

# plot grid
p25 <- ggplot() +
  geom_point(
    data = grid_dt_25patches,
    mapping = aes(x = x_class, y = y_class, color = factor(siteID)), alpha = 0.5, size = 5) +
  xlab("x") + ylab("y")+
  theme_classic()+
  theme(
    panel.grid.minor = element_line(color = "grey90", linetype = "solid"),
    panel.grid.major = element_line(color = "grey10", linetype = "solid"),
    legend.position = "bottom"
  )+
  scale_x_continuous(breaks = seq(0, 1000, x_length), minor_breaks = seq(0, 1000, x_length)) +
  scale_y_continuous(breaks = seq(0, 500, y_length), minor_breaks = seq(0, 500, y_length))+
  coord_fixed()

# plot a map of all trees for each census
p_trees <- ggplot() +
  geom_point(
    data = all_trees[!is.na(status2) & status != "D" | (status == "D" & status2 == "died")],
    mapping = aes(x = gx, y = gy, color = status2), alpha = 0.05, size = 0.05) +
  # theme(legend.position = "none")+
  xlab("x") + ylab("y")+
  facet_grid(census~status2)+
  #maintain correct scaling with fixed proportions
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))+
  theme_classic()+
  theme(
    panel.grid.minor = element_line(color = "grey90", linetype = "solid"),
    panel.grid.major = element_line(color = "grey10", linetype = "solid"),
    legend.position = "top"
  )+
  scale_x_continuous(breaks = seq(0, 1000, 250), minor_breaks = seq(0, 1000, x_length)) +
  scale_y_continuous(breaks = seq(0, 500, 250), minor_breaks = seq(0, 500, y_length))+
  # use color blind friendly colors
  scale_color_manual(name = "Status",values = c("alive" = "black", "regeneration" = "blue", "died" = "red"))+
  coord_fixed()
# Save as compressed TIFF
ggsave(
  filename = paste0(clean_site_dir, "/all_trees_map.tiff"),  # Output file name
  plot = p_trees,                   # ggplot object
  device = "tiff",            # File format
  width = 8.27,               # A4 width in inches (210 mm)
  height = 10,             # A4 height in inches (297 mm)
  units = "in",               # Specify units as inches
  dpi = 300,                  # High resolution
  compression = "lzw"         # Use LZW compression
)

ggsave(
  filename = paste0(clean_site_dir, "/grid_map25.tiff"),  # Output file name
  plot = p25,                   # ggplot object
  device = "tiff",            # File format
  width = 10,               # A4 width in inches (210 mm)
  height = 8.27,             # A4 height in inches (297 mm)
  units = "in",               # Specify units as inches
  dpi = 300,                  # High resolution
  compression = "lzw"         # Use LZW compression
)

all_trees_pft = copy(all_trees_grid)
all_trees_pft[,species := as.integer(as.factor(PFT_2axes)),]

all_trees_genus = copy(all_trees_grid)
all_trees_genus[,species := as.integer(as.factor(Genus)),]

# Stand variables ####
calc_stand_out_pft = calculate_stand_vars(
  all_trees = all_trees_pft,
  out_dir = paste0(out_dir0,"/data-cleaning/pft"),
  area_of_square_ha = area_of_square_ha
  )
all_trees_pft <- calc_stand_out_pft$all_trees
stand_dt_pft <- calc_stand_out_pft$stand_dt[,-"m_old"]

calc_stand_out_genus = calculate_stand_vars(
  all_trees = all_trees_genus,
  out_dir = paste0(out_dir0,"/data-cleaning/genus"),
  area_of_square_ha = area_of_square_ha
  )
all_trees_genus <- calc_stand_out_genus$all_trees
stand_dt_genus <- calc_stand_out_genus$stand_dt[,-"m_old"]

for(i_patches in c("1patch", "25patches")){
  if(i_patches == "1patch"){
    temp_all_trees_pft <- merge(all_trees_pft, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_pft <- merge(stand_dt_pft, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_all_trees_genus <- merge(all_trees_genus, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_genus <- merge(stand_dt_genus, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
  }
  if(i_patches == "25patches"){
    temp_all_trees_pft <- merge(all_trees_pft, grid_dt_25patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_pft <- merge(stand_dt_pft, grid_dt_25patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_all_trees_genus <- merge(all_trees_genus, grid_dt_25patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_genus <- merge(stand_dt_genus, grid_dt_25patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
  }
  # initial cohorts ####
  create_init_cohorts(
    all_trees = temp_all_trees_pft,
    stand_dt = temp_stand_dt_pft,
    out_dir = paste0(out_dir0,"/data-cleaning/pft/",i_patches),
    dbh_cm_class_size = 0.1,
    init_years = census_years[-1]
    )
  create_init_cohorts(
    all_trees = temp_all_trees_genus,
    stand_dt = temp_stand_dt_genus,
    out_dir = paste0(out_dir0,"/data-cleaning/genus/",i_patches),
    dbh_cm_class_size = 0.1,
    init_years = census_years[-1]
    )
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# meteorological data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# downloaded from
# https://www.osti.gov/biblio/1771850
meteo_dt <- fread("data/raw-data/BCI_1985_2018c_mod_2018substituted_20210129025453.csv")

# check correlation between variables
cor(as.matrix(meteo_dt[,.(Temp_deg_C, RH_prc, SR_W_m2, Rainfall_mm_hr)]), method = "spearman")

# convert 12/31/84 20:00 to year, month, day, hour columns
meteo_dt[, c("date", "hour") := tstrsplit(DateTime, " ", fixed = TRUE)]
meteo_dt[, c("month", "day", "year") := lapply(tstrsplit(date, "/", fixed = TRUE), as.integer)]
meteo_dt[year > 80, year := as.integer(paste0(19,year)),]
meteo_dt[year < 80 & year >= 10, year := as.integer(paste0(20,year)),]
meteo_dt[year < 80 & year < 10, year := as.integer(paste0(200,year)),]
meteo_dt[, c("hour") := as.integer(gsub(":00", "", hour))]

meteo_dt = meteo_dt[,.(
  Prec = sum(Rainfall_mm_hr),
  SR_W_m2 = sum(SR_W_m2),
  Temp = mean(Temp_deg_C),
  RH_prc = mean(RH_prc)
), by = .(year, month, day)]

meteo_dt <- meteo_dt[,.(
  Prec = sum(Prec),
  SR_kW_m2 = sum(SR_W_m2)/1000,
  RH_prc = mean(RH_prc),
  T_mean = mean(Temp),
  T_max = max(Temp),
  T_min = min(Temp)
), by = .(year)]

## save full meteo_dt ####
meteo_dt <- meteo_dt[year >= 1985, .(year,Prec, SR_kW_m2, RH_prc, T_max, T_min)]
meteo_dt35 <- copy(meteo_dt)
meto_dt35_out = paste0(out_dir0,"/data-cleaning/site/meteo_dt_period35.csv")
if(!dir.exists(dirname(meto_dt35_out))) dir.create(dirname(meto_dt35_out), recursive = T)
meteo_dt35[, year := year - 1984,]
fwrite(meteo_dt35, file = meto_dt35_out)

actual_census_years =  census_years[-1]
meteo_dt[,year2 :=
           as.integer(factor(
             year,
             levels = min(actual_census_years):(max(actual_census_years)),
             labels = c(1,rep(2:(length(actual_census_years)),each = 5)))
           ),]

## save meteo_dt_period ####
meteo_dt_period =
  meteo_dt[!is.na(year2),.(
    Prec = mean(Prec),
    SR_kW_m2 = mean(SR_kW_m2),
    RH_prc = mean(RH_prc),
    # T_mean = mean(T_mean),
    T_max = mean(T_max),
    T_min = mean(T_min)
  ), by = .(year = year2)]

meto_dt7_out = paste0(out_dir0,"/data-cleaning/site/meteo_dt_period7.csv")
if(!dir.exists(dirname(meto_dt7_out))) dir.create(dirname(meto_dt7_out), recursive = T)
fwrite(meteo_dt_period, file = meto_dt7_out)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# soil water Potential ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

early_dt <- fread("data/raw-data/Kupers_et_al/Output/BCI_SWP_map_early_dry_season_regular.txt")
mid_dt <- fread("data/raw-data/Kupers_et_al/Output/BCI_SWP_map_mid_dry_season_regular.txt")
late_dt <- fread("data/raw-data/Kupers_et_al/Output/BCI_SWP_map_late_dry_season_regular.txt")
late_drougth_dt <- fread("data/raw-data/Kupers_et_al/Output/BCI_SWP_map_mid_dry_season_drought.txt")
names(early_dt)[3] <- c("SWP_early")
names(mid_dt)[3] <- c("SWP_mid")
names(late_dt)[3] <- c("SWP_late")
names(late_drougth_dt)[3] <- c("SWP_late_drought")

# merge all files
all_swp_dt <- merge(early_dt, mid_dt, by = c("x", "y"))
all_swp_dt <- merge(all_swp_dt, late_dt, by = c("x", "y"))
all_swp_dt <- merge(all_swp_dt, late_drougth_dt, by = c("x", "y"))
var(as.matrix(all_swp_dt[,.(
  SWP_early,
  SWP_mid,
  SWP_late,
  SWP_late_drought
)]))
cor(as.matrix(all_swp_dt[,.(
  SWP_early,
  SWP_mid,
  SWP_late,
  SWP_late_drought
)]))

# only use late swp because it has the highest spatial variation
late_r <- raster("data/raw-data/Kupers_et_al/Output/BCI_SWP_map_late_dry_season_regular.tif")
swp_dt_1patch <- grid_dt_1patch
swp_dt_1patch$swp <- raster::extract(late_r, cbind(swp_dt_1patch$x_class, swp_dt_1patch$y_class))
fwrite(swp_dt_1patch,  paste0(out_dir0, "/data-cleaning/site/swp_dt_1patch.csv"))

swp_dt_25patches <- merge(swp_dt_1patch[,.(uniquePatch, swp)], grid_dt_25patches, by = "uniquePatch")
swp_dt_25patches <- swp_dt_25patches[,.(
  swp = mean(swp),
  x_class = mean(x_class),
  y_class = mean(y_class)
  ), by = siteID]
fwrite(swp_dt_25patches,  paste0(out_dir0, "/data-cleaning/site/swp_dt_25patches.csv"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# final preparation ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# i_patches = "1patch"
i_period = "period7"
i_species = "genus"
i_patches = "25patches"
for(i_species in c("pft","genus")){
  for(i_period in c("period7", "period35")){
    for(i_patches in c("1patch", "25patches")){
      out_dir = paste0(out_dir0,"/noSplits/",i_species,"-", i_period,"-",i_patches)
      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
      swp_dt = fread(paste0(out_dir0,"/data-cleaning/site/swp_dt_",i_patches,".csv"))
      grid_dt = fread(paste0(out_dir0,"/data-cleaning/site/grid_dt_",i_patches,".csv"))
      meteo_dt = fread(paste0(out_dir0,"/data-cleaning/site/meteo_dt_",i_period, ".csv"))

      if(i_species == "pft") stand_dt = copy(stand_dt_pft)
      if(i_species == "genus") stand_dt = copy(stand_dt_genus)

      stand_dt = merge(stand_dt, grid_dt[,.(siteID, patch, uniquePatch)], by = "uniquePatch")

      ## obs_dt.csv ####
      if(i_period == "period7"){
        obs_dt1 = stand_dt[,.(
          year = as.integer(as.factor(census)),
          siteID,
          patch,
          species,
          ba = ba,
          dbh = dbh_mean,
          trees = trees,
          growth = g_5yr,
          mort = m_5yr,
          reg = r,
          period_length = fifelse(census == 1985, NA_integer_, 1)
          )]
      }
      if(i_period == "period35"){
        obs_dt1 = stand_dt[,.(
          year = census - min(stand_dt$census)+1,
          siteID,
          patch,
          species,
          ba = ba,
          dbh = dbh_mean,
          trees = trees,
          growth = g,
          mort = m,
          reg = r,
          period_length
        ), ]
      }

      if(i_patches == "1patch") Npatches = 1
      if(i_patches == "25patches") Npatches = 25

      obs_dt <- obs_dt1[,.(
        ba = sum(ba, na.rm = T)/Npatches,
        dbh = sum(dbh, na.rm = T)/Npatches,
        trees = sum(trees, na.rm = T)/Npatches,
        growth = sum(growth, na.rm = T)/Npatches,
        mort = sum(mort, na.rm = T)/Npatches,
        reg = sum(reg, na.rm = T)/Npatches
      ),by = .(siteID, year, species, period_length)]

      # par(mfrow = c(3,2))
      # hist(obs_dt$ba, main = paste0(i_species, " ", i_period, " ", i_patches))
      # hist(obs_dt$dbh, main = paste0(i_species, " ", i_period, " ", i_patches))
      # hist(obs_dt$trees, main = paste0(i_species, " ", i_period, " ", i_patches))
      # hist(obs_dt$growth, main = paste0(i_species, " ", i_period, " ", i_patches))
      # hist(obs_dt$mort, main = paste0(i_species, " ", i_period, " ", i_patches))
      # hist(obs_dt$reg, main = paste0(i_species, " ", i_period, " ", i_patches))

      ## expand missing obs_dt ####
      # fill sites
      # Jedes jahr alle sites!
      empty_res =
        lapply(unique(obs_dt$species), function(sp) {
          empty_sites =
            lapply(unique(obs_dt$year), function(yy) {
              expand.grid(siteID = setdiff(obs_dt$siteID |> unique(), obs_dt[species==sp & year == yy]$siteID),
                          year = yy,
                          species = sp
              )
            })
          rbindlist(empty_sites, fill = TRUE)
        })
      obs_dt = rbindlist(list(obs_dt, rbindlist(empty_res, fill = TRUE)), fill = TRUE)
      obs_dt[is.na(ba), ba := 0,]
      obs_dt[is.na(trees), trees := 0,]
      obs_dt[is.na(reg) & year !=1 , reg := 0,]
      obs_dt$period_length <- unique(obs_dt[!is.na(period_length)]$period_length)

      if(!all(table(obs_dt$species) == unique(table(obs_dt$species)))) stop("species entries not complete")
      if(!all(table(obs_dt$siteID) == unique(table(obs_dt$siteID)))) stop("site entries not complete")
      if(!all(table(obs_dt$year) == unique(table(obs_dt$year)))) stop("year entries not complete")

      if(uniqueN(obs_dt$species) != max(obs_dt$species)) stop("species not continuous")

      ## env_dt.csv ####
      ### expand missing env_dt ####
      env_dt =
        lapply(unique(obs_dt$siteID), function(site) {
          tmp = meteo_dt
          tmp$siteID = site
          return(tmp)
        }) |> rbindlist()
      fwrite(swp_dt, paste0(out_dir,"/swp_dt.csv"))
      env_dt <- merge(env_dt, swp_dt[,.(siteID, swp)], by = "siteID")
      env_dt <- env_dt[year <= max(obs_dt$year)]
      for(i in names(env_dt)[!(names(env_dt) %in% c("siteID", "year"))]) env_dt[[i]] = scale(env_dt[[i]])

      # adjust years/periods for correct initialization
      # initial cohort will be the only information of the initial inventory year
      obs_dt <- obs_dt[year != 1]
      obs_dt[,year := year-1,]
      env_dt <- env_dt[year != 1]
      env_dt[,year := as.integer(as.factor(year)),]

      fwrite(obs_dt, paste0(out_dir,"/obs_dt.csv"))
      fwrite(env_dt, paste0(out_dir,"/env_dt.csv"))
      #copy init_cohorts
      # init_cohort_path =list.files(paste0(out_dir0, "/data-cleaning/", i_species), pattern = "initial_cohorts", full.names = T)
      init_cohort_path =list.files(paste0(out_dir0, "/data-cleaning/", i_species,"/",i_patches), pattern = "initial_cohorts", full.names = T)
      file.copy(init_cohort_path, out_dir, overwrite = T)
    }
  }
}

