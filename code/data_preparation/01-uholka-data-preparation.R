library(data.table)
library(ggplot2)
library(raster)
library(sf)

out_dir0 = "data/Uholka"

if(!dir.exists(out_dir0)) {
  dir.create(out_dir0,recursive = T)
  }else{
    unlink(out_dir0, recursive = T)
    dir.create(out_dir0,recursive = T)
  }

source("code/data_preparation/uholka-data-prep-functions.R")
dbh_cmTOba_m2 <- function(dbh) {
  dbh = dbh/100
  return(pi*dbh^2/4)
  }
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# raw data cleaning ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
all_trees <- fread("raw-data/ecy2845-sup-0001-datas1/uholka_trees.csv")
# remove very rare species
# rename species to simple scientific names
all_trees[grepl("Acer platanoides", species), scientific_name := "Acer platanoides"]
all_trees[grepl("Acer pseudoplatanus", species), scientific_name := "Acer pseudoplatanus"]
all_trees[grepl("Fagus sylvatica", species), scientific_name := "Fagus sylvatica"]
all_trees[grepl("Ulmus glabra", species), scientific_name := "Ulmus glabra"]


### define species ####
# remove all other species
all_trees <- all_trees[
  scientific_name %in% c("Acer platanoides", "Acer pseudoplatanus", "Fagus sylvatica", "Ulmus glabra")
]
# assign speciesID
all_trees[, species := as.integer(as.factor(species))]

species_assigned <- unique(
  all_trees[,.(scientific_name, species)][order(species)]
)

# add entry for each year and tree_nr
all_years <- all_trees[, .(year = c(2000,2005,2010,2015)), by = .(tree_nr, species, x_local, y_local)]
all_trees <- merge(all_years, all_trees, by = c( "tree_nr", "year", "species", "x_local", "y_local"), all.x = T)

# define status
# status 1 Tree found (1) or not found (0) during inventory
# status 2 Tree alive (1) or dead (0) 
# status 3 Tree standing (1) or lying (0) 
# status 4 Only for dead trees: entire tree (1) or broken tree (snag), (0)
all_trees[status_2 == 1, status := "A"]
all_trees[status_2 == 0, status := "D"]

# calculate mean dbh in cm for each tree from dbh_1 and dbh_2
all_trees[, dbh := ((dbh_1 + dbh_2)/2)]
all_trees[,treeID := tree_nr]
all_trees[,stemID := 1]
all_trees[,census := year]
all_trees[,hom := 1.3]

all_trees_clean <- clean_tree_level(all_trees, out_dir = paste0(out_dir0, "/data-cleaning"))
all_trees_clean[status=="A", status2 := "alive",]
all_trees_clean[year != 2000 & status=="A" & is.na(status_before), status2 := "regeneration",]
all_trees_clean[status=="D" & status_before=="A", status2 := "died",]

all_trees_clean[,gx := x_local,]
all_trees_clean[,gy := y_local,]

x_length = 31.62278
y_length = 31.62278
x_range = c(725, round(930/x_length)*x_length)
x_breaks = seq(x_range[1], x_range[2], x_length)
y_range = c(120, round(600/y_length)*y_length)
y_breaks = seq(y_range[1], y_range[2], y_length)
plot(all_trees_clean$x_local, all_trees_clean$y_local, pch = 1,asp = 1)
abline(h = y_breaks, v = x_breaks, col = "red")

grid_list = create_grid(
  all_trees_cleaned = all_trees_clean, 
  x_length = x_length, 
  y_length = y_length, 
  x_range = x_range, 
  y_range = y_range, 
  patches2sites = c(3,3))

grid_dt_1patch = grid_list$grid_dt_1patch
grid_dt_9patches = grid_list$grid_dt_coarse
all_trees_grid = grid_list$all_trees_grid

clean_site_dir = paste0(out_dir0, "/data-cleaning/site")
if(!dir.exists(clean_site_dir)) dir.create(clean_site_dir, recursive = T)
fwrite(grid_dt_1patch, paste0(clean_site_dir, "/grid_dt_1patch.csv"))
fwrite(grid_dt_9patches, paste0(clean_site_dir, "/grid_dt_9patches.csv"))

# plot grid
p9 <- ggplot() +
  geom_point(
    data = all_trees_grid,
    mapping = aes(x = gx, y = gy), alpha = 0.5, size = 0.5) +
  geom_point(
    data = grid_dt_9patches,
    mapping = aes(x = x_class, y = y_class, color = factor(siteID)), alpha = 0.7, size = 5) +
  xlab("x") + ylab("y")+
  theme_classic()+
  theme(
    panel.grid.minor = element_line(color = "grey90", linetype = "solid"),
    panel.grid.major = element_line(color = "grey10", linetype = "solid"),
    legend.position = "top"
  )+
  scale_x_continuous(breaks = round(x_breaks)) +
  scale_y_continuous(breaks = round(y_breaks))+
  coord_fixed()
p9
# plot a map of all trees for each census
p_trees <- ggplot() +
  geom_point(
    data = all_trees[!is.na(status2) & status != "D" | (status == "D" & status2 == "died")],
    mapping = aes(y = gx, x = gy, color = status2), alpha = 0.5, size = 0.5) +
  # theme(legend.position = "none")+
  xlab("y") + ylab("x")+
  facet_grid(census~status2)+
  #maintain correct scaling with fixed proportions
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))+
  theme_classic()+
  theme(
    panel.grid.minor = element_line(color = "grey90", linetype = "solid"),
    panel.grid.major = element_line(color = "grey10", linetype = "solid"),
    legend.position = "top"
  )+
  scale_y_continuous(breaks = round(x_breaks)) +
  scale_x_continuous(breaks = round(y_breaks))+
  # use color blind friendly colors
  scale_color_manual(name = "Status",values = c("alive" = "black", "regeneration" = "blue", "died" = "red"))+
  coord_fixed()
p_trees
# Save as compressed TIFF
ggsave(
  filename = paste0(clean_site_dir, "/all_trees_map.tiff"),  # Output file name
  plot = p_trees,                   # ggplot object
  device = "tiff",            # File format
  height = 8.27,               # A4 width in inches (210 mm)
  width = 10,             # A4 height in inches (297 mm)
  units = "in",               # Specify units as inches
  dpi = 300,                  # High resolution
  compression = "lzw"         # Use LZW compression
)

ggsave(
  filename = paste0(clean_site_dir, "/grid_map9.tiff"),  # Output file name
  plot = p9,                   # ggplot object
  device = "tiff",            # File format
  width = 10,               # A4 width in inches (210 mm)
  height = 8.27,             # A4 height in inches (297 mm)
  units = "in",               # Specify units as inches
  dpi = 300,                  # High resolution
  compression = "lzw"         # Use LZW compression
)

all_trees_species = copy(all_trees_grid)

# Stand variables ####
all_trees_species[,nostems := 1,]
calc_stand_out_species = calculate_stand_vars(
  all_trees = all_trees_species,
  out_dir = paste0(out_dir0,"/data-cleaning/species"),
  area_of_square_ha = 0.1
  )
all_trees_species <- calc_stand_out_species$all_trees
stand_dt_species <- calc_stand_out_species$stand_dt[,-"m_old"]

i_patches = "1patch"
for(i_patches in c("1patch", "9patches")){
  if(i_patches == "1patch"){
    temp_all_trees_species <- merge(all_trees_species, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_species <- merge(stand_dt_species, grid_dt_1patch[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
  }
  if(i_patches == "9patches"){
    temp_all_trees_species <- merge(all_trees_species, grid_dt_9patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
    temp_stand_dt_species <- merge(stand_dt_species, grid_dt_9patches[,.(siteID, uniquePatch, patch)], by = c("uniquePatch"))
  }
  # initial cohorts ####
  create_init_cohorts(
    all_trees = temp_all_trees_species,
    stand_dt = temp_stand_dt_species,
    out_dir = paste0(out_dir0,"/data-cleaning/species/",i_patches),
    dbh_cm_class_size = 0.1,
    init_years = unique(all_trees_species$year)
    )
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# environment ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## meteo data ####
# get lon and lat for site from centroid of all patches
points_dt <- fread("raw-data/ecy2845-sup-0001-datas1/uholka_points.csv")
centroid_dt <- points_dt[,.(
  x = mean(x_utm34n),
  y = mean(y_utm34n)
),]
# transform from utm34 to lat and lon
centroid_dt <- st_as_sf(centroid_dt, coords = c("x", "y"), crs = 32634)
centroid_dt <- st_transform(centroid_dt, crs = 4326)
uholka_coordinates <- data.table(st_coordinates(centroid_dt))

library(nasapower)
# check https://power.larc.nasa.gov/beta/parameters/ for details
nasapower_pars = c("RH2M", "T2M", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN")

# Example: get daily max/min temperature for 5 years at lat=50, lon=8
# Fetch single point climatology for air temperature
ag_c_point <- get_power(
  community = "AG",
  pars = nasapower_pars,
  lonlat =  c(uholka_coordinates$X, uholka_coordinates$Y),
  dates = c(2000,2015),
  temporal_api = "monthly"
)

pars_metadata_dt <- data.table()
for(i_par in nasapower_pars){
  tmp_dt <- as.data.table(
    nasapower::query_parameters(
      community = "ag",
      pars = i_par,
      temporal_api = "monthly"
    )[[1]])
  tmp_dt$par = i_par
  setcolorder(tmp_dt, neworder = names(tmp_dt)[c(length(names(tmp_dt)),1:(length(names(tmp_dt))-1))])
  pars_metadata_dt <- rbind(
    pars_metadata_dt, tmp_dt)
}

meteo_dt1 <- data.table(ag_c_point)[,-c("LON","LAT")]
meteo_dt1[, year := YEAR,]
meteo_dt1 <- meteo_dt1[, -c("YEAR")]
meteo_dt2 <- melt(meteo_dt1, id.vars = c("year", "PARAMETER"), variable.name = "month")
meteo_dt3 <- dcast(meteo_dt2, year + month ~ PARAMETER, value.var = "value")
# transform month JAN, FEB etc. to integers
# make lowercase
meteo_dt3[, month := match(month, toupper(month.abb))]
meteo_dt3 <- meteo_dt3[!is.na(month)]
# add actual number of days to each month
meteo_dt3[, days := c(31,28,31,30,31,30,31,31,30,31,30,31)[month]]
str(meteo_dt3)
meteo_dt4 <- meteo_dt3[,.(
  RH2M = mean(RH2M, na.rm = T),
  RH2M_sd = sd(RH2M, na.rm = T),
  T2M_mean = mean(T2M, na.rm = T),
  T2M_sd = sd(T2M, na.rm = T),
  T2M_max = max(T2M, na.rm = T),
  T2M_min = min(T2M, na.rm = T),
  PRECTOTCORR_sum = sum(PRECTOTCORR*days, na.rm = T),
  PRECTOTCORR_sd = sd(PRECTOTCORR*days, na.rm = T),
  ALLSKY_SFC_SW_DWN = sum(ALLSKY_SFC_SW_DWN*days, na.rm = T)
), by = .(year)]

# correlation matrix of meteo_dt4
cor(meteo_dt4[, -1], use = "pairwise.complete.obs")
# only keep T2M_max T2M_min PRECTOTCORR due to correlations
meteo_dt <- meteo_dt4[,.(
  Tmax = T2M_max, 
  Tmin = T2M_min, 
  Psum = PRECTOTCORR_sum,
  Psd = PRECTOTCORR_sd,
  year)]

## awc data ####
awc_raster <- raster("raw-data/awc/nfk.eps1.tif")
grid_dt_1patch[,x_center := as.numeric(as.character(x_class)),]
grid_dt_1patch[,y_center := as.numeric(as.character(y_class)),]

awc_dt <- unique(grid_dt_1patch[,.(x_center, y_center, uniquePatch)])
# create grid_shp from x_class and y_class in grid_dt
grid_shp <- st_as_sf(awc_dt, coords = c("x_center", "y_center"), crs = NA)

awc_dt$awc = extract(awc_raster, grid_shp)

# make plot margins smaller to plot
par(mar=c(2,2,2,2))
plot(awc_raster)

p_awc <- 
  ggplot(awc_dt)+
  geom_tile(aes(x = x_center, y = y_center, fill = awc))+
  # add heat colors 
  scale_fill_viridis_c()+
  coord_fixed()+
  theme_classic()+
  theme(legend.position = "right")+
  labs(title = "Available water capacity", x = "x", y = "y")
p_awc
#save p_awc as tiff
ggsave(
  filename = paste0(out_dir0, "/data-cleaning/site/awc_map.tiff"),  # Output file name
  plot = p_awc,                   # ggplot object
  device = "tiff",            # File format
  width = 210,               # A4 width in inches (210 mm)
  height = 297,             # A4 height in inches (297 mm)
  units = "mm",               # Specify units as inches
  dpi = 300,                  # High resolution
  compression = "lzw"         # Use LZW compression
)

### awc_dt.csv ####
awc_dt_1patch <- awc_dt[,.(siteID = uniquePatch, awc)]
fwrite(awc_dt_1patch, paste0(out_dir0, "/data-cleaning/site/awc_dt_1patch.csv"))
awc_dt_9patches <- merge(
  awc_dt[,.(uniquePatch, awc, x_center, y_center)], 
  grid_dt_9patches[,.(siteID, patch, uniquePatch)], by = "uniquePatch"
  )
awc_dt_9patches <- awc_dt_9patches[,.(
  awc = mean(awc, na.rm = T)
), by = .(siteID)]
fwrite(awc_dt_9patches, paste0(out_dir0, "/data-cleaning/site/awc_dt_9patches.csv"))

## save full meteo_dt ####
### meteo_dt.csv ####
meto_dt15_out = paste0(out_dir0,"/data-cleaning/site/meteo_dt_period15.csv")
meteo_dt15 <- copy(meteo_dt[year>2000])
if(!dir.exists(dirname(meto_dt15_out))) dir.create(dirname(meto_dt15_out), recursive = T)
meteo_dt15[, year := year - 2000,]
fwrite(meteo_dt15, meto_dt15_out)

actual_census_years =  c(2000,2005,2010,2015)
meteo_dt[,year2 :=
           as.integer(as.character(factor(
             year,
             levels = min(actual_census_years):(max(actual_census_years)),
             labels = c(0,rep(1:(length(actual_census_years)-1),each = 5))))
           ),]

## save meteo_dt_period ####
meteo_dt_period =
  meteo_dt[!is.na(year2),.(
    Tmax = mean(Tmax),
    Tmin = mean(Tmin),
    Psum = mean(Psum),
    Psd = mean(Psd)
  ), by = .(year = year2)]

meto_dt3_out = paste0(out_dir0,"/data-cleaning/site/meteo_dt_period3.csv")
if(!dir.exists(dirname(meto_dt3_out))) dir.create(dirname(meto_dt3_out), recursive = T)
fwrite(meteo_dt_period, file = meto_dt3_out)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# final preparation ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

i_patches = "1patch"
i_period = "period15"
i_species = "species"
i_patches = "9patches"
for(i_species in c("species")){
  for(i_period in c("period3", "period15")){
    for(i_patches in c("1patch", "9patches")){
      out_dir = paste0(out_dir0,"/noSplits/",i_species,"-", i_period,"-",i_patches)
      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
      awc_dt = fread(paste0(out_dir0,"/data-cleaning/site/awc_dt_",i_patches,".csv"))
      grid_dt = fread(paste0(out_dir0,"/data-cleaning/site/grid_dt_",i_patches,".csv"))
      meteo_dt = fread(paste0(out_dir0,"/data-cleaning/site/meteo_dt_",i_period, ".csv"))

      if(i_species == "species") stand_dt = copy(stand_dt_species)
      if(i_species == "genus") stand_dt = copy(stand_dt_genus)

      stand_dt = merge(stand_dt, grid_dt[,.(siteID, patch, uniquePatch)], by = "uniquePatch")

      ## obs_dt.csv ####
      if(i_period == "period3"){
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
          period_length = fifelse(census == 2000, NA_integer_, 1)
          )]
      }
      if(i_period == "period15"){
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
          period_length = fifelse(census == 2000, NA_integer_, 5)
        ), ]
      }

      if(i_patches == "1patch") Npatches = 1
      if(i_patches == "9patches") Npatches = 9

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
              expand.grid(
                siteID = setdiff(obs_dt$siteID |> unique(), obs_dt[species==sp & year == yy]$siteID),
                year = yy, species = sp
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
      fwrite(awc_dt, paste0(out_dir,"/awc_dt.csv"))
      env_dt <- merge(env_dt, awc_dt[,.(siteID, awc)], by = "siteID")
      env_dt <- env_dt[year <= max(obs_dt$year)]
      for(i in names(env_dt)[!(names(env_dt) %in% c("siteID", "year"))]) env_dt[[i]] = as.numeric(scale(env_dt[[i]]))

      # adjust years/periods for correct initialization
      # initial cohort will be the only information of the initial inventory year
      obs_dt <- obs_dt[year != 1]
      obs_dt[,year := year-1,]
      env_dt <- env_dt[year != 0]
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

