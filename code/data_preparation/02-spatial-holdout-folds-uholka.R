library(data.table)
library(blockCV)
library(sf)
library(gstat)
library(ggplot2)

for(i_patches in c("1patch", "9patches")){
  for(i_species in c("species")){
    for(i_period in c("period3", "period15")){
      dir0 = "data/Uholka/noSplits/"
      input_dir = paste0(dir0,i_species,"-",i_period,"-",i_patches)
      out_dir = input_dir

      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

      obs_dt <- fread(paste0(input_dir,"/obs_dt.csv"))
      env_dt <- fread(paste0(input_dir,"/env_dt.csv"))
      awc_dt <- fread(paste0(input_dir,"/awc_dt.csv"))
      grid_dt <- fread(paste0(dirname(dir0),"/data-cleaning/site/grid_dt_",i_patches,".csv"))
      awc_dt <- merge(awc_dt, grid_dt[,.(
        x_class = mean(x_class),
        y_class = mean(y_class)
      ), by = siteID], by = "siteID")

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## spatial holdout ####
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

      # Convert data.table to a spatial object (sf format)
      coordinates <- st_as_sf(awc_dt, coords = c("x_class", "y_class"), crs = NA) # Relative system has no CRS

      # Step 1: Calculate the spatial autocorrelation range using a variogram
      # Convert sf object to SpatialPointsDataFrame for gstat compatibility
      sp_data <- as(coordinates, "Spatial")
      variogram_model <- variogram(awc ~ 1, data = sp_data) # Compute experimental variogram
      # Step 2: Fit the variogram model with adjusted initial values
      # Adjust initial values for psill, range, and nugget based on variogram plot
      initial_model <- vgm(psill = var(awc_dt$awc, na.rm = TRUE), model = "Gau", range = 30, nugget = 0.1)

      # Fit the variogram model
      fitted_model <- fit.variogram(variogram_model, model = initial_model)

      # Check the fitted model parameters and plot the result
      print(fitted_model)
      plot(variogram_model, fitted_model, main = "Fitted Variogram")

      # Step 3: Use the fitted range parameter for block size in cv_spatial()
      autocorr_range <- fitted_model$range[2] # Extract the spatial range
      print(paste("Estimated spatial autocorrelation range:", autocorr_range))

      # Step 2: Use the estimated range as the block size in cv_spatial()
      set.seed(123) # Set a seed for reproducibility
      folds <- cv_spatial(
        x = coordinates,           # The spatial data
        # column = "awc",            # The variable of interest
        size = round(autocorr_range/10)*10,     # Use the autocorrelation range as block size
        k = 5,                     # Number of folds
        selection = "random",      # Randomly assign folds
        progress = TRUE            # Show progress
      )

      awc_dt$spatial_fold <- folds$folds_ids
      awc_dt$spatial_holdout <- factor(awc_dt$spatial_fold, levels = c(1:5), labels = c(rep("train",2),"test",rep("train",2)))

      p_spatial = ggplot(awc_dt, aes(x = x_class, y = y_class, fill = factor(awc, ordered = T))) +
        geom_tile()+
        theme_minimal() +
        labs(title = "Spatial Folds", color = "Fold")+
        # facet_wrap(~factor(fold)) +
        facet_wrap(~spatial_fold, ncol = 5) +
        theme(legend.position = "none")+
        xlab("X")+
        ylab("Y")+
        coord_fixed()

      obs_dt2 <- merge(obs_dt, awc_dt[,.(siteID, spatial_holdout)], by = "siteID")
      Nobs_dt1 <-
        obs_dt2[,.(
          growth_N = sum(!is.na(growth)),
          mort_N = sum(!is.na(mort)),
          reg_N = sum(!is.na(reg)),
          growth_mort_N = sum(!is.na(growth) & !is.na(mort) & reg > 0)
        ), by = .(species,spatial_holdout)]
      uniqueN(Nobs_dt1[spatial_holdout == "train" & growth_mort_N > 0]$species)
      uniqueN(Nobs_dt1[spatial_holdout == "test" & growth_mort_N > 0]$species)

      # Save as compressed TIFF
      ggsave(
        filename = paste0(out_dir,"/spatial_block.tiff"), # File name
        plot = p_spatial,                         # ggplot object
        device = "tiff",                  # File format
        compression = "lzw",              # Compression type
        dpi = 300,                        # Resolution
        width = 8,                        # Width in inches
        height = 3                        # Height in inches
      )

      spatial_folds_dt <- awc_dt[,.(siteID, spatial_fold, spatial_holdout)]
      fwrite(spatial_folds_dt, paste0(out_dir,"/spatial_folds_dt.csv"))
    }
  }
}
