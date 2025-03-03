plot_simulated_data <- function(models,
                                seed     = 123,
                                pdf_path,
                                dataset) {

  # Load required packages
  library(data.table)
  library(ggplot2)
  library(FINN)
  library(torch)

  # If output directory doesn't exist, create it
  if (!dir.exists(dirname(pdf_path))) {
    dir.create(dirname(pdf_path), recursive = TRUE)
  }

  # Define variables of interest
  vars <- c("growth", "reg", "mort", "ba", "trees", "dbh")

  # Helper: Spearman correlation for each species (returns long data.table)
  calc_species_corr <- function(dt, vars) {
    # dt has columns: growth.obs, growth.pred, reg.obs, reg.pred, ...
    # Return a data.table with columns: species, variable, spearman_corr
    wide <- dt[,
               {
                 corr_vals <- sapply(vars, function(v) {
                   obs_col  <- paste0(v, ".obs")
                   pred_col <- paste0(v, ".pred")
                   cor(.SD[[obs_col]], .SD[[pred_col]],
                       use = "complete.obs", method = "spearman")
                 })
                 as.list(corr_vals)  # convert named vector to list
               },
               by = species
    ]
    # Reshape from wide to long
    melt(wide,
         id.vars        = "species",
         variable.name  = "variable",
         value.name     = "spearman_corr")
  }

  # Helper: Spearman correlation overall (all species combined)
  calc_overall_corr <- function(dt, vars) {
    # dt has columns: growth.obs, growth.pred, ...
    # Return a data.table with columns: variable, spearman_corr
    corr_vals <- sapply(vars, function(v) {
      obs_col  <- paste0(v, ".obs")
      pred_col <- paste0(v, ".pred")
      cor(dt[[obs_col]], dt[[pred_col]],
          use = "complete.obs", method = "spearman")
    })
    data.table(variable = vars, spearman_corr = corr_vals)
  }

  # Load all models
  m_list  <- lapply(models, torch_load)

  # Extract model name prefixes
  m_names <- sapply(strsplit(models, "/"), tail, 1)
  names(m_list) <- sapply(strsplit(m_names, "_"), function(x) x[1])

  # Prepare to store correlation results across all models
  species_cors_list <- list()
  overall_cors_list <- list()

  # Open PDF device
  pdf(pdf_path, width = 10, height = 7)

  # Process each model
  i=1
  for (i in seq_along(m_list)) {
    FINN.seed(seed)

    i_name <- names(m_list)[i]
    m      <- m_list[[i]]

    # extract number from 1patch or 4patches or 3patch
    Npatches = as.integer(gsub("\\D", "", tstrsplit(i_name, "-")[[3]]))

    # Read data
    obs_dt     <- fread(paste0("data/",dataset,"/noSplits/", i_name, "/obs_dt.csv"))
    env_dt     <- fread(paste0("data/",dataset,"/noSplits/", i_name, "/env_dt.csv"))
    inicohort_files = list.files(paste0("data/",dataset,"/noSplits/", i_name), pattern = "initial_cohorts")
    initcohort_year = min(as.integer(gsub("\\D", "", inicohort_files)))
    cohorts_dt <- fread(paste0("data/",dataset,"/noSplits/", i_name, "/initial_cohorts",initcohort_year,".csv"))

    # Predict using the model
    cohort1 <- FINN::CohortMat(obs_df = cohorts_dt, sp = uniqueN(obs_dt$species))
    pred    <- m$simulate(env = env_dt, init_cohort = cohort1, patches = Npatches)
    pred_dt <- pred$wide$site
    pred_dt$pred <- "pred"

    # Correct the regeneration response
    pred_dt[, reg := reg / 0.1]
    obs_dt$pred <- "obs"

    # Combine observed & predicted
    comp_dt1 <- rbind(pred_dt, obs_dt, fill = TRUE)
    # Add "all" species row
    comp_dt1 <- rbind(
      comp_dt1,
      comp_dt1[
        ,
        .(
          species = "all",
          dbh     = mean(dbh, na.rm = TRUE),
          ba      = sum(ba),
          trees   = sum(trees),
          reg     = mean(reg, na.rm = TRUE),
          mort    = mean(mort, na.rm = TRUE),
          growth  = mean(growth, na.rm = TRUE)
        ),
        by = .(siteID, year, pred, period_length)
      ],
      fill = TRUE
    )

    # Pick 10 most abundant species
    top10_species <- obs_dt[, .(ba = sum(ba)), by = species][
      order(ba, decreasing = TRUE)
    ]$species[1:min(c(10, uniqueN(obs_dt$species)))]
    comp_dt1 <- comp_dt1[species %in% c(top10_species, "all")]

    # Melt for ggplot
    melt_dt <- melt(comp_dt1,
                    id.vars = c("siteID", "year", "species", "pred", "period_length"),
                    variable.name = "variable")
    p_dat <- melt_dt[, .(value = mean(value)), by = .(species, year, pred, variable)]

    # Plot with ggplot
    p1 <- ggplot(p_dat[variable != "r_mean_ha"],
                 aes(x = year, y = value, color = factor(species))) +
      geom_point(aes(shape = pred, size = pred)) +
      geom_line(aes(linetype = pred)) +
      scale_size_manual(values = c(pred = 2, obs = 3)) +
      facet_wrap(~variable, scales = "free_y") +
      ggtitle(i_name)
    print(p1)

    # Merge observed & predicted
    df <- merge.data.table(obs_dt, pred_dt,
                           by = c("siteID", "year", "species"),
                           suffixes = c(".obs", ".pred"))

    # Summarize across siteIDs if needed
    if (grepl("25patches", i_name)) {
      comp_allspecies_dt <- df[
        ,
        .(
          ba.obs     = sum(ba.obs) / uniqueN(siteID),
          ba.pred    = sum(ba.pred) / uniqueN(siteID),
          trees.obs  = sum(trees.obs) / uniqueN(siteID),
          trees.pred = sum(trees.pred) / uniqueN(siteID),
          dbh.obs    = mean(dbh.obs, na.rm = TRUE),
          dbh.pred   = mean(dbh.pred, na.rm = TRUE),
          reg.obs    = mean(reg.obs, na.rm = TRUE),
          reg.pred   = mean(r_mean_ha, na.rm = TRUE),
          mort.obs   = mean(mort.obs, na.rm = TRUE),
          mort.pred  = mean(mort.pred, na.rm = TRUE),
          growth.obs = mean(growth.obs, na.rm = TRUE),
          growth.pred= mean(growth.pred, na.rm = TRUE)
        ),
        by = .(siteID, year, species)
      ]
    } else {
      # e.g. "1patch"
      comp_allspecies_dt <- df[
        ,
        .(
          ba.obs     = sum(ba.obs) / uniqueN(siteID),
          ba.pred    = sum(ba.pred) / uniqueN(siteID),
          trees.obs  = sum(trees.obs) / uniqueN(siteID),
          trees.pred = sum(trees.pred) / uniqueN(siteID),
          dbh.obs    = mean(dbh.obs, na.rm = TRUE),
          dbh.pred   = mean(dbh.pred, na.rm = TRUE),
          reg.obs    = mean(reg.obs, na.rm = TRUE),
          reg.pred   = mean(r_mean_ha, na.rm = TRUE),
          mort.obs   = mean(mort.obs, na.rm = TRUE),
          mort.pred  = mean(mort.pred, na.rm = TRUE),
          growth.obs = mean(growth.obs, na.rm = TRUE),
          growth.pred= mean(growth.pred, na.rm = TRUE)
        ),
        by = .(year, species)
      ]
    }

    #### 1) Calculate Correlations ####
    # (a) By species
    species_cors_long <- calc_species_corr(comp_allspecies_dt, vars)
    # (b) Overall correlation across all species
    overall_cors_dt   <- calc_overall_corr(comp_allspecies_dt, vars)

    # Add "model" column so we can track which model each row is from
    species_cors_long[, model := i_name]
    overall_cors_dt[,   model := i_name]

    # Store in lists
    species_cors_list[[i]] <- species_cors_long
    overall_cors_list[[i]] <- overall_cors_dt

    #### 2) Barplots of correlations by species ####
    par(mfrow = c(2, 3)) # 6 variables => 2x3 layout
    for (v in vars) {
      sub_cors <- species_cors_long[variable == v]
      barplot(sub_cors$spearman_corr,
              names.arg = sub_cors$species,
              main      = paste0(i_name, "\n", v),
              xlab      = "species",
              ylab      = "spearman_cor",
              ylim      = c(-1, 1))
    }

    #### 3) Scatter plots (obs vs. pred) ####
    # (a) same x-/y-range for each variable
    par(mfrow = c(2, 3), pty = "s")
    for (v in vars) {
      var_obs  <- paste0(v, ".obs")
      var_pred <- paste0(v, ".pred")
      cor_val  <- round(
        cor(comp_allspecies_dt[[var_obs]],
            comp_allspecies_dt[[var_pred]],
            use = "complete.obs",
            method = "spearman"), 2
      )

      plot(comp_allspecies_dt[[var_obs]] ~ comp_allspecies_dt[[var_pred]],
           col  = comp_allspecies_dt$species,
           main = paste0(i_name, "\nCorr(", v, "): ", cor_val),
           xlim = range(c(comp_allspecies_dt[[var_obs]],
                          comp_allspecies_dt[[var_pred]]), na.rm = TRUE),
           ylim = range(c(comp_allspecies_dt[[var_obs]],
                          comp_allspecies_dt[[var_pred]]), na.rm = TRUE),
           asp  = 1)
      abline(0, 1)
    }

    # (b) Observed-based range
    par(mfrow = c(2, 3))
    for (v in vars) {
      var_obs  <- paste0(v, ".obs")
      var_pred <- paste0(v, ".pred")
      cor_val  <- round(
        cor(comp_allspecies_dt[[var_obs]],
            comp_allspecies_dt[[var_pred]],
            use = "complete.obs",
            method = "spearman"), 2
      )

      plot(comp_allspecies_dt[[var_obs]] ~ comp_allspecies_dt[[var_pred]],
           col  = comp_allspecies_dt$species,
           main = paste0(i_name, "\nCorr(", v, "): ", cor_val),
           xlim = range(comp_allspecies_dt[[var_obs]], na.rm = TRUE),
           ylim = range(comp_allspecies_dt[[var_obs]], na.rm = TRUE))
      abline(0, 1)
    }
  } # End per-model loop

  #### Combine correlation data for all models ####
  species_cors_all <- rbindlist(species_cors_list)
  overall_cors_all <- rbindlist(overall_cors_list)

  #### 4) Additional Summary Plots (in the same PDF) ####
  # (a) Heatmap: Overall vs. within-species (mean across species)
  species_cors_dt <- species_cors_all[,
                                      .(spearman_corr = mean(spearman_corr, na.rm = TRUE)),
                                      by = .(variable, model)
  ]
  p_dat <- rbind(
    data.table(overall_cors_all, aggregation = "overall R"),
    data.table(species_cors_dt,  aggregation = "within species R")
  )

  p_all <- ggplot(p_dat, aes(y = model, x = variable, fill = spearman_corr)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1, 1),
      name = "spearman's R"
    ) +
    theme_minimal() +
    geom_text(aes(label = round(spearman_corr, 2)), size = 3, color = "black") +
    theme(legend.position = "top") +
    facet_wrap(~aggregation, scales = "free")

  print(p_all)

  # (b) Boxplot: within-species R by variable, split by genus/period/patches
  # First, parse the model name "genus-period-patches" into three columns
  species_cors_all[, genus_pft := tstrsplit(model, "-")[[1]]]
  species_cors_all[, period    := tstrsplit(model, "-")[[2]]]
  species_cors_all[, patches   := tstrsplit(model, "-")[[3]]]

  p1 <- ggplot(species_cors_all,
               aes(x = variable, y = spearman_corr, fill = genus_pft)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position = "top") +
    coord_cartesian(ylim = c(-1, 1)) +
    facet_grid(period ~ patches, scales = "free") +
    ylab("within species\nspearman's R")

  print(p1)

  # Close PDF device
  dev.off()

  # Return the final correlation tables
  list(
    species_cors = species_cors_all,
    overall_cors = overall_cors_all
  )
}
