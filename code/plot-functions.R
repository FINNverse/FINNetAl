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
    if (grepl("25patches", i_name) | grepl("species", i_name)) {
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


library(data.table)
library(torch)
library(FINN)

library(data.table)
library(torch)
library(FINN)

build_model_dt <- function(pt_file) {
  # -------------------------------------------------------------------------
  # 1) Load the .pt model
  # -------------------------------------------------------------------------
  cat("\nLoading model:", pt_file, "\n")
  m <- torch_load(pt_file)

  # -------------------------------------------------------------------------
  # 2) Infer where to read CSV data from
  # -------------------------------------------------------------------------
  model_dir       <- dirname(pt_file)       # e.g. "results/02_simulated"
  i_result_folder <- basename(model_dir)    # e.g. "02_simulated"

  # Decide how to handle real vs sim vs hybrid
  if (i_result_folder == "02_realdata") {
    split_location <- "CVsplits-realdata/"
  } else if (i_result_folder == "02_simulated") {
    split_location <- "CVsplits-simdata/"
  } else if (i_result_folder == "02_realdata_hybrid") {
    split_location <- "CVsplits-realdata/"
  } else {
    stop("Unrecognized folder name: '", i_result_folder,
         "'. Please handle logic for new folder names.")
  }

  # Parse the filename (e.g. "sp1_t1_ba.pt")
  model_name <- gsub("\\.pt$", "", basename(pt_file))
  parts <- tstrsplit(model_name, "_", fixed = TRUE)

  folder  <- parts[[1]][1]             # e.g. "sp1"
  dataset <- if (!grepl("species", folder)) "BCI" else "Uholka"

  data_dir <- file.path("data", dataset, split_location, folder)

  # Suppose you use 4 parts: <folder>_<cv_S>_<cv_T>_<response>
  cv_S <- parts[[2]][1]
  cv_T <- parts[[3]][1]

  # -------------------------------------------------------------------------
  # 3) Read the CSV files (wide obs, environment, cohorts)
  # -------------------------------------------------------------------------
  cat("Reading observation and environment data...\n")
  obs_dt_train <- fread(file.path(data_dir, paste0("obs_dt_", cv_S, "_", cv_T, "_train.csv")))
  obs_dt_test  <- fread(file.path(data_dir, paste0("obs_dt_", cv_S, "_", cv_T, "_test.csv")))

  env_dt_train <- fread(file.path(data_dir, paste0("env_dt_", cv_S, "_", cv_T, "_train.csv")))
  env_dt_train <- env_dt_train[, !c("splitID", "holdout", "siteID_holdout"), with = FALSE]
  env_dt_test  <- fread(file.path(data_dir, paste0("env_dt_", cv_S, "_", cv_T, "_test.csv")))
  env_dt_test  <- env_dt_test[, !c("splitID", "holdout", "siteID_holdout"), with = FALSE]

  cohorts_dt_train <- fread(file.path(data_dir, paste0("initial_cohorts_", cv_S, "_", cv_T, "_train.csv")))
  cohorts_dt_test  <- fread(file.path(data_dir, paste0("initial_cohorts_", cv_S, "_", cv_T, "_test.csv")))

  Nspecies <- max(obs_dt_train$species)
  Npatches <- max(cohorts_dt_train$patchID)

  cohorts_train <- FINN::CohortMat(cohorts_dt_train, sp = Nspecies)
  cohorts_test  <- FINN::CohortMat(cohorts_dt_test,  sp = Nspecies)

  # -------------------------------------------------------------------------
  # 4) Simulate predictions (FINN automatically provides a wide version)
  # -------------------------------------------------------------------------
  cat("Simulating predictions...\n")
  pred_train <- m$simulate(env = env_dt_train, init_cohort = cohorts_train, patches = Npatches)
  pred_test  <- m$simulate(env = env_dt_test,  init_cohort = cohorts_test,  patches = Npatches)

  # `pred_train$wide$site` and `pred_test$wide$site` are already wide
  pred_train_wide <- pred_train$wide$site
  pred_test_wide  <- pred_test$wide$site

  # -------------------------------------------------------------------------
  # 5) Merge wide predictions with wide observations (train)
  # -------------------------------------------------------------------------
  # We'll rename predicted columns with ".pred" suffix, observed columns with ".obs".
  # We assume obs_dt_train is also wide: siteID, year, species, ba, dbh, growth, mort, trees, reg, etc.

  # Copy so we can rename
  pred_train_wide <- copy(pred_train_wide)
  obs_train_wide  <- copy(obs_dt_train)

  # Identify measurement columns to rename
  measure_cols <- c("ba", "dbh", "growth", "mort", "trees", "reg")
  # Some FINN sims have exactly these columns – adapt if yours differ

  # Rename predicted columns => e.g. "ba" => "ba.pred"
  old_names_pred <- measure_cols
  new_names_pred <- paste0(measure_cols, ".pred")
  setnames(pred_train_wide, old_names_pred, new_names_pred, skip_absent = TRUE)

  # Rename observed columns => e.g. "ba" => "ba.obs"
  old_names_obs <- measure_cols
  new_names_obs <- paste0(measure_cols, ".obs")
  setnames(obs_train_wide, old_names_obs, new_names_obs, skip_absent = TRUE)

  # Merge on (siteID, year, species) – adapt if your data uses other keys
  train_merged <- merge(
    pred_train_wide,
    obs_train_wide,
    by = c("siteID", "year", "species"),
    all = TRUE
  )
  train_merged[, test_train := "train"]

  # If predicted reg is NA => set to 0
  if ("reg.pred" %in% names(train_merged)) {
    train_merged[is.na(reg.pred), reg.pred := 0]
  }

  # -------------------------------------------------------------------------
  # 6) Merge wide predictions with wide observations (test)
  # -------------------------------------------------------------------------
  pred_test_wide <- copy(pred_test_wide)
  obs_test_wide  <- copy(obs_dt_test)

  setnames(pred_test_wide, old_names_pred, new_names_pred, skip_absent = TRUE)
  setnames(obs_test_wide,  old_names_obs,  new_names_obs,  skip_absent = TRUE)

  test_merged <- merge(
    pred_test_wide,
    obs_test_wide,
    by = c("siteID", "year", "species"),
    all = TRUE
  )
  test_merged[, test_train := "test"]

  if ("reg.pred" %in% names(test_merged)) {
    test_merged[is.na(reg.pred), reg.pred := 0]
  }

  # -------------------------------------------------------------------------
  # 7) Combine train + test
  # -------------------------------------------------------------------------
  cat("Combining train + test sets...\n")
  final_dt <- rbindlist(list(train_merged, test_merged), use.names = TRUE)

  final_dt[, trees.obs_before := shift(trees.obs, 1, type = "lag"), by = .(siteID,species,test_train)]
  final_dt[, trees.pred_before := shift(trees.pred, 1, type = "lag"), by = .(siteID,species,test_train)]
  final_dt[trees.obs_before == 0 & trees.obs == 0, ":="(mort.obs = NA, growth.obs = NA),]
  final_dt[trees.pred_before == 0 & trees.pred == 0, ":="(mort.pred = NA, growth.pred = NA),]
  final_dt[,reg.pred := r_mean_ha,]
  final_dt <- final_dt[,-c("holdout", "splitID", "trees.obs_before", "trees.pred_before","siteID_holdout","r_mean_ha", "period_length")]

  # Clean up large objects
  rm(m, pred_train, pred_test,
     pred_train_wide, obs_train_wide,
     pred_test_wide,  obs_test_wide,
     train_merged, test_merged)
  gc()

  # -------------------------------------------------------------------------
  # 8) Return a named list with exactly one entry: { pt_file = final_dt }
  # -------------------------------------------------------------------------
  ret_list <- list()
  ret_list[[pt_file]] <- final_dt
  return(ret_list)
}
