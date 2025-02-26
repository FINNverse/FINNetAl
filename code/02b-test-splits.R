library(data.table)

#' Validate the cross-validation split files for a given i_folder.
#'
#' This function reads the six CSV files produced by your splitting code
#' (obs_dt_{S,T}_train.csv, obs_dt_{S,T}_test.csv, env_dt_{S,T}_train.csv,
#'  env_dt_{S,T}_test.csv, initial_cohorts_{S,T}_train.csv,
#'  initial_cohorts_{S,T}_test.csv) for each (S,T) combination,
#' and performs a series of checks:
#'   1) obs_dt:
#'       - siteIDs should start at 1 and be consecutive.
#'       - years should start at 1 if "period7" in the folder name, else 5 if "period35".
#'       - number of siteIDs should match the subset in spatial_folds_dt.
#'       - number of years should match the subset in temporal_folds_dt.
#'   2) env_dt:
#'       - years should start at 1.
#'       - siteIDs count should match spatial_folds_dt.
#'       - number of years should match temporal_folds_dt.
#'   3) initial_cohorts_{train,test}:
#'       - The siteIDs must match those in the obs/env “train” or “test” sets.
#'       - The census (or year) in the cohort file must match init_year from temporal_folds_dt.
#'
#' The function returns a data.table summarizing pass/fail for each fold.
#'
#' @param i_folder     (char) e.g. "pft-period7-25patches"
#' @param out_dir0     (char) path to the folder with the final CV-split CSVs
#' @param noSplitsDir  (char) path to the folder with spatial_folds_dt_applied.csv and temporal_folds_dt_applied.csv
#'
#' @return data.table with one row per (S,T) fold, containing logical checks
validate_CV_files <- function(i_folder,
                              out_dir0    = "data/BCI/CVsplits-realdata",
                              noSplitsDir = "data/BCI/noSplits") {

  #---------------------------------------------------------------------------
  # 1) Read the "applied" folds you saved
  #---------------------------------------------------------------------------
  spatial_folds_path  <- file.path(noSplitsDir, i_folder, "spatial_folds_dt.csv")
  temporal_folds_path <- file.path(noSplitsDir, i_folder, "temporal_folds_dt_applied.csv")

  if(!file.exists(spatial_folds_path) || !file.exists(temporal_folds_path)) {
    stop("Cannot find spatial_folds_dt_applied.csv or temporal_folds_dt_applied.csv for folder: ", i_folder)
  }

  spatial_folds_dt  <- fread(spatial_folds_path)
  temporal_folds_dt <- fread(temporal_folds_path)

  # Unique fold IDs
  spatial_splits  <- 1:5
  temporal_splits <- sort(unique(temporal_folds_dt$splitID))

  # In your creation code, you also used "fold=0" for all sites,
  # and "splitID=0" for all years:
  spatial_splits  <- c(0, spatial_splits)
  temporal_splits <- c(0, temporal_splits)

  #---------------------------------------------------------------------------
  # 2) Determine the correct expected min(year) for obs_dt
  #    depending on "period7" vs. "period35"
  #---------------------------------------------------------------------------
  # We check the folder name to see if it has "period7" or "period35"
  # If "period7", obs_dt$year should start == 1
  # If "period35", obs_dt$year should start == 5
  check_obs_year_start <- if (grepl("period7", i_folder)) 1 else 5

  # For environment data we keep the old rule: it always starts at 1
  # because env_dt$year = as.integer(as.factor(.)) was done.

  # If you are 100% sure it must be only "period7" or "period35" folder names,
  # the above line is fine. Otherwise you might add a fallback or extra checks.

  #---------------------------------------------------------------------------
  # 3) Prepare a table to store results
  #---------------------------------------------------------------------------
  results <- data.table()

  #---------------------------------------------------------------------------
  # 4) Loop over each combination of (spatial fold, temporal fold)
  #---------------------------------------------------------------------------
  for(i_S in spatial_splits) {
    for(i_T in temporal_splits) {
      cv_label <- paste0("S", i_S, "_T", i_T)

      #-----------------------------------------------------------------------
      # Identify the six CSV files we want to check
      #-----------------------------------------------------------------------
      obs_train_file  <- file.path(out_dir0, i_folder, paste0("obs_dt_", cv_label, "_train.csv"))
      obs_test_file   <- file.path(out_dir0, i_folder, paste0("obs_dt_", cv_label, "_test.csv"))
      env_train_file  <- file.path(out_dir0, i_folder, paste0("env_dt_", cv_label, "_train.csv"))
      env_test_file   <- file.path(out_dir0, i_folder, paste0("env_dt_", cv_label, "_test.csv"))
      init_train_file <- file.path(out_dir0, i_folder, paste0("initial_cohorts_", cv_label, "_train.csv"))
      init_test_file  <- file.path(out_dir0, i_folder, paste0("initial_cohorts_", cv_label, "_test.csv"))

      # Check existence:
      files_exist <- file.exists(obs_train_file)  &
        file.exists(obs_test_file)   &
        file.exists(env_train_file)  &
        file.exists(env_test_file)   &
        file.exists(init_train_file) &
        file.exists(init_test_file)

      # If any file is missing, log failure and continue
      if(!files_exist) {
        results <- rbind(
          results,
          data.table(
            i_folder = i_folder,
            cv_label = cv_label,
            pass_all = FALSE,
            message  = "One or more files are missing."
          ),
          fill = TRUE
        )
        next
      }

      #-----------------------------------------------------------------------
      # Read in the data
      #-----------------------------------------------------------------------
      obs_train_dt  <- fread(obs_train_file)
      obs_test_dt   <- fread(obs_test_file)
      env_train_dt  <- fread(env_train_file)
      env_test_dt   <- fread(env_test_file)
      init_train_dt <- fread(init_train_file)
      init_test_dt  <- fread(init_test_file)

      #-----------------------------------------------------------------------
      # A) obs_train_dt / obs_test_dt checks
      #-----------------------------------------------------------------------

      # 1) siteIDs start from 1 and form a consecutive sequence
      check_obs_train_siteID_min  <- (min(obs_train_dt$siteID) == 1)
      check_obs_train_siteID_seq  <- all(sort(unique(obs_train_dt$siteID)) ==
                                           seq_len(max(obs_train_dt$siteID)))

      check_obs_test_siteID_min   <- (min(obs_test_dt$siteID) == 1)
      check_obs_test_siteID_seq   <- all(sort(unique(obs_test_dt$siteID)) ==
                                           seq_len(max(obs_test_dt$siteID)))

      # 2) Years: for "period7" => min(year)==1, for "period35" => min(year)==5
      check_obs_train_year_ok <- (min(obs_train_dt$year) == check_obs_year_start)
      check_obs_test_year_ok  <- (min(obs_test_dt$year) == check_obs_year_start)

      # 3) Number of unique siteIDs matches the # of sites in train/test (spatial folds)
      if(i_S == 0) {
        # "fold 0" => all sites in both train and test
        n_train_sites <- length(unique(spatial_folds_dt$siteID))
        n_test_sites  <- length(unique(spatial_folds_dt$siteID))
      } else {
        n_train_sites <- spatial_folds_dt[spatial_fold != i_S, uniqueN(siteID)]
        n_test_sites  <- spatial_folds_dt[spatial_fold == i_S, uniqueN(siteID)]
      }
      check_obs_train_nSites <- (uniqueN(obs_train_dt$siteID) == n_train_sites)
      check_obs_test_nSites  <- (uniqueN(obs_test_dt$siteID)  == n_test_sites)

      # 4) Number of years matches temporal folds.
      #    If i_T=0 => entire range; otherwise check the relevant range in temporal_folds_dt.
      if(i_T == 0) {
        n_years_obs <- length(seq(min(temporal_folds_dt$obs_start),
                                  max(temporal_folds_dt$obs_end),
                                  by = 1))  # period_length=1 step in the original code for the "obs" index
        # But double check your original code. If you used period_length=5 for "period35"
        # you'd do `by=5`. For demonstration, we do `by=1` here or adapt to your code logic.

        # If in your creation script you used `seq(..., by=period_length)`,
        # then we need to detect period_length from the folder or from the data.

        # We'll do something minimal here:
        if (grepl("period7", i_folder)) {
          # "period7" => obs taken every 5 years => step=5 ??? or step=1 ???
          # Actually, the code for obs says:
          #   obs_train_dt[, year := year-min(year)+period_length]
          # If period_length=1 for period7 => increments by 1 for each obs step
          # So we do `by=1` for counting.
          n_years_obs <- length(seq(min(temporal_folds_dt$obs_start),
                                    max(temporal_folds_dt$obs_end),
                                    by = 1))
        } else {
          # "period35" => period_length=5
          n_years_obs <- length(seq(min(temporal_folds_dt$obs_start),
                                    max(temporal_folds_dt$obs_end),
                                    by = 5))
        }

        check_obs_train_nYears <- (uniqueN(obs_train_dt$year) == n_years_obs)
        check_obs_test_nYears  <- (uniqueN(obs_test_dt$year)  == n_years_obs)
      } else {
        Tfolds <- temporal_folds_dt[splitID == i_T]

        # train portion
        train_part <- Tfolds[holdout == "train"]
        # in your script:
        # train_years_obs = seq(train_part$obs_start, train_part$obs_end, period_length)

        if(nrow(train_part) > 0) {
          if(grepl("period7", i_folder)) {
            train_seq <- seq(train_part$obs_start, train_part$obs_end, by = 1)
          } else {
            train_seq <- seq(train_part$obs_start, train_part$obs_end, by = 5)
          }
          check_obs_train_nYears <- (uniqueN(obs_train_dt$year) == length(train_seq))
        } else {
          check_obs_train_nYears <- (nrow(obs_train_dt) == 0)
        }

        # test portion
        test_part <- Tfolds[holdout == "test"]
        if(nrow(test_part) > 0) {
          if(grepl("period7", i_folder)) {
            test_seq <- seq(test_part$obs_start, test_part$obs_end, by = 1)
          } else {
            test_seq <- seq(test_part$obs_start, test_part$obs_end, by = 5)
          }
          check_obs_test_nYears <- (uniqueN(obs_test_dt$year) == length(test_seq))
        } else {
          check_obs_test_nYears <- (nrow(obs_test_dt) == 0)
        }
      }

      #-----------------------------------------------------------------------
      # B) env_train_dt / env_test_dt checks
      #-----------------------------------------------------------------------
      # 1) Years start from 1
      check_env_train_year_start1 <- (min(env_train_dt$year) == 1)
      check_env_test_year_start1  <- (min(env_test_dt$year) == 1)

      # 2) Number of siteIDs
      check_env_train_nSites <- (uniqueN(env_train_dt$siteID) == n_train_sites)
      check_env_test_nSites  <- (uniqueN(env_test_dt$siteID)  == n_test_sites)

      # 3) Number of years (env has yearly increments in your script)
      if(i_T == 0) {
        n_years_env <- length(seq(min(temporal_folds_dt$env_start),
                                  max(temporal_folds_dt$env_end), by = 1))
        check_env_train_nYears <- (uniqueN(env_train_dt$year) == n_years_env)
        check_env_test_nYears  <- (uniqueN(env_test_dt$year)  == n_years_env)
      } else {
        Tfolds <- temporal_folds_dt[splitID == i_T]

        train_part <- Tfolds[holdout == "train"]
        if(nrow(train_part) > 0) {
          train_seq_env <- seq(train_part$env_start, train_part$env_end, by = 1)
          check_env_train_nYears <- (uniqueN(env_train_dt$year) == length(train_seq_env))
        } else {
          check_env_train_nYears <- (nrow(env_train_dt) == 0)
        }

        test_part <- Tfolds[holdout == "test"]
        if(nrow(test_part) > 0) {
          test_seq_env <- seq(test_part$env_start, test_part$env_end, by = 1)
          check_env_test_nYears <- (uniqueN(env_test_dt$year) == length(test_seq_env))
        } else {
          check_env_test_nYears <- (nrow(env_test_dt) == 0)
        }
      }

      #-----------------------------------------------------------------------
      # C) initial_cohorts_ checks
      #-----------------------------------------------------------------------
      # 1) siteIDs match obs/env
      train_sites_cohort <- unique(init_train_dt$siteID)
      train_sites_obs    <- unique(obs_train_dt$siteID)
      train_sites_env    <- unique(env_train_dt$siteID)
      check_init_train_sites <- setequal(train_sites_cohort,
                                         union(train_sites_obs, train_sites_env))

      test_sites_cohort <- unique(init_test_dt$siteID)
      test_sites_obs    <- unique(obs_test_dt$siteID)
      test_sites_env    <- unique(env_test_dt$siteID)
      check_init_test_sites <- setequal(test_sites_cohort,
                                        union(test_sites_obs, test_sites_env))

      # 2) census year matches init_year from temporal_folds_dt
      check_init_train_censusyear <- TRUE
      check_init_test_censusyear  <- TRUE

      if("census" %in% names(init_train_dt)) {
        train_census_years <- unique(init_train_dt$census)
      } else {
        train_census_years <- NA
      }
      if("census" %in% names(init_test_dt)) {
        test_census_years <- unique(init_test_dt$census)
      } else {
        test_census_years <- NA
      }

      if(i_T == 0) {
        # entire range => your code picks init_year=1985 for both
        check_init_train_censusyear <- all(train_census_years == 1985)
        check_init_test_censusyear  <- all(test_census_years  == 1985)
      } else {
        Tfolds <- temporal_folds_dt[splitID == i_T]
        tr_init_year <- unique(Tfolds[holdout == "train"]$init_year)
        te_init_year <- unique(Tfolds[holdout == "test"]$init_year)

        if(length(tr_init_year) == 1 && !is.na(train_census_years[1])) {
          check_init_train_censusyear <- all(train_census_years == tr_init_year)
        }
        if(length(te_init_year) == 1 && !is.na(test_census_years[1])) {
          check_init_test_censusyear <- all(test_census_years == te_init_year)
        }
      }

      #-----------------------------------------------------------------------
      # Combine all checks into a single "pass_all" flag
      #-----------------------------------------------------------------------
      pass_all <- all(
        check_obs_train_siteID_min,
        check_obs_train_siteID_seq,
        check_obs_test_siteID_min,
        check_obs_test_siteID_seq,
        check_obs_train_year_ok,
        check_obs_test_year_ok,
        check_obs_train_nSites,
        check_obs_test_nSites,
        check_obs_train_nYears,
        check_obs_test_nYears,

        check_env_train_year_start1,
        check_env_test_year_start1,
        check_env_train_nSites,
        check_env_test_nSites,
        check_env_train_nYears,
        check_env_test_nYears,

        check_init_train_sites,
        check_init_test_sites,
        check_init_train_censusyear,
        check_init_test_censusyear
      )

      #-----------------------------------------------------------------------
      # 5) Append this row to "results"
      #-----------------------------------------------------------------------
      results <- rbind(
        results,
        data.table(
          i_folder = i_folder,
          cv_label = cv_label,
          pass_all = pass_all,

          check_obs_train_siteID_min,
          check_obs_train_siteID_seq,
          check_obs_test_siteID_min,
          check_obs_test_siteID_seq,
          check_obs_train_year_ok,
          check_obs_test_year_ok,
          check_obs_train_nSites,
          check_obs_test_nSites,
          check_obs_train_nYears,
          check_obs_test_nYears,

          check_env_train_year_start1,
          check_env_test_year_start1,
          check_env_train_nSites,
          check_env_test_nSites,
          check_env_train_nYears,
          check_env_test_nYears,

          check_init_train_sites,
          check_init_test_sites,
          check_init_train_censusyear,
          check_init_test_censusyear
        ),
        fill=TRUE
      )
    }
  }

  return(results)
}

# -------------------------------------------------------------------------
# 2. Define the folders vector you want to check
# -------------------------------------------------------------------------
folders <- c(
  "pft-period7-25patches", "pft-period7-1patch",
  "pft-period35-1patch",   "pft-period35-25patches",
  "genus-period7-25patches", "genus-period7-1patch",
  "genus-period35-1patch", "genus-period35-25patches"
)

# -------------------------------------------------------------------------
# 3. We will validate both "realdata" and "simdata" directories
# -------------------------------------------------------------------------
# out_dirs <- c("data/BCI/CVsplits-realdata", "data/BCI/CVsplits-simdata")
out_dirs <- c("data/BCI/CVsplits-realdata")

# -------------------------------------------------------------------------
# 4. Run validation for each combination of folder and out_dir
# -------------------------------------------------------------------------
# We store results in a list for convenience (optional).
all_results <- list()

for (folder in folders) {
  for (odir in out_dirs) {
    message("\n--- Validating folder: '", folder,
            "' in directory: '", odir, "' ---")

    # Call your validation function
    val_dt <- validate_CV_files(
      i_folder    = folder,
      out_dir0    = odir,
      noSplitsDir = "data/BCI/noSplits"
    )

    # Summarize
    n_total <- nrow(val_dt)
    n_pass  <- sum(val_dt$pass_all)
    n_fail  <- n_total - n_pass

    message("  --> Out of ", n_total, " folds, ", n_pass, " passed, ", n_fail, " failed.")

    # Print details if there are any failures
    if(n_fail > 0) {
      message("  --> The following folds did NOT pass:")
      print(val_dt[pass_all == FALSE])
    }

    # Store results in our list (optional)
    key_name <- paste(folder, basename(odir), sep = "_")
    all_results[[key_name]] <- val_dt
  }
}
