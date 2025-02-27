
create_grid <- function(
    all_trees_cleaned,
    x_length = 40, 
    y_length = 25,
    x_range = c(0, 500), 
    y_range = c(0, 1000),
    patches2sites = c(2, 2)
) {
  # 1) Filter out any rows with missing gx/gy
  all_trees <- all_trees_cleaned[!is.na(gx) & !is.na(gy)]
  
  # 2) Create the "fine" grid classes for x and y
  all_trees[, c("x_class", "y_class") := .(
    cut(
      gx,
      breaks = seq(x_range[1], x_range[2], x_length),
      labels = seq(x_range[1] + x_length, x_range[2], x_length) - x_length / 2,
      include.lowest = TRUE
    ),
    cut(
      gy,
      breaks = seq(y_range[1], y_range[2], y_length),
      labels = seq(y_range[1] + y_length, y_range[2], y_length) - y_length / 2,
      include.lowest = TRUE
    )
  )]
  
  # Convert the factors to numeric for convenience
  all_trees[, x_class := as.numeric(as.character(x_class))]
  all_trees[, y_class := as.numeric(as.character(y_class))]
  # remove entries with x and y NA
  all_trees <- all_trees[!is.na(x_class) & !is.na(y_class)]
  
  # 3) Create a unique patch ID for every fine-grid cell
  #    ("1-patch" scenario, each cell is its own patch)
  all_trees[, uniquePatch := .GRP, by = .(x_class, y_class)]
  
  # 4) Make a small table of these unique patches
  grid_dt_1patch <- unique(
    all_trees[, .(
      siteID      = uniquePatch,
      uniquePatch = uniquePatch,
      patch       = 1L,
      x_class,
      y_class
    )]
  )
  
  # 5) Create a coarser patching, based on patches2sites = c(row_bins, col_bins).
  #    This lumps the fine grid cells into fewer patches.
  #    For example, patches2sites = c(2,2) => 4 patches total
  #                 patches2sites = c(4,4) => 16 patches total
  #    The logic here is:
  #    - Convert x_class, y_class to “index” form (1..N)
  #    - Cut that index range into the desired number of rowwise/colwise bins
  #    - Combine those bins into a “coarse” patch ID
  
  # First, get the integer indices for each unique x_class, y_class
  all_trees[, x_class_idx := as.integer(as.factor(x_class))]
  all_trees[, y_class_idx := as.integer(as.factor(y_class))]
  
  # Figure out how many distinct x- and y- bins we have in the fine grid
  Nx <- max(all_trees$x_class_idx)
  Ny <- max(all_trees$y_class_idx)
  
  # Define rowwise/columnwise cut breaks for coarser grouping
  x_breaks <- seq(0.5, Nx+0.5, patches2sites[1])
  y_breaks <- seq(0.5, Ny+0.5, patches2sites[2])
  
  # Cut the indices into coarser bins
  all_trees[, x_class_group := cut(
    x_class_idx, breaks = x_breaks, 
    labels = FALSE, include.lowest = TRUE
  )]
  all_trees[, y_class_group := cut(
    y_class_idx, breaks = y_breaks, 
    labels = FALSE, include.lowest = TRUE
  )]
  
  # Combine rowwise and columnwise bins into a single coarse patch ID
  # We can do it via .GRP or by a formula:
  #   patch_coarse = (x_class_group - 1) * patches2sites[2] + y_class_group
  # For uniqueness, let's just use .GRP for clarity:
  all_trees[, patch_coarse := .GRP, by = .(x_class_group, y_class_group)]
  
  # Create a data.table summarizing those coarse patches
  grid_dt_coarse <- unique(
    all_trees[,
              .(siteID      = patch_coarse,
                uniquePatch = uniquePatch,  # fine patch
                x_class,
                y_class
              )
    ]
  )
  plot(all_trees$gx, all_trees$gy, pch = 1, col = "grey90", asp = 1)
  points(grid_dt_coarse$x_class, grid_dt_coarse$y_class, pch = 19, col = factor(grid_dt_coarse$siteID))
  # If you also want an index within each coarse site:
  grid_dt_coarse[, patch := seq_len(.N), by = siteID]
  # If you also want an index within each coarse site:
  
  # 6) Return a list that includes:
  #    - The original data with all columns added
  #    - The fine-grid "1-patch" table
  #    - The coarser patch table
  return(list(
    all_trees_grid  = all_trees,
    grid_dt_1patch  = grid_dt_1patch,
    grid_dt_coarse  = grid_dt_coarse
  ))
}


clean_tree_level = function(all_trees, out_dir){
  if(!dir.exists(out_dir)) {
    dir.create(out_dir,recursive = T)
  }
  #### status ####
  all_trees[, status_orig := status,]
  all_trees[status %in% c("A", "AD", "AR"), status := "A",]
  all_trees[is.na(dbh) & status != "P", status := "D",]
  all_trees[, gdbh_cm := dbh/10] # for growth
  all_trees[, dbh_cm := dbh/10]
  all_trees[,NstemIDs := uniqueN(stemID), by = treeID]


  # define variables of entries before and after
  all_trees[, ":="(
    status_before = data.table::shift(status, 1, type = "lag"),
    status_after = data.table::shift(status, 1, type = "lead"),
    # gdbh_cm_before = data.table::shift(gdbh_cm, 1, type = "lag"),
    # gdbh_cm_after = data.table::shift(gdbh_cm, 1, type = "lead"),
    period_length = census - data.table::shift(census, 1, type = "lag"),
    period_length_after = data.table::shift(census, 1, type = "lead") - census,
    hom_before = data.table::shift(hom, 1, type = "lag"),
    hom_after = data.table::shift(hom, 1, type = "lead")
  ),
  by = .(treeID)]

  #### growth plausibility ####
  cat("\niterating through time series to check growth plausibility")
  growth_checked = F
  gdbh_cm_interpol_sum = 0
  all_trees[,dbh_interpolated := F,]
  while(!growth_checked){
    check_all_trees <- copy(all_trees)
    # update gdbh_cm_before and gdbh_cm_after
    all_trees[, ":="(
      gdbh_cm_before = data.table::shift(gdbh_cm, 1, type = "lag"),
      gdbh_cm_after = data.table::shift(gdbh_cm, 1, type = "lead")
    ),
    by = .(treeID)]

    # individual tree growth calculation
    all_trees[, ":="(
      absolute_growth_cm = (gdbh_cm-gdbh_cm_before)/period_length,
      absolute_growth_cm_after = (gdbh_cm_after-gdbh_cm)/period_length_after,
      relative_growth = (gdbh_cm/gdbh_cm_before)^(1/period_length)-1,
      relative_growth_5yr = (gdbh_cm/gdbh_cm_before)^(1/1)-1
    ),]
    # cat("NA in relative growth: ",sum(is.na(all_trees$relative_growth))/nrow(all_trees), "\n")

    all_trees[census != 1 & hom != hom_before | status != "A" | NstemIDs > 1,
              ":="(
                absolute_growth_cm = NA,
                relative_growth = NA,
                relative_growth_5yr = NA
              ), ]

    # define ranges for growth that are considered unreasonable
    quantile(all_trees$absolute_growth_cm, c(seq(0, 0.005,0.001), seq(0.99,1,0.001)), na.rm = T)
    absolute_growth_range = c(-0.2, 4.5)
    quantile(all_trees$relative_growth, c(seq(0, 0.005,0.001), seq(0.999,1,0.0001)), na.rm = T)
    relative_growth_range = c(-0.01, 5)
    all_trees[absolute_growth_cm > absolute_growth_range[2] | absolute_growth_cm < absolute_growth_range[1],
              ":="(absolute_growth_cm = NA, relative_growth = NA), ]
    all_trees[relative_growth > relative_growth_range[2] | relative_growth < relative_growth_range[1],
              ":="(absolute_growth_cm = NA, relative_growth = NA),]
    all_trees[absolute_growth_cm_after < absolute_growth_range[1],
              ":="(gdbh_cm = NA), ]

    # interpolate growth
    all_trees[census != 1 & is.na(gdbh_cm) & status == "A" & !is.na(gdbh_cm_before) & !is.na(gdbh_cm_after) & NstemIDs == 1,
              ":="(
                gdbh_cm = mean(c(gdbh_cm_before, gdbh_cm_after), na.rm = F),
                dbh_interpolated = T
              )
              , by = treeID]
    # all_trees[
    #   !is.na(gdbh_cm) & status == "A" & !is.na(gdbh_cm_before) & !is.na(gdbh_cm_after),
    #   gdbh_cm2 := (gdbh_cm_before+gdbh_cm_after/2), by = treeID]
    if(identical(check_all_trees, all_trees)) growth_checked = T
    gdbh_cm_interpol_sum = nrow(all_trees[dbh_interpolated == T])
    Ngdbh_cm = length(all_trees[status == "A"]$gdbh_cm)
    # # cat("interpolated ", gdbh_cm_interpol, "of", Ngdbh_cm, round(gdbh_cm_interpol/Ngdbh_cm,4)  ," %" ,"in this iterations\n")
    # cat("interpolated a total of", gdbh_cm_interpol_sum, "of", Ngdbh_cm, round(gdbh_cm_interpol_sum/Ngdbh_cm,4)  ," %" ,"\n")
    # cat("Growth checked: ", growth_checked, "\n")
  }
  cat("\nGrowth plausibility checked\n")


  all_trees[status=="A", status2 := "alive",]
  all_trees[census != 1985 & status=="A" & status_before=="P", status2 := "regeneration",]
  all_trees[status=="D" & status_before=="A", status2 := "died",]

  return(all_trees)
}

calculate_stand_vars = function(all_trees, out_dir, area_of_square_ha){
  if(!dir.exists(out_dir)) {
    dir.create(out_dir,recursive = T)
  }

  ## stand variables ####
  stand_dt <- all_trees[,.(
    ba = sum(dbh_cmTOba_m2(dbh_cm)*nostems*(status=="A"), na.rm = T),
    trees = sum((status=="A" & !is.na(dbh_cm))*nostems, na.rm = T),
    dbh_mean = sum(dbh_cm*(status=="A" & !is.na(dbh_cm))*nostems, na.rm = T)/sum((status=="A" & !is.na(dbh_cm))*nostems, na.rm = T),
    g = mean(relative_growth, na.rm = T),
    g_5yr = mean(relative_growth_5yr, na.rm = T),
    fresh_dead = sum(status=="D" & status_before=="A", na.rm = T),
    r = sum(year != 2000 & status=="A" & is.na(status_before), na.rm = T)/area_of_square_ha
  ), by=.(census, uniquePatch = uniquePatch, species, period_length)]

  stand_dt[census == 1985, ":="(
    r = NA,
    m = NA,
    m_5yr = NA,
    g = NA,
    g_5yr = NA
  ),]

  # define trees_before in stand_dt
  stand_dt[, trees_before := data.table::shift(trees, 1, type = "lag"), by=.(uniquePatch, species)]
  stand_dt[, m_old := fresh_dead/trees_before, ]
  stand_dt[, m := 1-((trees_before-fresh_dead)/trees_before)^(1/period_length), ]
  stand_dt[is.infinite(m), m := 1]
  stand_dt[, m_5yr := 1-((trees_before-fresh_dead)/trees_before)^(1/1), ]
  stand_dt[is.infinite(m_5yr), m_5yr := 1]

  # calculate .N for growth, mort, and reg for each species
  Nobs_dt2 <-
    stand_dt[,.(
      growth_N = sum(!is.na(g)),
      mort_N = sum(!is.na(m)),
      reg_N = sum(!is.na(r)),
      growth_mort_N = sum(!is.na(g) & !is.na(m))
    ), by = species]
  Nobs_dt2[order(growth_N)]

  ### plotting stand data ####
  pdat_census_site_all <-  all_trees[,.(
    ba = sum(dbh_cmTOba_m2(dbh_cm)*(status=="A"), na.rm = T),
    trees = sum((status=="A" & !is.na(dbh_cm))*nostems, na.rm = T),
    dbh_mean = mean(gdbh_cm, na.rm = T),
    g = mean(relative_growth, na.rm = T),
    g_5yr = mean(relative_growth_5yr, na.rm = T),
    fresh_dead = sum(status=="D" & status_before=="A", na.rm = T),
    r = sum(census != 1985 & status=="A" & status_before=="P", na.rm = T)/area_of_square_ha
  ), by=.(census, uniquePatch = uniquePatch, period_length)]
  pdat_census_site_all[, trees_before := data.table::shift(trees, 1, type = "lag"), by=.(uniquePatch)]
  pdat_census_site_all[, m := 1-((trees_before-fresh_dead)/trees_before)^(1/period_length), ]
  pdat_census_site_all[is.infinite(m), m := 1]
  pdat_census_site_all[, m_5yr := 1-((trees_before-fresh_dead)/trees_before)^(1/1), ]
  pdat_census_site_all[is.infinite(m_5yr), m_5yr := 1]
  pdat_census_site_all[census == 1985, ":="(r = NA, m = NA, m_5yr = NA,g = NA, g_5yr = NA),]
  pdat_census_site_all[, species := "all species",]

  Nsites = uniqueN(stand_dt$uniquePatch)
  pdat_census_species <- stand_dt[, .(
    ba = sum(ba, na.rm = T)/Nsites,
    trees = sum(trees, na.rm = T)/Nsites,
    dbh_mean = mean(dbh_mean, na.rm = T),
    g = mean(g, na.rm = T),
    g_5yr = mean(g_5yr, na.rm = T),
    r = sum(r, na.rm = T)/Nsites,
    fresh_dead = sum(fresh_dead, na.rm = T)/Nsites,
    trees_before = sum(trees_before, na.rm = T)/Nsites,
    m_old = sum(fresh_dead, na.rm = T)/sum(trees_before, na.rm = T),
    m = 1-((trees_before-fresh_dead)/trees_before)^(1/period_length),
    m_5yr = 1-((trees_before-fresh_dead)/trees_before)^(1/1)
  ), by = .(census, species, period_length)]
  pdat_census_species[census == 1985, ":="(r = NA, m = NA, m_5yr = NA, g = NA, g_5yr = NA),]
  pdat_census_species[, uniquePatch := "all sites",]

  pdat_census_all <- pdat_census_species[, .(
    ba = sum(ba, na.rm = T),
    trees = sum(trees, na.rm = T),
    dbh_mean = mean(dbh_mean, na.rm = T),
    g = mean(g, na.rm = T),
    g_5yr = mean(g_5yr, na.rm = T),
    r = sum(r, na.rm = T),
    # m = sum(fresh_dead, na.rm = T)/sum(trees_before, na.rm = T),
    m = 1-((trees_before-fresh_dead)/trees_before)^(1/period_length),
    m_5yr = 1-((trees_before-fresh_dead)/trees_before)^(1/1)
  ), by = .(census)]
  pdat_census_all[census == 1985, ":="(r = NA, m = NA, m_5yr = NA, g = NA, g_5yr = NA, species = "all species"),]
  pdat_census_all[, species := "all species",]
  pdat_census_all[, uniquePatch := "all sites",]

  library(ggplot2)
  pdat_list = list(
    "patch-species" = stand_dt,
    "patch-all" = pdat_census_site_all,
    "full-species" = pdat_census_species,
    "full-all" = pdat_census_all)
  pdf(paste0(out_dir, "/pdat.pdf"), width = 10, height = 5)
  i=1
  for(i in 1:length(pdat_list)){
    pdat_i = pdat_list[[i]]
    name_i = names(pdat_list)[i]
    alpha = ifelse(grepl("all", name_i), 0.5, 0.2)
    p_g <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = g, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = g, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("growth rate")+
      xlab("census")+
      ggtitle(paste0(name_i, " growth rate"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_g)

    p_r <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = r, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = r, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("regeneration rate")+
      xlab("census")+
      ggtitle(paste0(name_i, " regeneration rate"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_r)

    p_m <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = m, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = m, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("mortality rate")+
      xlab("census")+
      ggtitle(paste0(name_i, " mortality rate"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_m)

    p_ba <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = ba, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = ba, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("basal area")+
      xlab("census")+
      ggtitle(paste0(name_i, " basal area"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_ba)

    p_trees <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = trees, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = trees, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("trees")+
      xlab("census")+
      ggtitle(paste0(name_i, " trees"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_trees)

    p_dbh <- ggplot() +
      geom_line(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = dbh_mean, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha) +
      geom_jitter(
        data = pdat_i,
        mapping = aes(
          x = factor(census), y = dbh_mean, color = factor(species),
          group = interaction(factor(species), uniquePatch)), alpha = alpha, width = 0.1) +
      theme(legend.position = "none")+
      ylab("dbh mean")+
      xlab("census")+
      ggtitle(paste0(name_i, " dbh mean"))+
      coord_cartesian(ylim = c(0, NA))
    print(p_dbh)

    # save pdat
    fwrite(pdat_i, paste0(out_dir, "/", name_i, ".csv"))

    # save summary statistic of each rate for each species
    padt_summary <- pdat_i[, .(
      g_mean = mean(g, na.rm = T),
      g_sd = sd(g, na.rm = T),
      min_g = min(g, na.rm = T),
      max_g = max(g, na.rm = T),
      r_mean = mean(r, na.rm = T),
      r_sd = sd(r, na.rm = T),
      min_r = min(r, na.rm = T),
      max_r = max(r, na.rm = T),
      m_mean = mean(m, na.rm = T),
      m_sd = sd(m, na.rm = T),
      min_m = min(m, na.rm = T),
      max_m = max(m, na.rm = T)
    ), by = .(species)]
    fwrite(padt_summary, paste0(out_dir, "/", name_i, "_summary.csv"))
    # save quantiles of each rate over all species with rate type as rows and quantiles as columns
    sum_quantiles = seq(0,1,0.05)
    padt_summary_quantiles <- t(pdat_i[, .(
      g = round(quantile(g, sum_quantiles, na.rm = T),3),
      r = round(quantile(r, sum_quantiles, na.rm = T),3),
      m = round(quantile(m, sum_quantiles, na.rm = T),3)
    ),])
    colnames(padt_summary_quantiles) <- paste0(sum_quantiles*100,"%")
    fwrite(padt_summary_quantiles, paste0(out_dir, "/", name_i, "_summary_quantiles.csv"))
    cat("Summary statistics and quantiles of", name_i, "saved\n")
  }
  dev.off()


  # make histogram with 100 breaks with base r for each rate and save it to file
  pdf(paste0(out_dir, "/histograms.pdf"), width = 5, height = 5)
  for(i in 1:length(pdat_list)){
    pdat_i = pdat_list[[i]]
    name_i = names(pdat_list)[i]
    hist(pdat_i$g, breaks = 100, main = paste0(name_i ," Growth rate"), xlab = "Growth rate")
    hist(pdat_i$r, breaks = 100, main = paste0(name_i ," Regeneration rate"), xlab = "Regeneration rate")
    hist(pdat_i$m, breaks = 100, main = paste0(name_i ," Mortality rate"), xlab = "Mortality rate")
    hist(pdat_i$ba, breaks = 100, main = paste0(name_i ," Basal area"), xlab = "Basal area")
    hist(pdat_i$trees, breaks = 100, main = paste0(name_i ," Trees"), xlab = "Trees")
  }
  dev.off()

  ## stand_dt.csv ####
  fwrite(stand_dt, paste0(out_dir,"/stand_dt.csv"))

  return(list(stand_dt = stand_dt, all_trees = all_trees))
}


create_init_cohorts = function(all_trees, stand_dt, out_dir, dbh_cm_class_size = 0.1, init_years){
  if(!dir.exists(out_dir)) {
    dir.create(out_dir,recursive = T)
  }

  for(i_year in init_years){
    # split cohorts by dbh size class for efficiency
    dbh_cm_class_size = dbh_cm_class_size
    initial_trees <- all_trees[census == i_year & status=="A" & !is.na(dbh_cm), .(
      trees = sum((status == "A")*nostems)
    ), by = .(dbh = cut(
      dbh_cm,
      breaks = seq(6, 351, dbh_cm_class_size),
      labels = seq(6, 351-dbh_cm_class_size, dbh_cm_class_size),
      include.lowest = T), species, siteID, patchID = patch, census)
    ]
    initial_trees <- initial_trees[order(siteID, patchID)]
    initial_trees[,cohortID := 1:.N,by = .(siteID, patchID, census)]
    initial_trees[, dbh := as.numeric(as.character(dbh)),]
    
    # check if basic stand properties in initial_trees are similar to all_trees and stand_dt
    comp_stand_dt_from_initial_trees <- initial_trees[, .(
      ba = sum(dbh_cmTOba_m2(dbh)*trees, na.rm = T),
      trees = sum(trees, na.rm = T),
      dbh_mean = sum(dbh*trees, na.rm = T)/sum(trees, na.rm = T)
    ), by = .(census, siteID, patchID, species)]

    stand_dt[, patchID := patch,]
    compdt1 <- merge(
      comp_stand_dt_from_initial_trees, stand_dt[census == i_year & (m == 0 | is.na(m)) & !is.na(dbh_mean)],
      by = c("census", "siteID", "patchID", "species"), suffixes = c(".init", ".stand_dt")
    )

    compdt1[, ba.diff := abs(ba.init - ba.stand_dt),]
    compdt1[, dbh.diff := abs(dbh_mean.init - dbh_mean.stand_dt),]
    compdt1[, trees.diff := abs(trees.init - trees.stand_dt),]

    if(any(compdt1$ba.diff > 0.1)) {
      plot(compdt1$ba.init, compdt1$ba.stand_dt)
      stop("ba difference is larger than 0.01")
    }

    # hist(all_trees[census == 2000 & siteID == 3 & patch == 1 & species == 3]$dbh)
    # mean(all_trees[census == 2000 & siteID == 3 & patch == 1 & species == 3]$dbh, na.rm = T)
    # mean(all_trees[census == i_year & status=="A" & !is.na(dbh_cm) & siteID == 3 & patch == 1 & species == 3, .(tree_nr, status, nostems, dbh)]$dbh)
    # initial_trees[census == 2000 & siteID == 3 & patchID == 1 & species == 3]
    # comp_stand_dt_from_initial_trees[census == 2000 & siteID == 3 & patchID == 1 & species == 3]
    
    compdt1[dbh.diff > 3]
    if(any(compdt1$dbh.diff > dbh_cm_class_size+0.001)) {
      par(pty = "s")
      plot(compdt1$dbh_mean.init, compdt1$dbh_mean.stand_dt)
      abline(0,1)
      stop(paste("dbh difference is larger than", dbh_cm_class_size,"+0.001"))
    }

    if(any(compdt1$trees.diff > 1)) {
      plot(compdt1$trees.init, compdt1$trees.stand_dt)
      stop("trees difference is larger than 1")
    }

    ## initial_cohorts.csv ####
    fwrite(initial_trees, paste0(out_dir,"/initial_cohorts",i_year,".csv"))
  }
}
