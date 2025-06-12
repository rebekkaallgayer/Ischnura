library(blockCV)

# =============================================================

# 2. Test for autocorrealation (within each fold)

# =============================================================

# 2.1 Open a connection to write logs

output_file <- file(paste0(outDir, "m_rf_spatblock.txt"), "w")



for (i in 1:10){
  
  # 2.2 Select all data for the current fold (both presences and pseudo-absences)
  
  ones <- occ[occ$pa == 1, ]         # Filter rows where 'pa' is 1
  
  print(nrow(ones)) # 17575
  
  print(head(ones))
  
  print(summary(ones))
  
  
  
  fold_zeros <- occ[!is.na(occ$subset) & occ$subset == i, ] # Filter rows where 'subset' matches 'i'
  
  print(nrow(fold_zeros)) # 35150
  
  print(head(fold_zeros))
  
  print(summary(fold_zeros))
  
  
  
  fold_data <- rbind(ones, fold_zeros)
  
  print(nrow(fold_data)) # 52725
  
  print(head(fold_data))
  
  print(summary(fold_data))
  
  
  
  # remove subset col
  
  fold_data = subset(fold_data, select = -subset)
  
  
  
  print("head of fold_data")
  
  print(head(fold_data))
  
  print(summary(fold_data))
  
  
  
  # 2.3 Test for spatial autocorrelation
  
  # 2.3.1 Converting to km
  
  dat.sf_env <- st_as_sf(fold_data, coords = c("X", "Y"), crs = 4326)
  
  dat.sf_lambert_env <- st_transform(dat.sf_env, 3035) # EPSG:3035 is the ETRS89-extended / LAEA Europe coordinate reference system = PANEUROPEAN
  
  
  
  # 2.3.2 Spatial autocorrelation of presence/pseudo-absence
  
  sac1 <- cv_spatial_autocor(x = dat.sf_lambert_env,
                             
                             column = "pa", # binary
                             
                             plot = TRUE)
  
  
  
  # Extract suggeste distance
  
  range_dist <- sac1$range
  
  
  
  # =============================================================
  
  # 3. Creating spatial blocks
  
  # =============================================================
  
  # 3.1 Spatial blocking
  
  sb2 <- cv_spatial(x = dat.sf_lambert_env,
                    
                    column = "pa",
                    
                    size = range_dist,
                    
                    k = 5, # 5 blocks
                    
                    hexagon = FALSE,
                    
                    selection = "systematic")
  
  str(sb2)
  
  
  
  # =============================================================
  
  # 4. Cross-validation with spatial blcoks
  
  # =============================================================
  
  dat.sf_lambert_env$block <- sb2$folds_ids  # Assign fold IDs to each record
  
  # class(dat.sf_lambert)
  
  
  
  coords_env <- st_coordinates(dat.sf_lambert_env)  # Extract X and Y coordinates
  
  dat_df_env <- cbind(st_drop_geometry(dat.sf_lambert_env), coords_env)
  
  
  
  # Take unique blocks to loop over
  
  unique_blocks <- unique(dat_df_env$block)
  
  
  
  for (block in unique_blocks){
    
    # Split into training (9 block) and testing (1 block)
    
    train_data <- dat_df_env[dat_df_env$block != block, ]
    
    print("No of train data")
    
    print(nrow(train_data))
    
    print(head(train_data))
    
    test_data <- dat_df_env[dat_df_env$block == block, ]
    
    print("No of test data")
    
    print(nrow(test_data))
    
    print(head(test_data))
    
    
    
    # Write to the output file
    
    cat(paste("Iteration/block", block, "(block) in", i, "(fold)"), "\n", file = output_file, append = TRUE)
    
    
    
    d <- sdmData(formula_d, train = train_data, test = test_data)
    
    
    
    m <- sdm(formula_m,
             
             data = d,
             
             methods = "rf", # RF, BRT, GAM
             
             parallelSetting = list(ncore = 5, method = 'parallel')) # No internal resampling
    
    
    
    print(m)
    
    
    
    print(paste("Fold", i))
    
    cat("Model summary:\n", file = output_file, append = TRUE)
    
    capture.output(print(m), file = output_file, append = TRUE)
    
    cat("\n\n", file = output_file, append = TRUE)
    
    
    
    remove(train_data)
    
    remove(test_data)
    
    remove(m)
    
    remove(d)
    
    
    
  }
  
  
  
}