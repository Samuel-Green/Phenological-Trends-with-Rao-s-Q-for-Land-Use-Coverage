############################################################################ ###
### Elliot Samuel Shayle - University of Marburg                               #
### 04.3_Kilimanjaro_PPI_a_Grid-search.R                                       #
### TWDTW a (Steepness) Optimization via Transect Grid Search                  #
############################################################################ ###

library(rasterdiv)
library(twdtw)
library(vegan)
library(terra)
library(dplyr)
library(here) # Ensures dynamic, drive-agnostic pathing
library(parallel)

# 1. Get every single folder path inside your main project directory
all_dirs <- list.dirs(here(), recursive = TRUE)

# 2. Define a regex pattern of dangerous Linux characters (Spaces, brackets, shell triggers, and non-ASCII/emojis)
danger_pattern <- "[ \\(\\)\\[\\]\\{\\}&\\|\\*\\?\\$\\<\\>\\'\"]|[^\x01-\x7F]"

# 3. Filter for any paths that match the danger pattern
bad_dirs <- all_dirs[grepl(danger_pattern, all_dirs)]

# 4. Print the results
if(length(bad_dirs) > 0) {
  message("?? DANGER: Found the following problematic paths:")
  print(bad_dirs)
} else {
  message("? ALL CLEAR: Your folder paths are perfectly Linux-safe!")
}

### 1. Path definitions ####

InputData <- file.path("~/TWDTW_Paper/Data/Input_Data")
ProcessedData <- file.path("~/TWDTW_Paper/Data/Processed_Data")
Results <- file.path("~/TWDTW_Paper/Results")

## KiliNP-specific paths

KiliNP_Input <- file.path(InputData, "Kilimanjaro")
KiliNP_Processed <- file.path(ProcessedData, "Kilimanjaro")
KiliNP_Results <- file.path(Results, "Kilimanjaro")
dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

message("Importing Spatial Objects for Transect Subsetting...")

# 2. Import Spatial Vectors
kili.tiling.grid <- vect(file.path(KiliNP_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON")) # Load it back in
transect_line <- vect(file.path(KiliNP_Processed, "KiliNP_Transect_Line.geojson"))
KiliNP_LandCover_Vector <- vect(file.path(KiliNP_Input, "Kili_Ground_Truthing_Land_Cover_Classifications", "VegAug1_KILI_SES_withnewcof.shp"))

# 3. Import Full Cleaned PPI Timeseries
KiliNP_Timeseries_Clean <- rast(file.path(KiliNP_Processed, "KiliNP_PPI_2017-2021_Timeseries_SG-filtered.tif"))

# Ensure the transect line shares the exact CRS of the tiling grid
if (crs(transect_line) != crs(kili.tiling.grid)) {
  transect_line <- project(transect_line, crs(kili.tiling.grid))
}

# 4. Extract the Time Vector dynamically from layer names
Kili.dates <- as.Date(
  paste0(
    sub(" - .*", "", names(KiliNP_Timeseries_Clean)), "-",
    match(sub(".* - ", "", names(KiliNP_Timeseries_Clean)), month.name),
    "-01"
  ),
  format = "%Y-%m-%d"
)

### Transect Subsetting ####
message("Intersecting transect line with tiling grid...")

# Keep only the tiles that physically intersect the transect line
transect_tiles <- kili.tiling.grid[transect_line, ]

# Dissolve the selected tiles into a single continuous polygon for easy cropping
transect_poly <- aggregate(transect_tiles)

message("Cropping and masking PPI timeseries to transect footprint...")
transect_ts <- crop(KiliNP_Timeseries_Clean, transect_poly)
transect_ts <- mask(transect_ts, transect_poly)
transect_ts <- trim(transect_ts)

# ---> NEW: The Array Severance Trick (NO BINNING) <---
message("Severing lazy C++ pointers and mapping strict numeric memory...")

# 1. Extract the lazy raster into a raw base R array
transect_array <- terra::as.array(transect_ts)

# 2. Multiply by 1.0 to forcefully block ALTREP from converting NA windows to logicals
transect_array <- transect_array * 1.0 

# 3. Rebuild the raster directly into physical RAM
transect_ts_RAM <- terra::rast(
  transect_array, 
  crs = crs(transect_ts), 
  ext = ext(transect_ts)
)

# 4. BAKE the strictly numeric raster to the hard drive to protect parallel workers
message("Caching full-resolution numeric raster to disk for safe parallel transmission...")
tmp.tif_path <- file.path(KiliNP_Processed, "Kili_PPI_Grid_Search_Transect_FullRes.tif")
writeRaster(transect_ts_RAM, tmp.tif_path, datatype = "FLT8S", overwrite = TRUE)

# 5. Free up the fragile temporary RAM objects
rm(transect_array, transect_ts_RAM)
gc()

# 6. Load the file-backed raster (This safely passes the file path to the workers!)
transect_ts <- terra::rast(tmp.tif_path)
message("Done ??")

# Create a spatial template for the Land Cover rasterisation
transect_mean <- app(transect_ts, fun = mean, na.rm = TRUE)

### Pre-Loop Land Cover Rasterisation ####
message("Preparing Land Cover Vector and Raster Template...")

# Align CRS
if (crs(transect_mean) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(transect_mean))
}
message("CRSs aligned")

# Rasterise using the transect mean as the precise spatial template
KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  transect_mean,
  field = "grid_code"
)
message("Land cover rasterised successfully ??")

### Grid Search for Optimal TWDTW a ####
message("Starting two-stage grid search for optimal TWDTW a (steepness)...")

# Define the log file path
log_csv <- file.path(KiliNP_Results, "Kili_PPI_Alpha_GridSearch_NoBin_Log.csv")
message("CSV log file initialised ??")

alpha_coarse_grid <- c(-0.1, -0.3, -0.5, -0.7, -0.9)
alpha_results <- data.frame(Alpha = numeric(), R2 = numeric(), p_value = numeric())

message("Commensing coarse grid search ????")

## STAGE 1: Coarse Grid
for (a in alpha_coarse_grid) {
  message(paste("Coarse testing \U03B1 =", a))
  
  tmp.RaoQ <- paRao(
    x = transect_ts,
    time_vector = Kili.dates,
    window = 3,
    alpha = 2,
    na.tolerance = 0,
    simplify = 2,
    np = max(1, detectCores()), 
    progBar = FALSE,
    method = "multidimension",
    dist_m = "twdtw",
    midpoint = 6, 
    stepness = a, 
    cycle_length = "year",
    time_scale = "month"
  )
  
  message(paste0(Sys.time()," The paRao function ran successfully")) 

  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack <- c(tmp.raster, KiliNP_LandCover_Raster) 
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  
  # Protect against PERMANOVA memory crashes on large transects
  if(nrow(tmp.df) > 10000) tmp.df <- tmp.df[sample(nrow(tmp.df), 10000), ]
  
  message(paste0(Sys.time()," Commensing PERMANOVA"))
  
  tmp.permanova <- adonis2(tmp.df$RaosQ ~ tmp.df$Veg_GroundTruth, permutations = 999, parallel = max(1, detectCores()))
  
  # Log result to dataframe
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  
  # Overwrite the entire CSV with the updated dataframe as a live checkpoint
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  
  message("PERMANOVA result successfully logged ????")
  
  # Force R to dump the C++ matrices from RAM before starting the next alpha value
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df)
  gc()
  
  message("Commensing next loop iteration ?")
}

best.coarse.alpha <- alpha_results$Alpha[which.max(alpha_results$R2)]
message(paste("Best coarse a is:", best.coarse.alpha))

## STAGE 2: Fine Grid
alpha_fine_grid <- seq(best.coarse.alpha + 0.1, best.coarse.alpha - 0.1, by = -0.05)
alpha_fine_grid <- alpha_fine_grid[alpha_fine_grid < 0] 

message("Commensing fine grid search ??")

for (a in alpha_fine_grid) {
  if (a %in% alpha_coarse_grid) next 
  
  message(paste("Fine testing \U03B1 =", a))
  
  tmp.RaoQ <- paRao(
    x = transect_ts, time_vector = Kili.dates, window = 3, alpha = 2,
    na.tolerance = 0, simplify = 2, np = max(1, detectCores()), progBar = FALSE, 
    method = "multidimension", dist_m = "twdtw", midpoint = 6, stepness = a, 
    cycle_length = "year", time_scale = "month"
  )
  
  message(paste0(Sys.time()," The paRao function ran successfully")) 

  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack <- c(tmp.raster, KiliNP_LandCover_Raster) 
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  if(nrow(tmp.df) > 10000) tmp.df <- tmp.df[sample(nrow(tmp.df), 10000), ]
  
  message(paste0(Sys.time()," Commensing PERMANOVA")) 
  
  tmp.permanova <- adonis2(tmp.df$RaosQ ~ tmp.df$Veg_GroundTruth, permutations = 999, parallel = max(1, detectCores()))
  
  # Log result to dataframe
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  
  # Overwrite the entire CSV with the updated dataframe as a live checkpoint
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  
  message("PERMANOVA result successfully logged ????")
  
  # Force R to dump the C++ matrices from RAM before starting the next alpha value
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df)
  gc()
}

message("Grid search completed successfully ??")

# Sort final results to show best at the top
alpha_results <- alpha_results[order(-alpha_results$R2), ]
print("Full Grid Search Complete. Results:")
print(alpha_results)

Kili_PPI_Optimal_Alpha <- alpha_results$Alpha[1]
message(paste("The absolute optimal a value is:", Kili_PPI_Optimal_Alpha))

### Exports ####

# Save the numerical optimal alpha as an RDS object
saveRDS(Kili_PPI_Optimal_Alpha, file.path(KiliNP_Results, "Kili_PPI_Optimal_Alpha_NoBin.rds"))
message(paste("RDS results file saved to", file.path(KiliNP_Results, "Kili_PPI_Optimal_Alpha_NoBin.rds")))

# Export the final, sorted dataframe to guarantee the log is complete and ordered
write.csv(alpha_results, file = log_csv, row.names = FALSE)
message("Final CSV log file exported successfully")

message("All grid search results were exported successfully ??.")