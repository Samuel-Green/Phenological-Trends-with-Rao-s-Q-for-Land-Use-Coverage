############################################################################ ###
### Elliot Samuel Shayle - University of Marburg                               #
### 04.3_Kilimanjaro_PPI_α_Grid-search.R                                       #
### TWDTW α (Steepness) Optimization via Transect Grid Search                  #
############################################################################ ###

library(rasterdiv)
library(twdtw)
library(vegan)
library(terra)
library(dplyr)
library(here) # Ensures dynamic, drive-agnostic pathing
library(parallel)

### 1. Path definitions ####

InputData <- here::here("Data/Input Data")
ProcessedData <- here::here("Data/Processed Data")
Results <- here::here("Results 📈📉")

## KiliNP-specific paths

KiliNP_Input <- file.path(InputData, "Kilimanjaro")
KiliNP_Processed <- file.path(ProcessedData, "Kilimanjaro")
KiliNP_Results <- file.path(Results, "Kilimanjaro")
dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

message("Importing Spatial Objects for Transect Subsetting...")

# 2. Import Spatial Vectors
kili.tiling.grid <- vect(file.path(KiliNP_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON"))
transect_line <- vect(file.path(KiliNP_Processed, "KiliNP_Transect_Line.geojson"))
KiliNP_LandCover_Vector <- vect(file.path(KiliNP_Input, "Kili Ground Truthing Land Cover Classifications", "VegAug1_KILI_SES_withnewcof.shp"))

# 3. Import Full Cleaned PPI Timeseries
KiliNP_PPI_Timeseries_SG <- rast(file.path(KiliNP_Processed, "KiliNP_PPI_2017-2021_Timeseries_SG-filtered.tif"))

# Ensure the transect line shares the exact CRS of the tiling grid
if (crs(transect_line) != crs(kili.tiling.grid)) {
  transect_line <- project(transect_line, crs(kili.tiling.grid))
}

# 4. Extract the Time Vector dynamically from layer names
Kili.dates <- as.Date(
  paste0(
    sub(" - .*", "", names(KiliNP_PPI_Timeseries_SG)), "-",
    match(sub(".* - ", "", names(KiliNP_PPI_Timeseries_SG)), month.name),
    "-15"
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
transect_ts <- crop(KiliNP_PPI_Timeseries_SG, transect_poly)
transect_ts <- mask(transect_ts, transect_poly)
transect_ts <- trim(transect_ts)

# Create a spatial template for the Land Cover rasterisation
transect_mean <- app(transect_ts, fun = mean, na.rm = TRUE)

### Pre-Loop Land Cover Rasterisation ####
message("Preparing Land Cover Vector and Raster Template...")

# Align CRS
if (crs(transect_mean) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(transect_mean))
}

# Rasterise using the transect mean as the precise spatial template
# (Retains raw original 1-19 numerical grid codes)
KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  transect_mean,
  field = "grid_code"
)

### Grid Search for Optimal TWDTW α ####
message("Starting two-stage grid search for optimal TWDTW \U03B1 (steepness)...")

# Define the log file path
log_csv <- file.path(KiliNP_Results, "Kili_PPI_Alpha_GridSearch_Log.csv")

alpha_coarse_grid <- c(-0.1, -0.3, -0.5, -0.7, -0.9)
alpha_results <- data.frame(Alpha = numeric(), R2 = numeric(), p_value = numeric())

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
    np = max(1, detectCores() / 2), 
    progBar = FALSE,
    method = "multidimension",
    dist_m = "twdtw",
    midpoint = 6, 
    stepness = a, 
    cycle_length = "year",
    time_scale = "month"
  )
  
  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack <- c(tmp.raster, KiliNP_LandCover_Raster) 
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  
  # Protect against PERMANOVA memory crashes on large transects
  if(nrow(tmp.df) > 10000) tmp.df <- tmp.df[sample(nrow(tmp.df), 10000), ]
  
  tmp.permanova <- adonis2(tmp.df$RaosQ ~ tmp.df$Veg_GroundTruth, permutations = 999)
  
  # Log result to dataframe
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  
  # Overwrite the entire CSV with the updated dataframe as a live checkpoint
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  
  # Force R to dump the C++ matrices from RAM before starting the next alpha value
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df)
  gc()
}

best.coarse.alpha <- alpha_results$Alpha[which.max(alpha_results$R2)]
message(paste("Best coarse \U03B1 is:", best.coarse.alpha))

## STAGE 2: Fine Grid
alpha_fine_grid <- seq(best.coarse.alpha + 0.1, best.coarse.alpha - 0.1, by = -0.05)
alpha_fine_grid <- alpha_fine_grid[alpha_fine_grid < 0] 

for (a in alpha_fine_grid) {
  if (a %in% alpha_coarse_grid) next 
  
  message(paste("Fine testing \U03B1 =", a))
  
  tmp.RaoQ <- paRao(
    x = transect_ts, time_vector = Kili.dates, window = 3, alpha = 2,
    na.tolerance = 0, simplify = 2, np = max(1, detectCores() / 2), progBar = FALSE, 
    method = "multidimension", dist_m = "twdtw", midpoint = 6, stepness = a, 
    cycle_length = "year", time_scale = "month"
  )
  
  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack <- c(tmp.raster, KiliNP_LandCover_Raster) 
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  if(nrow(tmp.df) > 10000) tmp.df <- tmp.df[sample(nrow(tmp.df), 10000), ]
  
  tmp.permanova <- adonis2(tmp.df$RaosQ ~ tmp.df$Veg_GroundTruth, permutations = 999)
  
  # Log result to dataframe
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  
  # Overwrite the entire CSV with the updated dataframe as a live checkpoint
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  
  # Force R to dump the C++ matrices from RAM before starting the next alpha value
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df)
  gc()
}

# Sort final results to show best at the top
alpha_results <- alpha_results[order(-alpha_results$R2), ]
print("Full Grid Search Complete. Results:")
print(alpha_results)

# ---> CHANGED: Object renamed for PPI
Kili_PPI_Optimal_Alpha <- alpha_results$Alpha[1]
message(paste("The absolute optimal \U03B1 value is:", Kili_PPI_Optimal_Alpha))

### Exports ####

# Save the numerical optimal alpha as an RDS object
saveRDS(Kili_PPI_Optimal_Alpha, file.path(KiliNP_Results, "Kili_PPI_Optimal_Alpha.rds"))

# Export the final, sorted dataframe to guarantee the log is complete and ordered
write.csv(alpha_results, file = log_csv, row.names = FALSE)

message("PPI Grid search complete and exported successfully.")