################################################################################
### Elliot Samuel Shayle - University of Marburg
### 03.3_Kilimanjaro_NDVI_α_Grid-search.R
### TWDTW α (Steepness) Optimization via Transect Grid Search
################################################################################

library(rasterdiv)
library(twdtw)
library(vegan)
library(terra)
library(dplyr)
library(here) # Ensures dynamic, drive-agnostic pathing

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
kili.tiling.grid <- vect(file.path(KiliNP_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON")) # Load it back in
transect_line <- vect(file.path(KiliNP_Processed, "KiliNP_Transect_Line.geojson"))
KiliNP_LandCover_Vector <- vect(file.path(KiliNP_Input, "Kili Ground Truthing Land Cover Classifications", "VegAug1_KILI_SES_withnewcof.shp"))

# 3. Import Full Cleaned NDVI Timeseries
KiliNP_Timeseries_Clean <- rast(file.path(KiliNP_Processed, "KiliNP_2017-2021_Cropped_&_Masked.tif"))

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

message("Cropping and masking NDVI timeseries to transect footprint...")
transect_ts <- crop(KiliNP_Timeseries_Clean, transect_poly)
transect_ts <- mask(transect_ts, transect_poly)

# Create a spatial template for the Land Cover rasterisation
transect_mean <- app(transect_ts, fun = mean, na.rm = TRUE)

### Pre-Loop Land Cover Rasterisation ####
message("Preparing Land Cover Vector and Raster Template...")

# Align CRS
if (crs(transect_mean) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(transect_mean))
}

# Rasterise using the transect mean as the precise spatial template

KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  transect_mean,
  field = "grid_code"
)

### Grid Search for Optimal TWDTW α ####
message("Starting two-stage grid search for optimal TWDTW \U03B1 (steepness)...")

# Initialize Active Log File
log_csv <- file.path(KiliNP_Results, "Kili_NDVI_Alpha_GridSearch_Log.csv")
write.csv(data.frame(Alpha = numeric(), R2 = numeric(), p_value = numeric()), log_csv, row.names = FALSE)

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
    np = max(1, detectCores() - 2), 
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
  
  # Log result to dataframe and immediately write to the CSV checkpoint
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  write.table(step_result, file = log_csv, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
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
    na.tolerance = 0, simplify = 2, np = max(1, detectCores() - 2), progBar = FALSE, 
    method = "multidimension", dist_m = "twdtw", midpoint = 6, stepness = a, 
    cycle_length = "year", time_scale = "month"
  )
  
  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack <- c(tmp.raster, KiliNP_LandCover_Raster) 
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  if(nrow(tmp.df) > 10000) tmp.df <- tmp.df[sample(nrow(tmp.df), 10000), ]
  
  tmp.permanova <- adonis2(tmp.df$RaosQ ~ tmp.df$Veg_GroundTruth, permutations = 999)
  
  # Log result to dataframe and immediately write to the CSV checkpoint
  step_result <- data.frame(Alpha = a, R2 = tmp.permanova$R2[1], p_value = tmp.permanova$`Pr(>F)`[1])
  alpha_results <- rbind(alpha_results, step_result)
  write.table(step_result, file = log_csv, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
}

# Sort final results to show best at the top
alpha_results <- alpha_results[order(-alpha_results$R2), ]
print("Full Grid Search Complete. Results:")
print(alpha_results)

Kili_NDVI_Optimal_Alpha <- alpha_results$Alpha[1]
message(paste("The absolute optimal \U03B1 value is:", Kili_NDVI_Optimal_Alpha))

### Exports ####

# Save the numerical optimal alpha as an RDS object
saveRDS(Kili_NDVI_Optimal_Alpha, file.path(KiliNP_Results, "Kili_NDVI_Optimal_Alpha.rds"))

message("Grid search complete and exported successfully.")