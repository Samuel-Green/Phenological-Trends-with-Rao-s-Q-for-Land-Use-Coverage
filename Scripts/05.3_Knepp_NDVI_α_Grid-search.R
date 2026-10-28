############################################################################ ###
### Elliot Samuel Shayle - University of Marburg - 27/10/2026                  #
### 05.3_Knepp_NDVI_α_Grid-search.R                                            #
### TWDTW α (Steepness) Optimization via Consensus UKCEH Ground Truth (No Bin) #
############################################################################ ###

library(rasterdiv)
library(twdtw)
library(vegan)
library(terra)
library(dplyr)
library(here)     # Ensures dynamic, drive-agnostic pathing
library(parallel)

### 0. Linux Path Safety Check ####
all_dirs       <- list.dirs(here(), recursive = TRUE)
danger_pattern <- "[ \\(\\)\\[\\]\\{\\}&\\|\\*\\?\\$\\<\\>\\'\"]|[^\x01-\x7F]"
bad_dirs       <- all_dirs[grepl(danger_pattern, all_dirs)]

if (length(bad_dirs) > 0) {
  message("!! WARNING: Found paths with spaces/special characters:")
  print(bad_dirs)
} else {
  message("ALL CLEAR: Your folder paths are perfectly Linux-safe!")
}

### Define the file paths and import the site data ####
## 1. Internal Project Paths (Works on Keller Computer, or uncomment for MaRC3a)
# Keller Computer paths
InputData     <- here("Data/Input Data")
ProcessedData <- here("Data/Processed Data")
Results       <- here("Results 📈📉")

# MaRC3a server paths
# InputData     <- here("Data/Input_Data")
# ProcessedData <- here("Data/Processed_Data")
# Results       <- here("Results")

# Keller Computer paths
Knepp_Input     <- file.path(InputData, "Knepp Estate")
Knepp_Processed <- file.path(ProcessedData, "Knepp Estate")
Knepp_Results   <- file.path(Results, "Knepp Estate")

# MaRC3a server paths
# Knepp_Input     <- file.path(InputData, "Knepp_Estate")
# Knepp_Processed <- file.path(ProcessedData, "Knepp_Estate")
# Knepp_Results   <- file.path(Results, "Knepp_Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

## 2. External Drive Paths (for the massive raw Sentinel-2 inputs)
# If your external drive sits outside the R project on the Keller Computer, 
# you can use a simple if() check to point to it automatically:
if (dir.exists("E:/Elliot Shayle/Knepp Estate Geodata")) {
  # We are on the Keller Computer (Windows E: Drive)
  Knepp_Input_Ext     <- "E:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Input Data"
  Knepp_Processed_Ext <- "E:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Processed Data"
} else {
  # We are on MaRC3a (Linux), so just use the standard project folders
  Knepp_Input_Ext     <- Knepp_Input
  Knepp_Processed_Ext <- Knepp_Processed
}
dir.create(Knepp_Processed_Ext, showWarnings = FALSE, recursive = TRUE)

# Cap parallel workers to a safe limit (e.g., 12 cores) to prevent RAM duplication OOM crashes
knepp.cores <- max(1, detectCores() - 2)
message("Using ", knepp.cores, " parallel CPU cores for paRao and PERMANOVA.")

### 2. Import Spatial Vectors & Cleaned 2019-2024 NDVI Timeseries ####
message("Importing spatial vectors and 2019-2024 SG-filtered NDVI timeseries...")

# Load the 750-tile grid created in Step 2 of 05_Analyse_Knepp_Estate_Sentinel_NDVI.R
knepp.tiling.grid <- vect(file.path(Knepp_Processed, "Knepp_Tiling_Grid_Polygons.geojson"))

# Load the inner Knepp Estate boundary and outer 5 km buffered boundary
KneppEstate_Boundaries          <- vect(file.path(Knepp_Input, "Schulte-to-Bühne_Knepp_Site_Data", "knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed, "Knepp_Estate_Buffered_Boundary.shp"))

# Load the full 72-layer (2019-2024) SG-smoothed float NDVI raster
Knepp_Timeseries_Clean <- rast(file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Cleaned_SG-method.tif"))

# Ensure all vectors share the exact CRS of the NDVI raster stack
if (crs(knepp.tiling.grid) != crs(Knepp_Timeseries_Clean)) {
  knepp.tiling.grid <- project(knepp.tiling.grid, crs(Knepp_Timeseries_Clean))
}
if (crs(KneppEstate_Boundaries) != crs(Knepp_Timeseries_Clean)) {
  KneppEstate_Boundaries <- project(KneppEstate_Boundaries, crs(Knepp_Timeseries_Clean))
}
if (crs(KneppEstate_Buffered_Boundaries) != crs(Knepp_Timeseries_Clean)) {
  KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(Knepp_Timeseries_Clean))
}

# Dynamically extract the 72-month Date vector from the layer names (e.g., "2019-01_NDVI")
current.knepp.names <- names(Knepp_Timeseries_Clean)
knepp.dates.string  <- paste0(sub("_NDVI$", "", current.knepp.names), "-15")
Knepp.dates         <- as.Date(knepp.dates.string, format = "%Y-%m-%d")

stopifnot(length(Knepp.dates) == nlyr(Knepp_Timeseries_Clean))
message("Loaded ", length(Knepp.dates), " monthly layers spanning ", min(Knepp.dates), " to ", max(Knepp.dates), ".")

### 3. Stratified Random Tile Subsampling (No Transect Required) ####
message("Selecting a stratified random subsample of 25 valid tiles (10 Inside Estate, 15 Outside Buffer)...")

# Add a unique Tile_ID attribute so we can distinguish inner vs. outer tiles
knepp.tiling.grid$Tile_ID <- seq_len(nrow(knepp.tiling.grid))

# Step 3a: Keep only tiles that fall within the 5 km buffered study area (excludes empty NA corners)
valid_tiles <- knepp.tiling.grid[KneppEstate_Buffered_Boundaries, ]

# Step 3b: Identify tiles that intersect the rewilded Knepp Estate core vs. surrounding farmland
inner_tiles <- valid_tiles[KneppEstate_Boundaries, ]
outer_tiles <- valid_tiles[!valid_tiles$Tile_ID %in% inner_tiles$Tile_ID, ]

# Step 3c: Randomly sample 10 inner tiles and 15 outer tiles for a balanced 25-tile tuning set
set.seed(42)
n_inner <- min(10, nrow(inner_tiles))
n_outer <- min(15, nrow(outer_tiles))

sampled_inner <- inner_tiles[sample(seq_len(nrow(inner_tiles)), n_inner), ]
sampled_outer <- outer_tiles[sample(seq_len(nrow(outer_tiles)), n_outer), ]

sample_tiles <- rbind(sampled_inner, sampled_outer)

# Dissolve selected tiles into a single multi-polygon mask
sample_poly <- aggregate(sample_tiles)

message("Cropping, masking, and trimming NDVI timeseries to the 25 sampled tiles...")
sample_ts <- crop(Knepp_Timeseries_Clean, sample_poly)
sample_ts <- mask(sample_ts, sample_poly)
sample_ts <- trim(sample_ts)

### 4. The Array Severance Trick (Full 10m Resolution - NO BINNING) ####
message("Severing lazy C++ pointers and mapping strict numeric memory (No Binning)...")

# 1. Extract the lazy raster into a raw base R array
sample_array <- terra::as.array(sample_ts)

# 2. Multiply by 1.0 to block ALTREP from converting NA windows to logicals
sample_array <- sample_array * 1.0

# 3. Rebuild the raster directly in physical RAM
sample_ts_RAM <- terra::rast(
  sample_array,
  crs = crs(sample_ts),
  ext = ext(sample_ts)
)
names(sample_ts_RAM) <- names(sample_ts)

# 4. Bake the strictly numeric raster to disk so parallel workers safely read from a file path
message("Caching full-resolution numeric subset raster to disk for safe parallel transmission...")
tmp.tif_path <- file.path(Knepp_Processed, "Knepp_NDVI_Grid_Search_Sample_FullRes.tif")
writeRaster(sample_ts_RAM, tmp.tif_path, datatype = "FLT8S", overwrite = TRUE)

# 5. Free temporary RAM objects
rm(sample_array, sample_ts_RAM)
gc()

# 6. Reload the file-backed raster
sample_ts <- terra::rast(tmp.tif_path)

# Create a single-layer spatial template for aligning all UKCEH Land Cover maps
sample_mean <- app(sample_ts, fun = mean, na.rm = TRUE)

### 5. Pre-Loop Strict Consensus UKCEH Land Cover Preparation (2019-2024) ####
message("Loading and aligning all available UKCEH Land Cover rasters (2019-2024)...")

# Define the 21 UKCEH Land Cover Classes dictionary
UKCEH_LC_Classes <- c(
  "1"  = "Broadleaved woodland",
  "2"  = "Coniferous woodland",
  "3"  = "Arable and horticulture",
  "4"  = "Improved grassland",
  "5"  = "Neutral grassland",
  "6"  = "Calcareous grassland",
  "7"  = "Acid grassland",
  "8"  = "Fen, marsh and swamp",
  "9"  = "Heather",
  "10" = "Heather grassland",
  "11" = "Bog",
  "12" = "Inland rock",
  "13" = "Saltwater",
  "14" = "Freshwater",
  "15" = "Supralittoral rock",
  "16" = "Supralittoral sediment",
  "17" = "Littoral rock",
  "18" = "Littoral sediment",
  "19" = "Saltmarsh",
  "20" = "Urban",
  "21" = "Suburban"
)

# Define the dominant, ecologically meaningful land cover classes for stratified PERMANOVA
target_classes <- c(
  "Broadleaved woodland",
  "Coniferous woodland",
  "Arable and horticulture",
  "Improved grassland",
  "Neutral grassland",
  "Acid grassland",
  "Urban",
  "Suburban"
)

# Dynamically discover all UKCEH LCM .tif files inside the CEH Land-cover data directory
ceh_base_dir <- file.path(Knepp_Input, "CEH Land-cover data")
if (!dir.exists(ceh_base_dir)) {
  ceh_base_dir <- file.path(Knepp_Input, "CEH_Land-cover_data")
}

ceh_tif_files <- list.files(
  ceh_base_dir,
  pattern    = "^LCM.*\\.tif$",
  recursive  = TRUE,
  full.names = TRUE
)

# Extract the 4-digit year from each path and keep only years within our 2019-2024 study window
ceh_years     <- sub(".*(20\\d{2}).*", "\\1", ceh_tif_files)
valid_ceh_idx <- which(ceh_years %in% as.character(2019:2024))

ceh_tif_files <- ceh_tif_files[valid_ceh_idx]
ceh_years     <- ceh_years[valid_ceh_idx]
names(ceh_tif_files) <- ceh_years

if (length(ceh_tif_files) == 0) {
  stop("No UKCEH LCM .tif files found for 2019-2024 in: ", ceh_base_dir)
}

message("Found UKCEH Land Cover maps for years: ", paste(ceh_years, collapse = ", "))

# Project and snap each UKCEH year's Band 1 (integer class 1-21) to the 10m sample_mean template
Knepp_LC_Aligned_List <- lapply(ceh_years, function(y) {
  message("  -> Aligning UKCEH Land Cover for year: ", y)
  lc_raw     <- rast(ceh_tif_files[[y]])[[1]]
  lc_aligned <- project(lc_raw, sample_mean, method = "near") # Nearest-neighbour preserves integers 1:21
  lc_aligned <- mask(lc_aligned, sample_mean)
  return(lc_aligned)
})

# Stack all aligned integer UKCEH layers
Knepp_LC_Stack <- rast(Knepp_LC_Aligned_List)

# STRICT CONSENSUS MASK (Option A):
# A pixel's land cover is invariant across 2019-2024 if and only if its minimum class ID
# equals its maximum class ID across all years. Any pixel that changed class is masked to NA.
message("Computing strict multi-year UKCEH consensus mask (retaining only invariant pixels)...")
lc_min         <- min(Knepp_LC_Stack, na.rm = FALSE)
lc_max         <- max(Knepp_LC_Stack, na.rm = FALSE)
consensus_mask <- (lc_min == lc_max)

Knepp_LandCover_Consensus <- mask(Knepp_LC_Stack[[1]], consensus_mask, maskvalues = 0)
names(Knepp_LandCover_Consensus) <- "Veg_GroundTruth"

### 6. Helper Function: Run paRao & Stratified PERMANOVA for a Single Alpha ####
run_knepp_alpha_step <- function(a_val, sample_ts, dates_vec, lc_consensus, ukceh_lookup, target_classes, n_cores) {
  
  message(paste0(Sys.time(), " | Testing \U03B1 = ", a_val, "..."))
  
  tmp.RaoQ <- paRao(
    x            = sample_ts,
    time_vector  = dates_vec,
    window       = 3,
    alpha        = 2,
    na.tolerance = 0,
    simplify     = 2,
    np           = n_cores,
    progBar      = FALSE,
    method       = "multidimension",
    dist_m       = "twdtw",
    midpoint     = 2,     # 2 months (or change to 3 for 3-month temperate season inflection)
    stepness     = a_val,
    cycle_length = "year",
    time_scale   = "month"
  )
  
  message(paste0(Sys.time(), " | paRao completed for \U03B1 = ", a_val))
  
  tmp.raster <- tmp.RaoQ$window.3$alpha.2
  tmp.stack  <- c(tmp.raster, lc_consensus)
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  
  # Map raw integer UKCEH codes (1:21) to descriptive strings in base R (bypasses terra factor bugs)
  tmp.df$Veg_GroundTruth <- ukceh_lookup[as.character(tmp.df$Veg_GroundTruth)]
  
  # Filter out 0/NA artifacts and subset to target classes
  tmp.df_filtered <- tmp.df %>%
    filter(!is.na(Veg_GroundTruth) & Veg_GroundTruth %in% target_classes) %>%
    filter(!is.na(RaosQ) & RaosQ != 0)
  
  # Stratified 10,000-pixel subsample across the present classes
  present_classes  <- unique(tmp.df_filtered$Veg_GroundTruth)
  pixels_per_class <- floor(10000 / length(present_classes))
  
  set.seed(42)
  tmp.df_sub <- tmp.df_filtered %>%
    group_by(Veg_GroundTruth) %>%
    slice_sample(n = pixels_per_class, replace = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  
  message(paste0(Sys.time(), " | Running PERMANOVA on ", nrow(tmp.df_sub), " stratified consensus pixels..."))
  
  tmp.permanova <- adonis2(
    tmp.df_sub$RaosQ ~ tmp.df_sub$Veg_GroundTruth,
    permutations = 999,
    parallel     = n_cores
  )
  
  step_result <- data.frame(
    Alpha   = a_val,
    R2      = tmp.permanova$R2[1],
    F_stat  = tmp.permanova$F[1],
    p_value = tmp.permanova$`Pr(>F)`[1]
  )
  
  # Free memory immediately before the next alpha iteration
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df, tmp.df_filtered, tmp.df_sub, tmp.permanova)
  gc()
  
  return(step_result)
}

### 7. Two-Stage Grid Search for Optimal TWDTW α ####
message("Starting two-stage grid search for optimal TWDTW \U03B1 (steepness)...")

log_csv <- file.path(Knepp_Results, "Knepp_NDVI_Alpha_GridSearch_NoBin_Log.csv")

alpha_coarse_grid <- c(-0.1, -0.3, -0.5, -0.7, -0.9)
alpha_results     <- data.frame(Alpha = numeric(), R2 = numeric(), F_stat = numeric(), p_value = numeric())

## STAGE 1: Coarse Grid
for (a in alpha_coarse_grid) {
  step_result <- run_knepp_alpha_step(
    a_val          = a,
    sample_ts      = sample_ts,
    dates_vec      = Knepp.dates,
    lc_consensus   = Knepp_LandCover_Consensus,
    ukceh_lookup   = UKCEH_LC_Classes,
    target_classes = target_classes,
    n_cores        = knepp.cores
  )
  
  alpha_results <- rbind(alpha_results, step_result)
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  message(paste0("Logged \U03B1 = ", a, " | R2 = ", round(step_result$R2, 4)))
}

best.coarse.alpha <- alpha_results$Alpha[which.max(alpha_results$R2)]
message(paste("Best coarse \U03B1 is:", best.coarse.alpha))

## STAGE 2: Fine Grid
alpha_fine_grid <- seq(best.coarse.alpha + 0.1, best.coarse.alpha - 0.1, by = -0.05)
alpha_fine_grid <- round(alpha_fine_grid[alpha_fine_grid < 0], 2)

for (a in alpha_fine_grid) {
  if (a %in% round(alpha_coarse_grid, 2)) next
  
  step_result <- run_knepp_alpha_step(
    a_val          = a,
    sample_ts      = sample_ts,
    dates_vec      = Knepp.dates,
    lc_consensus   = Knepp_LandCover_Consensus,
    ukceh_lookup   = UKCEH_LC_Classes,
    target_classes = target_classes,
    n_cores        = knepp.cores
  )
  
  alpha_results <- rbind(alpha_results, step_result)
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  message(paste0("Logged \U03B1 = ", a, " | R2 = ", round(step_result$R2, 4)))
}

### 8. Sort and Export Final Results ####
alpha_results <- alpha_results[order(-alpha_results$R2), ]
print("Full Knepp NDVI Grid Search Complete. Results:")
print(alpha_results)

Knepp_NDVI_Optimal_Alpha <- alpha_results$Alpha[1]
message(paste("The absolute optimal \U03B1 value is:", Knepp_NDVI_Optimal_Alpha))

# Save optimal alpha to the exact RDS filename expected by 05.2A_Knepp_NDVI_TWDTW-RaoQ_SG-Filtered.R
saveRDS(Knepp_NDVI_Optimal_Alpha, file.path(Knepp_Results, "Knepp_NDVI_Optimal_Alpha.rds"))
saveRDS(Knepp_NDVI_Optimal_Alpha, file.path(Knepp_Results, "Knepp_NDVI_Optimal_Alpha_NoBin.rds"))

write.csv(alpha_results, file = log_csv, row.names = FALSE)
message("All Knepp NDVI grid search results exported successfully!")