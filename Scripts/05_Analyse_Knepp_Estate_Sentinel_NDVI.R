############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 23/10/2026                    #
# 05_Analyse_Knepp_Estate_Sentinel_NDVI.R                                      #
# Comparative NDVI analysis of TWDTW Rao's Q & classic Rao's Q in Knepp Estate #
############################################################################ ###

### Install and load the necessary packages ####
## This should already be done from the setup file
# rasterdiv now contains the TWDTW-enabled paRao()

library(rasterdiv)
library(twdtw)
library(vegan)
library(pROC)
library(terra)
library(parallel)

# Core spatial stack already loaded in 00_setup.R:
# terra, sf, here, dplyr, stringr

### Define the file paths and import the site data ####
## Output directory for this script

Knepp_Results <- file.path(Results, "Knepp Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Define CPU cores for parallel operations (consistent with kili.cores in 03_Analyse_Kilimanjaro_NDVI.R)

knepp.cores <- max(1, detectCores() - 2)

# Explicitly set the input/processed directories on the external drive to manage storage

Knepp_Input_Ext     <- "K:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Input Data"
Knepp_Processed_Ext <- "K:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Processed Data"
dir.create(Knepp_Processed_Ext, showWarnings = FALSE, recursive = TRUE)

# Load Knepp Estate boundaries (core estate and 5 km buffered boundary)

KneppEstate_Boundaries          <- vect(file.path(Knepp_Input, "Schulte-to-Bühne_Knepp_Site_Data", "knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed, "Knepp_Estate_Buffered_Boundary.shp"))

### Step 1A: Crop, mask, clean NoData, and scale raw annual Sentinel-2 NDVI rasters ####
## Only run this loop if you haven't already generated the cropped annual NDVI files

message("Beginning cropping, masking, and scaling of raw Sentinel-2 NDVI rasters...")

tmp.files <- list.files(
  Knepp_Input_Ext,
  pattern = ".*_HL_TSA_SEN2L_NDVI_FBM\\.tif$",
  full.names = TRUE
)

for (f in tmp.files) {
  
  tmp.Knepp.raster <- rast(f)
  
  # Ensure the buffered boundary shapefile shares the exact CRS of the raw raster
  
  if (crs(tmp.Knepp.raster) != crs(KneppEstate_Buffered_Boundaries)) {
    KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(tmp.Knepp.raster))
  }
  
  # Attempt to crop and mask to the 5 km buffered boundary.
  # If the extents do not overlap (e.g., the corrupted 2018 file), tryCatch safely skips it.
  
  tmp.Knepp.raster.masked <- tryCatch({
    crop(tmp.Knepp.raster, KneppEstate_Buffered_Boundaries, mask = TRUE)
  }, error = function(e) {
    message(paste("  -> Skipping non-overlapping raster:", basename(f)))
    return(NULL)
  })
  
  # If tryCatch returned NULL, clean up and move to the next file
  
  if (is.null(tmp.Knepp.raster.masked)) {
    rm(tmp.Knepp.raster)
    next
  }
  
  # CRITICAL: Remove Sentinel-2 integer NoData flags (-9999 / -10000) and scale to [-1, 1] floats.
  # Doing this now ensures SG smoothing, Shannon's H, and paRao all receive true float NDVI values.
  
  tmp.Knepp.raster.masked[tmp.Knepp.raster.masked <= -9000] <- NA_real_
  tmp.Knepp.raster.masked <- tmp.Knepp.raster.masked / 10000
  
  # Extract the 4-digit year from the start of the filename
  
  tmp.year <- sub("^(\\d{4}).*", "\\1", basename(f))
  
  # Assign clean monthly layer names for this year (e.g., "2019-01_NDVI" ... "2019-12_NDVI")
  
  names(tmp.Knepp.raster.masked) <- paste0(tmp.year, "-", sprintf("%02d", 1:12), "_NDVI")
  
  # Save the cropped, scaled 12-month raster to the processed folder
  
  tmp.new.name <- paste0("Knepp_", tmp.year, "_NDVI_Cropped.tif")
  
  writeRaster(
    tmp.Knepp.raster.masked,
    filename = file.path(Knepp_Processed_Ext, tmp.new.name),
    overwrite = TRUE
  )
  message(paste0("Finished writing scaled raster: ", tmp.new.name))
  
  # Free RAM before next year [pun intended]
  
  rm(tmp.Knepp.raster, tmp.Knepp.raster.masked)
  gc()
}

### Step 1B: Import the cropped annual rasters into a single multi-year stack ####

message("Importing cropped annual NDVI rasters into a continuous multi-year stack...")

Knepp_Cropped_NDVI_Files <- list.files(
  Knepp_Processed_Ext,
  pattern = "^Knepp_\\d{4}_NDVI_Cropped\\.tif$",
  full.names = TRUE
)

# Load and stack all valid annual files chronologically

Knepp_Buffered_Timeseries.NDVI <- rast(Knepp_Cropped_NDVI_Files)

# Extract the years from the filenames and assign standardized layer names ("YYYY-MM_NDVI")

Knepp.years.loaded <- sub("^Knepp_(\\d{4})_NDVI_Cropped\\.tif$", "\\1", basename(Knepp_Cropped_NDVI_Files))
message("Years found on disk: ", paste(Knepp.years.loaded, collapse = ", "))

Knepp.layer.names <- unlist(lapply(Knepp.years.loaded, function(y) {
  paste0(y, "-", sprintf("%02d", 1:12), "_NDVI")
}))
names(Knepp_Buffered_Timeseries.NDVI) <- Knepp.layer.names

## Exclude edge years with incomplete seasonal data (as done in 06_Analyse_Knepp_Estate_Sentinel_PPI.R)
# Note: Since the 2017 year is disconnected from the rest of the timeseries, I will exclude it

years_to_exclude <- c("2017")
exclude_regex    <- paste(years_to_exclude, collapse = "|")

if (any(grepl(exclude_regex, names(Knepp_Buffered_Timeseries.NDVI)))) {
  message("Excluding incomplete edge years (", paste(years_to_exclude, collapse = ", "), ") from continuous stack...")
  Knepp_Buffered_Timeseries.NDVI <- Knepp_Buffered_Timeseries.NDVI[[grep(exclude_regex, names(Knepp_Buffered_Timeseries.NDVI), invert = TRUE)]]
}

### Inspect temporal structure ####
## rasterdiv::paRao() requires an explicit time_vector matching the exact layer order

message("Constructing time vector for continuous Knepp NDVI time series...")

# 1. Grab the remaining layer names (e.g., "2019-01_NDVI")

current.knepp.names <- names(Knepp_Buffered_Timeseries.NDVI)

# 2. Strip "_NDVI" and append "-15" to represent mid-month dates

knepp.dates.string <- paste0(sub("_NDVI$", "", current.knepp.names), "-15")

# 3. Convert to Date objects

Knepp.full.dates <- as.Date(knepp.dates.string, format = "%Y-%m-%d")

# 4. Safety checks: verify no NAs and confirm exact match with number of raster layers

if (any(is.na(Knepp.full.dates))) stop("Date parsing failed! Check layer names.")
stopifnot(length(Knepp.full.dates) == nlyr(Knepp_Buffered_Timeseries.NDVI))

# 5. Check whether the multi-year monthly sequence is strictly contiguous

expected_dates <- seq(min(Knepp.full.dates), max(Knepp.full.dates), by = "month")
if (length(expected_dates) != length(Knepp.full.dates)) {
  warning("Watch out: There is a missing year/month gap in your sequence! Missing dates: ",
          paste(expected_dates[!expected_dates %in% Knepp.full.dates], collapse = ", "))
} else {
  message("Temporal sequence is 100% contiguous from ", min(Knepp.full.dates), " to ", max(Knepp.full.dates),
          " (", length(Knepp.full.dates), " monthly layers).")
}

# Export the unsmoothed continuous stack for safekeeping

writeRaster(
  Knepp_Buffered_Timeseries.NDVI,
  filename = file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Continuous_Raw_Timeseries.tif"),
  overwrite = TRUE
)

Knepp_Buffered_Timeseries.NDVI <- rast(
  file.path(Knepp_Processed,"Knepp_Buffered_NDVI_Continuous_Raw_Timeseries.tif")) # And load it back in as necessary!

### Savitzky-Golay Gap-Filling and Smoothing ####
## Define the pixel-wise gap-filling function (identical to Kilimanjaro case study)

sg_gapfill <- function(x) {
  
  # If a pixel is NA across all layers (outside the 5 km boundary), return NAs immediately
  
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  
  # 1. Linear interpolation across cloud gaps (rule = 2 pads leading/trailing NAs with nearest valid value)
  
  x_interp <- zoo::na.approx(x, na.rm = FALSE, rule = 2)
  
  # 2. Savitzky-Golay smoothing (fl = 5 preserves seasonal amplitude; forder = 2 is standard quadratic)
  
  x_smoothed <- pracma::savgol(
    x_interp,
    fl = 5,
    forder = 2
  )
  
  return(x_smoothed)
}

message("Gap-filling and smoothing the continuous multi-year NDVI dataset with Savitzky-Golay...")

Knepp_Buffered_Timeseries.NDVI.Cleaned <- app(
  Knepp_Buffered_Timeseries.NDVI,
  fun = sg_gapfill,
  cores = knepp.cores
)

# Restore layer names stripped by terra::app()

names(Knepp_Buffered_Timeseries.NDVI.Cleaned) <- names(Knepp_Buffered_Timeseries.NDVI)

# Re-mask against the buffered boundary to ensure background pixels remain strict NAs

if (crs(KneppEstate_Buffered_Boundaries) != crs(Knepp_Buffered_Timeseries.NDVI.Cleaned)) {
  KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(Knepp_Buffered_Timeseries.NDVI.Cleaned))
}
Knepp_Buffered_Timeseries.NDVI.Cleaned <- mask(Knepp_Buffered_Timeseries.NDVI.Cleaned, KneppEstate_Buffered_Boundaries)

message("Gap-filling complete! Writing cleaned multi-year NDVI stack to disk...")

writeRaster(
  Knepp_Buffered_Timeseries.NDVI.Cleaned,
  filename = file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Cleaned_SG-method.tif"),
  overwrite = TRUE
)

# Reload from disk so the SpatRaster points cleanly to the file

Knepp_Buffered_Timeseries.NDVI.Cleaned <- rast(
  file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Cleaned_SG-method.tif")
) # Load it back in!

print(Knepp_Buffered_Timeseries.NDVI.Cleaned)

### Calculate and plot the mean seasonal trajectory across the site ####
# 1. Crop and mask to the inner Knepp Estate boundary

Knepp_Timeseries.NDVI <- crop(
  Knepp_Buffered_Timeseries.NDVI.Cleaned,
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.NDVI.Cleaned)),
  mask = TRUE
)

Knepp_Timeseries_Mean.NDVI <- global(Knepp_Timeseries.NDVI, fun = "mean", na.rm = TRUE)

# 2. Mask out the inside of the Knepp Estate (inverse = TRUE) to isolate the surrounding 5 km farmland

Knepp_Environs_Timeseries.NDVI <- mask(
  Knepp_Buffered_Timeseries.NDVI.Cleaned,
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.NDVI.Cleaned)),
  inverse = TRUE
)

Knepp_Environs_Timeseries.NDVI <- global(Knepp_Environs_Timeseries.NDVI, fun = "mean", na.rm = TRUE)

# 3. Export the 4K comparison plot of the continuous multi-year trajectories

png(file.path(Knepp_Results, "Knepp_Estate_NDVI_Mean_Timeseries.png"),
    width = 3840, height = 2160, res = 150)

plot(Knepp.full.dates, Knepp_Timeseries_Mean.NDVI[, 1],
     type = "l",
     lwd = 2,
     col = "forestgreen",
     ylim = range(c(Knepp_Timeseries_Mean.NDVI[, 1], Knepp_Environs_Timeseries.NDVI[, 1]), na.rm = TRUE),
     xlab = "Date",
     ylab = "Mean NDVI",
     main = "Knepp Estate vs. Surrounding Farmland - Continuous Mean NDVI",
     xaxt = "n")

lines(Knepp.full.dates, Knepp_Environs_Timeseries.NDVI[, 1],
      lwd = 2,
      col = "chocolate")

grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# Minor monthly ticks and major yearly ticks

axis.Date(1, at = Knepp.full.dates, labels = FALSE, tcl = -0.25)
axis.Date(1,
          at = as.Date(paste0(unique(format(Knepp.full.dates, "%Y")), "-01-15")),
          format = "%Y",
          tcl = -0.5)

legend("topright",
       legend = c("Knepp Estate", "Environs surrounding Knepp"),
       col = c("forestgreen", "chocolate"),
       lwd = 2,
       bty = "n")

dev.off()

### 1. Shannon-Wiener Index ####

message("Calculating Shannon-Wiener diversity index across the 2019-2024 SG-filtered stack...")

# Collapse the 72-layer time series into a single multi-year mean NDVI raster

Knepp_Mean_SG <- app(Knepp_Buffered_Timeseries.NDVI.Cleaned, fun = mean, na.rm = TRUE)
writeRaster(Knepp_Mean_SG, file.path(Knepp_Processed, "Knepp_2019-2024_MeanNDVI_SG.tif"), overwrite = TRUE)
Knepp_Mean_SG <- rast(file.path(Knepp_Processed, "Knepp_2019-2024_MeanNDVI_SG.tif")) # Load back in from disk

# Simplify the raster to 2 decimal places to prevent numerical oversaturation in Shannon's H
# (Note: Rasters were already scaled by / 10000 in Step 1, so no division is needed here!)

Knepp_Mean_SG2dec <- trim(round(Knepp_Mean_SG, 2))

# Convert the SpatRaster to a wide matrix (required by ShannonS) and compute Shannon's H (3x3 window)

Knepp.ShannonH.SG.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(Knepp_Mean_SG2dec, wide = TRUE),
  window = 3,
  na.tolerance = 0
)

# Convert the resulting matrix back into a SpatRaster and restore spatial metadata

Knepp_ShannonH_SG <- rast(Knepp.ShannonH.SG.matrix)
ext(Knepp_ShannonH_SG) <- ext(Knepp_Mean_SG2dec)
crs(Knepp_ShannonH_SG) <- crs(Knepp_Mean_SG2dec)
names(Knepp_ShannonH_SG) <- "ShannonH_SG"

# Extend back to the untrimmed Knepp_Mean_SG template so extents match 100% when stacking later

Knepp_ShannonH_SG <- terra::extend(Knepp_ShannonH_SG, Knepp_Mean_SG)

# Export the final Shannon's H raster to Results and reload

writeRaster(
  Knepp_ShannonH_SG,
  filename = file.path(Knepp_Results, "Knepp_2019-2024_NDVI_ShannonH_SG.tif"),
  overwrite = TRUE
)
Knepp_ShannonH_SG <- rast( # Load it back in!
  file.path(Knepp_Results, "Knepp_2019-2024_NDVI_ShannonH_SG.tif"))

message("Shannon's H calculation complete!")

### 2. Classic Rao's Q ####
## Due to the size of the 10m Sentinel-2 site, we tile the dataset into 500 tiles
## for parallel processing on the MaRC3a HPC cluster, mirroring the Kilimanjaro pipeline.

message("Preparing 1000-tile grid for Knepp Estate Classic & TWDTW Rao's Q...")

### Step 1: Create a master grid to define zones for tiling

# Trim outer NA borders from the mean raster so tiles are only generated over valid data

trimmed.Knepp_Mean_Raster <- trim(Knepp_Mean_SG)

# Define target number of tiles (750 tiles for Knepp Estate) and raster aspect ratio

knepp.total.tiles  <- 750
knepp.aspect.ratio <- ncol(trimmed.Knepp_Mean_Raster) / nrow(trimmed.Knepp_Mean_Raster)

# Find all factor pairs that multiply to exactly 750

knepp.tiling.factors <- expand.grid(
  ncols = 1:knepp.total.tiles,
  nrows = 1:knepp.total.tiles
)

knepp.tiling.factors <- knepp.tiling.factors[
  knepp.tiling.factors$ncols * knepp.tiling.factors$nrows == knepp.total.tiles, 
]

# Choose the factor pair whose column-to-row ratio is closest to the raster's aspect ratio

knepp.tiling.factors$ratio_diff <- abs(
  (knepp.tiling.factors$ncols / knepp.tiling.factors$nrows) - knepp.aspect.ratio
)

knepp.best.tile.size <- knepp.tiling.factors[which.min(knepp.tiling.factors$ratio_diff), ]

knepp.cols <- knepp.best.tile.size$ncols
knepp.rows <- knepp.best.tile.size$nrows

# Generate the spatial polygon tiling grid

knepp.tiling.grid <- as.polygons(
  rast(
    ext(trimmed.Knepp_Mean_Raster),
    ncols = knepp.cols,
    nrows = knepp.rows,
    crs   = crs(trimmed.Knepp_Mean_Raster)
  )
)

# Export the tiling grid as a GeoJSON for reproducibility and use in the 05.3 alpha grid search

knepp.tiles.filepath <- file.path(Knepp_Processed, "Knepp_Tiling_Grid_Polygons.geojson")
writeVector(knepp.tiling.grid, knepp.tiles.filepath, filetype = "GeoJSON", overwrite = TRUE)
knepp.tiling.grid <- vect(knepp.tiles.filepath)

plot(knepp.tiling.grid, main = "Knepp Estate 750-Tile Grid (2019-2024)")

# Define moving window size and 1-pixel overlap buffer so Rao's Q has no tile-edge NAs

RaoQ.window.size   <- 3
knepp.tile.overlap <- floor(RaoQ.window.size / 2)

# Create directory for Classic Rao's Q input tiles (Mean NDVI)

knepp.tile.dir.sg <- file.path(Knepp_Processed, "Knepp_Tiles_NDVI_SG-Filtered")
dir.create(knepp.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

message("Generating 750 Mean NDVI tiles for Classic Rao's Q...")

makeTiles(
  trimmed.Knepp_Mean_Raster,
  y        = knepp.tiling.grid,
  buffer   = knepp.tile.overlap,
  filename = file.path(knepp.tile.dir.sg, "Knepp_MeanNDVI_SG_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Submit the Classic Rao's Q tiles for processing on MaRC3a
## Please see script "05.1A_Knepp_NDVI_Classic-RaoQ_SG-Filtered.R" and its SLURM wrapper
## for the HPC batch job execution. Run Step 3 below once all 750 tiles have finished on MaRC3a.

### Step 3: Demosaic the classical Rao's Q tiles

message("Demosaicing Classic Rao's Q tiles...")

knepp.rao.files.sg <- list.files(
  file.path(Knepp_Processed, "Knepp_Tiles_NDVI_SG-Filtered", "Rao-utputs"),
  pattern    = "Knepp_NDVI_Classic-RaoQ_SG_Tile-",
  full.names = TRUE
)

if (length(knepp.rao.files.sg) == 0) {
  stop("No Classic Rao's Q output tiles found! Ensure MaRC3a jobs have completed and synced.")
}

# Stitch all 750 tiles back together using a SpatRasterCollection (sprc)

Knepp_Classic_RaoQ_SG <- terra::mosaic(sprc(lapply(knepp.rao.files.sg, rast)))

# Extend to match the exact master extent of Knepp_Mean_SG

Knepp_Classic_RaoQ_SG <- terra::extend(Knepp_Classic_RaoQ_SG, Knepp_Mean_SG)
names(Knepp_Classic_RaoQ_SG) <- "RaosQ_Classic_SG"

writeRaster(
  Knepp_Classic_RaoQ_SG,
  filename = file.path(Knepp_Results, "Knepp_2019-2024_NDVI_Classic-RaoQ_SG.tif"),
  overwrite = TRUE
)
Knepp_Classic_RaoQ_SG <- rast(file.path(Knepp_Results, "Knepp_2019-2024_NDVI_Classic-RaoQ_SG.tif"))

plot(Knepp_Classic_RaoQ_SG, main = "Knepp Estate Classic Rao's Q (NDVI 2019-2024)")

### 3. Rao's Q with TWDTW ####

message("Preparing 72-layer timeseries tiles for TWDTW Rao's Q...")

### Step 1: Tile the continuous 2019-2024 NDVI timeseries stack

knepp.twdtw.tile.dir.sg <- file.path(Knepp_Processed, "Knepp_TS_Tiles_NDVI_SG-Filtered")
dir.create(knepp.twdtw.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

# Save the 72-month date vector to disk so the MaRC3a worker scripts can load it directly
# writeLines(
#   as.character(Knepp.full.dates),
#   file.path(knepp.twdtw.tile.dir.sg, "time_vector.txt")
# )

message("Generating 750 Timeseries NDVI tiles (2019-2024) for TWDTW Rao's Q...")

makeTiles(
  trim(Knepp_Buffered_Timeseries.NDVI.Cleaned),
  y        = knepp.tiling.grid,
  buffer   = knepp.tile.overlap,
  filename = file.path(knepp.twdtw.tile.dir.sg, "Knepp_TS_NDVI_SG_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Run the Alpha Grid Search (05.3) and submit TWDTW tiles on MaRC3a (05.2A)
## 1. Run "05.3_Knepp_NDVI_α_Grid-search.R" to find the optimal TWDTW stepness (α).
## 2. Submit "05.2A_Knepp_NDVI_TWDTW-RaoQ_SG-Filtered.R" via SLURM on MaRC3a.
## 3. Once the 750 TWDTW tiles are synced back, run Step 3 below.

### Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster

message("Demosaicing TWDTW Rao's Q tiles...")

knepp.twdtw.rao.files.sg <- list.files(
  file.path(Knepp_Processed, "Knepp_TS_Tiles_NDVI_SG-Filtered", "Rao-utputs"),
  pattern    = "Knepp_NDVI_TWDTW-RaoQ_SG_Tile-",
  full.names = TRUE
)

if (length(knepp.twdtw.rao.files.sg) == 0) {
  stop("No TWDTW Rao's Q output tiles found! Ensure MaRC3a jobs have completed and synced.")
}

Knepp_TWDTW_RaoQ_SG <- terra::mosaic(sprc(lapply(knepp.twdtw.rao.files.sg, rast)))

# Extend to match the exact master extent of Knepp_Mean_SG

Knepp_TWDTW_RaoQ_SG <- terra::extend(Knepp_TWDTW_RaoQ_SG, Knepp_Mean_SG)
names(Knepp_TWDTW_RaoQ_SG) <- "RaosQ_TWDTW_SG"

writeRaster(
  Knepp_TWDTW_RaoQ_SG,
  filename = file.path(Knepp_Results, "Knepp_2019-2024_NDVI_TWDTW-RaoQ_SG.tif"),
  overwrite = TRUE
)
Knepp_TWDTW_RaoQ_SG <- rast(file.path(Knepp_Results, "Knepp_2019-2024_NDVI_TWDTW-RaoQ_SG.tif"))

plot(Knepp_TWDTW_RaoQ_SG, main = "Knepp Estate TWDTW Rao's Q (NDVI 2019-2024)")

### Export all rasters for comparison ####

Knepp_Comparison_Rasters <- c(
  Knepp_Mean_SG,
  Knepp_ShannonH_SG,
  Knepp_Classic_RaoQ_SG,
  Knepp_TWDTW_RaoQ_SG
)

names(Knepp_Comparison_Rasters) <- c(
  "Sentinel-2 MeanNDVI SG Gap-Filled",
  "ShannonH SG (NDVI)",
  "RaosQ Classic SG (NDVI)",
  "RaosQ TWDTW SG (NDVI)"
)

writeRaster(
  Knepp_Comparison_Rasters,
  filename = file.path(Knepp_Results, "Knepp_2019-2024_NDVI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

Knepp_Comparison_Rasters <- rast(file.path(Knepp_Results, "Knepp_2019-2024_NDVI_Diversity_Comparison.tif"))

# Export a 2x2 multipanel comparison figure for the manuscript

png(file.path(Knepp_Results, "Knepp_NDVI-Based_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(Knepp_Comparison_Rasters)

dev.off()

plot(Knepp_Comparison_Rasters)