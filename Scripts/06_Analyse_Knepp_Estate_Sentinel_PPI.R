############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 03/11/2026                    #
# 06_Analyse_Knepp_Estate_Sentinel_PPI.R                                       #
# Comparative PPI analysis of TWDTW Rao's Q & classic Rao's Q in Knepp Estate  #
############################################################################ ###

### Install and load the necessary packages ####
## This should already be done from the setup file
# rasterdiv now contains the TWDTW-enabled paRao()

library(rasterdiv)
library(twdtw)
library(vegan)
library(pROC)
library(ppi)
library(terra)
library(parallel)

# Core spatial stack already loaded in 00_setup.R:
# terra, sf, here, dplyr, stringr

### Define the file paths and import the site data ####
## Output directory for this script

Knepp_Results <- file.path(Results, "Knepp Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Define CPU cores for parallel operations (consistent with 05_Analyse_Knepp_Estate_Sentinel_NDVI.R)

knepp.cores     <- max(1, detectCores() - 2)
permanova.cores <- min(8, knepp.cores) # Capped at 8 cores to prevent 17,500x17,500 distance matrix OOM crashes

# Explicitly set the input/processed directories on the external drive to manage storage

Knepp_Input_Ext     <- "K:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Input Data"
Knepp_Processed_Ext <- "K:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Processed Data"
dir.create(Knepp_Processed_Ext, showWarnings = FALSE, recursive = TRUE)

# Load Knepp Estate boundaries (core estate and 5 km buffered boundary)

KneppEstate_Boundaries          <- vect(file.path(Knepp_Input, "Schulte-to-Bühne_Knepp_Site_Data", "knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed, "Knepp_Estate_Buffered_Boundary.shp"))

### Step 1A: Crop and mask raw annual Sentinel-2 reflectance rasters (B02, B03, B04, B08) ####
## Only run this loop if you haven't already generated the cropped annual reflectance files

message("Beginning cropping and masking of raw Sentinel-2 reflectance rasters...")

## List raster files
# The regex explicitly looks for BLUE, GREEN, RED, or NIR
tmp.files <- list.files(
  Knepp_Input_Ext,
  pattern    = "SEN2L_(BLUE|GREEN|RED|NIR)_FBM\\.tif$",
  full.names = TRUE
)

# A quick mapping dictionary to convert the text colours to standard Sentinel-2 band numbers
band_map <- c("BLUE" = "B02", "GREEN" = "B03", "RED" = "B04", "NIR" = "B08")

for (f in tmp.files) {
  
  tmp.Knepp.raster <- rast(f)
  
  # Ensure the buffered boundary shapefile shares the exact CRS of the raw raster
  if (crs(tmp.Knepp.raster) != crs(KneppEstate_Buffered_Boundaries)) {
    KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(tmp.Knepp.raster))
  }
  
  # Attempt to crop and mask to the 5 km buffered boundary.
  # If the extents do not overlap (e.g., the corrupted 2018 files), tryCatch safely skips it.
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
  
  # Extract the 4-digit year from the start of the filename
  tmp.year <- sub("^(\\d{4}).*", "\\1", basename(f))
  
  # Extract the colour string (BLUE, GREEN, RED, or NIR) and map to standard band notation
  tmp.color <- sub(".*_SEN2L_([A-Z]+)_FBM\\.tif$", "\\1", basename(f))
  tmp.band  <- band_map[tmp.color]
  
  # Create new filename incorporating both year and band
  tmp.new.name <- paste0("Knepp_", tmp.year, "_", tmp.band, "_Cropped.tif")
  
  # Write to processed folder on the external drive
  writeRaster(
    tmp.Knepp.raster.masked,
    filename  = file.path(Knepp_Processed_Ext, tmp.new.name),
    overwrite = TRUE
  )
  message(paste0("Finished writing cropped reflectance raster: ", tmp.new.name))
  
  # Cleanup RAM
  rm(tmp.Knepp.raster, tmp.Knepp.raster.masked)
  gc()
}
rm(band_map) # Must sit outside the loop

### Step 1B: Import the cropped multi-band rasters into continuous multi-year stacks ####

message("Importing the cropped multi-band rasters...")

# 1. Blue Band (B02)
Knepp_B02_Files       <- list.files(Knepp_Processed_Ext, pattern = "^Knepp_\\d{4}_B02_Cropped\\.tif$", full.names = TRUE)
Knepp_Blue_Timeseries <- rast(Knepp_B02_Files)

# 2. Green Band (B03)
Knepp_B03_Files        <- list.files(Knepp_Processed_Ext, pattern = "^Knepp_\\d{4}_B03_Cropped\\.tif$", full.names = TRUE)
Knepp_Green_Timeseries <- rast(Knepp_B03_Files)

# 3. Red Band (B04)
Knepp_B04_Files      <- list.files(Knepp_Processed_Ext, pattern = "^Knepp_\\d{4}_B04_Cropped\\.tif$", full.names = TRUE)
Knepp_Red_Timeseries <- rast(Knepp_B04_Files)

# 4. Near-Infrared Band (B08)
Knepp_B08_Files      <- list.files(Knepp_Processed_Ext, pattern = "^Knepp_\\d{4}_B08_Cropped\\.tif$", full.names = TRUE)
Knepp_NIR_Timeseries <- rast(Knepp_B08_Files)

message("Applying standard layer names...")

# Extract years directly from the loaded Blue Timeseries file sources
Knepp.years.loaded <- sub(".*Knepp_(\\d{4})_B02_Cropped\\.tif$", "\\1", basename(Knepp_B02_Files))
message("Years found on disk: ", paste(Knepp.years.loaded, collapse = ", "))

# Create clean base layer names (e.g., "2019-01", "2019-02")
Knepp.layer.names.base <- unlist(lapply(Knepp.years.loaded, function(y) {
  paste0(y, "-", sprintf("%02d", 1:12))
}))

# Assign the generated names with explicit band suffixes to each object
names(Knepp_Blue_Timeseries)  <- paste0(Knepp.layer.names.base, "_B02")
names(Knepp_Green_Timeseries) <- paste0(Knepp.layer.names.base, "_B03")
names(Knepp_Red_Timeseries)   <- paste0(Knepp.layer.names.base, "_B04")
names(Knepp_NIR_Timeseries)   <- paste0(Knepp.layer.names.base, "_B08")

## Exclude disconnected/incomplete edge years (2017) prior to PPI calculation
# Filtering 2017 out now ensures M (maximum canopy DVI) is calibrated strictly on 2019-2024
years_to_exclude <- c("2017")
exclude_regex    <- paste(years_to_exclude, collapse = "|")

if (any(grepl(exclude_regex, names(Knepp_NIR_Timeseries)))) {
  message("Excluding disconnected edge years (", paste(years_to_exclude, collapse = ", "), ") from continuous band stacks...")
  keep_idx               <- grep(exclude_regex, names(Knepp_NIR_Timeseries), invert = TRUE)
  Knepp_Blue_Timeseries  <- Knepp_Blue_Timeseries[[keep_idx]]
  Knepp_Green_Timeseries <- Knepp_Green_Timeseries[[keep_idx]]
  Knepp_Red_Timeseries   <- Knepp_Red_Timeseries[[keep_idx]]
  Knepp_NIR_Timeseries   <- Knepp_NIR_Timeseries[[keep_idx]]
}

### Inspect temporal structure ####
## Construct and verify the 72-month time vector (2019-2024)

message("Constructing time vector for continuous Knepp PPI time series...")

current_names        <- names(Knepp_NIR_Timeseries)
dynamic_dates_string <- paste0(sub("_B08$", "", current_names), "-15")
Knepp.full.dates     <- as.Date(dynamic_dates_string, format = "%Y-%m-%d")

if (any(is.na(Knepp.full.dates))) stop("Date parsing failed! Check your layer names.")
stopifnot(length(Knepp.full.dates) == nlyr(Knepp_NIR_Timeseries))

expected_dates <- seq(min(Knepp.full.dates), max(Knepp.full.dates), by = "month")
if (length(expected_dates) != length(Knepp.full.dates)) {
  warning("Watch out: There is a missing year/month gap in your sequence! Missing dates: ",
          paste(expected_dates[!expected_dates %in% Knepp.full.dates], collapse = ", "))
} else {
  message("Temporal sequence is 100% contiguous from ", min(Knepp.full.dates), " to ", max(Knepp.full.dates),
          " (", length(Knepp.full.dates), " monthly layers).")
}

### Compute the Continuous 2019-2024 PPI Raster ####

message("Calculating Astronomical Solar Zenith Angles for Knepp Estate...")

## 1. Calculate Astronomical Solar Zenith Angles
knepp_ext     <- ext(Knepp_NIR_Timeseries)
centroid_geom <- vect(
  matrix(c(mean(c(knepp_ext$xmin, knepp_ext$xmax)),
           mean(c(knepp_ext$ymin, knepp_ext$ymax))), ncol = 2),
  crs = crs(Knepp_NIR_Timeseries)
)

# Project the centroid to WGS84 to get the latitude in degrees (~50.9 degrees)
centroid_latlon <- project(centroid_geom, "EPSG:4326")
knepp_lat_deg   <- geom(centroid_latlon)[, "y"]
knepp_lat_rad   <- knepp_lat_deg * (pi / 180)

day_of_year <- as.numeric(strftime(Knepp.full.dates, format = "%j"))

# Calculate solar declination angle (delta) in radians
declination_deg <- 23.45 * sin((2 * pi / 365) * (day_of_year - 81))
declination_rad <- declination_deg * (pi / 180)

# Calculate hour angle (h) in radians (Sentinel-2 descends at ~10:30 AM local solar time)
hour_angle_deg <- -1.5 * 15
hour_angle_rad <- hour_angle_deg * (pi / 180)

# Calculate final Solar Zenith Angle (theta_z) vector in radians
cos_theta_z    <- sin(knepp_lat_rad) * sin(declination_rad) +
                  cos(knepp_lat_rad) * cos(declination_rad) * cos(hour_angle_rad)
sza_rad_vector <- acos(cos_theta_z)

## 2. Calculate Difference Vegetation Index (DVI)
message("Calculating Difference Vegetation Index (DVI)...")

# Strip any Sentinel-2 integer NoData flags (<= -9000) before subtracting Red from NIR
Knepp_NIR_Timeseries[Knepp_NIR_Timeseries <= -9000] <- NA_real_
Knepp_Red_Timeseries[Knepp_Red_Timeseries <= -9000] <- NA_real_

Knepp_Buffered_Timeseries.DVI <- Knepp_NIR_Timeseries - Knepp_Red_Timeseries

## 3. Define the Safe Hybrid PPI Wrapper
message("Defining PPI Wrapper...")

calc_pixel_ppi <- function(dvi_ts, sza_vector) {
  
  # CRITICAL: Scale Sentinel-2 integers back to true float (0.0 to 1.0)
  dvi_ts_scaled <- dvi_ts / 10000
  
  # Mathematical safety clamps to remove implausible extreme values and NoData asymptotes:
  # 1. Neutralize the NoData Asymptote: Force DVI <= 0 to NA
  dvi_ts_scaled[dvi_ts_scaled <= 0] <- NA_real_
  
  # 2. Neutralize the M=1 Asymptote: Clamp extreme values to 0.95
  dvi_ts_scaled[dvi_ts_scaled > 0.95] <- 0.95
  
  n_dates   <- length(dvi_ts_scaled)
  valid_idx <- which(!is.na(dvi_ts_scaled))
  
  # If fewer than 2 valid points exist across the time series, return all NA
  if (length(valid_idx) < 2) {
    return(rep(NA_real_, n_dates))
  }
  
  dvi_valid <- dvi_ts_scaled[valid_idx]
  sza_valid <- sza_vector[valid_idx]
  
  # Run MitMat's PPI on the valid observations with tryCatch safety armor
  ppi_valid <- tryCatch({
    ppi::ppi(
      dvi          = dvi_valid,
      zenith.angle = sza_valid,
      G            = 0.2 # Lower Leaf Angle Distribution suited to Knepp's grassland/scrub mosaic
    )
  }, error = function(e) {
    return(rep(NA_real_, length(dvi_valid)))
  })
  
  # Reconstruct full-length output preserving temporal structure
  ppi_full            <- rep(NA_real_, n_dates)
  ppi_full[valid_idx] <- ppi_valid
  
  return(ppi_full)
}

## 4. Apply the PPI function pixel-by-pixel across the 2019-2024 stack
message("Calculating PPI timeseries across all pixels via MitMat's ppi package...")

Knepp_Buffered_Timeseries.PPI <- app(
  Knepp_Buffered_Timeseries.DVI,
  fun        = calc_pixel_ppi,
  sza_vector = sza_rad_vector,
  cores      = knepp.cores
)

# Clamp final PPI to ecologically valid bounds to prevent logarithmic asymptotes skewing TWDTW
Knepp_Buffered_Timeseries.PPI[Knepp_Buffered_Timeseries.PPI < -2] <- NA_real_
Knepp_Buffered_Timeseries.PPI[Knepp_Buffered_Timeseries.PPI > 10] <- NA_real_

# Assign clean monthly layer names (e.g., "2019-01_PPI")
names(Knepp_Buffered_Timeseries.PPI) <- paste0(sub("_B08$", "", names(Knepp_NIR_Timeseries)), "_PPI")

message("PPI calculation complete! Exporting unsmoothed continuous 2019-2024 PPI stack...")

writeRaster(
  Knepp_Buffered_Timeseries.PPI,
  filename  = file.path(Knepp_Processed, "Knepp_Buffered_PPI_Continuous_Raw_Timeseries.tif"),
  overwrite = TRUE
)

Knepp_Buffered_Timeseries.PPI <- rast(
  file.path(Knepp_Processed, "Knepp_Buffered_PPI_Continuous_Raw_Timeseries.tif")
) # And load it back in as necessary!

### Savitzky-Golay Gap-Filling and Smoothing ####
## Define the pixel-wise gap-filling function (identical to Kilimanjaro & Knepp NDVI)

sg_gapfill <- function(x) {
  
  # If a pixel is NA across all layers (outside the 5 km boundary or pure water), return NAs immediately
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  
  # 1. Linear interpolation across cloud gaps (rule = 2 pads leading/trailing NAs with nearest valid value)
  x_interp <- zoo::na.approx(x, na.rm = FALSE, rule = 2)
  
  # 2. Savitzky-Golay smoothing (fl = 5 preserves seasonal amplitude; forder = 2 is standard quadratic)
  x_smoothed <- pracma::savgol(
    x_interp,
    fl     = 5,
    forder = 2
  )
  
  return(x_smoothed)
}

message("Gap-filling and smoothing the continuous multi-year PPI dataset with Savitzky-Golay...")

Knepp_Buffered_Timeseries.PPI.Cleaned <- app(
  Knepp_Buffered_Timeseries.PPI,
  fun   = sg_gapfill,
  cores = knepp.cores
)

# Restore layer names stripped by terra::app()
names(Knepp_Buffered_Timeseries.PPI.Cleaned) <- names(Knepp_Buffered_Timeseries.PPI)

# Re-mask against the buffered boundary to ensure background pixels remain strict NAs
if (crs(KneppEstate_Buffered_Boundaries) != crs(Knepp_Buffered_Timeseries.PPI.Cleaned)) {
  KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(Knepp_Buffered_Timeseries.PPI.Cleaned))
}
Knepp_Buffered_Timeseries.PPI.Cleaned <- mask(Knepp_Buffered_Timeseries.PPI.Cleaned, KneppEstate_Buffered_Boundaries)

message("Gap-filling complete! Writing cleaned multi-year PPI stack to disk...")

writeRaster(
  Knepp_Buffered_Timeseries.PPI.Cleaned,
  filename  = file.path(Knepp_Processed, "Knepp_Buffered_PPI_Cleaned_SG-method.tif"),
  overwrite = TRUE
)

# Reload from disk so the SpatRaster points cleanly to the file
Knepp_Buffered_Timeseries.PPI.Cleaned <- rast(
  file.path(Knepp_Processed, "Knepp_Buffered_PPI_Cleaned_SG-method.tif")
) # Load it back in!

print(Knepp_Buffered_Timeseries.PPI.Cleaned)

### Calculate and plot the mean seasonal trajectory across the site ####
# 1. Crop and mask to the inner Knepp Estate boundary

Knepp_Timeseries.PPI <- crop(
  Knepp_Buffered_Timeseries.PPI.Cleaned,
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.PPI.Cleaned)),
  mask = TRUE
)

Knepp_Timeseries_Mean.PPI <- global(Knepp_Timeseries.PPI, fun = "mean", na.rm = TRUE)

# 2. Mask out the inside of the Knepp Estate (inverse = TRUE) to isolate the surrounding 5 km farmland

Knepp_Environs_Timeseries.PPI <- mask(
  Knepp_Buffered_Timeseries.PPI.Cleaned,
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.PPI.Cleaned)),
  inverse = TRUE
)

Knepp_Environs_Timeseries.PPI <- global(Knepp_Environs_Timeseries.PPI, fun = "mean", na.rm = TRUE)

# 3. Export the 4K comparison plot of the continuous multi-year PPI trajectories

png(file.path(Knepp_Results, "Knepp_Estate_PPI_Mean_Timeseries.png"),
    width = 3840, height = 2160, res = 150)

plot(Knepp.full.dates, Knepp_Timeseries_Mean.PPI[, 1],
     type = "l",
     lwd  = 2,
     col  = "forestgreen",
     ylim = range(c(Knepp_Timeseries_Mean.PPI[, 1], Knepp_Environs_Timeseries.PPI[, 1]), na.rm = TRUE),
     xlab = "Date",
     ylab = "Mean PPI",
     main = "Knepp Estate vs. Surrounding Farmland - Continuous Mean PPI",
     xaxt = "n")

lines(Knepp.full.dates, Knepp_Environs_Timeseries.PPI[, 1],
      lwd = 2,
      col = "chocolate")

grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# Minor monthly ticks and major yearly ticks
axis.Date(1, at = Knepp.full.dates, labels = FALSE, tcl = -0.25)
axis.Date(1,
          at     = as.Date(paste0(unique(format(Knepp.full.dates, "%Y")), "-01-15")),
          format = "%Y",
          tcl    = -0.5)

legend("topright",
       legend = c("Knepp Estate", "Environs surrounding Knepp"),
       col    = c("forestgreen", "chocolate"),
       lwd    = 2,
       bty    = "n")

dev.off()

### 1. Shannon-Wiener Index ####

message("Calculating Shannon-Wiener diversity index across the 2019-2024 SG-filtered PPI stack...")

# Collapse the 72-layer time series into a single multi-year mean PPI raster

Knepp_Mean_PPI_SG <- app(Knepp_Buffered_Timeseries.PPI.Cleaned, fun = mean, na.rm = TRUE)
writeRaster(Knepp_Mean_PPI_SG, file.path(Knepp_Processed, "Knepp_2019-2024_MeanPPI_SG.tif"), overwrite = TRUE)
Knepp_Mean_PPI_SG <- rast(file.path(Knepp_Processed, "Knepp_2019-2024_MeanPPI_SG.tif")) # Load back in from disk

# Simplify the raster to 2 decimal places to prevent numerical oversaturation in Shannon's H
# (PPI is already a float, so no / 10000 division is needed)

Knepp_Mean_PPI_SG2dec <- trim(round(Knepp_Mean_PPI_SG, 2))

# Convert the SpatRaster to a wide matrix (required by ShannonS) and compute Shannon's H (3x3 window)

Knepp.PPI.ShannonH.SG.matrix <- rasterdiv::ShannonS(
  x            = terra::as.matrix(Knepp_Mean_PPI_SG2dec, wide = TRUE),
  window       = 3,
  na.tolerance = 0
)

# Convert the resulting matrix back into a SpatRaster and restore spatial metadata

Knepp_PPI_ShannonH_SG        <- rast(Knepp.PPI.ShannonH.SG.matrix)
ext(Knepp_PPI_ShannonH_SG)   <- ext(Knepp_Mean_PPI_SG2dec)
crs(Knepp_PPI_ShannonH_SG)   <- crs(Knepp_Mean_PPI_SG2dec)
names(Knepp_PPI_ShannonH_SG) <- "ShannonH_SG"

# Extend back to the untrimmed Knepp_Mean_PPI_SG template so extents match 100% when stacking later

Knepp_PPI_ShannonH_SG <- terra::extend(Knepp_PPI_ShannonH_SG, Knepp_Mean_PPI_SG)

# Export the final Shannon's H raster to Results and reload

writeRaster(
  Knepp_PPI_ShannonH_SG,
  filename  = file.path(Knepp_Results, "Knepp_2019-2024_PPI_ShannonH_SG.tif"),
  overwrite = TRUE
)
Knepp_PPI_ShannonH_SG <- rast(
  file.path(Knepp_Results, "Knepp_2019-2024_PPI_ShannonH_SG.tif")
) # Load it back in!

message("PPI Shannon's H calculation complete!")

### 2. Classic Rao's Q ####
## Due to the size of the 10m Sentinel-2 site, we tile the dataset into 750 tiles
## for parallel processing on the MaRC3a HPC cluster, mirroring the NDVI pipeline.

message("Preparing 750-tile grid for Knepp Estate PPI Classic & TWDTW Rao's Q...")

### Step 1: Create or load the master 750-tile grid

trimmed.Knepp_Mean_PPI_Raster <- trim(Knepp_Mean_PPI_SG)

knepp.tiles.filepath <- file.path(Knepp_Processed, "Knepp_Tiling_Grid_Polygons.geojson")

if (file.exists(knepp.tiles.filepath)) {
  message("Loading existing 750-tile grid from Knepp_Tiling_Grid_Polygons.geojson...")
  knepp.tiling.grid <- vect(knepp.tiles.filepath)
} else {
  knepp.total.tiles  <- 750
  knepp.aspect.ratio <- ncol(trimmed.Knepp_Mean_PPI_Raster) / nrow(trimmed.Knepp_Mean_PPI_Raster)
  
  knepp.tiling.factors <- expand.grid(
    ncols = 1:knepp.total.tiles,
    nrows = 1:knepp.total.tiles
  )
  
  knepp.tiling.factors <- knepp.tiling.factors[
    knepp.tiling.factors$ncols * knepp.tiling.factors$nrows == knepp.total.tiles, 
  ]
  
  knepp.tiling.factors$ratio_diff <- abs(
    (knepp.tiling.factors$ncols / knepp.tiling.factors$nrows) - knepp.aspect.ratio
  )
  
  knepp.best.tile.size <- knepp.tiling.factors[which.min(knepp.tiling.factors$ratio_diff), ]
  
  knepp.cols <- knepp.best.tile.size$ncols
  knepp.rows <- knepp.best.tile.size$nrows
  
  knepp.tiling.grid <- as.polygons(
    rast(
      ext(trimmed.Knepp_Mean_PPI_Raster),
      ncols = knepp.cols,
      nrows = knepp.rows,
      crs   = crs(trimmed.Knepp_Mean_PPI_Raster)
    )
  )
  
  writeVector(knepp.tiling.grid, knepp.tiles.filepath, filetype = "GeoJSON", overwrite = TRUE)
  knepp.tiling.grid <- vect(knepp.tiles.filepath)
}

plot(knepp.tiling.grid, main = "Knepp Estate 750-Tile Grid (2019-2024)")

# Define moving window size and 1-pixel overlap buffer so Rao's Q has no tile-edge NAs

RaoQ.window.size   <- 3
knepp.tile.overlap <- floor(RaoQ.window.size / 2)

# Create directory for Classic Rao's Q input tiles (Mean PPI)

knepp.ppi.tile.dir.sg <- file.path(Knepp_Processed, "Knepp_Tiles_PPI_SG-Filtered")
dir.create(knepp.ppi.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

message("Generating 750 Mean PPI tiles for Classic Rao's Q...")

makeTiles(
  trimmed.Knepp_Mean_PPI_Raster,
  y         = knepp.tiling.grid,
  buffer    = knepp.tile.overlap,
  filename  = file.path(knepp.ppi.tile.dir.sg, "Knepp_MeanPPI_SG_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Submit the Classic Rao's Q tiles for processing on MaRC3a
## Please see script "06.1A_Knepp_PPI_Classic-RaoQ_SG-Filtered.R" and its SLURM wrapper
## for the HPC batch job execution. Run Step 3 below once all 750 tiles have finished on MaRC3a.

### Step 3: Demosaic the classical Rao's Q tiles

message("Demosaicing PPI Classic Rao's Q tiles...")

knepp.ppi.rao.files.sg <- list.files(
  file.path(Knepp_Processed, "Knepp_Tiles_PPI_SG-Filtered", "Rao-utputs"),
  pattern    = "Knepp_PPI_Classic-RaoQ_SG_Tile-",
  full.names = TRUE
)

if (length(knepp.ppi.rao.files.sg) == 0) {
  stop("No PPI Classic Rao's Q output tiles found! Ensure MaRC3a jobs have completed and synced.")
}

# Stitch all 750 tiles back together using a SpatRasterCollection (sprc)

Knepp_PPI_Classic_RaoQ_SG <- terra::mosaic(sprc(lapply(knepp.ppi.rao.files.sg, rast)))
Knepp_PPI_Classic_RaoQ_SG <- terra::extend(Knepp_PPI_Classic_RaoQ_SG, Knepp_Mean_PPI_SG)
names(Knepp_PPI_Classic_RaoQ_SG) <- "RaosQ_Classic_SG"

writeRaster(
  Knepp_PPI_Classic_RaoQ_SG,
  filename  = file.path(Knepp_Results, "Knepp_2019-2024_PPI_Classic-RaoQ_SG.tif"),
  overwrite = TRUE
)
Knepp_PPI_Classic_RaoQ_SG <- rast(file.path(Knepp_Results, "Knepp_2019-2024_PPI_Classic-RaoQ_SG.tif"))

plot(Knepp_PPI_Classic_RaoQ_SG, main = "Knepp Estate Classic Rao's Q (PPI 2019-2024)")

### 3. Rao's Q with TWDTW ####

message("Preparing 72-layer PPI timeseries tiles for TWDTW Rao's Q...")

### Step 1: Tile the continuous 2019-2024 PPI timeseries stack

knepp.ppi.twdtw.tile.dir.sg <- file.path(Knepp_Processed, "Knepp_TS_Tiles_PPI_SG-Filtered")
dir.create(knepp.ppi.twdtw.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

message("Generating 750 Timeseries PPI tiles (2019-2024) for TWDTW Rao's Q...")

makeTiles(
  trim(Knepp_Buffered_Timeseries.PPI.Cleaned),
  y         = knepp.tiling.grid,
  buffer    = knepp.tile.overlap,
  filename  = file.path(knepp.ppi.twdtw.tile.dir.sg, "Knepp_TS_PPI_SG_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Run the Alpha Grid Search (06.3) and submit TWDTW tiles on MaRC3a (06.2A)
## 1. Run "06.3_Knepp_PPI_α_Grid-search.R" to find the optimal TWDTW stepness (α).
## 2. Submit "06.2A_Knepp_PPI_TWDTW-RaoQ_SG-Filtered.R" via SLURM on MaRC3a.
## 3. Once the 750 TWDTW tiles are synced back, run Step 3 below.

### Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster

message("Demosaicing PPI TWDTW Rao's Q tiles...")

knepp.ppi.twdtw.rao.files.sg <- list.files(
  file.path(Knepp_Processed, "Knepp_TS_Tiles_PPI_SG-Filtered", "Rao-utputs"),
  pattern    = "Knepp_PPI_TWDTW-RaoQ_SG_Tile-",
  full.names = TRUE
)

if (length(knepp.ppi.twdtw.rao.files.sg) == 0) {
  stop("No PPI TWDTW Rao's Q output tiles found! Ensure MaRC3a jobs have completed and synced.")
}

Knepp_PPI_TWDTW_RaoQ_SG <- terra::mosaic(sprc(lapply(knepp.ppi.twdtw.rao.files.sg, rast)))
Knepp_PPI_TWDTW_RaoQ_SG <- terra::extend(Knepp_PPI_TWDTW_RaoQ_SG, Knepp_Mean_PPI_SG)
names(Knepp_PPI_TWDTW_RaoQ_SG) <- "RaosQ_TWDTW_SG"

writeRaster(
  Knepp_PPI_TWDTW_RaoQ_SG,
  filename  = file.path(Knepp_Results, "Knepp_2019-2024_PPI_TWDTW-RaoQ_SG.tif"),
  overwrite = TRUE
)
Knepp_PPI_TWDTW_RaoQ_SG <- rast(file.path(Knepp_Results, "Knepp_2019-2024_PPI_TWDTW-RaoQ_SG.tif"))

plot(Knepp_PPI_TWDTW_RaoQ_SG, main = "Knepp Estate TWDTW Rao's Q (PPI 2019-2024)")

### Export all rasters for comparison ####
## Stack the 2019-2024 Mean PPI and all three SG-filtered diversity indices

message("Stacking and exporting 2019-2024 Knepp PPI diversity rasters for comparison...")

Knepp_PPI_Comparison_Rasters <- c(
  Knepp_Mean_PPI_SG,
  Knepp_PPI_ShannonH_SG,
  Knepp_PPI_Classic_RaoQ_SG,
  Knepp_PPI_TWDTW_RaoQ_SG
)

names(Knepp_PPI_Comparison_Rasters) <- c(
  "Sentinel-2 MeanPPI SG Gap-Filled",
  "ShannonH SG (PPI)",
  "RaosQ Classic SG (PPI)",
  "RaosQ TWDTW SG (PPI)"
)

# Export the combined 4-layer comparison stack so it doesn't need to be rebuilt each session
writeRaster(
  Knepp_PPI_Comparison_Rasters,
  filename  = file.path(Knepp_Results, "Knepp_2019-2024_PPI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

Knepp_PPI_Comparison_Rasters <- rast(
  file.path(Knepp_Results, "Knepp_2019-2024_PPI_Diversity_Comparison.tif")
) # Load it back in!

# Export a 2x2 multipanel comparison plot for the manuscript
png(file.path(Knepp_Results, "Knepp_PPI-Based_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(Knepp_PPI_Comparison_Rasters)

dev.off()

plot(Knepp_PPI_Comparison_Rasters) # Plot in RStudio viewer to inspect visually

### Assess index performance using Strict Consensus UKCEH Land Cover (2019-2024) ####

message("Assessing 2019-2024 PPI diversity indices against Strict Consensus UKCEH Land Cover...")

# 1. Define the 21 UKCEH Land Cover Classes lookup dictionary

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

# 2. Load the pre-computed Strict Consensus UKCEH Land Cover raster (or build it if not yet saved)

consensus_lc_path <- file.path(Knepp_Processed, "Knepp_2019-2024_UKCEH_Consensus_LandCover.tif")

if (file.exists(consensus_lc_path)) {
  message("Loading pre-computed 2019-2024 UKCEH Consensus Land Cover raster from disk...")
  Knepp_LandCover_Consensus <- rast(consensus_lc_path)
  
  # Ensure exact grid alignment with Knepp_PPI_Comparison_Rasters
  if (!compareGeom(Knepp_PPI_Comparison_Rasters[[1]], Knepp_LandCover_Consensus, stopOnError = FALSE)) {
    Knepp_LandCover_Consensus <- project(Knepp_LandCover_Consensus, Knepp_PPI_Comparison_Rasters[[1]], method = "near")
  }
} else {
  message("Consensus raster not found on disk; building from raw 2019-2024 UKCEH files...")
  
  ceh_base_dir <- file.path(Knepp_Input, "CEH Land-cover data")
  if (!dir.exists(ceh_base_dir)) {
    ceh_base_dir <- file.path(Knepp_Input, "CEH_Land-cover_data")
  }
  
  ceh_tif_files <- list.files(ceh_base_dir, pattern = "^LCM.*\\.tif$", recursive = TRUE, full.names = TRUE)
  ceh_years     <- sub(".*(20\\d{2}).*", "\\1", ceh_tif_files)
  valid_ceh_idx <- which(ceh_years %in% as.character(2019:2024))
  
  ceh_tif_files <- ceh_tif_files[valid_ceh_idx]
  ceh_years     <- ceh_years[valid_ceh_idx]
  names(ceh_tif_files) <- ceh_years
  
  Knepp_LC_Aligned_List <- lapply(ceh_years, function(y) {
    message("  -> Projecting and snapping UKCEH Land Cover for year: ", y)
    lc_raw <- rast(ceh_tif_files[[y]])[[1]]
    project(lc_raw, Knepp_PPI_Comparison_Rasters[[1]], method = "near")
  })
  
  Knepp_LC_Stack <- rast(Knepp_LC_Aligned_List)
  lc_min         <- min(Knepp_LC_Stack, na.rm = FALSE)
  lc_max         <- max(Knepp_LC_Stack, na.rm = FALSE)
  consensus_mask <- (lc_min == lc_max)
  
  Knepp_LandCover_Consensus <- mask(Knepp_LC_Stack[[1]], consensus_mask, maskvalues = 0)
  names(Knepp_LandCover_Consensus) <- "Veg_GroundTruth"
  
  writeRaster(Knepp_LandCover_Consensus, filename = consensus_lc_path, overwrite = TRUE)
  Knepp_LandCover_Consensus <- rast(consensus_lc_path)
}

# 3. Bind the PPI comparison rasters with the consensus ground truth and mask to the 5 km buffered boundary

if (crs(KneppEstate_Buffered_Boundaries) != crs(Knepp_PPI_Comparison_Rasters)) {
  KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(Knepp_PPI_Comparison_Rasters))
}

Knepp_PPI_Indices_Comparison_Raster <- c(
  Knepp_PPI_Comparison_Rasters,
  Knepp_LandCover_Consensus
)

Knepp_PPI_Indices_Comparison_Raster <- mask(
  crop(Knepp_PPI_Indices_Comparison_Raster, KneppEstate_Buffered_Boundaries),
  KneppEstate_Buffered_Boundaries
)

# Export the aligned 5-layer PPI stack for safekeeping

writeRaster(
  Knepp_PPI_Indices_Comparison_Raster,
  filename  = file.path(Knepp_Results, "Knepp_2019-2024_Sentinel_Aligned_Indices_PPI.tif"),
  overwrite = TRUE
)

Knepp_PPI_Indices_Comparison_Raster <- rast(
  file.path(Knepp_Results, "Knepp_2019-2024_Sentinel_Aligned_Indices_PPI.tif")
) # Load it back in for conversion to dataframe

### Convert the index rasters to a dataframe for performance analysis ####

message("Converting aligned PPI raster stack to dataframe...")

Knepp_PPI_Indices_Comparison_DF <- as.data.frame(
  Knepp_PPI_Indices_Comparison_Raster,
  na.rm = TRUE
)

colnames(Knepp_PPI_Indices_Comparison_DF) <- c(
  "Sentinel-2 MeanPPI SG Gap-Filled",
  "ShannonH SG (PPI)",
  "RaosQ Classic SG (PPI)",
  "RaosQ TWDTW SG (PPI)",
  "Veg_GroundTruth"
)

# Translate the numeric UKCEH class codes (1-21) to descriptive UKCEH habitat names

Knepp_PPI_Indices_Comparison_DF$Veg_GroundTruth <- UKCEH_LC_Classes[
  as.character(Knepp_PPI_Indices_Comparison_DF$Veg_GroundTruth)
]

### PERMANOVA ####
## Filter out 0 and NA edge artifacts and execute a stratified 17,500-pixel subsample
## across the 7 UKCEH land cover classes, matching 05_Analyse_Knepp_Estate_Sentinel_NDVI.R.

# 1. Strip out NA and 0 value background/boundary pixels across all index columns

Knepp_PPI_Filtered <- Knepp_PPI_Indices_Comparison_DF %>%
  filter(!is.na(Veg_GroundTruth)) %>%
  filter(if_all(-Veg_GroundTruth, ~ . != 0 & !is.na(.)))

# 2. Specify the 7 ecologically significant UKCEH classes present at the Knepp Estate

target_classes <- c(
  "Broadleaved woodland",
  "Coniferous woodland",
  "Arable and horticulture",
  "Improved grassland",
  "Freshwater",
  "Urban",
  "Suburban"
)

Knepp_PPI_Filtered <- Knepp_PPI_Filtered %>%
  filter(Veg_GroundTruth %in% target_classes)

present_classes <- unique(Knepp_PPI_Filtered$Veg_GroundTruth)

# 3. Calculate the equal split for the 17,500 pixel target (2,500 pixels per class across 7 classes)

total_target_pixels <- 17500
n_classes           <- length(present_classes)
pixels_per_class    <- floor(total_target_pixels / n_classes)

message("Sampling up to ", pixels_per_class, " consensus pixels across ", n_classes, " UKCEH land cover classes.")

# 4. Execute the Stratified Subsample

set.seed(123)

subset.Knepp_PPI_Indices_Comparison_DF <- Knepp_PPI_Filtered %>%
  group_by(Veg_GroundTruth) %>%
  slice_sample(n = pixels_per_class, replace = FALSE) %>%
  ungroup() %>%
  as.data.frame()

# 5. Save the filtered, stratified dataframe for reproducibility and standalone server execution (06.4)

saveRDS(
  subset.Knepp_PPI_Indices_Comparison_DF,
  file = file.path(Knepp_Processed, "Knepp_PPI_Comparison_DF_Filtered.rds")
)

subset.Knepp_PPI_Indices_Comparison_DF <- readRDS(
  file.path(Knepp_Processed, "Knepp_PPI_Comparison_DF_Filtered.rds")
)

### Conduct a series of PERMANOVAs ####
## Create the PERMANOVAs subfolder inside the Knepp results directory

Knepp_PERMANOVA_Dir <- file.path(Knepp_Results, "PERMANOVAs")
dir.create(Knepp_PERMANOVA_Dir, showWarnings = FALSE, recursive = TRUE)

# Create an explicit PSOCK cluster (capped at 8 cores) to prevent mclapply SIGPIPE OOM crashes
perm_cluster <- makeCluster(permanova.cores)

message(paste0(Sys.time(), ": Starting SG-Filtered Knepp PPI raster PERMANOVAs..."))

# 1. SG-Filtered Shannon's H (PPI)

message(paste0(Sys.time(), ": [1/3] Running PERMANOVA_Knepp_PPI_ShannonsH_SG..."))
PERMANOVA_Knepp_PPI_ShannonsH_SG <- adonis2(
  subset.Knepp_PPI_Indices_Comparison_DF$`ShannonH SG (PPI)` ~ subset.Knepp_PPI_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = perm_cluster
)
saveRDS(PERMANOVA_Knepp_PPI_ShannonsH_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_PPI_ShannonsH_SG.rds"))
message(paste0(Sys.time(), ": [1/3] Completed and saved PERMANOVA_Knepp_PPI_ShannonsH_SG.rds"))
gc()

# 2. SG-Filtered Classic Rao's Q (PPI)

message(paste0(Sys.time(), ": [2/3] Running PERMANOVA_Knepp_PPI_RaosQ_Classic_SG..."))
PERMANOVA_Knepp_PPI_RaosQ_Classic_SG <- adonis2(
  subset.Knepp_PPI_Indices_Comparison_DF$`RaosQ Classic SG (PPI)` ~ subset.Knepp_PPI_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = perm_cluster
)
saveRDS(PERMANOVA_Knepp_PPI_RaosQ_Classic_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_PPI_RaosQ_Classic_SG.rds"))
message(paste0(Sys.time(), ": [2/3] Completed and saved PERMANOVA_Knepp_PPI_RaosQ_Classic_SG.rds"))
gc()

# 3. SG-Filtered TWDTW Rao's Q (PPI)

message(paste0(Sys.time(), ": [3/3] Running PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG..."))
PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG <- adonis2(
  subset.Knepp_PPI_Indices_Comparison_DF$`RaosQ TWDTW SG (PPI)` ~ subset.Knepp_PPI_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = perm_cluster
)
saveRDS(PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG.rds"))
message(paste0(Sys.time(), ": [3/3] Completed and saved PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG.rds"))

stopCluster(perm_cluster)
gc()

### Compile and export final PERMANOVA Summary Table ####

message(paste0(Sys.time(), ": Compiling Knepp PPI PERMANOVA results..."))

Knepp_PPI_PERMANOVA_Results <- data.frame(
  Pipeline = c("SG Filtered", "SG Filtered", "SG Filtered"),
  Index    = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2       = c(
    PERMANOVA_Knepp_PPI_ShannonsH_SG$R2[1],
    PERMANOVA_Knepp_PPI_RaosQ_Classic_SG$R2[1],
    PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG$R2[1]
  ),
  F        = c(
    PERMANOVA_Knepp_PPI_ShannonsH_SG$F[1],
    PERMANOVA_Knepp_PPI_RaosQ_Classic_SG$F[1],
    PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG$F[1]
  ),
  p_value  = c(
    PERMANOVA_Knepp_PPI_ShannonsH_SG$`Pr(>F)`[1],
    PERMANOVA_Knepp_PPI_RaosQ_Classic_SG$`Pr(>F)`[1],
    PERMANOVA_Knepp_PPI_RaosQ_TWDTW_SG$`Pr(>F)`[1]
  )
)

write.csv(
  Knepp_PPI_PERMANOVA_Results,
  file      = file.path(Knepp_Results, "Knepp_PPI_PERMANOVA_Summary.csv"),
  row.names = FALSE
)
saveRDS(
  Knepp_PPI_PERMANOVA_Results,
  file = file.path(Knepp_Results, "Knepp_PERMANOVA_Summary_PPI.rds")
)

print(Knepp_PPI_PERMANOVA_Results)
message("Knepp Estate 2019-2024 PPI analysis complete!") # End of script ####