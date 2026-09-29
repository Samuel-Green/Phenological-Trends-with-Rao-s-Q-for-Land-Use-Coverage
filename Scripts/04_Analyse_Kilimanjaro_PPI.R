############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 26/06/2026                    #
# 04_Analyse_Kilimanjaro_PPI.R                                                 #
# Conducting comparative analysis of TWDTW Rao's Q and classic Rao's Q in Kili #
############################################################################ ###

### Install and load the necessary packages ####
## This should already be done from the setup file
# rasterdiv now contains the TWDTW-enabled paRao()

library(rasterdiv)
library(twdtw)
library(vegan)
library(pROC)
library(terra)

# Core spatial stack already loaded in 00_setup.R:
# terra, sf, here, dplyr, stringr

### Define the  file paths and import the site data ####
# Output directory for this script

KiliNP_Results <- file.path(Results, "Kilimanjaro")
dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

kili.cores <- max(1, detectCores() - 2)

# Explicitly set the input directory to the external drive to manage storage

KiliNP_PPI_Input <- "D:/Elliot Shayle/Kilimanjaro Geodata/Surface Reflectance Rasters" # Temporary location because my computer is low on storage
KiliNP_PPI_Processed <- "D:/Elliot Shayle/Kilimanjaro Geodata/Processed" # Temporary location because my computer is low on storage

# Load KiliNP_LandCover_Vector boundary (this is our land cover ground truth data)

KiliNP_LandCover_Vector <- 
  vect(file.path(KiliNP_Input,
                 "/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp")) # Load in the ground truth data

## List raster files
# The pattern now looks for the Sentinel-2 BOA files specifically

tmp.files <- list.files(KiliNP_PPI_Input, pattern = "SEN2_L3_BOA_.*_FBM\\.tif$", full.names = TRUE)

for (f in tmp.files) {
  
  tmp.KiliNP.raster <- rast(f)
  
  # Ensure CRS matches
  # It is best practice to project the vector to the raster's CRS to avoid altering raster values
  
  if (crs(tmp.KiliNP.raster) != crs(KiliNP_LandCover_Vector)) {
    KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(tmp.KiliNP.raster))}
    
  # Crop to bounding box
  
  tmp.KiliNP.raster.cropped <- crop(tmp.KiliNP.raster, KiliNP_LandCover_Vector)
  
  # Mask to exact boundary
  
  tmp.KiliNP.raster.masked <- mask(tmp.KiliNP.raster.cropped, KiliNP_LandCover_Vector)
  
  
  # Extract the year and the band from the filename 
  # e.g., KM_EPSG32737_SEN2_L3_BOA_B02_2017_FBM.tif
  
  tmp.year <- sub(".*_(\\d{4})_.*", "\\1", basename(f))
  tmp.band <- sub(".*_(B\\d{2})_.*", "\\1", basename(f)) # Extracts B02, B03, B04, B08, etc.
  
  # Create new filename incorporating both year and band to prevent overwriting
  
  tmp.new.name <- paste0("KiliNP_", tmp.year, "_", tmp.band, "_Cropped.tif")
  
  # Write to processed folder
  
  writeRaster(
    tmp.KiliNP.raster.masked,
    filename = file.path(KiliNP_PPI_Processed, tmp.new.name),
    overwrite = TRUE
  )
  message(paste0("Writing raster: ", tmp.new.name))
  
  # To keep memory usage under control, remove the temporary raster objects
  
  rm(tmp.KiliNP.raster, tmp.KiliNP.raster.cropped, tmp.KiliNP.raster.masked)
  gc() # Explicit garbage collection to free up RAM before the next heavy file
}

## Import the cropped multi-band rasters
# Get the files's names and locations for each reflectance band

# 1. Blue Band (B02)

KiliNP_B02_Files <- list.files(
  KiliNP_PPI_Processed,
  pattern = "^KiliNP_\\d{4}_B02_Cropped\\.tif$",
  full.names = TRUE
)
KiliNP_c_Timeseries <- rast(KiliNP_B02_Files)

# 2. Green Band (B03)

KiliNP_B03_Files <- list.files(
  KiliNP_PPI_Processed,
  pattern = "^KiliNP_\\d{4}_B03_Cropped\\.tif$",
  full.names = TRUE
)
KiliNP_Green_Timeseries <- rast(KiliNP_B03_Files)

# 3. Red Band (B04)

KiliNP_B04_Files <- list.files(
  KiliNP_PPI_Processed,
  pattern = "^KiliNP_\\d{4}_B04_Cropped\\.tif$",
  full.names = TRUE
)
KiliNP_Red_Timeseries <- rast(KiliNP_B04_Files)

# 4. Near-Infrared Band (B08)

KiliNP_B08_Files <- list.files(
  KiliNP_PPI_Processed,
  pattern = "^KiliNP_\\d{4}_B08_Cropped\\.tif$",
  full.names = TRUE
)
KiliNP_NIR_Timeseries <- rast(KiliNP_B08_Files)

## Rename the layers to something a bit more readable
# This regex safely skips over the band designation (e.g., _B02_) to grab the year

Kili.years <- sub(".*_(\\d{4})_B\\d{2}_Cropped\\.tif", "\\1", basename(sources(KiliNP_Blue_Timeseries)))

# Create new layer names

Kili.layer.names <- unlist(lapply(Kili.years, function(y) {
  paste0(y, " - ", month.name)
}))

# Assign the generated names to all four surface reflectance timeseries objects

names(KiliNP_Blue_Timeseries) <- Kili.layer.names
names(KiliNP_Green_Timeseries) <- Kili.layer.names
names(KiliNP_Red_Timeseries) <- Kili.layer.names
names(KiliNP_NIR_Timeseries) <- Kili.layer.names

### Generate the PPI raster stack ####

message("Calculating theoretical SZA vector for uncleaned rasters...")

### 1. Calculate Astronomical Solar Zenith Angles
# Extract the centroid coordinates of the Kilimanjaro bounding box using the uncleaned raster extent

kili_ext <- ext(KiliNP_Blue_Timeseries)
centroid_geom <- vect(
  matrix(c(mean(c(kili_ext$xmin, kili_ext$xmax)), 
           mean(c(kili_ext$ymin, kili_ext$ymax))), ncol=2), 
  crs = crs(KiliNP_Blue_Timeseries)
)

# Project the centroid to WGS84 to get the latitude in degrees

centroid_latlon <- project(centroid_geom, "EPSG:4326")
kili_lat_deg <- geom(centroid_latlon)[,"y"]
kili_lat_rad <- kili_lat_deg * (pi / 180) # Convert to radians

# Create a sequence of dates representing the middle of each composite month

dates <- seq(as.Date("2017-01-15"), as.Date("2021-12-15"), by = "1 month")
day_of_year <- as.numeric(strftime(dates, format = "%j"))

# Calculate solar declination angle (delta) in radians

declination_deg <- 23.45 * sin((2 * pi / 365) * (day_of_year - 81))
declination_rad <- declination_deg * (pi / 180)

# Calculate hour angle (h) in radians
# Sentinel-2 descends at approx 10:30 AM local solar time. 
# Hour angle = 15 degrees * (Hours from solar noon). 10.5 - 12.0 = -1.5 hours.

hour_angle_deg <- -1.5 * 15
hour_angle_rad <- hour_angle_deg * (pi / 180)

# Calculate Solar Zenith Angle (theta_z) in radians

cos_theta_z <- sin(kili_lat_rad) * sin(declination_rad) + 
  cos(kili_lat_rad) * cos(declination_rad) * cos(hour_angle_rad)
sza_rad_vector <- acos(cos_theta_z)

### 2. Calculate DVI

message("Calculating uncleaned Difference Vegetation Index (DVI)...")

# DVI is strictly NIR - Red, using the raw, uncleaned stacks

KiliNP_DVI <- KiliNP_NIR_Timeseries - KiliNP_Red_Timeseries

### 3. Calculate PPI via terra::app()

message("Applying Plant Phenology Index (PPI) formula across uncleaned time series...")

# Define the PPI maths as a custom function to gracefully handle the ~5% NA gaps

calc_ppi <- function(dvi_vals, sza_vector) {
  
  # If the pixel is fully masked (e.g., an entirely cloud-covered pixel across all 5 years), return NAs
  
  if(all(is.na(dvi_vals))) return(rep(NA_real_, length(dvi_vals)))
  
  # Calculate M: Canopy maximum of DVI plus a small constant (0.005)
  
  M <- max(dvi_vals, na.rm = TRUE) + 0.005
  
  # Initialise the output vector
  
  ppi_out <- rep(NA_real_, length(dvi_vals))
  valid <- !is.na(dvi_vals)
  
  if(any(valid)) {
    v_dvi <- dvi_vals[valid]
    v_sza <- sza_vector[valid]
    
    # Radiative transfer equations from Jin & Eklundh (2014)
    # Assuming G = 0.5 (spherical leaf angle distribution)
    
    d_c <- 0.0336 + 0.0477 / cos(v_sza)
    Q_E <- d_c + (1 - d_c) * 0.5 / cos(v_sza)
    K <- 1 / (4 * Q_E) * (1 + M) / (1 - M)
    
    # Logarithmic transformation (assuming bare soil DVI = 0.09)
    
    log_arg <- (M - v_dvi) / (M - 0.09)
    
    # Prevent negative values inside the logarithm (can occur with extreme data anomalies)
    
    log_arg[log_arg <= 0] <- NA_real_
    
    ppi_out[valid] <- -K * log(log_arg)
  }
  
  return(ppi_out)
}

# Apply the function across the z-dimension

KiliNP_PPI_Timeseries <- app(KiliNP_DVI, fun = calc_ppi, sza_vector = sza_rad_vector)

# Carry over the layer names for consistency

names(KiliNP_PPI_Timeseries) <- names(KiliNP_DVI)

# Export the raw PPI stack to your external drive to manage local storage

writeRaster(
  KiliNP_PPI_Timeseries, 
  filename = file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries.tif"), 
  overwrite = TRUE
)

# Load the PPI raster back in!

KiliNP_PPI_Timeseries <- rast(file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries.tif"))

### Parallel Track: Masked (Non-Gap-Filled) Pipeline ####

### Mask pixels in the raster stack which don't have a complete timeseries of data
message("Creating parallel non-gap-filled (Masked) pipeline...")

# Create logical mask: TRUE only where ALL layers are non-NA
Kili.pixel.mask <- app(KiliNP_PPI_Timeseries, function(x) all(!is.na(x)))

# Mask out incomplete pixels (FALSE becomes NA)
KiliNP_PPI_Timeseries_Masked <- mask(KiliNP_PPI_Timeseries, Kili.pixel.mask, maskvalues = 0)

# Export and load raster so I don't have to calculate it every time
writeRaster(KiliNP_PPI_Timeseries_Masked, file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries_Masked.tif"), overwrite = TRUE)
KiliNP_PPI_Timeseries_Masked <- rast(file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries_Masked.tif"))

### Run Savitzky-Golay filtering ####
## 1. Define the gap-filling function

sg_gapfill <- function(x) {
  
  # Check if the pixel is blank across all 279 layers. If TRUE, return NAs to save time.
  
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  
  ## Linear Interpolation
  # zoo::na.approx draws a straight line between the data points before and after gap 
  # 'rule = 2' is for if the 1st or last raster layers are NA, they are filled using the nearest valid observation
  
  x_interp <- zoo::na.approx(x, na.rm = FALSE, rule = 2) # No NA values prevents crashes
  
  ## Savitzky-Golay Smoothing
  
  x_smoothed <- pracma::savgol(x_interp, 
                               fl = 5, # fl = Filter length, must be an odd number, and `fl = 5` provides a seasonal smoothing window, preserving the seasonality
                               forder = 2) # forder: Filter order (polynomial degree), and 2 or 3 is standard for NDVI
  
  return(x_smoothed)
}

## 2. Check for missing months in the temporal sequence

message("Checking for missing months in the temporal sequence...")

current.names <- names(KiliNP_PPI_Timeseries)

# Adapt the date parsing to match the "2017 - January" format

current.dates <- as.Date(
  paste0(
    sub(" - .*", "", current.names), "-",
    match(sub(".* - ", "", current.names), month.name),
    "-15"
  ),
  format = "%Y-%m-%d"
) 

# Generate a perfectly continuous monthly sequence

Kili.full.dates <- seq(min(current.dates), max(current.dates), by = "month")
missing.dates <- Kili.full.dates[!Kili.full.dates %in% current.dates]

if (length(missing.dates) > 0) {
  message(paste("Found", length(missing.dates), "missing months. Creating blank template layers..."))
  
  Empty_Kili_Raster <- terra::init(KiliNP_PPI_Timeseries[[1]], NA)
  missing.rasters <- terra::rast(replicate(length(missing.dates), Empty_Kili_Raster))
  
  # Format missing names to match the "YYYY - Month" convention
  
  names(missing.rasters) <- paste0(format(missing.dates, "%Y"), " - ", months(missing.dates))
  
  KiliNP_PPI_Timeseries <- c(KiliNP_PPI_Timeseries, missing.rasters)
  
  # Sort the entire stack chronologically so NAs are in the correct sequence
  
  sorted.names <- paste0(format(Kili.full.dates, "%Y"), " - ", months(Kili.full.dates))
  KiliNP_PPI_Timeseries <- KiliNP_PPI_Timeseries[[sorted.names]]
  
} else {
  message("No missing months found. Temporal sequence is already contiguous.")
}

## 3. Apply the gap filling with Savitzky-Golay filter

message("Gap filling and smoothing the dataset...")

kili.cores <- max(1, detectCores() - 2)

KiliNP_PPI_Timeseries.SGfilter <- app(
  KiliNP_PPI_Timeseries, 
  fun = sg_gapfill, 
  cores = kili.cores 
)

# Carry over the names to the smoothed raster

names(KiliNP_PPI_Timeseries.SGfilter) <- names(KiliNP_PPI_Timeseries)

# 4. Export to a file (because it's not saved into the R data environment)

message("Exporting smoothed PPI timeseries...")
writeRaster(
  KiliNP_PPI_Timeseries.SGfilter,
  filename = file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries_SG-filtered.tif"),
  overwrite = TRUE
)

KiliNP_PPI_Timeseries.SGfilter <- rast( # Load it back in
  file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries_SG-filtered.tif"))

### Inspect temporal structure ####
## Unlike Macchia Sacra NetCDF, this GeoTIFF stack does not contain explicit time metadata
## Therefore, we must construct the time vector manually from layer names

message("Constructing time vector for Kilimanjaro time series...")

# Extract year and month from layer names
# Layer format: "2017 - January"

Kili.dates <- as.Date(
  paste0(
    sub(" - .*", "", names(KiliNP_PPI_Timeseries.SGfilter)), "-",
    match(sub(".* - ", "", names(KiliNP_PPI_Timeseries.SGfilter)), month.name),
    "-15" # As SZA used the 15th of the month, I will set the date here to be the 15th
  ),
  format = "%Y-%m-%d"
)

stopifnot(length(Kili.dates) == nlyr(KiliNP_PPI_Timeseries.SGfilter))

message(paste("Temporal length:", length(Kili.dates), "layers"))

### 1. Shannon-Wiener Index ####

message("Calculating Shannon-Wiener diversity index for Masked and SG tracks...")

## For the masked raster

KiliNP_Mean_PPI_Masked <- app(KiliNP_PPI_Timeseries_Masked, fun = mean, na.rm = TRUE)
writeRaster(KiliNP_Mean_PPI_Masked, file.path(KiliNP_PPI_Processed, "KiliNP_MeanPPI_Masked.tif"), overwrite = TRUE)
KiliNP_Mean_PPI_Masked <- rast(file.path(KiliNP_PPI_Processed, "KiliNP_MeanPPI_Masked.tif"))

KiliNP_Mean_PPI_Masked2dec <- trim(round(KiliNP_Mean_PPI_Masked, 2))
KiliNP.PPI.ShannonH.Masked.matrix <- rasterdiv::ShannonS(terra::as.matrix(KiliNP_Mean_PPI_Masked2dec, wide = TRUE), window = 3, na.tolerance = 0)

KiliNP_PPI_ShannonH_Masked <- rast(KiliNP.PPI.ShannonH.Masked.matrix)
ext(KiliNP_PPI_ShannonH_Masked) <- ext(KiliNP_Mean_PPI_Masked2dec)
crs(KiliNP_PPI_ShannonH_Masked) <- crs(KiliNP_Mean_PPI_Masked2dec)
names(KiliNP_PPI_ShannonH_Masked) <- "ShannonH_Masked"

writeRaster(KiliNP_PPI_ShannonH_Masked, file.path(KiliNP_Results, "Kilimanjaro_PPI_ShannonH_Masked.tif"), overwrite = TRUE)
KiliNP_PPI_ShannonH_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_ShannonH_Masked.tif"))

## For the gap-filled raster

KiliNP_Mean_PPI_SG <- app(KiliNP_PPI_Timeseries.SGfilter, fun = mean, na.rm = TRUE)
writeRaster(KiliNP_Mean_PPI_SG, file.path(KiliNP_PPI_Processed, "KiliNP_MeanPPI_SG.tif"), overwrite = TRUE)
KiliNP_Mean_PPI_SG <- rast(file.path(KiliNP_PPI_Processed, "KiliNP_MeanPPI_SG.tif"))

KiliNP_Mean_PPI_SG2dec <- trim(round(KiliNP_Mean_PPI_SG, 2))
KiliNP.PPI.ShannonH.SG.matrix <- rasterdiv::ShannonS(terra::as.matrix(KiliNP_Mean_PPI_SG2dec, wide = TRUE), window = 3, na.tolerance = 0)

KiliNP_PPI_ShannonH_SG <- rast(KiliNP.PPI.ShannonH.SG.matrix)
ext(KiliNP_PPI_ShannonH_SG) <- ext(KiliNP_Mean_PPI_SG2dec)
crs(KiliNP_PPI_ShannonH_SG) <- crs(KiliNP_Mean_PPI_SG2dec)
names(KiliNP_PPI_ShannonH_SG) <- "ShannonH_SG"

writeRaster(KiliNP_PPI_ShannonH_SG, file.path(KiliNP_Results, "Kilimanjaro_PPI_ShannonH_SG.tif"), overwrite = TRUE)
KiliNP_PPI_ShannonH_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_ShannonH_SG.tif"))

### 2. Classic Rao's Q  ####
## Due to the large size of the raster, I need to tile it so that it can be run
## The tiles will be stitched back together once they're computed

message("Calculating classical Rao's Q for Kilimanjaro...")

### Step 1: Create a grid to define zones for tiling

# Optional but recommended: trim outer NA borders

trimmed.KiliNP_Mean_PPI_Raster <- trim(KiliNP_Mean_PPI_Raster)

# For effective parallelisation on MaRC3a, 2000 tiny tiles is effective

kili.total.tiles <- 2000
kili.aspect.ratio <- ncol(trimmed.KiliNP_Mean_PPI_Raster) / nrow(trimmed.KiliNP_Mean_PPI_Raster)

# find factor pairs of 72 

tiling.factors <- expand.grid(
  ncols = 1:kili.total.tiles,
  nrows = 1:kili.total.tiles
)

# Subset to just pairs which equal "kili.total.tiles" when multiplied

tiling.factors <- tiling.factors[tiling.factors$ncols * tiling.factors$nrows == kili.total.tiles, ]

## choose the pair closest to the raster's aspect ratio
# First, create a new column calculating the difference in aspect ratio

tiling.factors$ratio_diff <- abs((tiling.factors$ncols / tiling.factors$nrows) - kili.aspect.ratio)

# Now find the pair which has the lowest distance from a square 1:1 aspect ratio

best.tile.size <- tiling.factors[which.min(tiling.factors$ratio_diff), ]

# And save that for later usage

kili.cols <- best.tile.size$ncols
kili.rows <- best.tile.size$nrows

# Create spatial polygons to specify what the tile sizes should be

kili.tiling.grid <- as.polygons(
  rast(
    ext(trimmed.KiliNP_Mean_PPI_Raster),
    ncols = kili.cols,
    nrows = kili.rows,
    crs = crs(trimmed.KiliNP_Mean_PPI_Raster)
  )
)

writeVector(kili.tiling.grid, file.path(KiliNP_PPI_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON"), filetype = "GeoJSON" , overwrite = TRUE) # Export for later use

kili.tiling.grid <- vect(file.path(KiliNP_PPI_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON")) # Load it back in

plot(kili.tiling.grid) # Plot it to make sure that it's loaded in

# Window size (I think 3 is the default, but this can be changed as necessary)

RaoQ.window.size <- 3
kili.tile.overlap <- floor(RaoQ.window.size / 2)

## Create directories and tiles
# For the masked raster

kili.tile.dir.masked <- file.path(KiliNP_PPI_Processed, "Kili_Tiles_PPI_Masked")
dir.create(kili.tile.dir.masked, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Mean_PPI_Masked), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.tile.dir.masked, "KiliNP_MeanPPI_Masked_Tile-.tif"),
  overwrite = FALSE
)

# For the gap-filled raster

kili.tile.dir.sg <- file.path(KiliNP_PPI_Processed, "Kili_Tiles_PPI_SG-Filtered")
dir.create(kili.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Mean_PPI_SG), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.tile.dir.sg, "KiliNP_MeanPPI_SG_Tile-.tif"),
  overwrite = FALSE
)

# ### Step 2: Compute classic Rao's Q for each tile
# ## Firstly, I need to setup the environment for parallelisation
# # Create a subfolder to store the classic Rao's Q output tiles
 
kili.rao.dir  <- file.path(KiliNP_Processed,"MeanPPI_Rao-utputs") 
dir.create(kili.rao.dir, recursive = TRUE, showWarnings = FALSE)
 
### Step 3: Demosaic the classical Rao's Q tiles
# For the masked raster

kili.rao.files.masked <- list.files(
  file.path(KiliNP_PPI_Processed, "Kili_Tiles_PPI_Masked"), 
  pattern = "KiliNP_MeanPPI_Masked_Tile-",
  full.names = TRUE
)
KiliNP_PPI_Classic_RaoQ_Masked <- terra::mosaic(sprc(lapply(kili.rao.files.masked, rast)))
writeRaster(KiliNP_PPI_Classic_RaoQ_Masked, file.path(KiliNP_Results, "Kilimanjaro_PPI_Classic-RaoQ_Masked.tif"), overwrite = TRUE)
KiliNP_PPI_Classic_RaoQ_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_Classic-RaoQ_Masked.tif"))

# For the gap-filled raster

kili.rao.files.sg <- list.files(
  file.path(KiliNP_PPI_Processed, "Kili_Tiles_PPI_SG-Filtered"), 
  pattern = "KiliNP_MeanPPI_SG_Tile-",
  full.names = TRUE
)
KiliNP_PPI_Classic_RaoQ_SG <- terra::mosaic(sprc(lapply(kili.rao.files.sg, rast)))
writeRaster(KiliNP_PPI_Classic_RaoQ_SG, file.path(KiliNP_Results, "Kilimanjaro_PPI_Classic-RaoQ_SG.tif"), overwrite = TRUE)
KiliNP_PPI_Classic_RaoQ_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_Classic-RaoQ_SG.tif"))

### 3. Rao's Q with TWDTW ####

message("Calculating Rao's Q with TWDTW distance for Kilimanjaro...")

## Step 1: I'll have to tile this as well because it is too large to compute as a single object
# This tiling script is copied from Step 1 of the classical Rao's Q analysis
# Some objects like "kili.tiling.grid" are assumed to be loaded, and "kili.tile.dir" is overwritten

# For the masked raster

kili.twdtw.tile.dir.masked <- file.path(KiliNP_PPI_Processed, "Kili_TS_Tiles_PPI_Masked")
dir.create(kili.twdtw.tile.dir.masked, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_PPI_Timeseries_Masked), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.twdtw.tile.dir.masked, "KiliNP_TS_PPI_Masked_Tile-.tif"),
  overwrite = FALSE
)

# For the gap-filled raster

kili.twdtw.tile.dir.sg <- file.path(KiliNP_PPI_Processed, "Kili_TS_Tiles_PPI_SG-Filtered")
dir.create(kili.twdtw.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_PPI_Timeseries.SGfilter), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.twdtw.tile.dir.sg, "KiliNP_TS_PPI_SG_Tile-.tif"),
  overwrite = FALSE
)

## Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster
# For the masked raster

kili.twdtw.rao.files.masked <- list.files(
  file.path(KiliNP_PPI_Processed, "Kili_TS_Tiles_PPI_Masked"), 
  pattern = "KiliNP_TS_PPI_Masked_Tile-",
  full.names = TRUE
)
KiliNP_PPI_TWDTW_RaoQ_Masked <- terra::mosaic(sprc(lapply(kili.twdtw.rao.files.masked, rast)))
writeRaster(KiliNP_PPI_TWDTW_RaoQ_Masked, file.path(KiliNP_Results, "Kilimanjaro_PPI_TWDTW-RaoQ_Masked.tif"), overwrite = TRUE)
KiliNP_PPI_TWDTW_RaoQ_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_TWDTW-RaoQ_Masked.tif"))

# For the gap-filled raster

kili.twdtw.rao.files.sg <- list.files(
  file.path(KiliNP_PPI_Processed, "Kili_TS_Tiles_PPI_SG-Filtered"), 
  pattern = "KiliNP_TS_PPI_SG_Tile-",
  full.names = TRUE
)
KiliNP_PPI_TWDTW_RaoQ_SG <- terra::mosaic(sprc(lapply(kili.twdtw.rao.files.sg, rast)))
writeRaster(KiliNP_PPI_TWDTW_RaoQ_SG, file.path(KiliNP_Results, "Kilimanjaro_PPI_TWDTW-RaoQ_SG.tif"), overwrite = TRUE)
KiliNP_PPI_TWDTW_RaoQ_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_PPI_TWDTW-RaoQ_SG.tif"))

### Export all rasters for comparison ####

KiliNP_PPI_Comparison_Rasters <- c(
  trim(KiliNP_Mean_PPI_Masked), 
  KiliNP_PPI_ShannonH_Masked,
  KiliNP_PPI_Classic_RaoQ_Masked,
  KiliNP_PPI_TWDTW_RaoQ_Masked,
  KiliNP_PPI_ShannonH_SG,
  KiliNP_PPI_Classic_RaoQ_SG,
  KiliNP_PPI_TWDTW_RaoQ_SG
)

names(KiliNP_PPI_Comparison_Rasters) <- c(
  "Sentinel-2_MeanPPI",
  "ShannonH_Masked",
  "RaosQ_Classic_Masked",
  "RaosQ_TWDTW_Masked",
  "ShannonH_SG",
  "RaosQ_Classic_SG",
  "RaosQ_TWDTW_SG"
)

writeRaster(
  KiliNP_PPI_Comparison_Rasters,
  filename = file.path(KiliNP_Results, "KiliNP_PPI_Diversity_Comparison.tif"),
  overwrite = TRUE
)
KiliNP_PPI_Comparison_Rasters <- rast(file.path(KiliNP_Results, "KiliNP_PPI_Diversity_Comparison.tif"))

png(file.path(KiliNP_Results, "KiliNP_PPI_Indices_Comparison.png"), width = 2560, height = 1440, res = 150)
plot(KiliNP_PPI_Comparison_Rasters)
dev.off()
plot(KiliNP_PPI_Comparison_Rasters)

### Assess index performance using vegetation ground truth ####

message("Assessing diversity indices against vegetation ground truth...")

KiliNP_LandCover_Vector <- vect(file.path(KiliNP_Input, "/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp"))
if (crs(KiliNP_PPI_Comparison_Rasters) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(KiliNP_PPI_Comparison_Rasters))
}

# Crop and mask the ENTIRE comparison stack simultaneously

masked.KiliNP_PPI_Comparison_Rasters <- mask(
  crop(KiliNP_PPI_Comparison_Rasters, KiliNP_LandCover_Vector), 
  KiliNP_LandCover_Vector
)

kili.land.cover.lookup <- c( 
  "0"  = NA, "1"  = "Anthropogenic", "2"  = "Dry Natural Vegetation (Savannah)", "3"  = NA, 
  "4"  = "Anthropogenic", "7"  = "Anthropogenic", "9"  = NA, "10" = "Moist/Evergreen Cloud Forest", 
  "11" = "Alpine/Subalpine Shrub & Grassland", "12" = "Moist/Evergreen Cloud Forest", 
  "13" = "Alpine/Subalpine Shrub & Grassland", "14" = "Anthropogenic", 
  "15" = "Alpine/Subalpine Shrub & Grassland", "16" = "Moist/Evergreen Cloud Forest", 
  "17" = NA, "18" = "Moist/Evergreen Cloud Forest", "19" = "Anthropogenic" 
) 

KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  KiliNP_PPI_Comparison_Rasters,
  field = "grid_code"
)

### Convert the index rasters to a dataframe for performance analysis ####

KiliNP_Indices_Comparison_Raster <- c(
  masked.KiliNP_PPI_Comparison_Rasters,
  KiliNP_LandCover_Raster
)

names(KiliNP_Indices_Comparison_Raster) <- c(
  names(masked.KiliNP_PPI_Comparison_Rasters),
  "Veg_GroundTruth_Numeric" 
)

KiliNP_Indices_Comparison_DF <- as.data.frame(
  KiliNP_Indices_Comparison_Raster,
  na.rm = TRUE
)

KiliNP_Indices_Comparison_DF.CombinedCoverClasses <- KiliNP_Indices_Comparison_DF
KiliNP_Indices_Comparison_DF$Veg_GroundTruth <- kili.land.cover.lookup[as.character(KiliNP_Indices_Comparison_DF$Veg_GroundTruth_Numeric)]
KiliNP_Indices_Comparison_DF.CombinedCoverClasses$Veg_GroundTruth <- kili.land.cover.lookup[as.character(KiliNP_Indices_Comparison_DF$Veg_GroundTruth_Numeric)]
KiliNP_Indices_Comparison_DF$Veg_GroundTruth_Numeric <- NULL

### Intra class variation exploration:

message("Calculating intra-class variance for original 19 classes...")
KiliNP_Stats_19Classes <- KiliNP_Indices_Comparison_DF %>%
  filter(!is.na(Veg_GroundTruth)) %>%
  group_by(Veg_GroundTruth) %>%
  summarise(across(c(ShannonH_Masked, RaosQ_Classic_Masked, RaosQ_TWDTW_Masked, ShannonH_SG, RaosQ_Classic_SG, RaosQ_TWDTW_SG), \(x) var(x, na.rm = TRUE)))
print(KiliNP_Stats_19Classes)

message("Calculating intra-class variance for combined classes...")
KiliNP_Stats_Combined <- KiliNP_Indices_Comparison_DF.CombinedCoverClasses %>%
  filter(!is.na(Veg_GroundTruth)) %>%
  group_by(Veg_GroundTruth) %>%
  summarise(across(c(ShannonH_Masked, RaosQ_Classic_Masked, RaosQ_TWDTW_Masked, ShannonH_SG, RaosQ_Classic_SG, RaosQ_TWDTW_SG), \(x) var(x, na.rm = TRUE)))
print(KiliNP_Stats_Combined)

message("Plotting box and whisker plots for combined classes...")
subset.KiliNP_Combined <- KiliNP_Indices_Comparison_DF.CombinedCoverClasses[sample(which(!is.na(KiliNP_Indices_Comparison_DF.CombinedCoverClasses$Veg_GroundTruth)), 10000), ]

KiliNP_Boxplot_Combined <- subset.KiliNP_Combined %>%
  pivot_longer(
    cols = c(ShannonH_Masked, RaosQ_Classic_Masked, RaosQ_TWDTW_Masked, ShannonH_SG, RaosQ_Classic_SG, RaosQ_TWDTW_SG), 
    names_to = "Index", 
    values_to = "Value"
  ) %>%
  ggplot(aes(x = Veg_GroundTruth, y = Value, fill = Veg_GroundTruth)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.2, outlier.size = 0.5) +
  facet_wrap(~Index, scales = "free_y") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" 
  ) +
  labs(
    title = "Intra-class Variance: Combined Land Cover Categories",
    x = "Land Cover Type",
    y = "Index Value"
  )
print(KiliNP_Boxplot_Combined)

### PERMANOVA ####

subset.KiliNP_Indices_Comparison_DF <- KiliNP_Indices_Comparison_DF[
  sample(which(!is.na(KiliNP_Indices_Comparison_DF$Veg_GroundTruth)), 10000), ]

# --- MASKED PERMANOVAs ---
PERMANOVA_ShannonsH_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
PERMANOVA_RaosQ_Classic_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
PERMANOVA_RaosQ_TWDTW_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

# --- SG FILTERED PERMANOVAs ---
PERMANOVA_ShannonsH_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
PERMANOVA_RaosQ_Classic_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
PERMANOVA_RaosQ_TWDTW_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

# Put the PERMANOVA results into a dataframe for effective presentation

KiliNP_PERMANOVA_Results <- data.frame(
  Pipeline = c("Masked", "Masked", "Masked", "SG Filtered", "SG Filtered", "SG Filtered"),
  Index = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q", "Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2 = c(
    PERMANOVA_ShannonsH_Masked$R2[1],
    PERMANOVA_RaosQ_Classic_Masked$R2[1],
    PERMANOVA_RaosQ_TWDTW_Masked$R2[1],
    PERMANOVA_ShannonsH_SG$R2[1],
    PERMANOVA_RaosQ_Classic_SG$R2[1],
    PERMANOVA_RaosQ_TWDTW_SG$R2[1]),
  F = c(
    PERMANOVA_ShannonsH_Masked$F[1],
    PERMANOVA_RaosQ_Classic_Masked$F[1],
    PERMANOVA_RaosQ_TWDTW_Masked$F[1],
    PERMANOVA_ShannonsH_SG$F[1],
    PERMANOVA_RaosQ_Classic_SG$F[1],
    PERMANOVA_RaosQ_TWDTW_SG$F[1]),
  p_value = c(
    PERMANOVA_ShannonsH_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_Masked$`Pr(>F)`[1],
    PERMANOVA_ShannonsH_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_SG$`Pr(>F)`[1])
)

print(KiliNP_PERMANOVA_Results)
message("Kilimanjaro PPI analyses complete!")

saveRDS(KiliNP_PERMANOVA_Results, file.path(KiliNP_Results, "KiliNP_PPI_PERMANOVA_Results.rds")) # Save it so I don't have to recalculate repeatedly
KiliNP_PERMANOVA_Results <- readRDS(file.path(KiliNP_Results, "KiliNP_PPI_PERMANOVA_Results.rds")) # And load it back in if necessary

message("Kilimanjaro analysis complete.")