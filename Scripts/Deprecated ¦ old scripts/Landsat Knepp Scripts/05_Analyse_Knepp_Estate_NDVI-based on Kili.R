############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 31/07/2026                    #
# 05_Analyse_Knepp_Estate_NDVI.R                                               #
# Comparative analysis of TWDTW Rao's Q and classic Rao's Q in the Knepp Estate#
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

Knepp_Results <- file.path(Results, "Knepp Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Explicitly set the input directory to the external drive to manage storage

Knepp_Buffered_NDVI_Input <- "D:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Input Data" # Temporary location because my computer is low on storage
Knepp_Buffered_NDVI_Processed <- "D:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Processed Data" # Temporary location because my computer is low on storage

# Load KiliNP_LandCover_Vector boundary (this is our land cover ground truth data)

Knepp_Buffered_LandCover_Vector <- 
  vect(file.path(KiliNP_Input,
                 "/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp")) # Load in the ground truth data

## List raster files
# The pattern now looks for the Sentinel-2 BOA files specifically

tmp.files <- list.files(Knepp_Buffered_NDVI_Input, pattern = "SEN2_L3_BOA_.*_FBM\\.tif$", full.names = TRUE)

for (f in tmp.files) {
  
  tmp.Knepp_Buffered.raster <- rast(f)
  
  # Ensure CRS matches
  # It is best practice to project the vector to the raster's CRS to avoid altering raster values
  
  if (crs(tmp.Knepp_Buffered.raster) != crs(Knepp_Buffered_LandCover_Vector)) {
    Knepp_Buffered_LandCover_Vector <- project(Knepp_Buffered_LandCover_Vector, crs(tmp.Knepp_Buffered.raster))}
    
  # Crop to bounding box
  
  tmp.Knepp_Buffered.raster.cropped <- crop(tmp.Knepp_Buffered.raster, Knepp_Buffered_LandCover_Vector)
  
  # Mask to exact boundary
  
  tmp.Knepp_Buffered.raster.masked <- mask(tmp.Knepp_Buffered.raster.cropped, Knepp_Buffered_LandCover_Vector)
  
  
  # Extract the year and the band from the filename 
  # e.g., KM_EPSG32737_SEN2_L3_BOA_B02_2017_FBM.tif
  
  tmp.year <- sub(".*_(\\d{4})_.*", "\\1", basename(f))
  tmp.band <- sub(".*_(B\\d{2})_.*", "\\1", basename(f)) # Extracts B02, B03, B04, B08, etc.
  
  # Create new filename incorporating both year and band to prevent overwriting
  
  tmp.new.name <- paste0("Knepp_Buffered_", tmp.year, "_", tmp.band, "_Cropped.tif")
  
  # Write to processed folder
  
  writeRaster(
    tmp.Knepp_Buffered.raster.masked,
    filename = file.path(Knepp_Buffered_NDVI_Processed, tmp.new.name),
    overwrite = TRUE
  )
  message(paste0("Writing raster: ", tmp.new.name))
  
  # To keep memory usage under control, remove the temporary raster objects
  
  rm(tmp.Knepp_Buffered.raster, tmp.Knepp_Buffered.raster.cropped, tmp.Knepp_Buffered.raster.masked)
  gc() # Explicit garbage collection to free up RAM before the next heavy file
}

# ## Import the cropped multi-band rasters
# # Get the files's names and locations for each reflectance band
# 
# # 1. Blue Band (B02)
# 
# Knepp_Buffered_B02_Files <- list.files(
#   Knepp_Buffered_NDVI_Processed,
#   pattern = "^KiliNP_\\d{4}_B02_Cropped\\.tif$",
#   full.names = TRUE
# )
# KiliNP_c_Timeseries <- rast(KiliNP_B02_Files)
# 
# # 2. Green Band (B03)
# 
# KiliNP_B03_Files <- list.files(
#   KiliNP_NDVI_Processed,
#   pattern = "^KiliNP_\\d{4}_B03_Cropped\\.tif$",
#   full.names = TRUE
# )
# KiliNP_Green_Timeseries <- rast(KiliNP_B03_Files)
# 
# # 3. Red Band (B04)
# 
# KiliNP_B04_Files <- list.files(
#   KiliNP_NDVI_Processed,
#   pattern = "^KiliNP_\\d{4}_B04_Cropped\\.tif$",
#   full.names = TRUE
# )
# KiliNP_Red_Timeseries <- rast(KiliNP_B04_Files)
# 
# # 4. Near-Infrared Band (B08)
# 
# KiliNP_B08_Files <- list.files(
#   KiliNP_NDVI_Processed,
#   pattern = "^KiliNP_\\d{4}_B08_Cropped\\.tif$",
#   full.names = TRUE
# )
# KiliNP_NIR_Timeseries <- rast(KiliNP_B08_Files)
# 
# ## Rename the layers to something a bit more readable
# # This regex safely skips over the band designation (e.g., _B02_) to grab the year
# 
# Kili.years <- sub(".*_(\\d{4})_B\\d{2}_Cropped\\.tif", "\\1", basename(sources(KiliNP_Blue_Timeseries)))
# 
# # Create new layer names
# 
# Kili.layer.names <- unlist(lapply(Kili.years, function(y) {
#   paste0(y, " - ", month.name)
# }))
# 
# # Assign the generated names to all four surface reflectance timeseries objects
# 
# names(KiliNP_Blue_Timeseries) <- Kili.layer.names
# names(KiliNP_Green_Timeseries) <- Kili.layer.names
# names(KiliNP_Red_Timeseries) <- Kili.layer.names
# names(KiliNP_NIR_Timeseries) <- Kili.layer.names
# 
# ### Generate the PPI raster stack ####
# 
# message("Calculating theoretical SZA vector for uncleaned rasters...")
# 
# ### 1. Calculate Astronomical Solar Zenith Angles
# # Extract the centroid coordinates of the Knepp_Buffered bounding box using the uncleaned raster extent
# 
# kili_ext <- ext(KiliNP_Blue_Timeseries)
# centroid_geom <- vect(
#   matrix(c(mean(c(kili_ext$xmin, kili_ext$xmax)), 
#            mean(c(kili_ext$ymin, kili_ext$ymax))), ncol=2), 
#   crs = crs(KiliNP_Blue_Timeseries)
# )
# 
# # Project the centroid to WGS84 to get the latitude in degrees
# 
# centroid_latlon <- project(centroid_geom, "EPSG:4326")
# kili_lat_deg <- geom(centroid_latlon)[,"y"]
# kili_lat_rad <- kili_lat_deg * (pi / 180) # Convert to radians
# 
# # Create a sequence of dates representing the middle of each composite month
# 
# dates <- seq(as.Date("2017-01-15"), as.Date("2021-12-15"), by = "1 month")
# day_of_year <- as.numeric(strftime(dates, format = "%j"))
# 
# # Calculate solar declination angle (delta) in radians
# 
# declination_deg <- 23.45 * sin((2 * pi / 365) * (day_of_year - 81))
# declination_rad <- declination_deg * (pi / 180)
# 
# # Calculate hour angle (h) in radians
# # Sentinel-2 descends at approx 10:30 AM local solar time. 
# # Hour angle = 15 degrees * (Hours from solar noon). 10.5 - 12.0 = -1.5 hours.
# 
# hour_angle_deg <- -1.5 * 15
# hour_angle_rad <- hour_angle_deg * (pi / 180)
# 
# # Calculate Solar Zenith Angle (theta_z) in radians
# 
# cos_theta_z <- sin(kili_lat_rad) * sin(declination_rad) + 
#   cos(kili_lat_rad) * cos(declination_rad) * cos(hour_angle_rad)
# sza_rad_vector <- acos(cos_theta_z)
# 
# ### 2. Calculate DVI
# 
# message("Calculating uncleaned Difference Vegetation Index (DVI)...")
# 
# # DVI is strictly NIR - Red, using the raw, uncleaned stacks
# 
# KiliNP_DVI <- KiliNP_NIR_Timeseries - KiliNP_Red_Timeseries
# 
# ### 3. Calculate PPI via terra::app()
# 
# message("Applying Plant Phenology Index (PPI) formula across uncleaned time series...")
# 
# # Define the PPI maths as a custom function to gracefully handle the ~5% NA gaps
# 
# calc_ppi <- function(dvi_vals, sza_vector) {
#   
#   # If the pixel is fully masked (e.g., an entirely cloud-covered pixel across all 5 years), return NAs
#   
#   if(all(is.na(dvi_vals))) return(rep(NA_real_, length(dvi_vals)))
#   
#   # Calculate M: Canopy maximum of DVI plus a small constant (0.005)
#   
#   M <- max(dvi_vals, na.rm = TRUE) + 0.005
#   
#   # Initialise the output vector
#   
#   ppi_out <- rep(NA_real_, length(dvi_vals))
#   valid <- !is.na(dvi_vals)
#   
#   if(any(valid)) {
#     v_dvi <- dvi_vals[valid]
#     v_sza <- sza_vector[valid]
#     
#     # Radiative transfer equations from Jin & Eklundh (2014)
#     # Assuming G = 0.5 (spherical leaf angle distribution)
#     
#     d_c <- 0.0336 + 0.0477 / cos(v_sza)
#     Q_E <- d_c + (1 - d_c) * 0.5 / cos(v_sza)
#     K <- 1 / (4 * Q_E) * (1 + M) / (1 - M)
#     
#     # Logarithmic transformation (assuming bare soil DVI = 0.09)
#     
#     log_arg <- (M - v_dvi) / (M - 0.09)
#     
#     # Prevent negative values inside the logarithm (can occur with extreme data anomalies)
#     
#     log_arg[log_arg <= 0] <- NA_real_
#     
#     ppi_out[valid] <- -K * log(log_arg)
#   }
#   
#   return(ppi_out)
# }
# 
# # Apply the function across the z-dimension
# 
# KiliNP_PPI_Timeseries <- app(KiliNP_DVI, fun = calc_ppi, sza_vector = sza_rad_vector)
# 
# # Carry over the layer names for consistency
# 
# names(KiliNP_PPI_Timeseries) <- names(KiliNP_DVI)
# 
# # Export the raw PPI stack to your external drive to manage local storage
# 
# writeRaster(
#   KiliNP_PPI_Timeseries, 
#   filename = file.path(KiliNP_PPI_Processed, "KiliNP_PPI_2017-2021_Timeseries.tif"), 
#   overwrite = TRUE
# )
# 
# # Load the NDVI raster back in!
# 
# KiliNP_NDVI_Timeseries <- rast(file.path(KiliNP_NDVI_Processed, "KiliNP_NDVI_2017-2021_Timeseries.tif"))

# ### Mask pixels in the raster stack which don't have a complete timeseries of data
# ## Only pixels with a complete set of data for every layer are suitable for analysis
# # There are many pixels with NA values scattered throughout the raster stack
# # Create logical mask: TRUE only where ALL layers are non-NA
# # I'm making the mask from just the NIR band, but as all bands have the same coverage, that's okay
# 
# Kili.pixel.mask <- app(KiliNP_NIR_Timeseries, function(x) all(!is.na(x))) # Will consume lots of RAM
# 
# # Mask out incomplete pixels (FALSE becomes NA)
# 
# KiliNP_Blue_Timeseries_Clean <- mask(KiliNP_Blue_Timeseries, Kili.pixel.mask, maskvalues = 0) # This is very computationally challenging, run on HPC
# KiliNP_Green_Timeseries_Clean <- mask(KiliNP_Green_Timeseries, Kili.pixel.mask, maskvalues = 0)
# KiliNP_Red_Timeseries_Clean <- mask(KiliNP_Red_Timeseries, Kili.pixel.mask, maskvalues = 0)
# KiliNP_NIR_Timeseries_Clean <- mask(KiliNP_NIR_Timeseries, Kili.pixel.mask, maskvalues = 0)
# 
# # Export rasters so I don't have to calculate it every time
# 
# writeRaster(KiliNP_Blue_Timeseries_Clean, file.path(KiliNP_NDVI_Processed, "KiliNP_BlueBand_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE)
# writeRaster(KiliNP_Green_Timeseries_Clean, file.path(KiliNP_NDVI_Processed, "KiliNP_GreenBand_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE)
# writeRaster(KiliNP_Red_Timeseries_Clean, file.path(KiliNP_NDVI_Processed, "KiliNP_RedBand_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE)
# writeRaster(KiliNP_NIR_Timeseries_Clean, file.path(KiliNP_NDVI_Processed, "KiliNP_NIRBand_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE)
# 
# # And then load back in the rasters
# 
# KiliNP_Blue_Timeseries_Clean <- rast(file.path(KiliNP_NDVI_Processed, "KiliNP_BlueBand_2017-2021_Cropped_&_Masked.tif")) # Load it in
# KiliNP_Green_Timeseries_Clean <- rast(file.path(KiliNP_NDVI_Processed, "KiliNP_GreenBand_2017-2021_Cropped_&_Masked.tif")) # Load it in
# KiliNP_Red_Timeseries_Clean <- rast(file.path(KiliNP_NDVI_Processed, "KiliNP_RedBand_2017-2021_Cropped_&_Masked.tif")) # Load it in
# KiliNP_NIR_Timeseries_Clean <- rast(file.path(KiliNP_NDVI_Processed, "KiliNP_NIRBand_2017-2021_Cropped_&_Masked.tif")) # Load it in

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

current.names <- names(Knepp_Buffered_NDVI_Timeseries)

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

Knepp_Buffered.full.dates <- seq(min(current.dates), max(current.dates), by = "month")
missing.dates <- Knepp_Buffered.full.dates[!Knepp_Buffered.full.dates %in% current.dates]

if (length(missing.dates) > 0) {
  message(paste("Found", length(missing.dates), "missing months. Creating blank template layers..."))
  
  Empty_Knepp_Buffered_Raster <- terra::init(Knepp_Buffered_NDVI_Timeseries[[1]], NA)
  missing.rasters <- terra::rast(replicate(length(missing.dates), Empty_Knepp_Buffered_Raster))
  
  # Format missing names to match the "YYYY - Month" convention
  
  names(missing.rasters) <- paste0(format(missing.dates, "%Y"), " - ", months(missing.dates))
  
  Knepp_Buffered_NDVI_Timeseries <- c(Knepp_Buffered_NDVI_Timeseries, missing.rasters)
  
  # Sort the entire stack chronologically so NAs are in the correct sequence
  
  sorted.names <- paste0(format(Knepp_Buffered.full.dates, "%Y"), " - ", months(Knepp_Buffered.full.dates))
  Knepp_Buffered_NDVI_Timeseries <- Knepp_Buffered_NDVI_Timeseries[[sorted.names]]
  
} else {
  message("No missing months found. Temporal sequence is already contiguous.")
}

## 3. Apply the gap filling with Savitzky-Golay filter

message("Gap filling and smoothing the dataset...")

Knepp_Buffered.cores <- max(1, detectCores() - 2)

Knepp_Buffered_NDVI_Timeseries.SGfilter <- app(
  Knepp_Buffered_NDVI_Timeseries, 
  fun = sg_gapfill, 
  cores = Knepp_Buffered.cores 
)

# Carry over the names to the smoothed raster

names(Knepp_Buffered_NDVI_Timeseries.SGfilter) <- names(Knepp_Buffered_NDVI_Timeseries)

# 4. Export to a file (because it's not saved into the R data environment)

message("Exporting smoothed NDVI timeseries...")
writeRaster(
  Knepp_Buffered_NDVI_Timeseries.SGfilter,
  filename = file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_NDVI_2017-2021_Timeseries_SG-filtered.tif"),
  overwrite = TRUE
)

Knepp_Buffered_NDVI_Timeseries.SGfilter <- rast( # Load it back in
  file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_NDVI_2017-2021_Timeseries_SG-filtered.tif"))

### Inspect temporal structure ####
## Unlike Macchia Sacra NetCDF, this GeoTIFF stack does not contain explicit time metadata
## Therefore, we must construct the time vector manually from layer names

message("Constructing time vector for Knepp_Buffered time series...")

# Extract year and month from layer names
# Layer format: "2017 - January"

Knepp_Buffered.dates <- as.Date(
  paste0(
    sub(" - .*", "", names(Knepp_Buffered_NDVI_Timeseries.SGfilter)), "-",
    match(sub(".* - ", "", names(Knepp_Buffered_NDVI_Timeseries.SGfilter)), month.name),
    "-15" # As SZA used the 15th of the month, I will set the date here to be the 15th
  ),
  format = "%Y-%m-%d"
)

stopifnot(length(Knepp_Buffered.dates) == nlyr(Knepp_Buffered_NDVI_Timeseries.SGfilter))

message(paste("Temporal length:", length(Knepp_Buffered.dates), "layers"))

### 1. Shannon-Wiener Index ####
## As in Macchia Sacra, collapse the time series to mean annual trajectory first

message("Calculating Shannon-Wiener diversity index for Knepp_Buffered...")

Knepp_Buffered_Mean_NDVI_Raster <- app(Knepp_Buffered_NDVI_Timeseries.SGfilter, fun = mean, na.rm = TRUE)

# Export raster so I don't have to calculate it every time

writeRaster(Knepp_Buffered_Mean_NDVI_Raster, file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_MeanNDVI-SGfiltered.tif"), overwrite = TRUE)

# And then load back in the raster

Knepp_Buffered_Mean_NDVI_Raster <- rast(file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_MeanNDVI-SGfiltered.tif")) # Load it in

# Round to 2 decimals to avoid numerical saturation

Knepp_Buffered_Mean_NDVI_Raster2dec <- round(Knepp_Buffered_Mean_NDVI_Raster, 2)
Knepp_Buffered_Mean_NDVI_Raster2dec <- trim(Knepp_Buffered_Mean_NDVI_Raster2dec) # Trimming it to avoid unnecessary computation

# Export and reload the rounded raster so I don't have to calculate it every time

writeRaster(Knepp_Buffered_Mean_NDVI_Raster2dec, file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_MeanNDVI_Rounded2DP-SGfiltered.tif"), overwrite = TRUE)
Knepp_Buffered_Mean_NDVI_Raster2dec <- rast(file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_MeanNDVI_Rounded2DP-SGfiltered.tif")) # Load it back in

# Run ShannonS

Knepp_Buffered.NDVI.ShannonH.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(Knepp_Buffered_Mean_NDVI_Raster2dec, wide = TRUE), # The function here converts my raster to a matrix suitable for analysis
  window = 3,
  na.tolerance = 0
)

## Now turn the matrix of Shannon H values into a spatial raster
# Put it in a raster signature

Knepp_Buffered_NDVI_ShannonH_Raster <- rast(Knepp_Buffered.NDVI.ShannonH.matrix)

# Make the raster's extent and CRS match the original raster

ext(Knepp_Buffered_NDVI_ShannonH_Raster) <- ext(Knepp_Buffered_Mean_NDVI_Raster2dec)
crs(Knepp_Buffered_NDVI_ShannonH_Raster) <- crs(Knepp_Buffered_Mean_NDVI_Raster2dec)

names(Knepp_Buffered_NDVI_ShannonH_Raster) <- "Shannon's H"

# Export the raster or reload it if necessary

writeRaster(
  Knepp_Buffered_NDVI_ShannonH_Raster,
  filename = file.path(Knepp_Buffered_Results, "Knepp_Buffered_MeanNDVI_ShannonH-SGfiltered.tif"),
  overwrite = TRUE
)

Knepp_Buffered_NDVI_ShannonH_Raster <- rast(file.path(Knepp_Buffered_Results, "Knepp_Buffered_MeanNDVI_ShannonH-SGfiltered.tif")) # Load it back in

### 2. Classic Rao's Q  ####
## Due to the large size of the raster, I need to tile it so that it can be run
## The tiles will be stitched back together once they're computed

message("Calculating classical Rao's Q for Knepp_Buffered...")

### Step 1: Create a grid to define zones for tiling

# Optional but recommended: trim outer NA borders

trimmed.Knepp_Buffered_Mean_NDVI_Raster <- trim(Knepp_Buffered_Mean_NDVI_Raster)

# For effective parallelisation on MaRC3a, 2000 tiny tiles is effective

Knepp_Buffered.total.tiles <- 2000
Knepp_Buffered.aspect.ratio <- ncol(trimmed.Knepp_Buffered_Mean_NDVI_Raster) / nrow(trimmed.Knepp_Buffered_Mean_NDVI_Raster)

# find factor pairs of 72 

tiling.factors <- expand.grid(
  ncols = 1:Knepp_Buffered.total.tiles,
  nrows = 1:Knepp_Buffered.total.tiles
)

# Subset to just pairs which equal "Knepp_Buffered.total.tiles" when multiplied

tiling.factors <- tiling.factors[tiling.factors$ncols * tiling.factors$nrows == Knepp_Buffered.total.tiles, ]

## choose the pair closest to the raster's aspect ratio
# First, create a new column calculating the difference in aspect ratio

tiling.factors$ratio_diff <- abs((tiling.factors$ncols / tiling.factors$nrows) - Knepp_Buffered.aspect.ratio)

# Now find the pair which has the lowest distance from a square 1:1 aspect ratio

best.tile.size <- tiling.factors[which.min(tiling.factors$ratio_diff), ]

# And save that for later usage

Knepp_Buffered.cols <- best.tile.size$ncols
Knepp_Buffered.rows <- best.tile.size$nrows

# Create spatial polygons to specify what the tile sizes should be

Knepp_Buffered.tiling.grid <- as.polygons(
  rast(
    ext(trimmed.Knepp_Buffered_Mean_NDVI_Raster),
    ncols = Knepp_Buffered.cols,
    nrows = Knepp_Buffered.rows,
    crs = crs(trimmed.Knepp_Buffered_Mean_NDVI_Raster)
  )
)

writeVector(Knepp_Buffered.tiling.grid, file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_Tiling_Grid_Polygons.geoJSON"), filetype = "GeoJSON" , overwrite = TRUE) # Export for later use

Knepp_Buffered.tiling.grid <- vect(file.path(Knepp_Buffered_NDVI_Processed, "Knepp_Buffered_Tiling_Grid_Polygons.geoJSON")) # Load it back in

plot(Knepp_Buffered.tiling.grid) # Plot it to make sure that it's loaded in

# Window size (I think 3 is the default, but this can be changed as necessary)

RaoQ.window.size <- 3
Knepp_Buffered.tile.overlap <- floor(RaoQ.window.size / 2)

# Create a directory to put the tiles in

#Knepp_Buffered.tile.dir <- file.path(Knepp_Buffered_Processed,"Mean NDVI tiles") # I disabled this line so I don't overwrite my 72 larger tiles
Knepp_Buffered.tile.dir <- file.path(Knepp_Buffered_NDVI_Processed,"Mean_NDVI_Knepp_Buffered_SG_Tiles")
dir.create(Knepp_Buffered.tile.dir, recursive = TRUE, showWarnings = FALSE)

## Finally, create the tiles

Knepp_Buffered.tiles <- makeTiles(
  trimmed.Knepp_Buffered_Mean_NDVI_Raster, # Trimmed for easier computation
  y = Knepp_Buffered.tiling.grid, # Specifies how the tiles should be allocated
  buffer = Knepp_Buffered.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
  filename = file.path(Knepp_Buffered.tile.dir, "Knepp_Buffered-SG_MeanNDVI_Tile-.tif"),
  overwrite = FALSE
)

# ### Step 2: Compute classic Rao's Q for each tile
# ## Firstly, I need to setup the environment for parallelisation
# # Create a subfolder to store the classic Rao's Q output tiles
 
Knepp_Buffered.rao.dir  <- file.path(Knepp_Buffered_Processed,"MeanNDVI_Rao-utputs") 
dir.create(Knepp_Buffered.rao.dir, recursive = TRUE, showWarnings = FALSE)
 
# ## Create a computing cluster to parallelise the calculation at the tile level
# # Set the number of cores to be used by the cluster
# 
# Knepp_Buffered.cores <- max(1, detectCores() - 2)
# 
# # Initialise a log file so I can actually see what's going on
# 
# Knepp_Buffered.log.file <- file.path(Knepp_Buffered.rao.dir, "Knepp_Buffered_RaoQ_processing_log.txt")
# 
# # If the log file doesn't exist already, create one
# 
# if(!file.exists(Knepp_Buffered.log.file)) file.create(Knepp_Buffered.log.file)
# 
# # Create the cluster (alliterative and punny names are mandatory)
# 
# Knepp_Buffered.cluster <- makeCluster(Knepp_Buffered.cores)
# 
# clusterEvalQ(Knepp_Buffered.cluster, {
#   library(terra)
#   library(rasterdiv)
# })
# 
# clusterExport(Knepp_Buffered.cluster, c(
#   "Knepp_Buffered.tiles",
#   "Knepp_Buffered.rao.dir",
#   "RaoQ.window.size",
#   "Knepp_Buffered.log.file"
# ))
# 
# # Identify tiles still needing processing (so resources aren't wasted processing tiles already done)
# 
# tile.outputs <- file.path(
#   Knepp_Buffered.rao.dir,
#   paste0("Knepp_Buffered_Classic-RaoQ_Tile-", seq_along(Knepp_Buffered.tiles), ".tif")
# )
# 
# tiles.to.process <- which(!file.exists(tile.outputs))
# 
# cat(length(tiles.to.process), "tiles remaining.\n")

## Now actually run the code
# This version creates a process for each CPU core and runs each tile as a single process
## REVIEWERS: Due to the CPU overhead of this workload, I ultimately decided to run it on my university's supercomputer instead
# Please see "03.1C_Knepp_Buffered_Classical-RaoQ_MaRC3a.R" for the job file I submitted
# 
# Knepp_Buffered.classic.rao.results <- parLapply( # Function call
#   Knepp_Buffered.cluster,
#   tiles.to.process,
#   function(i){
#     
#     library(terra)
#     library(rasterdiv)
#     
#     log_file <- Knepp_Buffered.log.file
#     
#     log_msg <- function(msg){
#       cat(
#         paste0(Sys.time(), " | Worker ", Sys.getpid(), " | ", msg, "\n"),
#         file = log_file,
#         append = TRUE
#       )
#     }
#     
#     out.file <- file.path(
#       Knepp_Buffered.rao.dir,
#       paste0("Knepp_Buffered_Classic-RaoQ_Tile-", i, ".tif")
#     )
#     
#     if(file.exists(out.file)){
#       log_msg(paste0("Tile", i, "already exists — skipped"))
#       return(NULL)
#     }
#     
#     log_msg(paste("Tile", i, "STARTED"))
#     
#     tmp.tile <- rast(Knepp_Buffered.tiles[i]) # Load in the raster for processing
#     
#     tmp.result <- paRao(
#       tmp.tile,
#       window = RaoQ.window.size,
#       alpha = 2,
#       simplify = 2, # This is necessary to maintain consistency with the Shannon's H test (keeps just 2 decimal places)
#       method = "classic", # Because this is not looking at timeseries Rao's Q, just regular unidimensional Rao's Q
#       np = 1 # Explicitly prevents nested parallelisation (or set above 1 if you want to melt your CPU)
#     )
#     
#     tmp.rao_raster <- tmp.result[[1]][[1]] # Subsetting avoids hardcoding "$window.3$alpha.2"
#     
#     writeRaster(
#       tmp.rao_raster,
#       filename = out.file,
#       overwrite = TRUE
#     )
#     
#     rm(tmp.tile,tmp.result,tmp.rao_raster)
#     gc()
#     
#     log_msg(paste("Tile №", i, "'s classic Rao's Q calculated successfully."))
#     
#     return(NULL) # So that each worker doesn't fill up R's memory with bloat upon completion
#   }
# )
# 
# ## This for loop is an alternative computational approach which uses all cores to work sequentially over each tile
# # This version seems computationally safer because each tile outputted is like a mini-checkpoint in the event that computation is interrupted
# 
# for(i in seq_along(Knepp_Buffered.tiles)){
#   
#   log_file <- Knepp_Buffered.log.file
#   
#   log_msg <- function(msg){
#     cat(
#       paste0(Sys.time(), " | Tile ", i, " | ", msg, "\n"),
#       file = log_file,
#       append = TRUE
#     )
#   }
#   
#   out.file <- file.path(
#     Knepp_Buffered.rao.dir,
#     paste0("Knepp_Buffered_Classic-RaoQ_Tile-", i, ".tif")
#   )
#   
#   # Skip tiles which already exist (prevents recomputation)
#   
#   if(file.exists(out.file)){
#     log_msg("already exists — skipped")
#     next
#   }
#   
#   log_msg("STARTED")
#   
#   tmp.tile <- rast(Knepp_Buffered.tiles[i]) # Load in the raster for processing
#   
#   tmp.result <- paRao(
#     tmp.tile,
#     window = RaoQ.window.size,
#     alpha = 2,
#     simplify = 2, # This is necessary to maintain consistency with the Shannon's H test (keeps just 2 decimal places)
#     method = "classic", # Because this is not looking at timeseries Rao's Q, just regular unidimensional Rao's Q
#     np = Knepp_Buffered.cores # Parallelise INSIDE paRao for faster per-tile processing
#   )
#   
#   tmp.rao_raster <- tmp.result[[1]][[1]] # Subsetting avoids hardcoding "$window.3$alpha.2"
#   
#   writeRaster(
#     tmp.rao_raster,
#     filename = out.file,
#     overwrite = TRUE
#   )
#   
#   rm(tmp.tile,tmp.result,tmp.rao_raster)
#   gc()
#   
#   log_msg("classic Rao's Q calculated successfully.")
# }

### Step 3: Demosaic the classical Rao's Q tiles
## Gather up all the files

Knepp_Buffered.rao.files <- list.files(
  Knepp_Buffered.rao.dir,
  pattern = "Knepp_Buffered_NDVI_Classic-RaoQ_SG_Tile-",
  full.names = TRUE
)

# Tell R to apply the `rast` function to them 

Knepp_Buffered.rao.tiles <- lapply(Knepp_Buffered.rao.files, rast)

# Convert them into a spatial raster collection:

Knepp_Buffered.rao.tiles <- sprc(Knepp_Buffered.rao.tiles)

# Run the demosaic function (which is curiously called `mosaic`)

Knepp_Buffered_NDVI_Classic_RaoQ <- terra::mosaic(Knepp_Buffered.rao.tiles)

# Export the final raster, and load it back in if necessary

writeRaster(
  Knepp_Buffered_NDVI_Classic_RaoQ,
  file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_Classic-RaoQ-SGfiltered.tif"),
  overwrite = TRUE
)

Knepp_Buffered_NDVI_Classic_RaoQ <- rast(file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_Classic-RaoQ-SGfiltered.tif")) # Load in the raster

plot(Knepp_Buffered_NDVI_Classic_RaoQ) # Plot it! (Good for data exploration and checking that the raster loaded in as normal)

### 3. Rao's Q with TWDTW ####

message("Calculating Rao's Q with TWDTW distance for Knepp_Buffered...")

## Step 1: I'll have to tile this as well because it is too large to compute as a single object
# This tiling script is copied from Step 1 of the classical Rao's Q analysis
# Some objects like "Knepp_Buffered.tiling.grid" are assumed to be loaded, and "Knepp_Buffered.tile.dir" is overwritten

# Create a directory to put the tiles in

Knepp_Buffered.twdtw.rao.dir <- file.path(Knepp_Buffered_NDVI_Processed,"Timeseries_NDVI_Knepp_Buffered_SG_Tiles")
dir.create(Knepp_Buffered.twdtw.rao.dir, recursive = TRUE, showWarnings = FALSE)

# Create the timeseries tiles

Knepp_Buffered.twdtw.tiles <- makeTiles(
  trim(Knepp_Buffered_NDVI_Timeseries.SGfilter), # Trimmed for easier computation
  y = Knepp_Buffered.tiling.grid, # Make sure this is still loaded in from the previous step!
  buffer = Knepp_Buffered.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
  filename = file.path(Knepp_Buffered.twdtw.rao.dir, "Knepp_Buffered-SG_TimeseriesNDVI_Tile-.tif"),
  overwrite = FALSE
)

## Step 2: Submit the tiles for processing on University of Marburg's supercomputer
# Please see script 03.2A_Knepp_Buffered_TWDTW-RaoQ_MaRC3a.R for the actual job scripts used
# If you wish to attempt computation on your local machine, then please see the code below

######### Beginning of not actually used section
# Knepp_Buffered_Rao_TWDTW <- paRao(
#   x = KiliNP_Timeseries_Clean,
#   time_vector = Kili.dates,
#   window = 3,
#   alpha = 2,
#   na.tolerance = 0,
#   simplify = 2,
#   np = detectCores() -1,
#   progBar = TRUE,
#   method = "multidimension",
#   dist_m = "twdtw",
#   midpoint = 6,          # Midpoint of annual cycle (June)
#   stepness = -0.5,
#   cycle_length = "year",
#   time_scale = "month"   # Now explicitly monthly data
# )
# 
# writeRaster(
#   Kili_Rao_TWDTW$window.3$alpha.2,
#   filename = file.path(KiliNP_Results, "KiliNP_RaoQ_TWDTW.tif"),
#   overwrite = TRUE
# )
######### End of not actually used section

## Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster
## Gather up all the files

Knepp_Buffered.twdtw.rao.files <- list.files(
  file.path(Knepp_Buffered_Processed,"TWDTWNDVI_Rao-utputs") ,
  pattern = "Knepp_Buffered_NDVI_TWDTW-RaoQ_SG_Tile",
  full.names = TRUE
)

# Tell R to apply the `rast` function to them 

Knepp_Buffered.twdtw.rao.files <- lapply(Knepp_Buffered.twdtw.rao.files, rast)

# Convert them into a spatial raster collection:

Knepp_Buffered.twdtw.rao.files <- sprc(Knepp_Buffered.twdtw.rao.files)

# Run the demosaic function (which is curiously called `mosaic`)

Knepp_Buffered_NDVI_TWDTW_RaoQ <- terra::mosaic(Knepp_Buffered.twdtw.rao.files)

# Export the final raster, and load it back in if necessary

writeRaster(
  Knepp_Buffered_NDVI_TWDTW_RaoQ,
  file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_TWDTW-RaoQ-SGfiltered.tif"),
  overwrite = TRUE
)

Knepp_Buffered_NDVI_TWDTW_RaoQ <- rast(file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_TWDTW-RaoQ-SGfiltered.tif")) # Load in the raster

plot(Knepp_Buffered_NDVI_TWDTW_RaoQ)

### Export rasters for comparison ####

Knepp_Buffered_NDVI_Comparison_Rasters <- c(
  trim(Knepp_Buffered_Mean_NDVI_Raster), # Trimmed so that it matches the extent of the other rasters
  Knepp_Buffered_NDVI_ShannonH_Raster,
  Knepp_Buffered_NDVI_Classic_RaoQ,
  Knepp_Buffered_NDVI_TWDTW_RaoQ
)

names(Knepp_Buffered_NDVI_Comparison_Rasters) <- c( # This sets nice layer names for easier browsing
  "Sentinel-2 NDVI",
  "Shannon's H",
  "Classic Rao's Q",
  "TWDTW Rao's Q"
)

writeRaster( # So I don't have to compute it every time
  Knepp_Buffered_NDVI_Comparison_Rasters,
  filename = file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_Diversity_Comparison-SGfiltered.tif"),
  overwrite = TRUE
)

Knepp_Buffered_NDVI_Comparison_Rasters <- rast(file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_Diversity_Comparison-SGfiltered.tif")) # Load it back in

png(file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_Indices_Comparison.png"), # Exported for the paper
    width = 2560, height = 1440, res = 150)

plot(Knepp_Buffered_NDVI_Comparison_Rasters) # Plot it in a file for export

dev.off()

plot(Knepp_Buffered_NDVI_Comparison_Rasters) # Plot it to see how it looks

### Assess index performance using vegetation ground truth ####

message("Assessing diversity indices against vegetation ground truth...")

# Load Knepp_Buffered_LandCover_Vector boundary again (this is our land cover ground truth data)

Knepp_Buffered_LandCover_Vector <- 
  vect(file.path(Knepp_Buffered_Input,
                 "/Knepp_Buffered Ground Truthing Land Cover Classifications/VegAug1_Knepp_Buffered_SES_withnewcof.shp")) # Load in the ground truth data

# Ensure CRS matches

if (crs(Knepp_Buffered_NDVI_Comparison_Rasters) != crs(Knepp_Buffered_LandCover_Vector)){
  Knepp_Buffered_LandCover_Vector <- project(Knepp_Buffered_LandCover_Vector, crs(Knepp_Buffered_NDVI_Comparison_Rasters))
}

# Crop and mask diversity rasters to ground truth extent

masked.Knepp_Buffered_NDVI_ShannonH_Raster <- mask(crop(Knepp_Buffered_NDVI_ShannonH_Raster, Knepp_Buffered_LandCover_Vector), # Crop
                                      Knepp_Buffered_LandCover_Vector) # and mask
masked.Knepp_Buffered_NDVI_Classic_RaoQ <- mask(crop(Knepp_Buffered_NDVI_Classic_RaoQ, Knepp_Buffered_LandCover_Vector), 
                                   Knepp_Buffered_LandCover_Vector)
masked.Knepp_Buffered_NDVI_TWDTW_RaoQ <- mask(crop(Knepp_Buffered_NDVI_TWDTW_RaoQ, Knepp_Buffered_LandCover_Vector), 
                                 Knepp_Buffered_LandCover_Vector)

## Rasterise vegetation class
# First I need to update the ground truth vector to use the proper category names
# I'll use a lookup table

Knepp_Buffered.land.cover.lookup <- c(
  "0"  = "Snow/glacier",
  "1"  = "Agriculture (MAI)",
  "2"  = "Savannah (SAV)",
  "3"  = "Swamp",
  "4"  = "Overgrown clearing",
  "7"  = "Forest plantation",
  "9"  = "Riverine",
  "10" = "Upper montane Erica excelsa forest (FPO Podocarpus disturbed)",
  "11" = "Subalpine Erica trimera bushland (FED incl FER (Erica forest and bushland))",
  "12" = "Podocarpus forest (FPO)",
  "13" = "Subalpine tussock grassland",
  "14" = "Chagga homegardens (HOM)",
  "15" = "Alpine Helichrysum vegetation (HEL)",
  "16" = "Ocotea forest (FOC)",
  "17" = "Bare rock",
  "18" = "Sub/lower montane rainforest (FLM)",
  "19" = "Coffee plantations (COF)"
)

# Assign readable names to the grid code vector of the ground truth

Knepp_Buffered_LandCover_Vector$grid_code <- Knepp_Buffered.land.cover.lookup[as.character(Knepp_Buffered_LandCover_Vector$grid_code)]

# Rasterise the land cover vector

Knepp_Buffered_LandCover_Raster <- rasterize(
  Knepp_Buffered_LandCover_Vector,
  Knepp_Buffered_NDVI_Comparison_Rasters,
  field = "grid_code"
)

### Convert the index rasters to a dataframe for performance analysis ####

Knepp_Buffered_Indices_Comparison_Raster <- c(
  masked.Knepp_Buffered_NDVI_ShannonH_Raster,
  masked.Knepp_Buffered_NDVI_Classic_RaoQ,
  masked.Knepp_Buffered_NDVI_TWDTW_RaoQ,
  Knepp_Buffered_LandCover_Raster
)

names(Knepp_Buffered_Indices_Comparison_Raster) <- c(
  "ShannonsH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

Knepp_Buffered_Indices_Comparison_DF <- as.data.frame(
  Knepp_Buffered_Indices_Comparison_Raster,
  na.rm = TRUE
)

colnames(Knepp_Buffered_Indices_Comparison_DF) <- c(
  "ShannonsH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

### PERMANOVA ####
## These datasets are too large to conduct a PERMANOVA upon (37121.9GB RAM required)
## Instead, I will use a random representative subset of the data
# Subset the data

subset.Knepp_Buffered_Indices_Comparison_DF <- Knepp_Buffered_Indices_Comparison_DF[sample(nrow(Knepp_Buffered_Indices_Comparison_DF), 10000), ]

# Conduct a series of PERMANOVAs

PERMANOVA_ShannonsH <- adonis2(
  subset.Knepp_Buffered_Indices_Comparison_DF$ShannonsH ~ subset.Knepp_Buffered_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel = Knepp_Buffered.cores # I've set it to parallelise using the Knepp_Buffered.cores argument from before
)

PERMANOVA_RaosQ_Classic <- adonis2(
  subset.Knepp_Buffered_Indices_Comparison_DF$RaosQ_Classic ~ subset.Knepp_Buffered_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel = Knepp_Buffered.cores
)

PERMANOVA_RaosQ_TWDTW <- adonis2(
  subset.Knepp_Buffered_Indices_Comparison_DF$RaosQ_TWDTW ~ subset.Knepp_Buffered_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel = Knepp_Buffered.cores
)

# Put the PERMANOVA results into a dataframe for effective presentation

Knepp_Buffered_PERMANOVA_Results <- data.frame(
  Index = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2 = c(
    PERMANOVA_ShannonsH$R2[1],
    PERMANOVA_RaosQ_Classic$R2[1],
    PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  F = c(
    PERMANOVA_ShannonsH$F[1],
    PERMANOVA_RaosQ_Classic$F[1],
    PERMANOVA_RaosQ_TWDTW$F[1]
  ),
  p_value = c(
    PERMANOVA_ShannonsH$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
  )
)

print(Knepp_Buffered_PERMANOVA_Results)

saveRDS(Knepp_Buffered_PERMANOVA_Results, file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_PERMANOVA_Results.rds")) # Save it so I don't have to recalculate repeatedly
Knepp_Buffered_PERMANOVA_Results <- readRDS(file.path(Knepp_Buffered_Results, "Knepp_Buffered_NDVI_PERMANOVA_Results.rds")) # And load it back in if necessary

message("Knepp_Buffered analysis complete.")