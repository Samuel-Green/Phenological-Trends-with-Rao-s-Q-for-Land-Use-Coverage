############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 06/08/2026                    #
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

# Core spatial stack already loaded in 00_setup.R:
# terra, sf, here, dplyr, stringr

### Define the  file paths and import the site data ####
## Output directory for this script

Knepp_Results <- file.path(Results, "Knepp Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Explicitly set the input directory to the external drive to manage storage

Knepp_Input_Ext <- "E:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Input Data" # External location because my computer is low on storage
Knepp_Processed_Ext <- "E:/Elliot Shayle/Knepp Estate Geodata/Sentinel-2 Processed Data" # External location because my computer is low on storage

message("Importing Sentinel-2 PPI rasters...")

# Load Knepp Buffered boundaries so that I can crop the raster to the right size

KneppEstate_Boundaries <- vect(file.path(Knepp_Input,"Schulte-to-Bühne_Knepp_Site_Data","knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed,"Knepp_Estate_Buffered_Boundary.shp"))

message("Beginning cropping and masking of raw Sentinel-2 rasters...")

## List raster files
# The regex now explicitly looks for BLUE, GREEN, RED, or NIR
tmp.files <- list.files(
  Knepp_Input_Ext, 
  pattern = "SEN2L_(BLUE|GREEN|RED|NIR)_FBM\\.tif$", 
  full.names = TRUE
)

# A quick mapping dictionary to convert the text colors to standard Sentinel-2 band numbers
band_map <- c("BLUE" = "B02", "GREEN" = "B03", "RED" = "B04", "NIR" = "B08")

for (f in tmp.files) {
  
  tmp.Knepp.raster <- rast(f)
  
  # Ensure CRS matches
  if (crs(tmp.Knepp.raster) != crs(KneppEstate_Buffered_Boundaries)) {
    KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(tmp.Knepp.raster))
  }
  
  # Attempt to crop and mask. If the extents don't overlap (like your 2018 files), 
  # tryCatch catches the error, prints a message, and returns NULL.
  tmp.Knepp.raster.masked <- tryCatch({
    crop(tmp.Knepp.raster, KneppEstate_Buffered_Boundaries, mask = TRUE)
  }, error = function(e) {
    message(paste("  -> Skipping non-overlapping raster:", basename(f)))
    return(NULL)
  })
  
  # If tryCatch returned NULL (i.e., it failed to crop), skip the rest of this loop iteration
  if (is.null(tmp.Knepp.raster.masked)) {
    rm(tmp.Knepp.raster) # clean up the raw raster
    next
  }
  
  # Extract the year (first 4 digits)
  tmp.year <- sub("^(\\d{4}).*", "\\1", basename(f)) 
  
  # Extract the color string (BLUE, GREEN, RED, or NIR)
  tmp.color <- sub(".*_SEN2L_([A-Z]+)_FBM\\.tif$", "\\1", basename(f)) 
  
  # Map the color to standard band notation (e.g., "BLUE" becomes "B02")
  tmp.band <- band_map[tmp.color]
  
  # Create new filename incorporating both year and band
  tmp.new.name <- paste0("Knepp_", tmp.year, "_", tmp.band, "_Cropped.tif")
  
  # Write to processed folder
  writeRaster(
    tmp.Knepp.raster.masked,
    filename = file.path(Knepp_Processed_Ext, tmp.new.name),
    overwrite = TRUE
  )
  message(paste0("Writing raster: ", tmp.new.name))
  
  # Cleanup RAM
  rm(tmp.Knepp.raster, tmp.Knepp.raster.masked)
  gc()
}
rm(band_map) # Has to be outside the loop or the 2nd iteration of the loop will crash

message("Importing the cropped multi-band rasters...")

# 1. Blue Band (B02)
Knepp_B02_Files <- list.files(
  Knepp_Processed_Ext,
  pattern = "^Knepp_\\d{4}_B02_Cropped\\.tif$",
  full.names = TRUE
)
Knepp_Blue_Timeseries <- rast(Knepp_B02_Files)

# 2. Green Band (B03)
Knepp_B03_Files <- list.files(
  Knepp_Processed_Ext,
  pattern = "^Knepp_\\d{4}_B03_Cropped\\.tif$",
  full.names = TRUE
)
Knepp_Green_Timeseries <- rast(Knepp_B03_Files)

# 3. Red Band (B04)
Knepp_B04_Files <- list.files(
  Knepp_Processed_Ext,
  pattern = "^Knepp_\\d{4}_B04_Cropped\\.tif$",
  full.names = TRUE
)
Knepp_Red_Timeseries <- rast(Knepp_B04_Files)

# 4. Near-Infrared Band (B08)
Knepp_B08_Files <- list.files(
  Knepp_Processed_Ext,
  pattern = "^Knepp_\\d{4}_B08_Cropped\\.tif$",
  full.names = TRUE
)
Knepp_NIR_Timeseries <- rast(Knepp_B08_Files)

message("Applying standard layer names...")

## Rename the layers to something a bit more readable
# Extract years directly from the loaded Blue Timeseries file sources
Knepp.years <- sub(".*Knepp_(\\d{4})_B02_Cropped\\.tif$", "\\1", basename(sources(Knepp_Blue_Timeseries)))

# Create clean base layer names (e.g., "2017-01", "2017-02")
Knepp.layer.names.base <- unlist(lapply(Knepp.years, function(y) {
  paste0(y, "-", sprintf("%02d", 1:12))
}))

# Assign the generated names with explicit band suffixes to each object
names(Knepp_Blue_Timeseries)  <- paste0(Knepp.layer.names.base, "_B02")
names(Knepp_Green_Timeseries) <- paste0(Knepp.layer.names.base, "_B03")
names(Knepp_Red_Timeseries)   <- paste0(Knepp.layer.names.base, "_B04")
names(Knepp_NIR_Timeseries)   <- paste0(Knepp.layer.names.base, "_B08")

message("Data prep complete! Ready for PPI calculation.")

### Compute a PPI raster ####

message("Calculating Astronomical Solar Zenith Angles for Knepp Estate...")

## 1. Calculate Astronomical Solar Zenith Angles
# Extract the centroid coordinates of the Knepp bounding box

knepp_ext <- ext(Knepp_NIR_Timeseries)
centroid_geom <- vect(
  matrix(c(mean(c(knepp_ext$xmin, knepp_ext$xmax)), 
           mean(c(knepp_ext$ymin, knepp_ext$ymax))), ncol=2), 
  crs = crs(Knepp_NIR_Timeseries)
)

# Project the centroid to WGS84 to get the latitude in degrees
centroid_latlon <- project(centroid_geom, "EPSG:4326")
knepp_lat_deg <- geom(centroid_latlon)[,"y"] # Will be ~50.9 degrees
knepp_lat_rad <- knepp_lat_deg * (pi / 180)  # Convert to radians

# Dynamically extract the dates directly from your raster layers
# Assumes names are like "2017-01_B08"
current_names <- names(Knepp_NIR_Timeseries)
dynamic_dates_string <- paste0(sub("_B08$", "", current_names), "-15")
dates <- as.Date(dynamic_dates_string, format = "%Y-%m-%d")

# Stop if dates failed to parse
if(any(is.na(dates))) stop("Date parsing failed. Check your layer names!")

day_of_year <- as.numeric(strftime(dates, format = "%j"))

# Calculate solar declination angle (delta) in radians
declination_deg <- 23.45 * sin((2 * pi / 365) * (day_of_year - 81))
declination_rad <- declination_deg * (pi / 180)

# Calculate hour angle (h) in radians
# Sentinel-2 descends at approx 10:30 AM local solar time. 
hour_angle_deg <- -1.5 * 15
hour_angle_rad <- hour_angle_deg * (pi / 180)

# Calculate final Solar Zenith Angle (theta_z) vector in radians
cos_theta_z <- sin(knepp_lat_rad) * sin(declination_rad) + 
  cos(knepp_lat_rad) * cos(declination_rad) * cos(hour_angle_rad)
sza_rad_vector <- acos(cos_theta_z)


## 2. Calculate DVI
message("Calculating Difference Vegetation Index (DVI)...")

# DVI is strictly NIR - Red
Knepp_Buffered_Timeseries.DVI <- Knepp_NIR_Timeseries - Knepp_Red_Timeseries


## 3. Define the Safe Hybrid Wrapper
message("Defining PPI Wrapper...")

calc_pixel_ppi <- function(dvi_ts, sza_vector) {
  
  # CRITICAL: Scale Sentinel-2 integers back to true float (0.0 to 1.0)
  # If you don't do this, the ppi package will crash on the massive integer values
  dvi_ts_scaled <- dvi_ts / 10000
  
  # Add a mathematical safety break to remove ecologically implausible extreme values and NA values
  # 1. Neutralize the NoData Asymptote: Force DVI <= 0 to NA
  dvi_ts_scaled[dvi_ts_scaled <= 0] <- NA_real_
  
  # 2. Neutralize the M=1 Asymptote: Clamp extreme values to 0.95
  dvi_ts_scaled[dvi_ts_scaled > 0.95] <- 0.95
  
  n_dates <- length(dvi_ts_scaled)
  
  # Identify valid (non-NA) time points
  valid_idx <- which(!is.na(dvi_ts_scaled))
  
  # If too few valid points, return all NA
  if (length(valid_idx) < 2) {
    return(rep(NA_real_, n_dates))
  }
  
  # Subset to valid observations only
  dvi_valid <- dvi_ts_scaled[valid_idx]
  sza_valid <- sza_vector[valid_idx]
  
  # Run MitMat's PPI on cleaned time series
  # Because we wrap this in tryCatch, if MitMat's code hits an edge-case 
  # mathematical error, it safely returns NA instead of killing your HPC job
  ppi_valid <- tryCatch({
    ppi::ppi(
      dvi = dvi_valid,
      zenith.angle = sza_valid,
      G = 0.2 # I set a lower value of Leaf Angle Distribution because most of the estate is grassland
    )
  }, error = function(e) {
    return(rep(NA_real_, length(dvi_valid)))
  })
  
  # Reconstruct full-length output (preserve time structure)
  ppi_full <- rep(NA_real_, n_dates)
  ppi_full[valid_idx] <- ppi_valid
  
  return(ppi_full)
}

## 4. Apply the function pixel-by-pixel
message("Calculating PPI timeseries across all pixels via MitMat's ppi package...")

# Notice how we pass sza_rad_vector as an explicit argument here!
Knepp_Buffered_Timeseries.PPI <- app(
  Knepp_Buffered_Timeseries.DVI, 
  fun = calc_pixel_ppi, 
  sza_vector = sza_rad_vector
)

## Clamp final PPI to ecologically valid extremes to prevent algorithmic asymptotes 
# This is necessary to avoid skewing the TWDTW cost matrices
Knepp_Buffered_Timeseries.PPI[Knepp_Buffered_Timeseries.PPI < -2] <- NA_real_
Knepp_Buffered_Timeseries.PPI[Knepp_Buffered_Timeseries.PPI > 10] <- NA_real_

# Assign names to the new PPI layers
names(Knepp_Buffered_Timeseries.PPI) <- paste0(sub("_B08$", "", names(Knepp_NIR_Timeseries)), "_PPI")

message("PPI Calculation Complete!")

# Export the raster for safe keeping

writeRaster(
  Knepp_Buffered_Timeseries.PPI, 
  filename = file.path(Knepp_Processed, "Knepp-Buffered_PPI_2019_2020_Timeseries.tif"), 
  overwrite = TRUE
)

# Load the PPI raster back in!

# Knepp_Buffered_Timeseries.PPI <- rast(file.path(Knepp_Processed, "Knepp-Buffered_PPI_2017_2019_2020_Timeseries.tif")) # Old raster from when 2017 was included
Knepp_Buffered_Timeseries.PPI <- rast(file.path(Knepp_Processed, "Knepp-Buffered_PPI_2019_2020_Timeseries.tif"))

### OPTION 1: Export the dataset as a NetCDF file for Gaussian Processing ####
## Saverio Vicario has a Python script for the Gaussian processing
## Export the data in the NetCDF format for use in the Python script
# Define the export path for the raw NetCDF

Knepp_NCDF_Raw_PPI_Path <- file.path(Knepp_Processed, "Knepp_Buffered_PPI_Raw.nc")

# Set the CRS to be exported as part of the NetCDF

NetCDF.CRS <- crs(Knepp_Buffered_Timeseries.PPI)

# Export the terra SpatRaster as a NetCDF file

message("Exporting raw PPI timeseries to NetCDF for Gaussian Process gap filling...")
terra::writeCDF(
  Knepp_Buffered_Timeseries.PPI, 
  filename = Knepp_NCDF_Raw_PPI_Path, 
  overwrite = TRUE, 
  varname = "PPI", 
  longname = "Normalized Difference Vegetation Index",
  gridmap = NetCDF.CRS,
  unit = "unitless"
)

### OPTION 2: Gap-filling and smoothing the dataset with Savitzky-Golay filter ####
## I need to drop the 2017 and 2024 data as it is missing valid springtime values

message("Excluding 2017 layers from the time series due to 100% missing spring data...")

# This uses grep to find the 2017 layers, invert the selection, and subset to just that selection (everything but 2017)

Knepp_Buffered_Timeseries.PPI <- Knepp_Buffered_Timeseries.PPI[[grep("2017|2024", names(Knepp_Buffered_Timeseries.PPI), invert = TRUE)]]

# Also, save the layer dates for other computations

Knepp.full.dates <- as.Date(dynamic_dates_string[-grep("2017|2024", dynamic_dates_string)], format = "%Y-%m-%d")

# Print the remaining layer names to verify that 2017 is removed

print("Remaining layers:")
print(names(Knepp_Buffered_Timeseries.PPI))

## Define a function to conduct gap filling
## This function will be applied to the temporal vector of every single pixel

sg_gapfill <- function(x) {
  
  # Check if the pixel is blank across all layers. If TRUE, return NAs to save time.
  
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
                               forder = 2) # forder: Filter order (polynomial degree), and 2 or 3 is standard for PPI
  
  return(x_smoothed)
}

## Applying the gap filling with Savitzky-Golay filter

message("Gap filling the dataset...")

# Apply the function across the SpatRaster

Knepp_Buffered_Timeseries.PPI.Cleaned <- app( # terra::app() to push the function through the Z-axis (time) of the raster
  Knepp_Buffered_Timeseries.PPI, 
  fun = sg_gapfill, 
  cores = detectCores() - 2 # Set cores as available
)

## Restore the original layer names (e.g., "2000-01_PPI")
# terra::app() stripped the layer names, so they're reassigned from the original object

names(Knepp_Buffered_Timeseries.PPI.Cleaned) <- names(Knepp_Buffered_Timeseries.PPI) # This only works if missing temporal slices have been inserted

message("Gap-filling complete. Ready for TWDTW analysis.")

## Write the raster to disk (and load back in if necessary)

writeRaster( # Write to disk
  Knepp_Buffered_Timeseries.PPI.Cleaned,
  filename = file.path(Knepp_Processed, "Knepp_Buffered_PPI_Cleaned_SG-method.tif"),
  overwrite = TRUE
)

Knepp_Buffered_Timeseries.PPI.Cleaned <- rast( # Load back in if necessary
  file.path(Knepp_Processed, "Knepp_Buffered_PPI_Cleaned_SG-method.tif")) # R environment doesn't save external C++ objects like rasters

# OPTIONAL: View the gap-free raster

print(Knepp_Buffered_Timeseries.PPI.Cleaned) # Summary

# 1. Calculate the global min and max across all 3 years and 12 months

tmp.global_range <- range(minmax(Knepp_Buffered_Timeseries.PPI.Cleaned), na.rm = TRUE)

for (i in 1:12) {
  
  # 2. Find all layers that match the current month (e.g., "-01_PPI")
  
  tmp.month_regex <- sprintf("-%02d_PPI", i)
  tmp.month_stack <- Knepp_Buffered_Timeseries.PPI.Cleaned[[grep(tmp.month_regex, names(Knepp_Buffered_Timeseries.PPI.Cleaned))]]
  
  # 3. Plot all years for this month side-by-side, locked to the global colour scale
  
  plot(
    tmp.month_stack, 
    range = tmp.global_range,
    main = names(tmp.month_stack)
  )
  
  # 4. Pause the loop so you can actually inspect the plot before it moves to the next month
  
  invisible(readline(prompt = paste("Displaying Month", sprintf("%02d", i), "- Press [Enter] in the console to view the next month...")))
}

# 5. Wipe all temporary variables to keep the workspace completely uncluttered

rm(i, tmp.global_range, tmp.month_regex, tmp.month_stack)

## Calculate the mean seasonal trajectory across the site ####
# Crop to Knepp Estate site boundaries

Knepp_Timeseries.PPI <- crop(
  Knepp_Buffered_Timeseries.PPI.Cleaned, 
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.PPI.Cleaned)),
  mask = TRUE
)

# Calculate the site's mean value for each month

Knepp_Timeseries_Mean.PPI <- global(Knepp_Timeseries.PPI, fun = "mean", na.rm = TRUE) 

# Calculate the mean for values outside the Knepp Estate
# Use inverse = TRUE to blank out the inside of the estate, keeping the surrounding farmland

Knepp_Environs_Timeseries.PPI <- mask(
  Knepp_Buffered_Timeseries.PPI.Cleaned, 
  project(KneppEstate_Boundaries, crs(Knepp_Buffered_Timeseries.PPI.Cleaned)), 
  inverse = TRUE
)

# Calculate the mean site value for each month for the outside farmland

Knepp_Environs_Timeseries.PPI <- global(Knepp_Environs_Timeseries.PPI, fun = "mean", na.rm = TRUE)

## Plot the seasonal PPI trajectories
# 1. Initiate the plot with the inside data

png(file.path(Knepp_Results, "Knepp_Estate_PPI_Mean_Timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(Knepp.full.dates, Knepp_Timeseries_Mean.PPI[,1], # Starts the plot
     type = "l",
     lwd = 2,
     col = "forestgreen",
     ylim = range(c(Knepp_Timeseries_Mean.PPI[,1], Knepp_Environs_Timeseries.PPI[,1]), na.rm = TRUE),
     xlab = "Date",
     ylab = "Mean PPI",
     main = "Knepp Estate vs. Surrounding Farmland - Mean PPI",
     xaxt = "n") 

# 2. Add the Outside data line

lines(Knepp.full.dates, Knepp_Environs_Timeseries.PPI[,1], 
      lwd = 2, 
      col = "chocolate")

# 3. Add horizontal and vertical gridlines

grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# 4. Add custom X-axis minor ticks (every month, no labels, shorter tick length)

axis.Date(1, at = Knepp.full.dates, labels = FALSE, tcl = -0.25)

# 5. Add custom X-axis major ticks (every year, with labels, longer tick length)
# This dynamically finds the unique years in your dataset to place the major labels

axis.Date(1, 
          at = as.Date(paste0(unique(format(Knepp.full.dates, "%Y")), "-01-15")), 
          format = "%Y", 
          tcl = -0.5)

# 6. Add a legend

legend("topright", 
       legend = c("Knepp Estate", "Environs surrounding Knepp"), 
       col = c("forestgreen", "chocolate"), 
       lwd = 2, 
       bty = "n") # 'n' removes the ugly bounding box around the legend

# Save the file

dev.off() # This actually exports the plot to the file

### Single year diversity analyses ####
## Select the years and dynamically set the years of interest for this analysis
# Extract year from the time vector

Knepp.years <- format(Knepp.full.dates, "%Y")

# Define Years of Interest (Abbreviation: YoI)

YoI <- unique(Knepp.years)

# Create list of PPI stacks per year

Knepp_PPI_YoI <- lapply(YoI, function(y) {
  Knepp_Buffered_Timeseries.PPI.Cleaned[[Knepp.years == y]]
})

# Name the list elements

names(Knepp_PPI_YoI) <- YoI

### Compute the Shannon-Wiener Index for each year ####
## Define a Shannon's H function

message("Calculating Shannon-Wiener diversity index...")

compute_shannon_year <- function(r_stack, year) {
  
  message("Processing Shannon's H for year: ", year)
  
  # Mean PPI (collapse time)
  
  mean_PPI <- app(r_stack, mean, na.rm = TRUE)
  
  # 1. Avoid numerical saturation by rounding to 2 decimal places to create classes
  # (PPI is already a float, so I'll skip the /10000 division from the NDVI script!)
  
  mean_PPI_2dec <- round(mean_PPI, 2)
  
  # 2. Calculate the Shannon-Wiener index
  
  shannon_mat <- rasterdiv::ShannonS(
    x = terra::as.matrix(mean_PPI_2dec, wide = TRUE),
    window = 3,
    na.tolerance = 0
  )
  
  # 3. Convert the matrix back to a raster
  
  shannon_rast <- rast(shannon_mat)
  ext(shannon_rast) <- ext(mean_PPI_2dec)
  crs(shannon_rast) <- crs(mean_PPI_2dec)
  names(shannon_rast) <- paste0("ShannonH_", year)
  
  # 4. Define an output file path
  
  out_path <- file.path(
    Knepp_Results,
    paste0("Knepp_", year,"_PPI_ShannonH", ".tif")
  )
  
  # 5. Write the raster
  
  writeRaster(
    shannon_rast,
    filename = out_path,
    overwrite = TRUE
  )
  
  return(shannon_rast)
}

## Apply the Shannon's H function to my rasters

Knepp_PPI_ShannonH_YoI <- mapply(
  compute_shannon_year,
  Knepp_PPI_YoI,
  names(Knepp_PPI_YoI),
  SIMPLIFY = FALSE
)

# Load them back in again if necessary

message("Loading Shannon's H rasters from disk...")

# Loop over the Years of Interest (YoI) using 'y' as the iterator

Knepp_PPI_ShannonH_YoI <- lapply(YoI, function(y) {
  
  # Ensure we use 'y' to construct the file string, not 'year'
  
  tmp.file.path <- file.path(Knepp_Results, paste0("Knepp_", y, "_PPI_ShannonH.tif"))
  
  # Import the raster
  
  tmp.imported.ShannonH.raster <- rast(tmp.file.path)
  
  return(tmp.imported.ShannonH.raster)
})

# Name the elements of the list so you can call them by year (e.g., Knepp_PPI_ShannonH_YoI[["2000"]])

names(Knepp_PPI_ShannonH_YoI) <- YoI

message("Shannon's H rasters successfully loaded.")

### Mosaic tiles for export and parallelisation on HPC (MaRC3a) ####
## Due to the size of the site and the need to run concurrently on 4 discrete years of data
## I will mosaic the dataset into 500 tiles (similar to the Kili analysis) for later processing on the HPC

Knepp_PPI_tile_info <- list()  # store metadata for later

for (y in YoI) {
  
  message("Creating tiles for year: ", y)
  
  # Get rasters
  
  tmp.r_stack <- Knepp_PPI_YoI[[y]] # This is the cleaned raster data (I checked)
  tmp.mean_raster <- app(tmp.r_stack, mean, na.rm = TRUE)
  
  # Trim outer NA borders (important for efficiency)
  
  tmp.trimmed_mean <- trim(tmp.mean_raster)
  
  ## Create a tiling grid
  
  knepp.tile.count <- 500  # same philosophy as Kili, but fewer tiles
  
  # Determine dimensions for tiling grid based on the aspect ratio of the raster
  
  knepp.aspect.ratio <- ncol(tmp.trimmed_mean) / nrow(tmp.trimmed_mean)
  
  knepp.tiling.factors <- expand.grid(
    ncols = 1:knepp.tile.count,
    nrows = 1:knepp.tile.count
  )
  
  knepp.tiling.factors <- knepp.tiling.factors[
    knepp.tiling.factors$ncols * knepp.tiling.factors$nrows == knepp.tile.count, ]
  
  knepp.tiling.factors$ratio_diff <- abs(
    (knepp.tiling.factors$ncols / knepp.tiling.factors$nrows) - knepp.aspect.ratio
  )
  
  knepp.tile.target.size <- knepp.tiling.factors[which.min(knepp.tiling.factors$ratio_diff), ]
  
  knepp.cols <- knepp.tile.target.size$ncols
  knepp.rows <- knepp.tile.target.size$nrows
  
  # Generate tiling grid as polygons (important for later spatial operations)
  
  knepp.tiling.grid <- as.polygons(
    rast(
      ext(tmp.trimmed_mean),
      ncols = knepp.cols,
      nrows = knepp.rows,
      crs = crs(tmp.trimmed_mean)
    )
  )
  
  # Save the grids (important for reproducibility + debugging)
  
  knepp.tiles.filepath <- file.path(Knepp_Processed, paste0("Knepp_Tiling_Grid_", y, ".geojson"))
  
  writeVector(knepp.tiling.grid, knepp.tiles.filepath, filetype = "GeoJSON", overwrite = TRUE)
  
  knepp.tiling.grid <- vect(knepp.tiles.filepath)
  
  ## Create tile overplotting margin so that the moving window of Rao's Q can be computed without edge effects
  
  RaoQ.window.size <- 3 # In `paRao` function, ensure moving window size is equal to this value!
  tmp.tile.overlap <- floor(RaoQ.window.size / 2)
  
  ## Define output directories
  
  tmp.knepp.mean.tile.dir <- file.path(Knepp_Processed, paste0("Knepp_", y, "_MeanPPI_Tiles"))
  tmp.knepp.ts.tile.dir   <- file.path(Knepp_Processed, paste0("Knepp_", y, "_PPI_Timeseries_Tiles"))
  
  dir.create(tmp.knepp.mean.tile.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp.knepp.ts.tile.dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Save the time vector (for TWDTW)
  
  tmp.time.vector <- Knepp_Time[Knepp.years == y]
  
  writeLines(
    as.character(tmp.time.vector),
    file.path(tmp.knepp.ts.tile.dir, "time_vector.txt")
  )
  
  ## Write the tiles to disk
  
  message("Tiling PPI MEAN raster for year: ", y)
  
  makeTiles(
    tmp.trimmed_mean,
    y = knepp.tiling.grid,
    buffer = tmp.tile.overlap,
    filename = file.path(tmp.knepp.mean.tile.dir, paste0("Knepp_", y, "_Mean_Tile-.tif")),
    overwrite = TRUE
  )
  
  message("Tiling PPI TIMESERIES raster for year: ", y)
  
  makeTiles(
    tmp.r_stack,
    y = knepp.tiling.grid,
    buffer = tmp.tile.overlap,
    filename = file.path(tmp.knepp.ts.tile.dir, paste0("Knepp_", y, "_TS_Tile-.tif")),
    overwrite = TRUE
  )
  
  ## Store the metadata for later :-)
  
  Knepp_PPI_tile_info[[y]] <- list(
    cols = knepp.cols,
    rows = knepp.rows,
    n_tiles = knepp.tile.count,
    mean_dir = tmp.knepp.mean.tile.dir,
    ts_dir = tmp.knepp.ts.tile.dir
  )
}

### Demosaic Knepp Estate tiles for further analysis and visualisation ####
## First, I will define a function that can demosaic tiles

Knepp_Demosaic_RaoTiles <- function(tile_dir, output_file) {
  
  message("Reading tiles from: ", tile_dir)
  
  tile_files <- list.files(tile_dir, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tile_files) == 0) {
    stop("No tiles found in: ", tile_dir)
  }
  
  # Read all tiles
  
  tile_list <- lapply(tile_files, rast)
  
  message("Mosaicing ", length(tile_list), " tiles...")
  
  # Merge tiles (mean handles overlap nicely)
  
  mosaiced <- do.call(mosaic, c(tile_list, fun = "mean"))
  
  message("Trimming outer NA borders...")
  
  mosaiced <- trim(mosaiced)
  
  message("Writing output...")
  
  writeRaster(mosaiced, output_file, overwrite = TRUE)
  
  message("Done: ", output_file)
  
  return(mosaiced)
}

## Secondly, I will apply this function to the tiles for each year and Rao's Q method

Knepp_Rao_Outputs_Dir <- file.path(Knepp_Processed) # Tell R where to find the files

# Now demosaic the classical Rao's Q tiles

for (y in YoI) {
  
  message("Demosaicing Classical Rao for year: ", y)
  
  tile_dir <- file.path(
    Knepp_Rao_Outputs_Dir,
    paste0("Knepp_", y, "_MeanPPI_Tiles"),
    "Rao-utputs"
  )
  
  out_file <- file.path(
    Knepp_Results,
    paste0("Knepp_", y, "_PPI_ClassicRao.tif")
  )
  
  Knepp_Demosaic_RaoTiles(tile_dir, out_file)
}

# And now demosaic the TWDTW Rao's Q tiles

for (y in YoI) {
  
  message("Demosaicing TWDTW Rao for year: ", y)
  
  tile_dir <- file.path(
    Knepp_Rao_Outputs_Dir,
    paste0("Knepp_", y, "_PPI_Timeseries_Tiles"),
    "Rao-utputs"
  )
  
  out_file <- file.path(
    Knepp_Results,
    paste0("Knepp_", y, "_PPI_TWDTW_Rao.tif")
  )
  
  Knepp_Demosaic_RaoTiles(tile_dir, out_file)
}

## As both of these are saved externally, I now need to load them back in to my R environment
# Classic Rao's Q

Knepp_PPI_Classic_Rao_YoI <- lapply(YoI, function(y) {
  rast(file.path(Knepp_Results, paste0("Knepp_", y, "_PPI_ClassicRao.tif")))
})

names(Knepp_PPI_Classic_Rao_YoI) <- YoI # Rename each layer to the corresponding year

# TWDTW Rao's Q

Knepp_PPI_TWDTW_Rao_YoI <- lapply(YoI, function(y) {
  rast(file.path(Knepp_Results, paste0("Knepp_", y, "_PPI_TWDTW_Rao.tif")))
})

names(Knepp_PPI_TWDTW_Rao_YoI) <- YoI

### Comparison to land cover data ####
## Load and process the land cover maps from the UK's Centre for Ecology and Hydrology
# Load the spatial rasters

#Knepp_CEH_LandCover_2017 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2017/data/LCM_2017.tif"))
Knepp_CEH_LandCover_2019 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2019/data/LCM_2019.tif"))
Knepp_CEH_LandCover_2020 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2020/data/LCM2020.tif"))
#Knepp_CEH_LandCover_2024 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2024/data/LCM2024.tif"))

# Put them in a list for convenience

Knepp_CEH_LandCover_YoI <- list(
  #"2017" = Knepp_CEH_LandCover_2017,
  "2019" = Knepp_CEH_LandCover_2019,
  "2020" = Knepp_CEH_LandCover_2020#,
  #"2024" = Knepp_CEH_LandCover_2024
)

message("Converting numeric land cover classes to descriptive UKCEH names...")

# Define the UKCEH Land Cover Classes dictionary

UKCEH_LC_Classes <- c( # These are only applicable to 2017 - 2020!
  "1" = "Broadleaved woodland",
  "2" = "Coniferous woodland",
  "3" = "Arable and horticulture",
  "4" = "Improved grassland",
  "5" = "Neutral grassland",
  "6" = "Calcareous grassland",
  "7" = "Acid grassland",
  "8" = "Fen, marsh and swamp",
  "9" = "Heather",
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

## Combine all my indices into one object for convenience

Knepp_Indices_YoI.PPI_derived <- lapply(YoI, function(y) {
  
  message("Rebuilding index object for year: ", y)
  
  list(
    Knepp_Mean_PPI_Values = app(Knepp_PPI_YoI[[y]], mean, na.rm = TRUE),
    Knepp_PPI_ShannonsH = Knepp_PPI_ShannonH_YoI[[y]],
    Knepp_PPI_Classic_RaoQ = Knepp_PPI_Classic_Rao_YoI[[y]],
    Knepp_PPI_TWDTW_RaoQ = Knepp_PPI_TWDTW_Rao_YoI[[y]]
  )
})

names(Knepp_Indices_YoI.PPI_derived) <- YoI

# Make a handy little comparison plot!

par(mfrow=c(2,2)) # Make a 2*2 multipanel plot

# Plot them!

for (y in YoI) {
  plot(Knepp_Indices_YoI.PPI_derived[[y]][[1]], 
       main = paste(y, names(Knepp_Indices_YoI.PPI_derived[[y]][1])))
  plot(Knepp_Indices_YoI.PPI_derived[[y]][[2]], 
       main = paste(y, names(Knepp_Indices_YoI.PPI_derived[[y]][2])))
  plot(Knepp_Indices_YoI.PPI_derived[[y]][[3]], 
       main = paste(y, names(Knepp_Indices_YoI.PPI_derived[[y]][3])))
  plot(Knepp_Indices_YoI.PPI_derived[[y]][[4]], 
       main = paste(y, names(Knepp_Indices_YoI.PPI_derived[[y]][4])))
}

## Align the computed indices to the CEH land cover data, and name the layers accordingly
# First, a helper function to compare and reproject CRSs

align_to_lc <- function(r, lc) {
  project(r, lc)  # always reproject so the CRSs match
}

# This function aligns my indices to the land cover data and stacks them for analysis

Knepp_LandCover_AlignedIndices_YoI.PPI <- lapply(names(Knepp_Indices_YoI.PPI_derived), function(y) {
  
  message("Aligning year: ", y)
  
  idx_list <- Knepp_Indices_YoI.PPI_derived[[y]]
  lc <- Knepp_CEH_LandCover_YoI[[y]][[1]]
  
  # Force CRS values to match
  
  mean_ppi <- project(idx_list$Knepp_Mean_PPI_Values, lc)
  shannon  <- project(idx_list$Knepp_PPI_ShannonsH, lc)
  classic  <- project(idx_list$Knepp_PPI_Classic_RaoQ, lc)
  twdtw    <- project(idx_list$Knepp_PPI_TWDTW_RaoQ, lc)
  
  # Resample to match the land cover grid
  
  mean_ppi_r <- resample(mean_ppi, lc, method = "bilinear")[[1]]
  shannon_r  <- resample(shannon, lc, method = "bilinear")[[1]]
  classic_r  <- resample(classic, lc, method = "bilinear")[[1]]
  twdtw_r    <- resample(twdtw,   lc, method = "bilinear")[[1]]
  
  aligned_stack <- c(mean_ppi_r, shannon_r, classic_r, twdtw_r, lc)
  
  names(aligned_stack) <- c(
    "Mean_PPI",
    "ShannonH",
    "RaosQ_Classic",
    "RaosQ_TWDTW",
    "LandCover"
  )
  
  ## Assign the categorical levels to the LandCover layer
  # terra expects a data.frame where the first column is the ID and the second is the label
    
  levels(aligned_stack$LandCover) <- data.frame(
    id = as.numeric(names(UKCEH_LC_Classes)),
    cover = as.character(UKCEH_LC_Classes)
    )
  
  return(aligned_stack)
})

names(Knepp_LandCover_AlignedIndices_YoI.PPI) <- names(Knepp_Indices_YoI.PPI_derived)

## Now I will crop and mask the land cover rasters to just the areas in my analysis
# I'm using the buffered site so I can compare inside the estate and out

Knepp_LandCover_AlignedIndices_YoI.PPI <- lapply(
  Knepp_LandCover_AlignedIndices_YoI.PPI,
  function(stack) {
    
    # Align polygon CRS to raster
    
    boundary_aligned <- project(KneppEstate_Buffered_Boundaries, stack)
    
    cropped <- crop(stack, boundary_aligned)
    masked  <- mask(cropped, boundary_aligned)
    
    return(masked)
  }
)

## OPTIONAL: export these rasters so I don't have to recreate them every time I run the script
# My object is a list of rasters, but not itself a raster, so  I need to export each year separately

for (y in names(Knepp_LandCover_AlignedIndices_YoI.PPI)) {
  
  r <- Knepp_LandCover_AlignedIndices_YoI.PPI[[y]]
  
  out_file <- file.path(
    Knepp_Results,
    paste0("Knepp_", y, "_Sentinel_Aligned_Indices_PPI.tif")
  )
  
  writeRaster(
    r,
    filename = out_file,
    overwrite = TRUE
  )
}

## Load the raster back in (if necessary)

Knepp_LandCover_AlignedIndices_YoI.PPI <- lapply(YoI, function(y) {
  
  rast(
    file.path(Knepp_Results, paste0("Knepp_", y, "_Sentinel_Aligned_Indices_PPI.tif"))
  )
})

names(Knepp_LandCover_AlignedIndices_YoI.PPI) <- YoI

## Now I convert the raster stacks to a dataframe so I can run analyses on them
# The dataframe will have one row per pixel, and columns for each index and the land cover class

Knepp_DF_YoI.PPI <- lapply(
  Knepp_LandCover_AlignedIndices_YoI.PPI,
  function(stack) {
    
    df <- as.data.frame(stack, na.rm = TRUE)
    
    colnames(df) <- c(
      "PPI_Yearly_Mean",
      "ShannonH",
      "RaosQ_Classic",
      "RaosQ_TWDTW",
      "LandCover"
    )
    
    # Ensure land cover is treated as categorical
    
    #df$LandCover <- as.factor(df$LandCover) # Not required anymore as using named LC types
    
    return(df)
  }
)

## Investigate frequency distribution of land cover types
# I thought it'd be interesting to identify what land cover types are present
# How they change for each year
# And understand the variance of the land cover classes better

View(data.frame( # View it in R
  Habitat = names(table(Knepp_DF_YoI.PPI[["2019"]]$LandCover)),
  Pixels_2019 = as.numeric(table(Knepp_DF_YoI.PPI[["2019"]]$LandCover)),
  Pixels_2020 = as.numeric(table(Knepp_DF_YoI.PPI[["2020"]]$LandCover))
  ))

write.csv(data.frame( # Save it externally as a .CSV
  Habitat = names(table(Knepp_DF_YoI.PPI[["2019"]]$LandCover)),
  Pixels_2019 = as.numeric(table(Knepp_DF_YoI.PPI[["2019"]]$LandCover)),
  Pixels_2020 = as.numeric(table(Knepp_DF_YoI.PPI[["2020"]]$LandCover))
  ), file.path(Knepp_Results,"KneppBuffered_LandCover_FrequencyDistribution.csv"))

## Subsample the dataframe
# Take a sample of 10,000 rows to prevent distance-matrix memory crashes
# Many land cover classes are ecologically insignificant in this dataset and reduce model performance
# Therefore, I'm specifying the important land cover types and then subsetting to include those only

target_classes <- c(
  "Broadleaved woodland", 
  "Coniferous woodland",
  "Arable and horticulture", 
  "Improved grassland", 
  "Urban",
  "Suburban" 
) # Add or remove classes as required

# Dynamically calculate the sample size per class

total_target_pixels <- 10000
n_classes <- length(target_classes)

# Use floor() to round down so we don't accidentally ask for a fraction of a pixel

pixels_per_class <- floor(total_target_pixels / n_classes) 

message("Sampling up to ", pixels_per_class, " pixels across ", n_classes, " classes.")

## And now we loop over the dataframe to run the PERMANOVAe upon it

Knepp_PERMANOVA_YoI.PPI <- lapply(
  names(Knepp_DF_YoI.PPI),
  function(y) {
    
    df <- Knepp_DF_YoI.PPI[[y]]
    
    # 3. Filter out the fringe/empty classes
    
    df_filtered <- df[df$LandCover %in% target_classes, ]
    
    # 4. Stratified Subsampling using dynamic values
    
    set.seed(123)
    df_subset <- df_filtered %>%
      group_by(LandCover) %>%
      # slice_sample automatically handles groups that are smaller than pixels_per_class!
      slice_sample(n = pixels_per_class, replace = FALSE) %>% 
      ungroup() %>%
      as.data.frame()
    
    list(
      
      # These are structured in the same way as the PERMANOVAe from the other case studies
      
      Shannon = adonis2(
        df_subset$ShannonH ~ df_subset$LandCover,
        #data = df_subset, # IDK why, but R cannot find ShannonH unless defined by a $
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_Classic = adonis2(
        df_subset$RaosQ_Classic ~ df_subset$LandCover,
        #data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_TWDTW = adonis2(
        df_subset$RaosQ_TWDTW ~ df_subset$LandCover,
        #data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      )
    )
  }
)

saveRDS(Knepp_PERMANOVA_YoI.PPI, file = file.path(Knepp_Results, "Knepp_PERMANOVA_YoI_PPI.rds"))

# Finally, putting the results into their own dataframe

Knepp_PERMANOVA_Summary.PPI <- do.call(rbind, lapply(
  names(Knepp_PERMANOVA_YoI.PPI),
  function(y) {
    
    res <- Knepp_PERMANOVA_YoI.PPI[[y]]
    
    data.frame(
      Year = y,
      Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
      R2 = c(
        res$Shannon$R2[1],
        res$Rao_Classic$R2[1],
        res$Rao_TWDTW$R2[1]
      ),
      F = c(
        res$Shannon$F[1],
        res$Rao_Classic$F[1],
        res$Rao_TWDTW$F[1]
      ),
      p_value = c(
        res$Shannon$`Pr(>F)`[1],
        res$Rao_Classic$`Pr(>F)`[1],
        res$Rao_TWDTW$`Pr(>F)`[1]
      )
    )
  }
))

saveRDS(Knepp_PERMANOVA_Summary.PPI, file = file.path(Knepp_Results, "Knepp_PERMANOVA_Summary_PPI.rds"))

message("Knepp Estate PPI analysis complete.") # End of script ####
