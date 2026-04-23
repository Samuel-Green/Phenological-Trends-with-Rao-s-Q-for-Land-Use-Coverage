############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025                    #
# 05B_Analyse_Knepp_Estate_NDVI.R                                              #
# Comparative NDVI analysis of TWDTW Rao's Q & classic Rao's Q in Knepp Estate #
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

## Import and stack all of my Knepp Estate GeoTIFFs into one object
# use the `rast` function of terra

Knepp_Buffered_Timeseries <- rast(list.files(Knepp_Input, pattern = ".tif$", full.names = TRUE))

# Get file paths

Knepp.input.files <- list.files(Knepp_Input, pattern = ".tif$", full.names = TRUE)

# Extract YYYY_MM from filenames

Knepp.dates <- sub(".*_(\\d{4}_\\d{2})\\.tif$", "\\1", basename(Knepp.input.files))

# Convert to YYYY-MM format

Knepp.dates <- gsub("_", "-", Knepp.dates)

# If the spectral band order changes, then downstream code will break!

Knepp.spectral.bands <- c("Blue", "Green", "Red", "NIR", "NDVI", "Solar Zenith Angle", "Solar Azimuth Angle")

# Repeat each date for each band

knew.Knepp.knames <- as.vector(sapply(Knepp.dates, function(d) paste0(d, "_", Knepp.spectral.bands)))

# Assign names to each layer of the raster

names(Knepp_Buffered_Timeseries) <- knew.Knepp.knames

### Inspect temporal structure ####
## rasterdiv::paRao() requires an explicit time_vector
## This must correspond EXACTLY to the layer order

Knepp.full.dates <- as.Date(paste0(Knepp.dates, "-15"))
str(Knepp.full.dates) # WARNING: Vector is shorter than length as raster stack includes multiple bands for each date

## Subset and crop the raster
# Remove extraneous layers (as this is an NDVI analysis)

Knepp_Buffered_Timeseries.NDVI <- Knepp_Buffered_Timeseries[[grep("_NDVI$", names(Knepp_Buffered_Timeseries))]]
Knepp_Buffered_Timeseries.NDVI

# Import the shapefile with the Knepp Estate's borders

KneppEstate_Boundaries <- vect(file.path(Knepp_Input,"Schulte-to-Bühne_Knepp_Site_Data","knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed,"Knepp_Estate_Buffered_Boundary.shp"))

# Check that the CRS matches the timeseries raster

if (crs(Knepp_Buffered_Timeseries.NDVI) != crs(KneppEstate_Boundaries)){
  Knepp_Buffered_Timeseries.NDVI <- project(Knepp_Buffered_Timeseries.NDVI, crs(KneppEstate_Boundaries))
  message("The raster's CRS differs from the shape file's CRS")
} else {message("The raster's CRS already matches the shapefile's CRS")}

### OPTION 1: Export the dataset as a NetCDF file for Gaussian Processing ####
## Saverio Vicario has a Python script for the Gaussian processing
## Export the data in the NetCDF format for use in the Python script
# Define the export path for the raw NetCDF

Knepp_NCDF_Raw_NDVI_Path <- file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Raw.nc")

# Set the CRS to be exported as part of the NetCDF

NetCDF.CRS <- "EPSG:32630"

# Export the terra SpatRaster as a NetCDF file

message("Exporting raw NDVI timeseries to NetCDF for Gaussian Process gap filling...")
terra::writeCDF(
  Knepp_Buffered_Timeseries.NDVI, 
  filename = Knepp_NCDF_Raw_NDVI_Path, 
  overwrite = TRUE, 
  varname = "NDVI", 
  longname = "Normalized Difference Vegetation Index",
  gridmap = NetCDF.CRS,
  unit = "unitless"
)

### OPTION 2: Gap-filling and smoothing the dataset with Savitzky-Golay filter ####
## Define a function to conduct gap filling
## This function will be applied to the temporal vector of every single pixel

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
                               fl = 11, # fl = Filter length, must be an odd number, and `fl = 11` provides a ~yearly smoothing window, preserving the seasonality
                               forder = 2) # forder: Filter order (polynomial degree), and 2 or 3 is standard for NDVI
  
  return(x_smoothed)
}

message("Checking for missing months in the temporal sequence...")

## Repair completely missing months by inserting blank layers filled with NAs
# Extract dates from current layer names using "-15" as a mid-month placeholder

current.names <- names(Knepp_Buffered_Timeseries.NDVI)
current.dates <- as.Date(paste0(sub("_NDVI$", "", current.names), "-15")) 

# Generate a perfectly continuous monthly sequence from start to end

Knepp.full.dates <- seq(min(current.dates), max(current.dates), by = "month")

# Identify any dates that are completely missing from the raster stack

missing.dates <- Knepp.full.dates[!Knepp.full.dates %in% current.dates]

if (length(missing.dates) > 0) {
  message(paste("Found", length(missing.dates), "missing months. Creating blank template layers..."))
  
  # Create a single blank raster (filled with NAs) using the first layer as a spatial template
  
  Empty_Knepp_Buffered_Raster <- terra::init(Knepp_Buffered_Timeseries.NDVI[[1]], NA)
  
  # Duplicate this blank template for every missing month and assign correct names
  
  missing.rasters <- terra::rast(replicate(length(missing.dates), Empty_Knepp_Buffered_Raster))
  names(missing.rasters) <- paste0(format(missing.dates, "%Y-%m"), "_NDVI")
  
  # Append the blank layers to the original raster stack
  
  Knepp_Buffered_Timeseries.NDVI <- c(Knepp_Buffered_Timeseries.NDVI, missing.rasters)
  
  # Sort the entire stack chronologically so the NAs are in the right position for the filter
  
  sorted.names <- paste0(format(Knepp.full.dates, "%Y-%m"), "_NDVI")
  Knepp_Buffered_Timeseries.NDVI <- Knepp_Buffered_Timeseries.NDVI[[sorted.names]]
  
} else {
  message("No missing months found. Temporal sequence is already contiguous.")
}

## Applying the gap filling with Savitzky-Golay filter

message("Gap filling the dataset...")

# Apply the function across the SpatRaster

Knepp_Buffered_Timeseries.NDVI.Cleaned <- app( # terra::app() to push the function through the Z-axis (time) of the raster
  Knepp_Buffered_Timeseries.NDVI, 
  fun = sg_gapfill, 
  cores = detectCores() - 2 # Set cores as available
)

## Restore the original layer names (e.g., "2000-01_NDVI")
# terra::app() stripped the layer names, so they're reassigned from the original object

names(Knepp_Buffered_Timeseries.NDVI.Cleaned) <- names(Knepp_Buffered_Timeseries.NDVI)

message("Gap-filling complete. Ready for TWDTW analysis.")

## Write the raster to disk (and load back in if necessary)

writeRaster( # Write to disk
  Knepp_Buffered_Timeseries.NDVI.Cleaned,
  filename = file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Cleaned_SG-method.tif"),
  overwrite = TRUE
)
Knepp_Buffered_Timeseries.NDVI.Cleaned <- rast( # Load back in if necessary
  file.path(Knepp_Processed, "Knepp_Buffered_NDVI_Cleaned_SG-method.tif")) # R environment doesn't save external C++ objects like rasters

# OPTIONAL: View the gap-free raster

print(Knepp_Buffered_Timeseries.NDVI.Cleaned) # Summary

for (i in 1:nlyr(Knepp_Buffered_Timeseries.NDVI.Cleaned)) { # Sequentially plot each individual raster layer (this takes a while)
  plot(Knepp_Buffered_Timeseries.NDVI.Cleaned[[i]], main = names(Knepp_Buffered_Timeseries.NDVI.Cleaned)[i])
}

## Calculate the mean seasonal trajectory across the site ####
# Crop to Knepp Estate site boundaries

Knepp_Timeseries.NDVI <- crop(Knepp_Buffered_Timeseries.NDVI.Cleaned, KneppEstate_Boundaries, mask = TRUE)

# Calculate the site's mean value for each month

Knepp_Timeseries_Mean.NDVI <- global(Knepp_Timeseries.NDVI, fun = "mean", na.rm = TRUE) 

png(file.path(Knepp_Results, "Knepp_Estate_NDVI_Mean_Timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(Knepp.full.dates, Knepp_Timeseries_Mean.NDVI[,1],
     type = "l",
     lwd = 2,
     xlab = "Date",
     ylab = "Mean NDVI",
     main = "Knepp Estate – Mean NDVI Time Series")

dev.off() # This actually exports the plot to the file

### Single year diversity analyses ####
## I have CEH land-cover data for 2000, 2007, 2015, and 2020
## To fairly compare the indices, I'm subsetting the spatial raster to those years
# Extract year from the time vector

Knepp.years <- format(Knepp.full.dates, "%Y")

# Define Years of Interest (Abbreviation: YoI)

YoI <- c("2000", "2007", "2015", "2020")

# Create list of NDVI stacks per year

Knepp_NDVI_YoI <- lapply(YoI, function(y) {
  Knepp_Buffered_Timeseries.NDVI.Cleaned[[Knepp.years == y]]
})

# Name the list elements

names(Knepp_NDVI_YoI) <- YoI

### Compute the Shannon-Wiener Index for each year ####
## Define a Shannon's H function

message("Calculating Shannon-Wiener diversity index...")

compute_shannon_year <- function(r_stack, year) {
  
  ## Restructure the raster into a simple data matrix
  
  message("Processing Shannon's H for year: ", year)
  
  # Mean NDVI (collapse time)
  
  mean_ndvi <- app(r_stack, mean, na.rm = TRUE)
  
  # Avoid numerical saturation by rounding to 2 decimal places
  
  mean_ndvi_2dec <- round(mean_ndvi, 2)
  
  ## Calculate the Shannon-Wiener index
  
  shannon_mat <- rasterdiv::ShannonS(
    x = terra::as.matrix(mean_ndvi_2dec, wide = TRUE),
    window = 3,
    na.tolerance = 0
  )
  
  ## Convert the matrix back to a raster
  
  shannon_rast <- rast(shannon_mat)
  ext(shannon_rast) <- ext(mean_ndvi_2dec)
  crs(shannon_rast) <- crs(mean_ndvi_2dec)
  names(shannon_rast) <- paste0("ShannonH_", year)
  
  ## Save the raster
  # Define an output file path
  
  out_path <- file.path(
    Knepp_Results,
    paste0("Knepp_", year,"_NDVI_ShannonH", ".tif")
  )
  
  # Write the raster
  
  writeRaster(
    shannon_rast,
    filename = out_path,
    overwrite = TRUE
  )
  
  return(shannon_rast)
}

## Apply the Shannon's H function to my rasters

Knepp_NDVI_ShannonH_YoI <- mapply(
  compute_shannon_year,
  Knepp_NDVI_YoI,
  names(Knepp_NDVI_YoI),
  SIMPLIFY = FALSE
)

# Load them back in again if necessary

message("Loading Shannon's H rasters from disk...")

# Loop over the Years of Interest (YoI) using 'y' as the iterator

Knepp_NDVI_ShannonH_YoI <- lapply(YoI, function(y) {
  
  # Ensure we use 'y' to construct the file string, not 'year'
  
  tmp.file.path <- file.path(Knepp_Results, paste0("Knepp_", y, "_NDVI_ShannonH.tif"))
  
  # Import the raster
  
  tmp.imported.ShannonH.raster <- rast(tmp.file.path)
  
  return(tmp.imported.ShannonH.raster)
})

# Name the elements of the list so you can call them by year (e.g., Knepp_NDVI_ShannonH_YoI[["2000"]])

names(Knepp_NDVI_ShannonH_YoI) <- YoI

message("Shannon's H rasters successfully loaded.")

### Mosaic tiles for export and parallelisation on HPC (MaRC3a) ####
## Due to the size of the site and the need to run concurrently on 4 discrete years of data
## I will mosaic the dataset into 1000 tiles (similar to the Kili analysis) for later processing on the HPC

Knepp_NDVI_tile_info <- list()  # store metadata for later

for (y in YoI) {
  
  message("Creating tiles for year: ", y)
  
  # Get rasters
  tmp.r_stack <- Knepp_NDVI_YoI[[y]] # This is the original, uncleaned stack. It will be reactivated when I have gap-filled data
  # tmp.r_stack <- Knepp_NDVI_YoI_Clean[[y]] # This is the cleaned stack, with only pixels with complete data across all layers
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
  
  tmp.knepp.mean.tile.dir <- file.path(Knepp_Processed, paste0("Knepp_", y, "_MeanNDVI_Tiles"))
  tmp.knepp.ts.tile.dir   <- file.path(Knepp_Processed, paste0("Knepp_", y, "_NDVI_Timeseries_Tiles"))
  
  dir.create(tmp.knepp.mean.tile.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp.knepp.ts.tile.dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Save the time vector (for TWDTW)
  
  tmp.time.vector <- Knepp_Time[Knepp.years == y]
  
  writeLines(
    as.character(tmp.time.vector),
    file.path(tmp.knepp.ts.tile.dir, "time_vector.txt")
  )
  
  ## Write the tiles to disk
  
  message("Tiling NDVI MEAN raster for year: ", y)
  
  makeTiles(
    tmp.trimmed_mean,
    y = knepp.tiling.grid,
    buffer = tmp.tile.overlap,
    filename = file.path(tmp.knepp.mean.tile.dir, paste0("Knepp_", y, "_Mean_Tile-.tif")),
    overwrite = TRUE
  )
  
  message("Tiling NDVI TIMESERIES raster for year: ", y)
  
  makeTiles(
    tmp.r_stack,
    y = knepp.tiling.grid,
    buffer = tmp.tile.overlap,
    filename = file.path(tmp.knepp.ts.tile.dir, paste0("Knepp_", y, "_TS_Tile-.tif")),
    overwrite = TRUE
  )
  
  ## Store the metadata for later :-)
  
  Knepp_NDVI_tile_info[[y]] <- list(
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
    paste0("Knepp_", y, "_MeanNDVI_Tiles"),
    "Rao-utputs"
  )
  
  out_file <- file.path(
    Knepp_Results,
    paste0("Knepp_", y, "_NDVI_ClassicRao.tif")
  )
  
  Knepp_Demosaic_RaoTiles(tile_dir, out_file)
}

# And now demosaic the TWDTW Rao's Q tiles

for (y in YoI) {
  
  message("Demosaicing TWDTW Rao for year: ", y)
  
  tile_dir <- file.path(
    Knepp_Rao_Outputs_Dir,
    paste0("Knepp_", y, "_NDVI_Timeseries_Tiles"),
    "Rao-utputs"
  )
  
  out_file <- file.path(
    Knepp_Results,
    paste0("Knepp_", y, "_NDVI_TWDTW_Rao.tif")
  )
  
  Knepp_Demosaic_RaoTiles(tile_dir, out_file)
}

## As both of these are saved externally, I now need to load them back in to my R environment
# Classic Rao's Q

Knepp_NDVI_Classic_Rao_YoI <- lapply(YoI, function(y) {
  rast(file.path(Knepp_Results, paste0("Knepp_", y, "_NDVI_ClassicRao.tif")))
})

names(Knepp_NDVI_Classic_Rao_YoI) <- YoI # Rename each layer to the corresponding year

# TWDTW Rao's Q

Knepp_NDVI_TWDTW_Rao_YoI <- lapply(YoI, function(y) {
  rast(file.path(Knepp_Results, paste0("Knepp_", y, "_NDVI_TWDTW_Rao.tif")))
})

names(Knepp_NDVI_TWDTW_Rao_YoI) <- YoI

### Comparison to land cover data ####
## Load and process the land cover maps from the UK's Centre for Ecology and Hydrology
# Load the spatial rasters

Knepp_CEH_LandCover_2000 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2000/data/LCM2000.tif"))
Knepp_CEH_LandCover_2007 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2007/data/LCM2007.tif"))
Knepp_CEH_LandCover_2015 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2015/data/LCM2015.tif"))
Knepp_CEH_LandCover_2020 <- rast(file.path(Knepp_Input, "CEH Land-cover data/Year_2020/data/LCM2020.tif"))

# Put them in a list for convenience

Knepp_CEH_LandCover_YoI <- list(
  "2000" = Knepp_CEH_LandCover_2000,
  "2007" = Knepp_CEH_LandCover_2007,
  "2015" = Knepp_CEH_LandCover_2015,
  "2020" = Knepp_CEH_LandCover_2020
)

## Combine all my indices into one object for convenience

Knepp_Indices_YoI.NDVI_derived <- lapply(YoI, function(y) {
  
  message("Rebuilding index object for year: ", y)
  
  list(
    Knepp_NDVI_ShannonsH = Knepp_NDVI_ShannonH_YoI[[y]],
    Knepp_NDVI_Classic_RaoQ = Knepp_NDVI_Classic_Rao_YoI[[y]],
    Knepp_NDVI_TWDTW_RaoQ = Knepp_NDVI_TWDTW_Rao_YoI[[y]]
  )
})

names(Knepp_Indices_YoI.NDVI_derived) <- YoI

## Align the computed indices to the CEH land cover data, and name the layers accordingly
# First, a helper function to compare and reproject CRSs

align_to_lc <- function(r, lc) {
  project(r, lc)  # always reproject to the CRSs to match
}

# This function aligns my indices to the land cover data and stacks them for analysis

Knepp_LandCover_AlignedIndices_YoI.NDVI <- lapply(names(Knepp_Indices_YoI.NDVI_derived), function(y) {
  
  message("Aligning year: ", y)
  
  idx_list <- Knepp_Indices_YoI.NDVI_derived[[y]]
  lc <- Knepp_CEH_LandCover_YoI[[y]][[1]]
  
  # FORCE CRS match
  
  shannon <- project(idx_list$Knepp_NDVI_ShannonsH, lc)
  classic <- project(idx_list$Knepp_NDVI_Classic_RaoQ, lc)
  twdtw   <- project(idx_list$Knepp_NDVI_TWDTW_RaoQ, lc)
  
  # Resample to LC grid
  
  shannon_r <- resample(shannon, lc, method = "bilinear")[[1]]
  classic_r <- resample(classic, lc, method = "bilinear")[[1]]
  twdtw_r   <- resample(twdtw,   lc, method = "bilinear")[[1]]
  
  aligned_stack <- c(shannon_r, classic_r, twdtw_r, lc)
  
  names(aligned_stack) <- c(
    "ShannonH",
    "RaosQ_Classic",
    "RaosQ_TWDTW",
    "LandCover"
  )
  
  return(aligned_stack)
})

names(Knepp_LandCover_AlignedIndices_YoI.NDVI) <- names(Knepp_Indices_YoI.NDVI_derived)

## Now I will crop and mask the land cover rasters to just the areas in my analysis
# I'm using the buffered site so I can compare inside the estate and out

Knepp_LandCover_AlignedIndices_YoI.NDVI <- lapply(
  Knepp_LandCover_AlignedIndices_YoI.NDVI,
  function(stack) {
    
    # Align polygon CRS to raster
    
    boundary_aligned <- project(KneppEstate_Buffered_Boundaries, stack)
    
    cropped <- crop(stack, boundary_aligned)
    masked  <- mask(cropped, boundary_aligned)
    
    return(masked)
  }
)

## Now I convert the raster stacks to a dataframe so I can run analyses on them
# The dataframe will have one row per pixel, and columns for each index and the land cover class

Knepp_DF_YoI.NDVI <- lapply(
  Knepp_LandCover_AlignedIndices_YoI.NDVI,
  function(stack) {
    
    df <- as.data.frame(stack, na.rm = TRUE)
    
    colnames(df) <- c(
      "ShannonH",
      "RaosQ_Classic",
      "RaosQ_TWDTW",
      "LandCover"
    )
    
    # Ensure land cover is treated as categorical
    
    df$LandCover <- as.factor(df$LandCover)
    
    return(df)
  }
)

# And now we loop over the dataframe to run the PERMANOVAe upon it

Knepp_PERMANOVA_YoI.NDVI <- lapply(
  Knepp_DF_YoI.NDVI,
  function(df) {
    
    ##Subsample the dataframe
    # Take a random sample of 10,000 rows to prevent distance-matrix memory crashes
    # The min() function acts as a safety net just in case a year has fewer than 10k valid pixels
    
    sample_size <- min(10000, nrow(df))
    df_subset <- df[sample(nrow(df), sample_size), ]
    
    list(
      
      # These are structured in the same way as the PERMANOVAe from the other case studies
      
      Shannon = adonis2(
        df_subset$ShannonH ~ df_subset$LandCover,
        #data = df_subset, # IDK why, but R cannot find ShannonH unless defined by a $
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_Classic = adonis2(
        RaosQ_Classic ~ LandCover,
        data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_TWDTW = adonis2(
        RaosQ_TWDTW ~ LandCover,
        data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      )
    )
  }
)

# Finally, putting the results into their own dataframe

Knepp_PERMANOVA_Summary.NDVI <- do.call(rbind, lapply(
  names(Knepp_PERMANOVA_YoI.NDVI),
  function(y) {
    
    res <- Knepp_PERMANOVA_YoI.NDVI[[y]]
    
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

message("Knepp Estate NDVI analysis complete.") # End of script ####
