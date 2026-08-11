############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025                    #
# 06_Analyse_Knepp_Estate_PPI.R                                                #
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

Knepp_Time <- as.Date(paste0(Knepp.dates, "-15"))
str(Knepp_Time) # WARNING: Vector is shorter than length as raster stack includes multiple bands for each date

### Compute the PPI values for each raster layer
## We will use the `ppi` package and apply the function to each pixel, then each layer of the raster stack
## 1. Subset the required bands for the PPI calculation
# PPI requires NIR, Red, and Solar Zenith Angle (SZA)

Knepp_Buffered_Timeseries.NIR <- Knepp_Buffered_Timeseries[[grep("_NIR$", names(Knepp_Buffered_Timeseries))]]
Knepp_Buffered_Timeseries.Red <- Knepp_Buffered_Timeseries[[grep("_Red$", names(Knepp_Buffered_Timeseries))]]
Knepp_Buffered_Timeseries.SZA <- Knepp_Buffered_Timeseries[[grep("_Solar Zenith Angle$", names(Knepp_Buffered_Timeseries))]]

## 2. Calculate DVI (Difference Vegetation Index)
# DVI = NIR - Red

Knepp_Buffered_Timeseries.DVI <- Knepp_Buffered_Timeseries.NIR - Knepp_Buffered_Timeseries.Red

## 3. Stack DVI and SZA together
# We stack them so we can feed a single object into our terra::app() wrapper function
# The first half of the layers will be DVI, the second half will be SZA

Knepp_DVI_SZA_Stack <- c(Knepp_Buffered_Timeseries.DVI, Knepp_Buffered_Timeseries.SZA)

## 4. Define a safe wrapper for the ppi() function
# The ppi function breaks if it encounters NAs. 
# This wrapper handles NAs gracefully by checking the pixel timeseries first.

calc_pixel_ppi <- function(x) {
  
  n_dates <- length(x) / 2
  dvi_ts <- x[1:n_dates]
  sza_deg_ts <- x[(n_dates + 1):(2 * n_dates)]
  
  # Identify valid (non-NA) time points
  valid_idx <- which(!is.na(dvi_ts) & !is.na(sza_deg_ts))
  
  # If too few valid points, return all NA
  if (length(valid_idx) < 2) {
    return(rep(NA, n_dates))
  }
  
  # Subset to valid observations only
  dvi_valid <- dvi_ts[valid_idx]
  sza_valid <- sza_deg_ts[valid_idx]
  
  # Convert SZA to radians
  sza_rad_valid <- sza_valid * (pi / 180)
  
  # Run PPI on cleaned time series
  ppi_valid <- ppi::ppi(
    dvi = dvi_valid,
    zenith.angle = sza_rad_valid
  )
  
  # Reconstruct full-length output (preserve time structure)
  ppi_full <- rep(NA, n_dates)
  ppi_full[valid_idx] <- ppi_valid
  
  return(ppi_full)
}

## 5. Apply the function pixel-by-pixel across the timeseries

message("Calculating PPI timeseries across all pixels...")
Knepp_Buffered_Timeseries.PPI <- app(Knepp_DVI_SZA_Stack, fun = calc_pixel_ppi)

# Assign names to the new PPI layers

names(Knepp_Buffered_Timeseries.PPI) <- paste0(Knepp.dates, "_PPI")

# Import the shapefile with the Knepp Estate's borders

KneppEstate_Boundaries <- vect(file.path(Knepp_Input,"Schulte-to-Bühne_Knepp_Site_Data","knepp_boundary_lcc_pretty.shp"))

# Check that the CRS matches the timeseries raster

if (crs(Knepp_Buffered_Timeseries.PPI) != crs(KneppEstate_Boundaries)){
  Knepp_Buffered_Timeseries.PPI <- project(Knepp_Buffered_Timeseries.PPI, crs(KneppEstate_Boundaries))
  message("The raster's CRS differs from the shape file's CRS")
} else {message("The raster's CRS already matches the shapefile's CRS")}

### Export the dataset as a NetCDF file for Gaussian Processing ####
## The Gaussian Process takes place separately in Python
## Export the data in the NetCDD format for use in the Python script
# Define the export path for the raw NetCDF

Knepp_NCDF_Raw_PPI_Path <- file.path(Knepp_Processed, "Knepp_Buffered_PPI_Raw.nc")

# Set the CRS to be exported as part of the NetCDF

NetCDF.CRS <- "EPSG:32630"

# Export the terra SpatRaster as a NetCDF file

message("Exporting raw PPI timeseries to NetCDF for Gaussian Process gap filling...")
terra::writeCDF(
  Knepp_Buffered_Timeseries.PPI, 
  filename = Knepp_NCDF_Raw_PPI_Path, 
  overwrite = TRUE, 
  varname = "PPI", 
  longname = "Plant Phenology Index",
  gridmap = NetCDF.CRS,
  unit = "unitless"
)

## Calculate the mean seasonal trajectory across the site
# Calculate the site's mean value for each month
# Crop to Knepp Estate site boundaries

Knepp_Timeseries.PPI <- crop(Knepp_Buffered_Timeseries.PPI, KneppEstate_Boundaries, mask = TRUE)

Knepp_Timeseries_Mean_PPI <- global(Knepp_Timeseries.PPI, fun = "mean", na.rm = TRUE) 

png(file.path(Knepp_Results, "Knepp_Estate_PPI_Mean_Timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(Knepp_Time, Knepp_Timeseries_Mean_PPI[,1],
     type = "l",
     lwd = 2,
     xlab = "Date",
     ylab = "Mean PPI",
     main = "Knepp Estate – Mean PPI Time Series")

dev.off() # This actually exports the plot to the file

#### Diversity analyses: ####
### 1. Shannon-Wiener Index ####
# Shannon's H is first applied to the collapsed timeseries of data

message("Calculating Shannon-Wiener diversity index...")

Knepp_PPI_Mean_Raster <- app(Knepp_Timeseries.PPI, fun = mean, na.rm = TRUE) # This takes all raster layers, and computes the mean pixel value based on all available layers

# Round to 2 decimals to avoid numerical saturation

Knepp_PPI_Mean_Raster2dec <- round(Knepp_PPI_Mean_Raster, 2)

# Calculate Shannon with moving window = 3

Knepp.PPI.ShannonH.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(Knepp_PPI_Mean_Raster2dec, wide = TRUE),
  window = 3,
  na.tolerance = 0
)

## Now turn the matrix of Shannon H values into a spatial raster
# Put it in a raster signature

Knepp_PPI_ShannonH_Raster <- rast(Knepp.PPI.ShannonH.matrix)

# Make the raster's extent and CRS match the original raster

ext(Knepp_PPI_ShannonH_Raster) <- ext(Knepp_PPI_Mean_Raster2dec)
crs(Knepp_PPI_ShannonH_Raster) <- crs(Knepp_PPI_Mean_Raster2dec)

names(Knepp_PPI_ShannonH_Raster) <- "Shannon's H (Derived from PPI)"

# Export the raster

writeRaster(
  Knepp_PPI_ShannonH_Raster,
  filename = file.path(Knepp_Results, "Knepp_Estate_PPI_ShannonH_Raster.tif"),
  overwrite = TRUE
)

### 2. Classic Rao's Q  ####

message("Calculating classical Rao's Q...")

## This function calculates parametric Rao's Q
# This version calculates "classic" Rao's Q on only one matrix
# The subsequent multidimensional paRao will be used to calculate TWDTW Rao's Q

Knepp_PPI_Classic_RaoQ <- paRao(
  x = Knepp_PPI_Mean_Raster, # This also uses the mean of PPI throughout time (i.e., all 146 observations are averaged into 1)
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  method = "classic"
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  Knepp_PPI_Classic_RaoQ$window.3$alpha.2, # Subsets to just the spatial raster (as far as I can tell)
  filename = file.path(Knepp_Results, "Knepp_Estate_PPI_Classic-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### 3. Rao's Q with TWDTW ####
## Parameters are copied from the Hackathon preprint

message("Calculating Rao's Q with TWDTW distance...")

# If the function is taking too long to compute, make sure to enable parallelisation
# For running on an HPC with unknown cores, I can use the function `parallel::detectCores()`
# And set the "np" argument to detectCores() - 1
# Also requires the package 'snow' to parallelise, so I've disabled it for now

Knepp_PPI_TWDTW_RaoQ <- paRao(
  x = Knepp_Timeseries,
  time_vector = Knepp_Time,
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  np = 6, # Number of cores to use (I'd rather not wait forever)
  progBar = TRUE,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 6, # This is not the midpoint of the vector, rather the ecological midpoint of the cycle
  stepness = -0.5, # I just noticed, shouldn't the α value be called "steepness" instead of "stepness"?
  cycle_length = "year",
  time_scale = "month" # this specifies that our midpoint, 35, occurs after 35 days
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  Knepp_PPI_TWDTW_RaoQ$window.3$alpha.2,
  filename = file.path(Knepp_Results, "Knepp_Estate_PPI_TWDTW-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### Export rasters for comparison as a stacked GeoTIF ####
## Stack outputs for easy visual comparison (as in Figure 3 of the Hackathon preprint)
# Combine all my rasters into one list object

Knepp_PPI_Index_Comparison_Rasters <- c(
  Knepp_PPI_Mean_Raster, # The mean of per pixel PPI for straightforward visual analysis
  Knepp_PPI_ShannonH_Raster,
  Knepp_PPI_Classic_RaoQ$window.3$alpha.2,
  Knepp_PPI_TWDTW_RaoQ$window.3$alpha.2
) # This previously contained the PPI raster, but I've decided to recompute in a separate file

# Name each layer

names(Knepp_PPI_Index_Comparison_Rasters) <- c(
  "Landsat Mean PPI",
  "Shannon's H (PPI Derived)",
  "Classic Rao's Q (PPI Derived)",
  "TWDTW Rao's Q (PPI Derived)"
)

# Export it for later viewing and observation in QGIS

writeRaster(
  Knepp_PPI_Index_Comparison_Rasters,
  filename = file.path(Knepp_Results, "Knepp_Estate_PPI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

png(file.path(Knepp_Results, "Knepp_Estate_PPI_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(Knepp_PPI_Index_Comparison_Rasters) # A nice comparison plot

dev.off()

### Assess and compare each index's performance ####
## Firstly, import the shape file of land cover classification
# Path to the land cover classification shapefile from Matteo

Knepp.land.cover.file <- file.path(Knepp_Input, "Knepp_Data_From_Matteo/KneppEstate_Shapefile/KneppEstate_32632.shp")

# Read vector landcover

KneppEstate_LandCover_Vector <- vect(Knepp.land.cover.file)

# Ensure CRS matches diversity rasters

if (crs(Knepp_PPI_Index_Comparison_Rasters) != crs(KneppEstate_LandCover_Vector)){
  message("The shapefile's CRS differs from the diversity index raster's CRS")
  KneppEstate_LandCover_Vector <- project(KneppEstate_LandCover_Vector, crs(Knepp_PPI_Index_Comparison_Rasters)) # reprojects the shapefile if CRSs differ
} else {message("The shapefile's CRS matches the diversity index raster's CRS")}

## The rasters exceed the boundaries of the shapefile, so they will be cropped to size
# Crop to polygon extent

cropped.Knepp_PPI_ShannonH_Raster <- crop(Knepp_PPI_ShannonH_Raster, KneppEstate_LandCover_Vector)
cropped.Knepp_PPI_Classic_RaoQ <- crop(Knepp_PPI_Classic_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)
cropped.Knepp_PPI_TWDTW_RaoQ <- crop(Knepp_PPI_TWDTW_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)

# Mask outside polygon

masked.Knepp_PPI_ShannonH_Raster <- mask(cropped.Knepp_PPI_ShannonH_Raster, KneppEstate_LandCover_Vector)
masked.Knepp_PPI_Classic_RaoQ <- mask(cropped.Knepp_PPI_Classic_RaoQ, KneppEstate_LandCover_Vector)
masked.Knepp_PPI_TWDTW_RaoQ <- mask(cropped.Knepp_PPI_TWDTW_RaoQ, KneppEstate_LandCover_Vector)

## Make land cover classes rasterised instead of vectorised
# Firstly, I need to give each polygon/vector within the spatial vector a proper category label
# (currently, it's just a number corresponding to which description it is)

KneppEstate_LandCover_Vector$CAT <- KneppEstate_LandCover_Vector$decsr

# This defines a new raster object with the dominant land cover type as the dominant pixel value

KneppEstate_LandCover_Raster <- rasterize(
  KneppEstate_LandCover_Vector,
  masked.Knepp_PPI_ShannonH_Raster,
  field = "CAT"
)

plot(KneppEstate_LandCover_Raster) # Looks cool

## Put the pixel values into a dataframe for easier computation
# Firstly, stack the spatial rasters into one object with vegetation class

masked.Knepp_PPI_Index_Comparison_Rasters <- c(
  masked.Knepp_PPI_ShannonH_Raster,
  masked.Knepp_PPI_Classic_RaoQ,
  masked.Knepp_PPI_TWDTW_RaoQ,
  KneppEstate_LandCover_Raster
)

# Rename the different layers of the spatial raster to something more memorable

names(masked.Knepp_PPI_Index_Comparison_Rasters) <- c(
  "ShannonH (PPI Derived)",
  "Rao's Q Classic (PPI Derived)",
  "Rao's Q TWDTW (PPI Derived)",
  "Vegetation Ground Truth"
)

# Convert the spatial raster to dataframe

masked.Knepp_PPI_Index_Comparison_DF <- as.data.frame(masked.Knepp_PPI_Index_Comparison_Rasters, na.rm = TRUE)

# Ensure vegetation is treated as factor [It is so I've disabled this line]
# 
# masked.Knepp_PPI_Index_Comparison_DF$Vegetation <- as.factor(masked.Knepp_PPI_Index_Comparison_DF$Vegetation)

# Rename the column names because otherwise R gets fussy

colnames(masked.Knepp_PPI_Index_Comparison_DF) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

## Run the PERMANOVA
# How much variance in each index is explained by each vegetation class?

Knepp_PPI_PERMANOVA_Shannon <- adonis2(
  masked.Knepp_PPI_Index_Comparison_DF$ShannonH ~ 
    masked.Knepp_PPI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

Knepp_PPI_PERMANOVA_RaosQ_Classic <- adonis2(
  masked.Knepp_PPI_Index_Comparison_DF$RaosQ_Classic ~ 
    masked.Knepp_PPI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

Knepp_PPI_PERMANOVA_RaosQ_TWDTW <- adonis2(
  masked.Knepp_PPI_Index_Comparison_DF$RaosQ_TWDTW ~ 
    masked.Knepp_PPI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

## Extract R² and p-values

KneppEstate_PPI_PERMANOVA_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  R2 = c(
    Knepp_PPI_PERMANOVA_Shannon$R2[1],
    Knepp_PPI_PERMANOVA_RaosQ_Classic$R2[1],
    Knepp_PPI_PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  `F` = c(
    Knepp_PPI_PERMANOVA_Shannon$F[1],
    Knepp_PPI_PERMANOVA_RaosQ_Classic$F[1],
    Knepp_PPI_PERMANOVA_RaosQ_TWDTW$F[1]
  ),
  p_value = c(
    Knepp_PPI_PERMANOVA_Shannon$`Pr(>F)`[1],
    Knepp_PPI_PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
    Knepp_PPI_PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
  )
)

print(KneppEstate_PPI_PERMANOVA_Results)

message("Knepp Estate PPI analysis complete.")
