############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025                    #
# 05_Analyse_Knepp_Estate_NDVI.R                                               #
# Comparative NDVI analysis of TWDTW Rao's Q & classic Rao's Q in Knepp Estate #
############################################################################ ###

### Install and load the necessary packages ####
## This should already be done from the setup file
# rasterdiv now contains the TWDTW-enabled paRao()

library(rasterdiv)
library(twdtw)
library(vegan)
library(pROC)

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

Knepp.spectral.bands <- c("Blue", "Green", "Red", "NIR", "NDVI", "PPI")

# Repeat each date for each band

knew.Knepp.knames <- as.vector(sapply(Knepp.dates, function(d) paste0(d, "_", Knepp.spectral.bands)))

# Assign names to each layer of the raster

names(Knepp_Buffered_Timeseries) <- knew.Knepp.knames

### Inspect temporal structure ####
## rasterdiv::paRao() requires an explicit time_vector
## This must correspond EXACTLY to the layer order

Knepp_Time <- as.Date(paste0(Knepp.dates, "-15"))
str(Knepp_Time) # WARNING: Vector is shorter than length as raster stack includes multiple bands for each date

## Subset and crop the raster
# Remove extraneous layers (as this is an NDVI analysis)

Knepp_Buffered_Timeseries.NDVI <- Knepp_Buffered_Timeseries[[grep("_NDVI$", names(Knepp_Buffered_Timeseries))]]
Knepp_Buffered_Timeseries.NDVI

# Import the shapefile with the Knepp Estate's borders

KneppEstate_Boundaries <- vect(file.path(Knepp_Input,"Schulte-to-Bühne_Knepp_Site_Data","knepp_boundary_lcc_pretty.shp"))

# Check that the CRS matches the timeseries raster

if (crs(Knepp_Buffered_Timeseries.NDVI) != crs(KneppEstate_Boundaries)){
  Knepp_Buffered_Timeseries.NDVI <- project(Knepp_Buffered_Timeseries.NDVI, crs(KneppEstate_Boundaries))
  message("The raster's CRS differs from the shape file's CRS")
} else {message("The raster's CRS already matches the shapefile's CRS")}

# Crop to Knepp Estate site boundaries

Knepp_Timeseries.NDVI <- crop(Knepp_Buffered_Timeseries.NDVI, KneppEstate_Boundaries, mask = TRUE)

## Calculate the mean seasonal trajectory across the site
# Calculate the site's mean value for each month

Knepp_Timeseries_Mean_NDVI <- global(Knepp_Timeseries.NDVI, fun = "mean", na.rm = TRUE) 

png(file.path(Knepp_Results, "Knepp_Estate_NDVI_Mean_Timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(Knepp_Time, Knepp_Timeseries_Mean_NDVI[,1],
     type = "l",
     lwd = 2,
     xlab = "Date",
     ylab = "Mean NDVI",
     main = "Knepp Estate – Mean NDVI Time Series")

dev.off() # This actually exports the plot to the file

#### Diversity analyses: ####
### 1. Shannon-Wiener Index ####
# Shannon's H is first applied to the collapsed timeseries of data

message("Calculating Shannon-Wiener diversity index...")

Knepp_NDVI_Mean_Raster <- app(Knepp_Timeseries.NDVI, fun = mean, na.rm = TRUE) # This takes all raster layers, and computes the mean pixel value based on all available layers

# Round to 2 decimals to avoid numerical saturation

Knepp_NDVI_Mean_Raster2dec <- round(Knepp_NDVI_Mean_Raster, 2)

# Calculate Shannon with moving window = 3

Knepp.NDVI.ShannonH.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(Knepp_NDVI_Mean_Raster2dec, wide = TRUE),
  window = 3,
  na.tolerance = 0
)

## Now turn the matrix of Shannon H values into a spatial raster
# Put it in a raster signature

Knepp_NDVI_ShannonH_Raster <- rast(Knepp.NDVI.ShannonH.matrix)

# Make the raster's extent and CRS match the original raster

ext(Knepp_NDVI_ShannonH_Raster) <- ext(Knepp_NDVI_Mean_Raster2dec)
crs(Knepp_NDVI_ShannonH_Raster) <- crs(Knepp_NDVI_Mean_Raster2dec)

names(Knepp_NDVI_ShannonH_Raster) <- "Shannon's H (Derived from NDVI)"

# Export the raster

writeRaster(
  Knepp_NDVI_ShannonH_Raster,
  filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_ShannonH_Raster.tif"),
  overwrite = TRUE
)

### 2. Classic Rao's Q  ####

message("Calculating classical Rao's Q...")

## This function calculates parametric Rao's Q
# This version calculates "classic" Rao's Q on only one matrix
# The subsequent multidimensional paRao will be used to calculate TWDTW Rao's Q

Knepp_NDVI_Classic_RaoQ <- paRao(
  x = Knepp_NDVI_Mean_Raster, # This also uses the mean of NDVI throughout time (i.e., all 146 observations are averaged into 1)
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  method = "classic"
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  Knepp_NDVI_Classic_RaoQ$window.3$alpha.2, # Subsets to just the spatial raster (as far as I can tell)
  filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_Classic-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### 3. Rao's Q with TWDTW ####
## Parameters are copied from the Hackathon preprint

message("Calculating Rao's Q with TWDTW distance...")

# If the function is taking too long to compute, make sure to enable parallelisation
# For running on an HPC with unknown cores, I can use the function `parallel::detectCores()`
# And set the "np" argument to detectCores() - 1
# Also requires the package 'snow' to parallelise, so I've disabled it for now

Knepp_NDVI_TWDTW_RaoQ <- paRao(
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
  Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2,
  filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_TWDTW-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### Export rasters for comparison as a stacked GeoTIF ####
## Stack outputs for easy visual comparison (as in Figure 3 of the Hackathon preprint)
# Combine all my rasters into one list object

Knepp_NDVI_Index_Comparison_Rasters <- c(
  Knepp_NDVI_Mean_Raster, # The mean of per pixel NDVI for straightforward visual analysis
  Knepp_NDVI_ShannonH_Raster,
  Knepp_NDVI_Classic_RaoQ$window.3$alpha.2,
  Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2
) # This previously contained the NDVI raster, but I've decided to recompute in a separate file

# Name each layer

names(Knepp_NDVI_Index_Comparison_Rasters) <- c(
  "Landsat Mean NDVI",
  "Shannon's H (NDVI Derived)",
  "Classic Rao's Q (NDVI Derived)",
  "TWDTW Rao's Q (NDVI Derived)"
)

# Export it for later viewing and observation in QGIS

writeRaster(
  Knepp_NDVI_Index_Comparison_Rasters,
  filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

png(file.path(Knepp_Results, "Knepp_Estate_NDVI_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(Knepp_NDVI_Index_Comparison_Rasters) # A nice comparison plot

dev.off()

### Assess and compare each index's performance ####
## Firstly, import the shape file of land cover classification
# Path to the land cover classification shapefile from Matteo

Knepp.land.cover.file <- file.path(Knepp_Input, "Knepp_Data_From_Matteo/KneppEstate_Shapefile/KneppEstate_32632.shp")

# Read vector landcover

KneppEstate_LandCover_Vector <- vect(Knepp.land.cover.file)

# Ensure CRS matches diversity rasters

if (crs(Knepp_NDVI_Index_Comparison_Rasters) != crs(KneppEstate_LandCover_Vector)){
  message("The shapefile's CRS differs from the diversity index raster's CRS")
  KneppEstate_LandCover_Vector <- project(KneppEstate_LandCover_Vector, crs(Knepp_NDVI_Index_Comparison_Rasters)) # reprojects the shapefile if CRSs differ
} else {message("The shapefile's CRS matches the diversity index raster's CRS")}

## The rasters exceed the boundaries of the shapefile, so they will be cropped to size
# Crop to polygon extent

cropped.Knepp_NDVI_ShannonH_Raster <- crop(Knepp_NDVI_ShannonH_Raster, KneppEstate_LandCover_Vector)
cropped.Knepp_NDVI_Classic_RaoQ <- crop(Knepp_NDVI_Classic_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)
cropped.Knepp_NDVI_TWDTW_RaoQ <- crop(Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)

# Mask outside polygon

masked.Knepp_NDVI_ShannonH_Raster <- mask(cropped.Knepp_NDVI_ShannonH_Raster, KneppEstate_LandCover_Vector)
masked.Knepp_NDVI_Classic_RaoQ <- mask(cropped.Knepp_NDVI_Classic_RaoQ, KneppEstate_LandCover_Vector)
masked.Knepp_NDVI_TWDTW_RaoQ <- mask(cropped.Knepp_NDVI_TWDTW_RaoQ, KneppEstate_LandCover_Vector)

## Make land cover classes rasterised instead of vectorised
# Firstly, I need to give each polygon/vector within the spatial vector a proper category label
# (currently, it's just a number corresponding to which description it is)

KneppEstate_LandCover_Vector$CAT <- KneppEstate_LandCover_Vector$decsr

# This defines a new raster object with the dominant land cover type as the dominant pixel value

KneppEstate_LandCover_Raster <- rasterize(
  KneppEstate_LandCover_Vector,
  masked.Knepp_NDVI_ShannonH_Raster,
  field = "CAT"
)

plot(KneppEstate_LandCover_Raster) # Looks cool

## Put the pixel values into a dataframe for easier computation
# Firstly, stack the spatial rasters into one object with vegetation class

masked.Knepp_NDVI_Index_Comparison_Rasters <- c(
  masked.Knepp_NDVI_ShannonH_Raster,
  masked.Knepp_NDVI_Classic_RaoQ,
  masked.Knepp_NDVI_TWDTW_RaoQ,
  KneppEstate_LandCover_Raster
)

# Rename the different layers of the spatial raster to something more memorable

names(masked.Knepp_NDVI_Index_Comparison_Rasters) <- c(
  "ShannonH (NDVI Derived)",
  "Rao's Q Classic (NDVI Derived)",
  "Rao's Q TWDTW (NDVI Derived)",
  "Vegetation Ground Truth"
)

# Convert the spatial raster to dataframe

masked.Knepp_NDVI_Index_Comparison_DF <- as.data.frame(masked.Knepp_NDVI_Index_Comparison_Rasters, na.rm = TRUE)

# Ensure vegetation is treated as factor [It is so I've disabled this line]
# 
# masked.Knepp_NDVI_Index_Comparison_DF$Vegetation <- as.factor(masked.Knepp_NDVI_Index_Comparison_DF$Vegetation)

# Rename the column names because otherwise R gets fussy

colnames(masked.Knepp_NDVI_Index_Comparison_DF) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

## Run the PERMANOVA
# How much variance in each index is explained by each vegetation class?

Knepp_NDVI_PERMANOVA_Shannon <- adonis2(
  masked.Knepp_NDVI_Index_Comparison_DF$ShannonH ~ 
    masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

Knepp_NDVI_PERMANOVA_RaosQ_Classic <- adonis2(
  masked.Knepp_NDVI_Index_Comparison_DF$RaosQ_Classic ~ 
    masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

Knepp_NDVI_PERMANOVA_RaosQ_TWDTW <- adonis2(
  masked.Knepp_NDVI_Index_Comparison_DF$RaosQ_TWDTW ~ 
    masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

## Extract R² and p-values

KneppEstate_NDVI_PERMANOVA_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  R2 = c(
    Knepp_NDVI_PERMANOVA_Shannon$R2[1],
    Knepp_NDVI_PERMANOVA_RaosQ_Classic$R2[1],
    Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  `F` = c(
    Knepp_NDVI_PERMANOVA_Shannon$F[1],
    Knepp_NDVI_PERMANOVA_RaosQ_Classic$F[1],
    Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$F[1]
  ),
  p_value = c(
    Knepp_NDVI_PERMANOVA_Shannon$`Pr(>F)`[1],
    Knepp_NDVI_PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
    Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
  )
)

print(KneppEstate_NDVI_PERMANOVA_Results)

message("Knepp Estate NDVI analysis complete.")
