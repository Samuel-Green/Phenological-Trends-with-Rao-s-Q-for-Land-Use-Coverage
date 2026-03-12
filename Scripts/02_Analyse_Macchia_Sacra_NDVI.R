############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025                    #
# 02_Analyse_Macchia_Sacra_NDVI.R                                              #
# Reproducing Macchia Sacra analysis (from Hackathon Preprint) using rasterdiv #
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

# Input NetCDF files provided by collaborators

MS_NDVI_File <- file.path(Macchia_Input, "ndviMSExp.nc")

# Output directory for this script

Macchia_Results <- file.path(Results, "Macchia Sacra")
dir.create(Macchia_Results, showWarnings = FALSE, recursive = TRUE)

## Read in the Macchia Sacra site data
# These NetCDFs are expected to be SpatRaster objects 
# Time encoded along the layer (Z) dimension, with 1 layer per acquisition date

message("Reading NDVI NetCDF files...")

MS_NDVI_Timeseries <- rast(MS_NDVI_File)

# Basic sanity checks

stopifnot(inherits(MS_NDVI_Timeseries, "SpatRaster"))

### Inspect temporal structure ####
## rasterdiv::paRao() requires an explicit time_vector
## This must correspond EXACTLY to the layer order

# Attempt to extract time directly from NetCDF metadata

MS_Time <- time(MS_NDVI_Timeseries)

if (is.null(MS_Time)) {
  stop("No temporal metadata found in NDVI NetCDF. Time vector must be supplied manually.")
} # Checks that time data is in the object "MS_Time"

# Convert to Date class if needed

MS_Time <- as.Date(MS_Time)

# Confirm dimensional agreement

stopifnot(length(MS_Time) == nlyr(MS_NDVI_Timeseries))

message(paste("Temporal length:", length(MS_Time), "layers"))

# Calculate the mean seasonal trajectory across the site

MS_Mean_NDVI_Timeseries <- global(MS_NDVI_Timeseries[[1:146]], fun = "mean", na.rm = TRUE) # I've subsetted MS_NDVI_Timeseries to just the first 146 layers as the subsequent 146 layers are SD instead.

png(file.path(Macchia_Results, "Macchia_Sacra_NDVI_Mean_Timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(MS_Time[1:146], MS_Mean_NDVI_Timeseries[,1],
     type = "l",
     lwd = 2,
     xlab = "Date",
     ylab = "Mean NDVI",
     main = "Macchia Sacra – Mean NDVI Time Series")

dev.off() # This actually exports the plot to the file

#### Diversity analyses: ####
### 1. Shannon-Wiener Index ####
## In the Hackathon, we first calculated, Shannon's H.
## It's applied to the mean yearly trajectory (i.e., collapse the time series first)

message("Calculating Shannon-Wiener diversity index...")

MS_NDVI_Mean_Raster <- app(MS_NDVI_Timeseries[[1:146]], fun = mean, na.rm = TRUE) # Subset to just the layers with observations (not standard deviations)

# Round to 2 decimals to avoid numerical saturation

MS_NDVI_Mean_Raster2dec <- round(MS_NDVI_Mean_Raster, 2)

# Calculate Shannon with moving window = 3

MS.NDVI.ShannonH.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(MS_NDVI_Mean_Raster2dec, wide = TRUE),
  window = 3,
  na.tolerance = 0
)

## Now turn the matrix of Shannon H values into a spatial raster
# Put it in a raster signature

MS_NDVI_ShannonH_Raster <- rast(MS.NDVI.ShannonH.matrix)

# Make the raster's extent and CRS match the original raster

ext(MS_NDVI_ShannonH_Raster) <- ext(MS_NDVI_Mean_Raster2dec)
crs(MS_NDVI_ShannonH_Raster) <- crs(MS_NDVI_Mean_Raster2dec)

names(MS_NDVI_ShannonH_Raster) <- "Shannon's H (Derived from NDVI)"

# Export the raster

writeRaster(
  MS_NDVI_ShannonH_Raster,
  filename = file.path(Macchia_Results, "Macchia_Sacra_NDVI_ShannonH_Raster.tif"),
  overwrite = TRUE
)

### 2. Classic Rao's Q  ####

message("Calculating classical Rao's Q...")

## This function calculates parametric Rao's Q
# This version calculates "classic" Rao's Q on only one matrix
# The subsequent multidimensional paRao will be used to calculate TWDTW Rao's Q

MS_NDVI_Classic_RaoQ <- paRao(
  x = MS_NDVI_Mean_Raster, # This also uses the mean of NDVI throughout time (i.e., all 146 observations are averaged into 1)
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  method = "classic"
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  MS_NDVI_Classic_RaoQ$window.3$alpha.2, # Subsets to just the spatial raster (as far as I can tell)
  filename = file.path(Macchia_Results, "Macchia_Sacra_NDVI_Classic-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### 3. Rao's Q with TWDTW ####
## Parameters are copied from the Hackathon preprint

message("Calculating Rao's Q with TWDTW distance...")

# If the function is taking too long to compute, make sure to enable parallelisation
# For running on an HPC with unknown cores, I can use the function `parallel::detectCores()`
# And set the "np" argument to detectCores() - 1
# Also requires the package 'snow' to parallelise, so I've disabled it for now

MS_NDVI_TWDTW_RaoQ <- paRao(
  x = MS_NDVI_Timeseries[[1:146]],
  time_vector = MS_Time[1:146],
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  np = 6, # Number of cores to use (I'd rather not wait forever)
  progBar = TRUE,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 35, # This is not the midpoint of the vector, rather the ecological midpoint of the cycle
  stepness = -0.5, # I just noticed, shouldn't the α value be called "steepness" instead of "stepness"?
  cycle_length = "year",
  time_scale = "day" # this specifies that our midpoint, 35, occurs after 35 days
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  MS_NDVI_TWDTW_RaoQ$window.3$alpha.2,
  filename = file.path(Macchia_Results, "Macchia_Sacra_NDVI_TWDTW-RaoQ_Raster.tif"),
  overwrite = TRUE
)

### Export rasters for comparison as a stacked GeoTIF ####
## Stack outputs for easy visual comparison (as in Figure 3 of the Hackathon preprint)
# Combine all my rasters into one list object

MS_NDVI_Index_Comparison_Rasters <- c(
  MS_NDVI_Mean_Raster, # The mean of per pixel NDVI for straightforward visual analysis
  MS_NDVI_ShannonH_Raster,
  MS_NDVI_Classic_RaoQ$window.3$alpha.2,
  MS_NDVI_TWDTW_RaoQ$window.3$alpha.2
) # This previously contained the NDVI raster, but I've decided to recompute in a separate file

# Name each layer

names(MS_NDVI_Index_Comparison_Rasters) <- c(
  "Sentinel-2 Mean NDVI",
  "Shannon's H (NDVI Derived)",
  "Classic Rao's Q (NDVI Derived)",
  "TWDTW Rao's Q (NDVI Derived)"
)

# Export it for later viewing and observation in QGIS

writeRaster(
  MS_NDVI_Index_Comparison_Rasters,
  filename = file.path(Macchia_Results, "Macchia_Sacra_NDVI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

png(file.path(Macchia_Results, "Macchia_Sacra_NDVI_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(MS_NDVI_Index_Comparison_Rasters) # A nice comparison plot

dev.off()

### Assess and compare each index's performance ####
## Firstly, import the shape file of land cover classification
# Path to the land cover classification shapefile from Matteo

MS.land.cover.file <- file.path(Macchia_Input, "MS_Data_From_Matteo/MacchiaSacra_Shapefile/MacchiaSacra_32632.shp")

# Read vector landcover

MacchiaSacra_LandCover_Vector <- vect(MS.land.cover.file)

# Ensure CRS matches diversity rasters

if (crs(MS_NDVI_Index_Comparison_Rasters) != crs(MacchiaSacra_LandCover_Vector)){
  message("The shapefile's CRS differs from the diversity index raster's CRS")
  MacchiaSacra_LandCover_Vector <- project(MacchiaSacra_LandCover_Vector, crs(MS_NDVI_Index_Comparison_Rasters)) # reprojects the shapefile if CRSs differ
} else {message("The shapefile's CRS matches the diversity index raster's CRS")}

## The rasters exceed the boundaries of the shapefile, so they will be cropped to size
# Crop to polygon extent

cropped.MS_NDVI_ShannonH_Raster     <- crop(MS_NDVI_ShannonH_Raster, MacchiaSacra_LandCover_Vector)
cropped.MS_NDVI_Classic_RaoQ <- crop(MS_NDVI_Classic_RaoQ$window.3$alpha.2, MacchiaSacra_LandCover_Vector)
cropped.MS_NDVI_TWDTW_RaoQ   <- crop(MS_NDVI_TWDTW_RaoQ$window.3$alpha.2, MacchiaSacra_LandCover_Vector)

# Mask outside polygon

masked.MS_NDVI_ShannonH_Raster     <- mask(cropped.MS_NDVI_ShannonH_Raster, MacchiaSacra_LandCover_Vector)
masked.MS_NDVI_Classic_RaoQ <- mask(cropped.MS_NDVI_Classic_RaoQ, MacchiaSacra_LandCover_Vector)
masked.MS_NDVI_TWDTW_RaoQ   <- mask(cropped.MS_NDVI_TWDTW_RaoQ, MacchiaSacra_LandCover_Vector)

## Make land cover classes rasterised instead of vectorised
# Firstly, I need to give each polygon/vector within the spatial vector a proper category label
# (currently, it's just a number corresponding to which description it is)

MacchiaSacra_LandCover_Vector$CAT <- MacchiaSacra_LandCover_Vector$decsr

# This defines a new raster object with the dominant land cover type as the dominant pixel value

MacchiaSacra_LandCover_Raster <- rasterize(
  MacchiaSacra_LandCover_Vector,
  masked.MS_NDVI_ShannonH_Raster,
  field = "CAT"
)

plot(MacchiaSacra_LandCover_Raster) # Looks cool

## Put the pixel values into a dataframe for easier computation
# Firstly, stack the spatial rasters into one object with vegetation class

masked.MS_NDVI_Index_Comparison_Rasters <- c(
  masked.MS_NDVI_ShannonH_Raster,
  masked.MS_NDVI_Classic_RaoQ,
  masked.MS_NDVI_TWDTW_RaoQ,
  MacchiaSacra_LandCover_Raster
)

# Rename the different layers of the spatial raster to something more memorable

names(masked.MS_NDVI_Index_Comparison_Rasters) <- c(
  "ShannonH (NDVI Derived)",
  "Rao's Q Classic (NDVI Derived)",
  "Rao's Q TWDTW (NDVI Derived)",
  "Vegetation Ground Truth"
)

# Convert the spatial raster to dataframe

masked.MS_NDVI_Index_Comparison_DF <- as.data.frame(masked.MS_NDVI_Index_Comparison_Rasters, na.rm = TRUE)

# Ensure vegetation is treated as factor [It is so I've disabled this line]
# 
# masked.MS_NDVI_Index_Comparison_DF$Vegetation <- as.factor(masked.MS_NDVI_Index_Comparison_DF$Vegetation)

# Rename the column names because otherwise R gets fussy

colnames(masked.MS_NDVI_Index_Comparison_DF) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

## Run the PERMANOVA
# How much variance in each index is explained by each vegetation class?

MS_NDVI_PERMANOVA_Shannon <- adonis2(
  masked.MS_NDVI_Index_Comparison_DF$ShannonH ~ 
    masked.MS_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

MS_NDVI_PERMANOVA_RaosQ_Classic <- adonis2(
  masked.MS_NDVI_Index_Comparison_DF$RaosQ_Classic ~ 
    masked.MS_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

MS_NDVI_PERMANOVA_RaosQ_TWDTW <- adonis2(
  masked.MS_NDVI_Index_Comparison_DF$RaosQ_TWDTW ~ 
    masked.MS_NDVI_Index_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

## Extract R² and p-values

MacchiaSacra_NDVI_PERMANOVA_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  R2 = c(
    MS_NDVI_PERMANOVA_Shannon$R2[1],
    MS_NDVI_PERMANOVA_RaosQ_Classic$R2[1],
    MS_NDVI_PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  `F` = c(
    MS_NDVI_PERMANOVA_Shannon$F[1],
    MS_NDVI_PERMANOVA_RaosQ_Classic$F[1],
    MS_NDVI_PERMANOVA_RaosQ_TWDTW$F[1]
  ),
  p_value = c(
    MS_NDVI_PERMANOVA_Shannon$`Pr(>F)`[1],
    MS_NDVI_PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
    MS_NDVI_PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
  )
)

print(MacchiaSacra_NDVI_PERMANOVA_Results)

#### Assess artefact (road) influence on diversity indices ####
### This analysis will test the ability of TWDTW Rao's Q to remove artefacts and compare it to other indices
## This analysis is broken up into multiple stages

message("Assessing road artefact influence...")

## 1. Import buffered road vector

MacchiaSacra_Road_File <- file.path(Macchia_Processed, "Italian_Roads/Italian roads bisecting Macchia Sacra.geojson")

MacchiaSacra_Road_Vector <- vect(MacchiaSacra_Road_File)

# Ensure CRS matches

if (crs(MacchiaSacra_Road_Vector) != crs(masked.MS_NDVI_ShannonH_Raster)){
  MacchiaSacra_Road_Vector <- project(MacchiaSacra_Road_Vector, crs(masked.MS_NDVI_ShannonH_Raster))
  message("The road's CRS differs from the diversity index raster's CRS")
} else {message("The road's CRS already matches the diversity index raster's CRS")}

## 2. Rasterise road mask
# Road pixels = 1, non-road = NA initially

MacchiaSacra_Road_Raster <- rasterize(
  MacchiaSacra_Road_Vector,
  MS_NDVI_Mean_Raster,
  field = 1
)

# Convert NA to 0 (non-road)

values(MacchiaSacra_Road_Raster)[is.na(values(MacchiaSacra_Road_Raster))] <- 0

names(MacchiaSacra_Road_Raster) <- "Road"

## 3. Stack indices with road mask

MS_NDVI_Road_Comparison_Raster <- c(
  MS_NDVI_Mean_Raster, # The mean NDVI raster for straightforward visual analysis
  MS_NDVI_ShannonH_Raster,
  MS_NDVI_Classic_RaoQ$window.3$alpha.2,
  MS_NDVI_TWDTW_RaoQ$window.3$alpha.2,
  MacchiaSacra_Road_Raster
)

names(MS_NDVI_Road_Comparison_Raster) <- c(
  "NDVI",
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Road"
)

MS_NDVI_Road_Comparison_DF <- as.data.frame(
  MS_NDVI_Road_Comparison_Raster,
  na.rm = TRUE
)

MS_NDVI_Road_Comparison_DF$Road <- as.factor(
  MS_NDVI_Road_Comparison_DF$Road
)

## 4. Compare mean contrast (effect size proxy)

Road_Contrast_Results <- data.frame(
  Index = c("NDVI", "Shannon", "Classic Rao", "TWDTW Rao"),
  Mean_Road = c(
    mean(MS_NDVI_Road_Comparison_DF$NDVI[MS_NDVI_Road_Comparison_DF$Road == 1]),
    mean(MS_NDVI_Road_Comparison_DF$ShannonH[MS_NDVI_Road_Comparison_DF$Road == 1]),
    mean(MS_NDVI_Road_Comparison_DF$RaosQ_Classic[MS_NDVI_Road_Comparison_DF$Road == 1]),
    mean(MS_NDVI_Road_Comparison_DF$RaosQ_TWDTW[MS_NDVI_Road_Comparison_DF$Road == 1])
  ),
  Mean_NonRoad = c(
    mean(MS_NDVI_Road_Comparison_DF$NDVI[MS_NDVI_Road_Comparison_DF$Road == 0]),
    mean(MS_NDVI_Road_Comparison_DF$ShannonH[MS_NDVI_Road_Comparison_DF$Road == 0]),
    mean(MS_NDVI_Road_Comparison_DF$RaosQ_Classic[MS_NDVI_Road_Comparison_DF$Road == 0]),
    mean(MS_NDVI_Road_Comparison_DF$RaosQ_TWDTW[MS_NDVI_Road_Comparison_DF$Road == 0])
  )
)

Road_Contrast_Results$Absolute_Difference <-
  abs(Road_Contrast_Results$Mean_Road -
        Road_Contrast_Results$Mean_NonRoad)

print(Road_Contrast_Results)

## 5. Wilcoxon tests (non-parametric contrast test)

MS_NDVI_Wilcox_Shannon <- wilcox.test(ShannonH ~ Road,
                              data = MS_NDVI_Road_Comparison_DF)

MS_NDVI_Wilcox_Classic <- wilcox.test(RaosQ_Classic ~ Road,
                              data = MS_NDVI_Road_Comparison_DF)

MS_NDVI_Wilcox_TWDTW <- wilcox.test(RaosQ_TWDTW ~ Road,
                            data = MS_NDVI_Road_Comparison_DF)

## 6. ROC–AUC comparison
# Higher AUC = stronger ability to detect road
# More robustness = LOWER AUC

MS_NDVI_ROC_Shannon <- roc(
  MS_NDVI_Road_Comparison_DF$Road,
  MS_NDVI_Road_Comparison_DF$ShannonH
)

MS_NDVI_ROC_Classic <- roc(
  MS_NDVI_Road_Comparison_DF$Road,
  MS_NDVI_Road_Comparison_DF$RaosQ_Classic
)

MS_NDVI_ROC_TWDTW <- roc(
  MS_NDVI_Road_Comparison_DF$Road,
  MS_NDVI_Road_Comparison_DF$RaosQ_TWDTW
)

MS_NDVI_Road_ROC_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  AUC = c(
    auc(MS_NDVI_ROC_Shannon),
    auc(MS_NDVI_ROC_Classic),
    auc(MS_NDVI_ROC_TWDTW)
  )
)

print(MS_NDVI_Road_ROC_Results)

### End of analyses ####

message("Macchia Sacra NDVI analysis complete.")
