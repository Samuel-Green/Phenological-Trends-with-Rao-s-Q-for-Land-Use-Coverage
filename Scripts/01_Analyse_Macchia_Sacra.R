############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025                    #
# 01_Analyse_Macchia_Sacra.R                                                   #
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

PPI_file  <- file.path(Macchia_Raw, "ppiMSExp.nc")
NDVI_file <- file.path(Macchia_Raw, "ndviMSExp.nc")

# Output directory for this script

Macchia_Results <- file.path(Results, "Macchia Sacra")
dir.create(Macchia_Results, showWarnings = FALSE, recursive = TRUE)

## Read in the Macchia Sacra site data
# These NetCDFs are expected to be SpatRaster objects 
# Time encoded along the layer (Z) dimension, with 1 layer per acquisition date

message("Reading PPI and NDVI NetCDF files...")

PPI_ts  <- rast(PPI_file)
NDVI_ts <- rast(NDVI_file)

# Basic sanity checks

stopifnot(inherits(PPI_ts, "SpatRaster"))
stopifnot(inherits(NDVI_ts, "SpatRaster"))

### Inspect temporal structure ####
## rasterdiv::paRao() requires an explicit time_vector
## This must correspond EXACTLY to the layer order

# Attempt to extract time directly from NetCDF metadata

PPI_time <- time(PPI_ts)

if (is.null(PPI_time)) {
  stop("No temporal metadata found in PPI NetCDF. Time vector must be supplied manually.")
} # Checks that time data is in the object "PPI_time"

# Convert to Date class if needed

PPI_time <- as.Date(PPI_time)

# Confirm dimensional agreement

stopifnot(length(PPI_time) == nlyr(PPI_ts))

message(paste("Temporal length:", length(PPI_time), "layers"))

# Calculate the mean seasonal trajectory across the site

PPI_mean_ts <- global(PPI_ts[[1:146]], fun = "mean", na.rm = TRUE) # I've subsetted PPI_ts to just the first 146 layers as the subsequent 146 layers are SD instead.

png(file.path(Macchia_Results, "PPI_mean_timeseries.png"), # Specifies that I want a 4K resolution .png file
    width = 3840, height = 2160, res = 150)

plot(PPI_time[1:146], PPI_mean_ts[,1],
     type = "l",
     lwd = 2,
     xlab = "Date",
     ylab = "Mean PPI",
     main = "Macchia Sacra – Mean PPI Time Series")

dev.off() # This actually exports the plot to the file

#### Diversity analyses: ####
### 1. Shannon-Wiener Index ####
## In the Hackathon, we first calculated, Shannon's H.
## It's applied to the mean yearly trajectory (i.e., collapse the time series first)

message("Calculating Shannon diversity...")

PPI_mean_raster <- app(PPI_ts[[1:146]], fun = mean, na.rm = TRUE) # Subset to just the layers with observations (not standard deviations)

PPI_mean_matrix <- matrix(NA, nrow = nrow(PPI_mean_raster), ncol = ncol(PPI_mean_raster)) # Create an empty matrix to store cell values from the spatial raster

PPI_mean_matrix[1:dim(PPI_mean_matrix)[1], 1:dim(PPI_mean_matrix)[2]] <- values(PPI_mean_raster) # Fill the matrix with values from the spatial raster, maintaining the same dimensions as the original raster

ShannonH_matrix <- rasterdiv::ShannonS( # Calculate Shannon's H (Shannon-Wiener Index) for the matrix
  x = PPI_mean_matrix,
  window = 1,
  na.tolerance = 0
)

ShannonH_raster <- rasterize(ShannonH_matrix, PPI_mean_raster)

values(ShannonH_raster) <- ShannonH_matrix

writeRaster(
  ShannonH_raster,
  filename = file.path(Macchia_Results, "Macchia_Sacra_ShannonH_raster.tif"),
  overwrite = TRUE
) # Export the data as a GeoTIF

### 2. Classic Rao's Q  ####

message("Calculating classical Rao's Q...")

## This function calculates parametric Rao's Q
# This version calculates "classic" Rao's Q on only one matrix
# The subsequent multidimensional paRao will be used to calculate TWDTW Rao's Q

Rao_classic <- paRao(
  x = PPI_mean_raster, # This also uses the mean of PPI throughout time (i.e., all 146 observations are averaged into 1)
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 3,
  method = "classic"
)

# This function exports my object as a GeoTIF for viewing in QGIS etc.

writeRaster(
  Rao_classic$window.3$alpha.2, # Subsets to just the spatial raster (as far as I can tell)
  filename = file.path(Macchia_Results, "Macchia_Sacra_RaoQ_Classic_raster.tif"),
  overwrite = TRUE
)

### 3. Rao's Q with TWDTW ####
## Parameters are copied from the Hackathon preprint

message("Calculating Rao's Q with TWDTW distance...")

# If the function is taking too long to compute, make sure to enable parallelisation
# For running on an HPC with unknown cores, I can use the function `parallel::detectCores()`
# And set the "np" argument to detectCores() - 1
# Also requires the package 'snow' to parallelise, so I've disabled it for now

Rao_TWDTW <- paRao(
  x = PPI_ts[[1:146]],
  time_vector = PPI_time[1:146],
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 3,
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
  Rao_TWDTW$window.3$alpha.2,
  filename = file.path(Macchia_Results, "RaoQ_TWDTW_PPI.tif"),
  overwrite = TRUE
)

### Export rasters for comparison as a stacked GeoTIF ####
## Stack outputs for easy visual comparison (as in Figure 3 of the Hackathon preprint)

# Combine all my rasters into one list object

MacchiaSacra_Comparison_Rasters <- c(
  PPI_ts[[146]], # The last PPI image for straightforward visual analysis
  NDVI_ts[[146]], # The last NDVI image because everyone wants to see NDVI
  ShannonH_raster,
  Rao_classic$window.3$alpha.2,
  Rao_TWDTW$window.3$alpha.2
)

names(MacchiaSacra_Comparison_Rasters) <- c(
  "Sentinel-2 PPI",
  "Sentinel-2 NDVI",
  "Shannon's H",
  "Classic Rao's Q",
  "TWDTW Rao's Q"
)

writeRaster(
  MacchiaSacra_Comparison_Rasters,
  filename = file.path(Macchia_Results, "Macchia_Sacra_Diversity_Comparison.tif"),
  overwrite = TRUE
)

png(file.path(Macchia_Results, "Macchia_Sacra_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(MacchiaSacra_Comparison_Rasters) # A nice comparison plot

dev.off()

### Assess and compare each index's performance ####
## Firstly, import the shape file of land cover classification

# Path to the land cover classification shapefile from Matteo

MacchiaSacra_LandCover_File <- file.path(Macchia_Raw, "MS_Data_From_Matteo/MacchiaSacra_Shapefile/MacchiaSacra_32632.shp")

# Read vector landcover

MacchiaSacra_LandCover_Vector <- vect(MacchiaSacra_LandCover_File)

# Ensure CRS matches diversity rasters

if (crs(ShannonH_raster) != crs(MacchiaSacra_LandCover_Vector)){
  message("The shapefile's CRS differs from the diversity index raster's CRS")
  MacchiaSacra_LandCover_Vector <- project(MacchiaSacra_LandCover_Vector, crs(ShannonH_raster)) # reprojects the shapefile if CRSs differ
} else {message("The shapefile's CRS matches the diversity index raster's CRS")}

## The rasters exceed the boundaries of the shapefile, so they will be cropped to size
# Crop to polygon extent

ShannonH_crop     <- crop(ShannonH_raster, MacchiaSacra_LandCover_Vector)
Rao_classic_crop <- crop(Rao_classic$window.3$alpha.2, MacchiaSacra_LandCover_Vector)
Rao_TWDTW_crop   <- crop(Rao_TWDTW$window.3$alpha.2, MacchiaSacra_LandCover_Vector)

# Mask outside polygon

ShannonH_masked     <- mask(ShannonH_crop, MacchiaSacra_LandCover_Vector)
Rao_classic_masked <- mask(Rao_classic_crop, MacchiaSacra_LandCover_Vector)
Rao_TWDTW_masked   <- mask(Rao_TWDTW_crop, MacchiaSacra_LandCover_Vector)

## Make land cover classes rasterised instead of vectorised
# Firstly, I need to give each polygon/vector within the spatial vector a proper category label
# (currently, it's just a number corresponding to which description it is)

MacchiaSacra_LandCover_Vector$CAT <- MacchiaSacra_LandCover_Vector$decsr

# This defines a new raster object with the dominant land cover type as the dominant pixel value

MacchiaSacra_LandCover_Raster <- rasterize(
  MacchiaSacra_LandCover_Vector,
  ShannonH_masked,
  field = "CAT"
)

plot(MacchiaSacra_LandCover_Raster) # Looks cool

## Put the pixel values into a dataframe for easier computation
# Firstly, stack the spatial rasters into one object with vegetation class

MacchiaSacra_Indices_Comparison_Raster <- c(
  ShannonH_masked,
  Rao_classic_masked,
  Rao_TWDTW_masked,
  MacchiaSacra_LandCover_Raster
)

# Rename the different layers of the spatial raster to something more memorable

names(MacchiaSacra_Indices_Comparison_Raster) <- c(
  "ShannonH",
  "Rao's Q Classic",
  "Rao's Q TWDTW",
  "Vegetation Ground Truth"
)

# Convert the spatial raster to dataframe

MacchiaSacra_Indices_Comparison_DF <- as.data.frame(MacchiaSacra_Indices_Comparison_Raster, na.rm = TRUE)

# Ensure vegetation is treated as factor [It is so I've disabled this line]
# 
# MacchiaSacra_Indices_Comparison_DF$Vegetation <- as.factor(MacchiaSacra_Indices_Comparison_DF$Vegetation)

# Rename the column names because otherwise R gets fussy

colnames(MacchiaSacra_Indices_Comparison_DF) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

## Run the PERMANOVA
# How much variance in each index is explained by each vegetation class?

PERMANOVA_Shannon <- adonis2(
  MacchiaSacra_Indices_Comparison_DF$ShannonH ~ 
    MacchiaSacra_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

PERMANOVA_RaosQ_Classic <- adonis2(
  MacchiaSacra_Indices_Comparison_DF$RaosQ_Classic ~ 
    MacchiaSacra_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

PERMANOVA_RaosQ_TWDTW <- adonis2(
  MacchiaSacra_Indices_Comparison_DF$RaosQ_TWDTW ~ 
    MacchiaSacra_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 9999
)

## Extract R² and p-values

MacchiaSacra_PERMANOVA_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  R2 = c(
    PERMANOVA_Shannon$R2[1],
    PERMANOVA_RaosQ_Classic$R2[1],
    PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  `F` = c(
    PERMANOVA_Shannon$F[1],
    PERMANOVA_RaosQ_Classic$F[1],
    PERMANOVA_RaosQ_TWDTW$F[1]
  ),
  p_value = c(
    PERMANOVA_Shannon$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
  )
)

print(MacchiaSacra_PERMANOVA_Results)

#### Assess artefact (road) influence on diversity indices ####
### This analysis will test the ability of TWDTW Rao's Q to remove artefacts and compare it to other indices
## This analysis is broken up into multiple stages

message("Assessing road artefact influence...")

## 1. Import buffered road vector

MacchiaSacra_Road_File <- file.path(Macchia_Processed, "Italian_Roads/Italian roads bisecting Macchia Sacra.geojson")

MacchiaSacra_Road_Vector <- vect(MacchiaSacra_Road_File)

# Ensure CRS matches

if (crs(MacchiaSacra_Road_Vector) != crs(ShannonH_masked)){
  MacchiaSacra_Road_Vector <- project(MacchiaSacra_Road_Vector, crs(ShannonH_masked))
  message("The road's CRS differs from the diversity index raster's CRS")
} else {message("The road's CRS already matches the diversity index raster's CRS")}

## 2. Rasterise road mask
# Road pixels = 1, non-road = NA initially

MacchiaSacra_Road_Raster <- rasterize(
  MacchiaSacra_Road_Vector,
  PPI_ts[[146]],
  field = 1
)

# Convert NA to 0 (non-road)

values(MacchiaSacra_Road_Raster)[is.na(values(MacchiaSacra_Road_Raster))] <- 0

names(MacchiaSacra_Road_Raster) <- "Road"

## 3. Stack indices with road mask

MacchiaSacra_Road_Comparison_Raster <- c(
  PPI_ts[[146]], # The last PPI image for straightforward visual analysis
  NDVI_ts[[146]], # The last NDVI image because everyone wants to see NDVI
  ShannonH_raster,
  Rao_classic$window.3$alpha.2,
  Rao_TWDTW$window.3$alpha.2,
  MacchiaSacra_Road_Raster
)

names(MacchiaSacra_Road_Comparison_Raster) <- c(
  "PPI",
  "NDVI",
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Road"
)

MacchiaSacra_Road_Comparison_DF <- as.data.frame(
  MacchiaSacra_Road_Comparison_Raster,
  na.rm = TRUE
)

MacchiaSacra_Road_Comparison_DF$Road <- as.factor(
  MacchiaSacra_Road_Comparison_DF$Road
)

## 4. Compare mean contrast (effect size proxy)

Road_Contrast_Results <- data.frame(
  Index = c("PPI", "NDVI", "Shannon", "Classic Rao", "TWDTW Rao"),
  Mean_Road = c(
    mean(MacchiaSacra_Road_Comparison_DF$PPI[MacchiaSacra_Road_Comparison_DF$Road == 1]),
    mean(MacchiaSacra_Road_Comparison_DF$NDVI[MacchiaSacra_Road_Comparison_DF$Road == 1]),
    mean(MacchiaSacra_Road_Comparison_DF$ShannonH[MacchiaSacra_Road_Comparison_DF$Road == 1]),
    mean(MacchiaSacra_Road_Comparison_DF$RaosQ_Classic[MacchiaSacra_Road_Comparison_DF$Road == 1]),
    mean(MacchiaSacra_Road_Comparison_DF$RaosQ_TWDTW[MacchiaSacra_Road_Comparison_DF$Road == 1])
  ),
  Mean_NonRoad = c(
    mean(MacchiaSacra_Road_Comparison_DF$PPI[MacchiaSacra_Road_Comparison_DF$Road == 0]),
    mean(MacchiaSacra_Road_Comparison_DF$NDVI[MacchiaSacra_Road_Comparison_DF$Road == 0]),
    mean(MacchiaSacra_Road_Comparison_DF$ShannonH[MacchiaSacra_Road_Comparison_DF$Road == 0]),
    mean(MacchiaSacra_Road_Comparison_DF$RaosQ_Classic[MacchiaSacra_Road_Comparison_DF$Road == 0]),
    mean(MacchiaSacra_Road_Comparison_DF$RaosQ_TWDTW[MacchiaSacra_Road_Comparison_DF$Road == 0])
  )
)

Road_Contrast_Results$Absolute_Difference <-
  abs(Road_Contrast_Results$Mean_Road -
        Road_Contrast_Results$Mean_NonRoad)

print(Road_Contrast_Results)

## 5. Wilcoxon tests (non-parametric contrast test)

Wilcox_Shannon <- wilcox.test(ShannonH ~ Road,
                              data = MacchiaSacra_Road_Comparison_DF)

Wilcox_Classic <- wilcox.test(RaosQ_Classic ~ Road,
                              data = MacchiaSacra_Road_Comparison_DF)

Wilcox_TWDTW <- wilcox.test(RaosQ_TWDTW ~ Road,
                            data = MacchiaSacra_Road_Comparison_DF)

## 6. ROC–AUC comparison
# Higher AUC = stronger ability to detect road
# More robustness = LOWER AUC

ROC_Shannon <- roc(
  MacchiaSacra_Road_Comparison_DF$Road,
  MacchiaSacra_Road_Comparison_DF$ShannonH
)

ROC_Classic <- roc(
  MacchiaSacra_Road_Comparison_DF$Road,
  MacchiaSacra_Road_Comparison_DF$RaosQ_Classic
)

ROC_TWDTW <- roc(
  MacchiaSacra_Road_Comparison_DF$Road,
  MacchiaSacra_Road_Comparison_DF$RaosQ_TWDTW
)

Road_ROC_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  AUC = c(
    auc(ROC_Shannon),
    auc(ROC_Classic),
    auc(ROC_TWDTW)
  )
)

print(Road_ROC_Results)

### End of analyses ####

message("Macchia Sacra analysis complete.")