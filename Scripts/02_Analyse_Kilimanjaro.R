############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/02/2026                    #
# 02_Analyse_Kilimanjaro.R                                                     #
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
## NOTICE: Elliot's computer is low on storage, so some large GeoTIFFs may have also been loaded from external hard drives not listed below

# Output directory for this script

KiliNP_Results <- file.path(Results, "Kilimanjaro")
dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

# Load KiliNP_LandCover_Vector boundary

KiliNP_LandCover_Vector <- vect(paste0(KiliNP_Input,"/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp"))

# List raster files
# Only run these lines if I haven't already generated my cropped raster files (it takes forever)

tmp.files <- list.files(KiliNP_Input, pattern = "\\FBM.tif$", full.names = TRUE)

for (f in tmp.files) {
  
 tmp.KiliNP.raster<- rast(f)
  
  # Ensure CRS matches
  if (!crs(tmp.KiliNP.raster) == crs(KiliNP_LandCover_Vector)) {
    KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(tmp.KiliNP.raster))
  }
  
  # Crop to bounding box
  tmp.KiliNP.raster.cropped <- crop(tmp.KiliNP.raster, KiliNP_LandCover_Vector)
  
  # Mask to exact boundary
  tmp.KiliNP.raster.masked <- mask(tmp.KiliNP.raster.cropped, KiliNP_LandCover_Vector)
  
  # Extract the year from the filename (because the original filenames are a *mess*)
  tmp.year <- sub(".*_(\\d{4})_.*", "\\1", basename(f))
  
  # Create new filename
  tmp.new.name <- paste0("KiliNP_", tmp.year, "_Cropped.tif")
  
  # Write to processed folder
  writeRaster(
    tmp.KiliNP.raster.masked,
    filename = file.path(KiliNP_Processed, tmp.new.name),
    overwrite = TRUE
  )
}

# To keep the memory usage under control, remove the tmp. raster files

rm(tmp.KiliNP.raster, tmp.KiliNP.raster.cropped, tmp.KiliNP.raster.masked)

## Import the cropped rasters and combine them into one object for analysis
# Get the files's names and locations

KiliNP_Cropped_Files <- list.files(
  KiliNP_Processed,
  pattern = "^KiliNP_\\d{4}_Cropped\\.tif$",
  full.names = TRUE
)

# Import them as one raster

KiliNP_Timeseries <- rast(KiliNP_Cropped_Files)

## Rename the layers to something a bit more readable
# Extract years from filenames

Kili.years <- sub(".*_(\\d{4})_Cropped\\.tif", "\\1", basename(KiliNP_Cropped_Files))

# Create new layer names

Kili.layer.names <- unlist(lapply(Kili.years, function(y) {
  paste0(y, " - ", month.name)
}))

# Assign names

names(KiliNP_Timeseries) <- Kili.layer.names

### Mask pixels in the raster stack which don't have a complete timeseries of data

# Check which layers are completely NA with a quick visual inspection

blank_layers <- for(i in 1:nlyr(KiliNP_Timeseries)) {
  plot(KiliNP_Timeseries[[i]], main = names(KiliNP_Timeseries)[i])
  if(i < nlyr(KiliNP_Timeseries)) readline(prompt = "Press [enter] to continue")
}

blank_layers # This is just to inspect each layer for the user's interest

## Only pixels with a complete set of data for every layer are suitable for analysis
# There are many pixels with NA values scattered throughout the raster stack
# Create logical mask: TRUE only where ALL layers are non-NA

Kili.pixel.mask <- app(KiliNP_Timeseries, function(x) all(!is.na(x))) # Will consume lots of RAM

# Mask out incomplete pixels (FALSE becomes NA)

KiliNP_Timeseries_Clean <- mask(KiliNP_Timeseries, Kili.pixel.mask, maskvalues = 0) # This is very computationally challenging, run on HPC

# Export raster so I don't have to calculate it every time

writeRaster(KiliNP_Timeseries_Clean, file.path(KiliNP_Processed, "KiliNP_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE)

# And then load back in the raster

KiliNP_Timeseries_Clean <- rast(file.path(KiliNP_Processed, "KiliNP_2017-2021_Cropped_&_Masked.tif"))

### Inspect temporal structure ####

## Unlike Macchia Sacra NetCDF, this GeoTIFF stack does not contain explicit time metadata
## Therefore, we must construct the time vector manually from layer names

message("Constructing time vector for Kilimanjaro time series...")

# Extract year and month from layer names
# Layer format: "2017 - January"

#Kili.layer.names <- names(KiliNP_Timeseries_Clean)

Kili.dates <- as.Date(
  paste0(
    sub(" - .*", "", names(KiliNP_Timeseries_Clean)), "-",
    match(sub(".* - ", "", names(KiliNP_Timeseries_Clean)), month.name),
    "-01"
  ),
  format = "%Y-%m-%d"
)

stopifnot(length(Kili.dates) == nlyr(KiliNP_Timeseries_Clean))

message(paste("Temporal length:", length(Kili.dates), "layers"))

### 1. Shannon-Wiener Index ####

## As in Macchia Sacra, collapse the time series to mean annual trajectory first

message("Calculating Shannon-Wiener diversity index for Kilimanjaro...")

KiliNP_Mean_Raster <- app(KiliNP_Timeseries_Clean, fun = mean, na.rm = TRUE)

# Export raster so I don't have to calculate it every time

writeRaster(KiliNP_Mean_Raster, file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Cropped_&_Masked.tif"), overwrite = TRUE)

# And then load back in the raster

KiliNP_Mean_Raster <- rast(file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Cropped_&_Masked.tif"))

# Round to 2 decimals to avoid numerical saturation

KiliNP_Mean_Raster2dec <- round(KiliNP_Mean_Raster, 2)

# Run ShannonS

KiliNP.ShannonH.matrix <- rasterdiv::ShannonS(
  x = terra::as.matrix(KiliNP_Mean_Raster2dec, wide = TRUE), # The function here converts my raster to a matrix suitable for analysis
  window = 3,
  na.tolerance = 0
)

## Now turn the matrix of Shannon H values into a spatial raster
# Put it in a raster signature

KiliNP_ShannonH_Raster <- rast(KiliNP.ShannonH.matrix)

# Make the raster's extent and CRS match the original raster

ext(KiliNP_ShannonH_Raster) <- ext(KiliNP_Mean_Raster2dec)
crs(KiliNP_ShannonH_Raster) <- crs(KiliNP_Mean_Raster2dec)

names(KiliNP_ShannonH_Raster) <- "Shannon's H"

# Export the raster

writeRaster(
  KiliNP_ShannonH_Raster,
  filename = file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_raster.tif"),
  overwrite = TRUE
)

KiliNP_ShannonH_Raster <- rast(file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_raster.tif"))

### 2. Classic Rao's Q  ####

message("Calculating classical Rao's Q for Kilimanjaro...")

KiliNP_Classic_RaoQ <- paRao(
  x = KiliNP_Mean_Raster,
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  method = "classic",
  np = detectCores() -1
)

writeRaster(
  KiliNP_Classic_RaoQ$window.3$alpha.2,
  filename = file.path(KiliNP_Results, "KiliNP_RaoQ_Classic_raster.tif"),
  overwrite = TRUE
)

### 3. Rao's Q with TWDTW ####

message("Calculating Rao's Q with TWDTW distance for Kilimanjaro...")

Kili_Rao_TWDTW <- paRao(
  x = KiliNP_Timeseries_Clean,
  time_vector = Kili.dates,
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  np = detectCores() -1,
  progBar = TRUE,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 6,          # Midpoint of annual cycle (June)
  stepness = -0.5,
  cycle_length = "year",
  time_scale = "month"   # Now explicitly monthly data
)

writeRaster(
  Kili_Rao_TWDTW$window.3$alpha.2,
  filename = file.path(KiliNP_Results, "KiliNP_RaoQ_TWDTW.tif"),
  overwrite = TRUE
)

### Export rasters for comparison ####

Kili_Comparison_Rasters <- c(
  KiliNP_Timeseries_Clean[[nlyr(KiliNP_Timeseries_Clean)]],
  KiliNP_ShannonH_Raster,
  KiliNP_Classic_RaoQ$window.3$alpha.2,
  Kili_Rao_TWDTW$window.3$alpha.2
)

names(Kili_Comparison_Rasters) <- c(
  "Sentinel-2 NDVI",
  "Shannon's H",
  "Classic Rao's Q",
  "TWDTW Rao's Q"
)

writeRaster(
  Kili_Comparison_Rasters,
  filename = file.path(KiliNP_Results, "KiliNP_Diversity_Comparison.tif"),
  overwrite = TRUE
)

png(file.path(KiliNP_Results, "KiliNP_Indices_Comparison.png"),
    width = 2560, height = 1440, res = 150)

plot(Kili_Comparison_Rasters)

dev.off()

### Assess index performance using vegetation ground truth ####

message("Assessing diversity indices against vegetation ground truth...")

# Ensure CRS matches

if (crs(KiliNP_ShannonH_Raster) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(KiliNP_ShannonH_Raster))
}

# Crop and mask diversity rasters to ground truth extent


Kili_Shannon_crop     <- crop(KiliNP_ShannonH_Raster, KiliNP_LandCover_Vector)
KiliNP_Classic_RaoQ_crop <- crop(KiliNP_Classic_RaoQ$window.3$alpha.2, KiliNP_LandCover_Vector)
Kili_Rao_TWDTW_crop   <- crop(Kili_Rao_TWDTW$window.3$alpha.2, KiliNP_LandCover_Vector)

Kili_Shannon_masked     <- mask(Kili_Shannon_crop, KiliNP_LandCover_Vector)
KiliNP_Classic_RaoQ_Masked <- mask(KiliNP_Classic_RaoQ_crop, KiliNP_LandCover_Vector)
Kili_Rao_TWDTW_masked   <- mask(Kili_Rao_TWDTW_crop, KiliNP_LandCover_Vector)

# Rasterise vegetation class

KiliNP_LandCover_Vector$CAT <- KiliNP_LandCover_Vector$decsr

KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  Kili_Shannon_masked,
  field = "CAT"
)

### Convert to dataframe for PERMANOVA ####

Kili_Indices_Comparison_Raster <- c(
  Kili_Shannon_masked,
  KiliNP_Classic_RaoQ_Masked,
  Kili_Rao_TWDTW_masked,
  KiliNP_LandCover_Raster
)

names(Kili_Indices_Comparison_Raster) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

Kili_Indices_Comparison_DF <- as.data.frame(
  Kili_Indices_Comparison_Raster,
  na.rm = TRUE
)

colnames(Kili_Indices_Comparison_DF) <- c(
  "ShannonH",
  "RaosQ_Classic",
  "RaosQ_TWDTW",
  "Veg_GroundTruth"
)

### PERMANOVA ####

PERMANOVA_Shannon <- adonis2(
  ShannonH ~ Veg_GroundTruth,
  data = Kili_Indices_Comparison_DF,
  permutations = 9999
)

PERMANOVA_RaosQ_Classic <- adonis2(
  RaosQ_Classic ~ Veg_GroundTruth,
  data = Kili_Indices_Comparison_DF,
  permutations = 9999
)

PERMANOVA_RaosQ_TWDTW <- adonis2(
  RaosQ_TWDTW ~ Veg_GroundTruth,
  data = Kili_Indices_Comparison_DF,
  permutations = 9999
)

Kili_PERMANOVA_Results <- data.frame(
  Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
  R2 = c(
    PERMANOVA_Shannon$R2[1],
    PERMANOVA_RaosQ_Classic$R2[1],
    PERMANOVA_RaosQ_TWDTW$R2[1]
  ),
  F = c(
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

print(Kili_PERMANOVA_Results)

message("Kilimanjaro analysis complete.")