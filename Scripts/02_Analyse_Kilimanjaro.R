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

# Load KiliNP_LandCover_Vector boundary

KiliNP_LandCover_Vector <- vect(paste0(KiliNP_Input,"/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp"))

# List raster files

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

# Define month names explicitly

month.names <- month.name  # built-in R constant

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

blank_layers

## Only pixels with a complete set of data for every layer are suitable for analysis
# There are many pixels with NA values scattered throughout the raster stack
# Create logical mask: TRUE only where ALL layers are non-NA

Kili.pixel.mask <- app(KiliNP_Timeseries, function(x) all(!is.na(x)))

# Mask out incomplete pixels (FALSE becomes NA)

KiliNP_Timeseries_Clean <- mask(KiliNP_Timeseries, Kili.pixel.mask, maskvalues = 0)
