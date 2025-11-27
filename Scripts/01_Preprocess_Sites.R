############################################################## ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025   ####
# 01_Preprocess_Sites.R - Crop Imagery for All Study Sites    ####
############################################################## ###

### Load project setup (libraries, reproducibility, path objects) ####
source("scripts/00_setup.R") # Run this if I haven't already

### Knepp Estate Preprocessing ####
# Load the Knepp Estate Area of Interest (AoI)
Knepp_AoI <- vect(file.path(Knepp_Raw, "/Schulte-to-Bühne_Knepp_Site_Data/knepp boundary lcc_pretty.shp"))

# Identify all Landsat GeoTIFFs for the Knepp Estate
Knepp_Landsat_Files <- list.files(
  "F:/Elliot Shayle/Knepp Estate Geodata", # This may change depending on where my photos are being stored
  pattern = "(?i)\\.tif$",   # (?i) = case insensitive regex
  full.names = TRUE
)

# Loop through each Landsat scene and crop it to the Knepp AOI
for (i in seq_along(Knepp_Landsat_Files)) {
  
  message("Processing Knepp file ", i, " of ", length(Knepp_Landsat_Files), ": ",
          basename(Knepp_Landsat_Files[i]))
  
  # Load Landsat raster
  tmp.landsat_raster <- rast(Knepp_Landsat_Files[i])
  
  # Reproject AoI if CRS differs from the raster CRS
  if (!identical(crs(tmp.landsat_raster), crs(Knepp_AoI))) {
    tmp.knepp_AoI_projected <- project(Knepp_AoI, crs(tmp.landsat_raster))
  } else {
    tmp.knepp_AoI_projected <- Knepp_AoI
  }
  
  # Crop the raster to the AOI and mask out non-AOI areas
  tmp.knepp_raster_cropped <- tmp.landsat_raster |>
    mask(tmp.knepp_AoI_projected) |>
    crop(tmp.knepp_AoI_projected)
  
  # Construct output filename
  tmp.knepp_cropped_filename <- file.path(
    #Knepp_Processed, #Re-enable this once I've sorted out my file storage
    "F:/Elliot Shayle/Temp Knepp Processed",
    paste0(
      tools::file_path_sans_ext(basename(Knepp_Landsat_Files[i])),
      "_Knepp_Cropped.tif"
    )
  )
  
  # Save the processed raster
  writeRaster(tmp.knepp_raster_cropped,
              tmp.knepp_cropped_filename,
              overwrite = T)
}

### ----------------------------- ###
### 2) Macchia Sacra Preprocessing ###
### ----------------------------- ###
# NOTE: Adjust steps here as needed for Sentinel-2, HRVPP, PPI, etc.

Macchia_AOI <- vect(file.path(Macchia_Raw, "Macchia_Sacra_AOI.shp"))

Macchia_Files <- list.files(
  Macchia_Raw,
  pattern = "\\.tif$",
  full.names = TRUE
)

for (i in seq_along(Macchia_Files)) {
  
  message("Processing Macchia Sacra file ", i, " of ", length(Macchia_Files), ": ",
          basename(Macchia_Files[i]))
  
  Macchia_Raster <- rast(Macchia_Files[i])
  
  if (!identical(crs(Macchia_Raster), crs(Macchia_AOI))) {
    Macchia_AOI_Projected <- project(Macchia_AOI, crs(Macchia_Raster))
  } else {
    Macchia_AOI_Projected <- Macchia_AOI
  }
  
  Macchia_Raster_Cropped <- Macchia_Raster |>
    mask(Macchia_AOI_Projected) |>
    crop(Macchia_AOI_Projected)
  
  Output_Filename <- file.path(
    Macchia_Processed,
    paste0(
      tools::file_path_sans_ext(basename(Macchia_Files[i])),
      "_Macchia_Cropped.tif"
    )
  )
  
  writeRaster(Macchia_Raster_Cropped,
              Output_Filename,
              overwrite = TRUE)
}

### ------------------------------ ###
### 3) Kilimanjaro NP Preprocessing ###
### ------------------------------ ###
# NOTE: This may require cloud masks or DEM-based clipping.
#       The structure below is ready for customisation.

KiliNP_AOI <- vect(file.path(KiliNP_Raw, "KiliNP_AOI.shp"))

KiliNP_Files <- list.files(
  KiliNP_Raw,
  pattern = "\\.tif$",
  full.names = TRUE
)

for (i in seq_along(KiliNP_Files)) {
  
  message("Processing Kilimanjaro file ", i, " of ", length(KiliNP_Files), ": ",
          basename(KiliNP_Files[i]))
  
  KiliNP_Raster <- rast(KiliNP_Files[i])
  
  if (!identical(crs(KiliNP_Raster), crs(KiliNP_AOI))) {
    KiliNP_AOI_Projected <- project(KiliNP_AOI, crs(KiliNP_Raster))
  } else {
    KiliNP_AOI_Projected <- KiliNP_AOI
  }
  
  KiliNP_Raster_Cropped <- KiliNP_Raster |>
    mask(KiliNP_AOI_Projected) |>
    crop(KiliNP_AOI_Projected)
  
  Output_Filename <- file.path(
    KiliNP_Processed,
    paste0(
      tools::file_path_sans_ext(basename(KiliNP_Files[i])),
      "_KiliNP_Cropped.tif"
    )
  )
  
  writeRaster(KiliNP_Raster_Cropped,
              Output_Filename,
              overwrite = TRUE)
}

### Script complete ###
message("All study sites have been successfully preprocessed.")
