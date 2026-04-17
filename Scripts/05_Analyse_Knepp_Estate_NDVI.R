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

Knepp_Time <- as.Date(paste0(Knepp.dates, "-15"))
str(Knepp_Time) # WARNING: Vector is shorter than length as raster stack includes multiple bands for each date

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

### Export the dataset as a NetCDF file for Gaussian Processing ####
## The Gaussian Process takes place separately in Python
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

## Calculate the mean seasonal trajectory across the site
# Crop to Knepp Estate site boundaries

Knepp_Timeseries.NDVI <- crop(Knepp_Buffered_Timeseries.NDVI, KneppEstate_Boundaries, mask = TRUE)

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

# #### Diversity analyses: ####
# ### 1. Shannon-Wiener Index ####
# # Shannon's H is first applied to the collapsed timeseries of data
# 
# message("Calculating Shannon-Wiener diversity index...")
# 
# Knepp_NDVI_Mean_Raster <- app(Knepp_Timeseries.NDVI, fun = mean, na.rm = TRUE) # This takes all raster layers, and computes the mean pixel value based on all available layers
# 
# # Round to 2 decimals to avoid numerical saturation
# 
# Knepp_NDVI_Mean_Raster2dec <- round(Knepp_NDVI_Mean_Raster, 2)
# 
# # Calculate Shannon with moving window = 3
# 
# Knepp.NDVI.ShannonH.matrix <- rasterdiv::ShannonS(
#   x = terra::as.matrix(Knepp_NDVI_Mean_Raster2dec, wide = TRUE),
#   window = 3,
#   na.tolerance = 0
# )
# 
# ## Now turn the matrix of Shannon H values into a spatial raster
# # Put it in a raster signature
# 
# Knepp_NDVI_ShannonH_Raster <- rast(Knepp.NDVI.ShannonH.matrix)
# 
# # Make the raster's extent and CRS match the original raster
# 
# ext(Knepp_NDVI_ShannonH_Raster) <- ext(Knepp_NDVI_Mean_Raster2dec)
# crs(Knepp_NDVI_ShannonH_Raster) <- crs(Knepp_NDVI_Mean_Raster2dec)
# 
# names(Knepp_NDVI_ShannonH_Raster) <- "Shannon's H (Derived from NDVI)"
# 
# # Export the raster
# 
# writeRaster(
#   Knepp_NDVI_ShannonH_Raster,
#   filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_ShannonH_Raster.tif"),
#   overwrite = TRUE
# )
# 
# ### 2. Classic Rao's Q  ####
# 
# message("Calculating classical Rao's Q...")
# 
# ## This function calculates parametric Rao's Q
# # This version calculates "classic" Rao's Q on only one matrix
# # The subsequent multidimensional paRao will be used to calculate TWDTW Rao's Q
# 
# Knepp_NDVI_Classic_RaoQ <- paRao(
#   x = Knepp_NDVI_Mean_Raster, # This also uses the mean of NDVI throughout time (i.e., all 146 observations are averaged into 1)
#   window = 3,
#   alpha = 2,
#   na.tolerance = 0,
#   simplify = 2,
#   method = "classic"
# )
# 
# # This function exports my object as a GeoTIF for viewing in QGIS etc.
# 
# writeRaster(
#   Knepp_NDVI_Classic_RaoQ$window.3$alpha.2, # Subsets to just the spatial raster (as far as I can tell)
#   filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_Classic-RaoQ_Raster.tif"),
#   overwrite = TRUE
# )
# 
# ### 3. Rao's Q with TWDTW ####
# ## Parameters are copied from the Hackathon preprint
# 
# message("Calculating Rao's Q with TWDTW distance...")
# 
# # If the function is taking too long to compute, make sure to enable parallelisation
# # For running on an HPC with unknown cores, I can use the function `parallel::detectCores()`
# # And set the "np" argument to detectCores() - 1
# # Also requires the package 'snow' to parallelise, so I've disabled it for now
# 
# Knepp_NDVI_TWDTW_RaoQ <- paRao(
#   x = Knepp_Timeseries,
#   time_vector = Knepp_Time,
#   window = 3,
#   alpha = 2,
#   na.tolerance = 0,
#   simplify = 2,
#   np = 6, # Number of cores to use (I'd rather not wait forever)
#   progBar = TRUE,
#   method = "multidimension",
#   dist_m = "twdtw",
#   midpoint = 6, # This is not the midpoint of the vector, rather the ecological midpoint of the cycle
#   stepness = -0.5, # I just noticed, shouldn't the α value be called "steepness" instead of "stepness"?
#   cycle_length = "year",
#   time_scale = "month" # this specifies that our midpoint, 35, occurs after 35 days
# )
# 
# # This function exports my object as a GeoTIF for viewing in QGIS etc.
# 
# writeRaster(
#   Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2,
#   filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_TWDTW-RaoQ_Raster.tif"),
#   overwrite = TRUE
# )
# 
# ### Export rasters for comparison as a stacked GeoTIF ####
# ## Stack outputs for easy visual comparison (as in Figure 3 of the Hackathon preprint)
# # Combine all my rasters into one list object
# 
# Knepp_NDVI_Index_Comparison_Rasters <- c(
#   Knepp_NDVI_Mean_Raster, # The mean of per pixel NDVI for straightforward visual analysis
#   Knepp_NDVI_ShannonH_Raster,
#   Knepp_NDVI_Classic_RaoQ$window.3$alpha.2,
#   Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2
# ) # This previously contained the NDVI raster, but I've decided to recompute in a separate file
# 
# # Name each layer
# 
# names(Knepp_NDVI_Index_Comparison_Rasters) <- c(
#   "Landsat Mean NDVI",
#   "Shannon's H (NDVI Derived)",
#   "Classic Rao's Q (NDVI Derived)",
#   "TWDTW Rao's Q (NDVI Derived)"
# )
# 
# # Export it for later viewing and observation in QGIS
# 
# writeRaster(
#   Knepp_NDVI_Index_Comparison_Rasters,
#   filename = file.path(Knepp_Results, "Knepp_Estate_NDVI_Diversity_Comparison.tif"),
#   overwrite = TRUE
# )
# 
# png(file.path(Knepp_Results, "Knepp_Estate_NDVI_Indices_Comparison.png"),
#     width = 2560, height = 1440, res = 150)
# 
# plot(Knepp_NDVI_Index_Comparison_Rasters) # A nice comparison plot
# 
# dev.off()
# 
# ### Assess and compare each index's performance ####
# ## Firstly, import the shape file of land cover classification
# # Path to the land cover classification shapefile from Matteo
# 
# Knepp.land.cover.file <- file.path(Knepp_Input, "Knepp_Data_From_Matteo/KneppEstate_Shapefile/KneppEstate_32632.shp")
# 
# # Read vector landcover
# 
# KneppEstate_LandCover_Vector <- vect(Knepp.land.cover.file)
# 
# # Ensure CRS matches diversity rasters
# 
# if (crs(Knepp_NDVI_Index_Comparison_Rasters) != crs(KneppEstate_LandCover_Vector)){
#   message("The shapefile's CRS differs from the diversity index raster's CRS")
#   KneppEstate_LandCover_Vector <- project(KneppEstate_LandCover_Vector, crs(Knepp_NDVI_Index_Comparison_Rasters)) # reprojects the shapefile if CRSs differ
# } else {message("The shapefile's CRS matches the diversity index raster's CRS")}
# 
# ## The rasters exceed the boundaries of the shapefile, so they will be cropped to size
# # Crop to polygon extent
# 
# cropped.Knepp_NDVI_ShannonH_Raster <- crop(Knepp_NDVI_ShannonH_Raster, KneppEstate_LandCover_Vector)
# cropped.Knepp_NDVI_Classic_RaoQ <- crop(Knepp_NDVI_Classic_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)
# cropped.Knepp_NDVI_TWDTW_RaoQ <- crop(Knepp_NDVI_TWDTW_RaoQ$window.3$alpha.2, KneppEstate_LandCover_Vector)
# 
# # Mask outside polygon
# 
# masked.Knepp_NDVI_ShannonH_Raster <- mask(cropped.Knepp_NDVI_ShannonH_Raster, KneppEstate_LandCover_Vector)
# masked.Knepp_NDVI_Classic_RaoQ <- mask(cropped.Knepp_NDVI_Classic_RaoQ, KneppEstate_LandCover_Vector)
# masked.Knepp_NDVI_TWDTW_RaoQ <- mask(cropped.Knepp_NDVI_TWDTW_RaoQ, KneppEstate_LandCover_Vector)
# 
# ## Make land cover classes rasterised instead of vectorised
# # Firstly, I need to give each polygon/vector within the spatial vector a proper category label
# # (currently, it's just a number corresponding to which description it is)
# 
# KneppEstate_LandCover_Vector$CAT <- KneppEstate_LandCover_Vector$decsr
# 
# # This defines a new raster object with the dominant land cover type as the dominant pixel value
# 
# KneppEstate_LandCover_Raster <- rasterize(
#   KneppEstate_LandCover_Vector,
#   masked.Knepp_NDVI_ShannonH_Raster,
#   field = "CAT"
# )
# 
# plot(KneppEstate_LandCover_Raster) # Looks cool
# 
# ## Put the pixel values into a dataframe for easier computation
# # Firstly, stack the spatial rasters into one object with vegetation class
# 
# masked.Knepp_NDVI_Index_Comparison_Rasters <- c(
#   masked.Knepp_NDVI_ShannonH_Raster,
#   masked.Knepp_NDVI_Classic_RaoQ,
#   masked.Knepp_NDVI_TWDTW_RaoQ,
#   KneppEstate_LandCover_Raster
# )
# 
# # Rename the different layers of the spatial raster to something more memorable
# 
# names(masked.Knepp_NDVI_Index_Comparison_Rasters) <- c(
#   "ShannonH (NDVI Derived)",
#   "Rao's Q Classic (NDVI Derived)",
#   "Rao's Q TWDTW (NDVI Derived)",
#   "Vegetation Ground Truth"
# )
# 
# # Convert the spatial raster to dataframe
# 
# masked.Knepp_NDVI_Index_Comparison_DF <- as.data.frame(masked.Knepp_NDVI_Index_Comparison_Rasters, na.rm = TRUE)
# 
# # Ensure vegetation is treated as factor [It is so I've disabled this line]
# # 
# # masked.Knepp_NDVI_Index_Comparison_DF$Vegetation <- as.factor(masked.Knepp_NDVI_Index_Comparison_DF$Vegetation)
# 
# # Rename the column names because otherwise R gets fussy
# 
# colnames(masked.Knepp_NDVI_Index_Comparison_DF) <- c(
#   "ShannonH",
#   "RaosQ_Classic",
#   "RaosQ_TWDTW",
#   "Veg_GroundTruth"
# )
# 
# ## Run the PERMANOVA
# # How much variance in each index is explained by each vegetation class?
# 
# Knepp_NDVI_PERMANOVA_Shannon <- adonis2(
#   masked.Knepp_NDVI_Index_Comparison_DF$ShannonH ~ 
#     masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
#   permutations = 9999
# )
# 
# Knepp_NDVI_PERMANOVA_RaosQ_Classic <- adonis2(
#   masked.Knepp_NDVI_Index_Comparison_DF$RaosQ_Classic ~ 
#     masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
#   permutations = 9999
# )
# 
# Knepp_NDVI_PERMANOVA_RaosQ_TWDTW <- adonis2(
#   masked.Knepp_NDVI_Index_Comparison_DF$RaosQ_TWDTW ~ 
#     masked.Knepp_NDVI_Index_Comparison_DF$Veg_GroundTruth,
#   permutations = 9999
# )
# 
# ## Extract R² and p-values
# 
# KneppEstate_NDVI_PERMANOVA_Results <- data.frame(
#   Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
#   R2 = c(
#     Knepp_NDVI_PERMANOVA_Shannon$R2[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_Classic$R2[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$R2[1]
#   ),
#   `F` = c(
#     Knepp_NDVI_PERMANOVA_Shannon$F[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_Classic$F[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$F[1]
#   ),
#   p_value = c(
#     Knepp_NDVI_PERMANOVA_Shannon$`Pr(>F)`[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_Classic$`Pr(>F)`[1],
#     Knepp_NDVI_PERMANOVA_RaosQ_TWDTW$`Pr(>F)`[1]
#   )
# )
# 
# print(KneppEstate_NDVI_PERMANOVA_Results)

### Yearly comparisons to CEH land-cover maps ####
## I have CEH land-cover data for 2000, 2007, 2015, and 2020
## To fairly compare the indices, I'm subsetting the spatial raster to those years
# Extract year from the time vector

Knepp_years <- format(Knepp_Time, "%Y")

# Define Years of Interest (Abbreviation: YoI)

YoI <- c("2000", "2007", "2015", "2020")

# Create list of NDVI stacks per year

Knepp_NDVI_YoI <- lapply(YoI, function(y) {
  Knepp_Buffered_Timeseries.NDVI[[Knepp_years == y]]
})

# Name the list elements

names(Knepp_NDVI_YoI) <- YoI

### Clean the dataset so that only pixels with complete data across all layers are retained ####
## First, a brief visual overview of the data
# This is optional, but it can be helpful to understand the data before applying the cleaning function

for(i in 1:nlyr(Knepp_NDVI_YoI[[4]])) {
  plot(Knepp_NDVI_YoI[[4]][[i]], main = names(Knepp_NDVI_YoI[[4]])[i])
  if(i < nlyr(Knepp_NDVI_YoI[[4]])) readline(prompt = "Press [enter] to continue")
}

# TEMPORARY line to drop 2 dodgy layers from the 2015 stack (these layers have >90% NA values, which is a problem for the subsequent analyses)

Knepp_NDVI_YoI[[3]] <- Knepp_NDVI_YoI[[3]][[1:10]]

# Now a function to create a mask of pixels with complete data across all layers, and apply it to the stack 

Knepp_NDVI_YoI_Clean <- mapply(function(r_stack, y) {
  
  message("Processing year: ", y)
  
  message("Creating complete-case mask...")
  
  # TRUE where all layers are non-NA
  
  pixel_mask <- app(r_stack, function(x) all(!is.na(x)))
  
  message("Applying mask...")
  
  r_clean <- mask(r_stack, pixel_mask, maskvalues = 0)
  
  ## Save cleaned rasters
  
  out_path <- file.path(
    Knepp_Processed,
    paste0("Knepp_Buffered_NDVI_", y, "_Clean.tif")
  )
  
  message("Writing cleaned raster to disk: ", out_path)
  
  writeRaster(
    r_clean,
    filename = out_path,
    overwrite = TRUE
  )
  
  return(r_clean)
  
}, Knepp_NDVI_YoI, names(Knepp_NDVI_YoI), SIMPLIFY = FALSE)

## Load  back in these rasters instead, if necessary and previously computed

Knepp_NDVI_YoI_Clean <- lapply(YoI, function(y) {
  rast(file.path(Knepp_Processed, paste0("Knepp_Buffered_NDVI_", y, "_Clean.tif")))
})

### Mosaic tiles for export and parallelisation on HPC (MaRC3a) ####
## Due to the size of the site and the need to run concurrently on 4 discrete years of data
## I will mosaic the dataset into 1000 tiles (similar to the Kili analysis) for later processing on the HPC

Knepp_NDVI_tile_info <- list()  # store metadata for later

for (y in YoI) {
  
  message("Creating tiles for year: ", y)
  
  # Get rasters
  #tmp.r_stack <- Knepp_NDVI_YoI[[y]] # This is the original, uncleaned stack. It will be reactivated when I have gap-filled data
  tmp.r_stack <- Knepp_NDVI_YoI_Clean[[y]] # This is the cleaned stack, with only pixels with complete data across all layers
  tmp.mean_raster <- app(tmp.r_stack, mean, na.rm = TRUE)
  
  # Trim outer NA borders (important for efficiency)
  tmp.trimmed_mean <- trim(tmp.mean_raster)
  
  # -------------------------------
  # DEFINE TILE GRID
  # -------------------------------
  
  knepp.tile.count <- 500  # same philosophy as Kili, but half as many tiles
  
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
  
  # -------------------------------
  # CREATE GRID
  # -------------------------------
  
  knepp.tiling.grid <- as.polygons(
    rast(
      ext(tmp.trimmed_mean),
      ncols = knepp.cols,
      nrows = knepp.rows,
      crs = crs(tmp.trimmed_mean)
    )
  )
  
  # Save grid (important for reproducibility + debugging)
  
  knepp.tiles.filepath <- file.path(Knepp_Processed, paste0("Knepp_Tiling_Grid_", y, ".geojson"))
  
  writeVector(knepp.tiling.grid, knepp.tiles.filepath, filetype = "GeoJSON", overwrite = TRUE)
  
  knepp.tiling.grid <- vect(knepp.tiles.filepath)
  
  # -------------------------------
  # TILE SETTINGS
  # -------------------------------
  
  RaoQ.window.size <- 3
  tmp.tile.overlap <- floor(RaoQ.window.size / 2)
  
  # -------------------------------
  # OUTPUT DIRECTORIES
  # -------------------------------
  
  tmp.knepp.mean.tile.dir <- file.path(Knepp_Processed, paste0("Knepp_", y, "_MeanNDVI_Tiles"))
  tmp.knepp.ts.tile.dir   <- file.path(Knepp_Processed, paste0("Knepp_", y, "_NDVI_Timeseries_Tiles"))
  
  dir.create(tmp.knepp.mean.tile.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp.knepp.ts.tile.dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------
  # SAVE TIME VECTOR (CRITICAL FOR TWDTW)
  # -------------------------------
  
  tmp.time.vector <- Knepp_Time[Knepp_years == y]
  
  writeLines(
    as.character(tmp.time.vector),
    file.path(tmp.knepp.ts.tile.dir, "time_vector.txt")
  )
  
  # -------------------------------
  # CREATE TILES
  # -------------------------------
  
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
  
  # -------------------------------
  # STORE METADATA (CRUCIAL LATER)
  # -------------------------------
  
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

# Now load in the classical Rao's Q tiles

for (y in c("2000", "2020")) {
  
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

# And now load in the TWDTW Rao's Q tiles

for (y in c("2000", "2020")) {
  
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

### Now I must compute Shannon's H and both Rao's Qs for each of those yearly rasters
## To keep the code neat, here's a function which uses the `rasterdiv` functions to compute these
#
# Knepp_NDVI_YearlyDiversityIndices <- function(r_stack, time_vector) {
#   
#   ## Mean NDVI (collapse time)
#   
#   mean_ndvi <- app(r_stack, mean, na.rm = TRUE)
#   
#   # Round to avoid Shannon saturation
#   
#   mean_ndvi_2dec <- round(mean_ndvi, 2)
#   
#   ## Compute Shannon's H
#   # This follows the same structure as elsewhere in this script
#   
#   shannon_mat <- rasterdiv::ShannonS(
#     x = terra::as.matrix(mean_ndvi_2dec, wide = TRUE),
#     window = 3,
#     na.tolerance = 0
#   )
#   
#   shannon_rast <- rast(shannon_mat)
#   ext(shannon_rast) <- ext(mean_ndvi_2dec)
#   crs(shannon_rast) <- crs(mean_ndvi_2dec)
#   
#   ## Classic Rao's Q
#   # `simplify = 2`, equivalent to 2 decimal places
#   
#   rao_classic <- paRao(
#     x = mean_ndvi,
#     window = 3,
#     alpha = 2,
#     na.tolerance = 0,
#     simplify = 2,
#     method = "classic"
#   )$window.3$alpha.2
#   
#   ## TWDTW Rao's Q (uses full time series)
#   # `simplify = 2`, equivalent to 2 decimal places
#   
#   rao_twdtw <- paRao(
#     x = r_stack,
#     time_vector = time_vector,
#     window = 3,
#     alpha = 2,
#     na.tolerance = 0,
#     simplify = 2,
#     method = "multidimension",
#     dist_m = "twdtw",
#     midpoint = 6,
#     stepness = -0.5,
#     cycle_length = "year",
#     time_scale = "month"
#   )$window.3$alpha.2
#   
#   ## Return everything nicely bundled
#   
#   return(list(
#     Knepp_Mean_NDVI = mean_ndvi,
#     Knepp_NDVI_ShannonsH = shannon_rast,
#     Knepp_NDVI_Classic_RaoQ = rao_classic,
#     Knepp_NDVI_TWDTW_RaoQ = rao_twdtw
#   ))
# }
# 
# ## And now I apply this newly defined function to each year of Knepp Estate data
# # This function goes over each year in the raster stack and computes the indices
# 
# Knepp_Indices_YoI.NDVI_derived <- lapply(names(Knepp_NDVI_YoI), function(y) {
#   
#   r <- Knepp_NDVI_YoI[[y]]
#   t <- Knepp_Time[Knepp_years == y]
#   
#   Knepp_NDVI_YearlyDiversityIndices(r, t)
# })
# 
# names(Knepp_Indices_YoI.NDVI_derived) <- names(Knepp_NDVI_YoI)

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

# Align the computed indices to the CEH land cover data, and name the layers accordingly

Knepp_LandCover_AlignedIndices_YoI.NDVI <- lapply(names(Knepp_Indices_YoI.NDVI_derived), function(y) {
  
  idx_list <- Knepp_Indices_YoI.NDVI_derived[[y]]
  lc <- Knepp_CEH_LandCover_YoI[[y]]
  
  # Using the CEH data as a template because it is higher resolution and a wider spatial extent
  
  aligned_stack <- c(
    resample(idx_list$Knepp_NDVI_ShannonsH, lc, method = "bilinear"),
    resample(idx_list$Knepp_NDVI_Classic_RaoQ, lc, method = "bilinear"),
    resample(idx_list$Knepp_NDVI_TWDTW_RaoQ, lc, method = "bilinear"),
    lc
  )
  
  names(aligned_stack) <- c(
    "Shannon's H (NDVI derived)",
    "Classic Rao's Q (NDVI derived)",
    "TWDTW Rao's Q (NDVI derived)",
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
    
    cropped <- crop(stack, KneppEstate_Buffered_Boundaries)
    masked  <- mask(cropped, KneppEstate_Buffered_Boundaries)
    
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
    list(
      
      # These are structured in the same way as the PERMANOVAe from the other case studies
      
      Shannon = adonis2(
        ShannonH ~ LandCover,
        data = df,
        permutations = 9999
      ),
      
      Rao_Classic = adonis2(
        RaosQ_Classic ~ LandCover,
        data = df,
        permutations = 9999
      ),
      
      Rao_TWDTW = adonis2(
        RaosQ_TWDTW ~ LandCover,
        data = df,
        permutations = 9999
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
