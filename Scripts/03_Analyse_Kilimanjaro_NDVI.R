############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 24/02/2026                    #
# 03_Analyse_Kilimanjaro_NDVI.R                                                #
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

# This script often requires parallelised computation, so this defines cores for computation

kili.cores <- max(1, detectCores() - 2)

# Load KiliNP_LandCover_Vector boundary (this is our land cover ground truth data)

KiliNP_LandCover_Vector <- 
  vect(file.path(KiliNP_Input,
                 "/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp")) # Load in the ground truth data

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

### Savitzky-Golay Gap Filling function ####

message("Running Savitzky-Golay filter for parallel non-masked pipeline...")

sg_gapfill <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  x_interp <- zoo::na.approx(x, na.rm = FALSE, rule = 2)
  x_smoothed <- pracma::savgol(x_interp, fl = 5, forder = 2)
  return(x_smoothed)
}

# Apply to the raw, unmasked timeseries

KiliNP_Timeseries_SG <- app(KiliNP_Timeseries, fun = sg_gapfill, cores = max(1, detectCores() - 2))
names(KiliNP_Timeseries_SG) <- names(KiliNP_Timeseries)

writeRaster(KiliNP_Timeseries_SG, file.path(KiliNP_Processed, "KiliNP_NDVI_2017-2021_Cropped_&_SGfiltered.tif"), overwrite = TRUE) # Export raster
KiliNP_Timeseries_SG <- rast(file.path(KiliNP_Processed, "KiliNP_NDVI_2017-2021_Cropped_&_SGfiltered.tif")) # Load it back in

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

writeRaster(KiliNP_Timeseries_Clean, file.path(KiliNP_Processed, "KiliNP_2017-2021_Cropped_&_Masked.tif"), overwrite = TRUE) # Export raster

# And then load back in the raster

KiliNP_Timeseries_Clean <- rast(file.path(KiliNP_Processed, "KiliNP_2017-2021_Cropped_&_Masked.tif")) # Load it in

### Inspect temporal structure ####

## Unlike Macchia Sacra NetCDF, this GeoTIFF stack does not contain explicit time metadata
## Therefore, we must construct the time vector manually from layer names

message("Constructing time vector for Kilimanjaro time series...")

# Extract year and month from layer names
# Layer format: "2017 - January"

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

message("Calculating Shannon-Wiener diversity index for Masked and SG tracks...")

## For the masked raster

KiliNP_Mean_Masked <- app(KiliNP_Timeseries_Clean, fun = mean, na.rm = TRUE) # As with Macchia Sacra, collapse the raster to a single mean
writeRaster(KiliNP_Mean_Masked, file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Masked.tif"), overwrite = TRUE)
KiliNP_Mean_Masked <- rast(file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Masked.tif")) # Load it back in!

# Simplify the raster to just 2 decimal places to prevent numerical oversaturation

KiliNP_Mean_Masked2dec <- trim(round(KiliNP_Mean_Masked, 2))

# Reduce raster to a simple data matrix (required by `ShannonS`) and then apply the Shannon's H function to it

KiliNP.ShannonH.Masked.matrix <- rasterdiv::ShannonS(terra::as.matrix(KiliNP_Mean_Masked2dec, wide = TRUE), window = 3, na.tolerance = 0)

# Convert the Shannon's H matrix back into a raster and apply a CRS to it

KiliNP_ShannonH_Masked <- rast(KiliNP.ShannonH.Masked.matrix)
ext(KiliNP_ShannonH_Masked) <- ext(KiliNP_Mean_Masked2dec)
crs(KiliNP_ShannonH_Masked) <- crs(KiliNP_Mean_Masked2dec)
names(KiliNP_ShannonH_Masked) <- "ShannonH_Masked"
KiliNP_ShannonH_Masked <- terra::extend(KiliNP_ShannonH_Masked, KiliNP_Mean_Masked) # Extend the masked raster, filling the uncomputed boundary tiles with NAs

# Export the raster for safe keeping

writeRaster(KiliNP_ShannonH_Masked, file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_Masked.tif"), overwrite = TRUE)
KiliNP_ShannonH_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_Masked.tif")) # And load it back in so I don't have to recompute it each time

## For the gap-filled raster
# As above, so below

KiliNP_Mean_SG <- app(KiliNP_Timeseries_SG, fun = mean, na.rm = TRUE)
writeRaster(KiliNP_Mean_SG, file.path(KiliNP_Processed, "KiliNP_MeanNDVI_SG.tif"), overwrite = TRUE)
KiliNP_Mean_SG <- rast(file.path(KiliNP_Processed, "KiliNP_MeanNDVI_SG.tif")) # Load it back in!

KiliNP_Mean_SG2dec <- trim(round(KiliNP_Mean_SG, 2))
KiliNP.ShannonH.SG.matrix <- rasterdiv::ShannonS(terra::as.matrix(KiliNP_Mean_SG2dec, wide = TRUE), window = 3, na.tolerance = 0)

KiliNP_ShannonH_SG <- rast(KiliNP.ShannonH.SG.matrix)
ext(KiliNP_ShannonH_SG) <- ext(KiliNP_Mean_SG2dec)
crs(KiliNP_ShannonH_SG) <- crs(KiliNP_Mean_SG2dec)
names(KiliNP_ShannonH_SG) <- "ShannonH_SG"

writeRaster(KiliNP_ShannonH_SG, file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_SG.tif"), overwrite = TRUE)
KiliNP_ShannonH_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_SG.tif")) # Load it back in!

### 2. Classic Rao's Q  ####
## Due to the large size of the raster, I need to tile it so that it can be run
## The tiles will be stitched back together once they're computed

message("Calculating classical Rao's Q for Kilimanjaro...")

### Step 1: Create a grid to define zones for tiling

# Optional but recommended: trim outer NA borders

trimmed.KiliNP_Mean_Raster <- trim(KiliNP_Timeseries_SG)

# We want approximately 72 tiles (not too many, not too few)
# Actually, I think this is far too many, having spent the last days trying to compute them all
# For effective parallelisation on MaRC3a, I actually think 2000 tiny tiles would be more effective

#kili.total.tiles <- 72
kili.total.tiles <- 2000
kili.aspect.ratio <- ncol(trimmed.KiliNP_Mean_Raster) / nrow(trimmed.KiliNP_Mean_Raster)

# find factor pairs of 72 

tiling.factors <- expand.grid(
  ncols = 1:kili.total.tiles,
  nrows = 1:kili.total.tiles
)

# Subset to just pairs which equal "kili.total.tiles" when multiplied

tiling.factors <- tiling.factors[tiling.factors$ncols * tiling.factors$nrows == kili.total.tiles, ]

## choose the pair closest to the raster's aspect ratio
# First, create a new column calculating the difference in aspect ratio

tiling.factors$ratio_diff <- abs((tiling.factors$ncols / tiling.factors$nrows) - kili.aspect.ratio)

# Now find the pair which has the lowest distance from a square 1:1 aspect ratio

best.tile.size <- tiling.factors[which.min(tiling.factors$ratio_diff), ]

# And save that for later usage

kili.cols <- best.tile.size$ncols
kili.rows <- best.tile.size$nrows

# Create spatial polygons to specify what the tile sizes should be

kili.tiling.grid <- as.polygons(
  rast(
    ext(trimmed.KiliNP_Mean_Raster),
    ncols = kili.cols,
    nrows = kili.rows,
    crs = crs(trimmed.KiliNP_Mean_Raster)
  )
)

writeVector(kili.tiling.grid, file.path(KiliNP_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON"), filetype = "GeoJSON" , overwrite = TRUE) # Export for later use

kili.tiling.grid <- vect(file.path(KiliNP_Processed, "KiliNP_Tiling_Grid_Polygons.geoJSON")) # Load it back in

plot(kili.tiling.grid) # Plot it to make sure that it's loaded in

# Window size (I think 3 is the default, but this can be changed as necessary)

RaoQ.window.size <- 3
kili.tile.overlap <- floor(RaoQ.window.size / 2)

# 
# 
# #kili.tile.dir <- file.path(KiliNP_Processed,"Mean NDVI tiles") # I disabled this line so I don't overwrite my 72 larger tiles
# kili.tile.dir <- file.path(KiliNP_Processed,"Tiny tiles")
# dir.create(kili.tile.dir, recursive = TRUE, showWarnings = FALSE)
# 
# ## Finally, create the tiles
# 
# kili.tiles <- makeTiles(
#   trimmed.KiliNP_Mean_Raster, # Trimmed for easier computation
#   y = kili.tiling.grid, # Specifies how the tiles should be allocated
#   buffer = kili.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
#   filename = file.path(kili.tile.dir, "KiliNP_MeanNDVI_Tile-.tif"),
#   overwrite = TRUE
# )

## Create a directory to put the tiles in
# For the masked raster

kili.tile.dir.masked <- file.path(KiliNP_Processed, "Kili_Tiles_NDVI_Masked")
dir.create(kili.tile.dir.masked, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Mean_Masked), # Trim it for easier computation
  y = kili.tiling.grid, # Specifies how the tiles should be allocated
  buffer = kili.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
  filename = file.path(kili.tile.dir.masked, "KiliNP_MeanNDVI_Masked_Tile-.tif"),
  overwrite = TRUE
)

# For the gap-filled raster

kili.tile.dir.sg <- file.path(KiliNP_Processed, "Kili_Tiles_NDVI_SG-Filtered")
dir.create(kili.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Mean_SG), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.tile.dir.sg, "KiliNP_MeanNDVI_SG_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Compute classic Rao's Q for each tile
## Firstly, I need to setup the environment for parallelisation
# Create a subfolder to store the classic Rao's Q output tiles

# kili.rao.dir  <- file.path(kili.tile.dir, "rao-utputs") 
# dir.create(kili.rao.dir, recursive = TRUE, showWarnings = FALSE)
# 
# ## Create a computing cluster to parallelise the calculation at the tile level
# # Set the number of cores to be used by the cluster
# 
# kili.cores <- max(1, detectCores() - 2)
# 
# # Initialise a log file so I can actually see what's going on
# 
# kili.log.file <- file.path(kili.rao.dir, "KiliNP_RaoQ_processing_log.txt")
# 
# # If the log file doesn't exist already, create one
# 
# if(!file.exists(kili.log.file)) file.create(kili.log.file)
# 
# # Create the cluster (alliterative and punny names are mandatory)
# 
# kili.cluster <- makeCluster(kili.cores)
# 
# clusterEvalQ(kili.cluster, {
#   library(terra)
#   library(rasterdiv)
# })
# 
# clusterExport(kili.cluster, c(
#   "kili.tiles",
#   "kili.rao.dir",
#   "RaoQ.window.size",
#   "kili.log.file"
# ))
# 
# # Identify tiles still needing processing (so resources aren't wasted processing tiles already done)
# 
# tile.outputs <- file.path(
#   kili.rao.dir,
#   paste0("KiliNP_Classic-RaoQ_Tile-", seq_along(kili.tiles), ".tif")
# )
# 
# tiles.to.process <- which(!file.exists(tile.outputs))
# 
# cat(length(tiles.to.process), "tiles remaining.\n")
# 
# ## Now actually run the code
# # This version creates a process for each CPU core and runs each tile as a single process
# ## REVIEWERS: Due to the CPU overhead of this workload, I ultimately decided to run it on my university's supercomputer instead
# # Please see "03.1C_Kilimanjaro_Classical-RaoQ_MaRC3a.R" for the job file I submitted
# 
# kili.classic.rao.results <- parLapply( # Function call
#   kili.cluster,
#   tiles.to.process,
#   function(i){
#     
#     library(terra)
#     library(rasterdiv)
#     
#     log_file <- kili.log.file
#     
#     log_msg <- function(msg){
#       cat(
#         paste0(Sys.time(), " | Worker ", Sys.getpid(), " | ", msg, "\n"),
#         file = log_file,
#         append = TRUE
#       )
#     }
#     
#     out.file <- file.path(
#       kili.rao.dir,
#       paste0("KiliNP_Classic-RaoQ_Tile-", i, ".tif")
#     )
#     
#     if(file.exists(out.file)){
#       log_msg(paste0("Tile", i, "already exists — skipped"))
#       return(NULL)
#     }
#     
#     log_msg(paste("Tile", i, "STARTED"))
#     
#     tmp.tile <- rast(kili.tiles[i]) # Load in the raster for processing
#     
#     tmp.result <- paRao(
#       tmp.tile,
#       window = RaoQ.window.size,
#       alpha = 2,
#       simplify = 2, # This is necessary to maintain consistency with the Shannon's H test (keeps just 2 decimal places)
#       method = "classic", # Because this is not looking at timeseries Rao's Q, just regular unidimensional Rao's Q
#       np = 1 # Explicitly prevents nested parallelisation (or set above 1 if you want to melt your CPU)
#     )
#     
#     tmp.rao_raster <- tmp.result[[1]][[1]] # Subsetting avoids hardcoding "$window.3$alpha.2"
#     
#     writeRaster(
#       tmp.rao_raster,
#       filename = out.file,
#       overwrite = TRUE
#     )
#     
#     rm(tmp.tile,tmp.result,tmp.rao_raster)
#     gc()
#     
#     log_msg(paste("Tile №", i, "'s classic Rao's Q calculated successfully."))
#     
#     return(NULL) # So that each worker doesn't fill up R's memory with bloat upon completion
#   }
# )
# 
# ## This for loop is an alternative computational approach which uses all cores to work sequentially over each tile
# # This version seems computationally safer because each tile outputted is like a mini-checkpoint in the event that computation is interrupted
# 
# for(i in seq_along(kili.tiles)){
#   
#   log_file <- kili.log.file
#   
#   log_msg <- function(msg){
#     cat(
#       paste0(Sys.time(), " | Tile ", i, " | ", msg, "\n"),
#       file = log_file,
#       append = TRUE
#     )
#   }
#   
#   out.file <- file.path(
#     kili.rao.dir,
#     paste0("KiliNP_Classic-RaoQ_Tile-", i, ".tif")
#   )
#   
#   # Skip tiles which already exist (prevents recomputation)
#   
#   if(file.exists(out.file)){
#     log_msg("already exists — skipped")
#     next
#   }
#   
#   log_msg("STARTED")
#   
#   tmp.tile <- rast(kili.tiles[i]) # Load in the raster for processing
#   
#   tmp.result <- paRao(
#     tmp.tile,
#     window = RaoQ.window.size,
#     alpha = 2,
#     simplify = 2, # This is necessary to maintain consistency with the Shannon's H test (keeps just 2 decimal places)
#     method = "classic", # Because this is not looking at timeseries Rao's Q, just regular unidimensional Rao's Q
#     np = kili.cores # Parallelise INSIDE paRao for faster per-tile processing
#   )
#   
#   tmp.rao_raster <- tmp.result[[1]][[1]] # Subsetting avoids hardcoding "$window.3$alpha.2"
#   
#   writeRaster(
#     tmp.rao_raster,
#     filename = out.file,
#     overwrite = TRUE
#   )
#   
#   rm(tmp.tile,tmp.result,tmp.rao_raster)
#   gc()
#   
#   log_msg("classic Rao's Q calculated successfully.")
# }

### Step 3: Demosaic the classical Rao's Q tiles
## Gather up all the files

# kili.rao.files <- list.files(
#   kili.rao.dir,
#   pattern = "Classic-RaoQ",
#   full.names = TRUE
# )
# 
# # Tell R to apply the `rast` function to them 
# 
# kili.rao.tiles <- lapply(kili.rao.files, rast)
# 
# # Convert them into a spatial raster collection:
# 
# kili.rao.tiles <- sprc(kili.rao.tiles)
# 
# # Run the demosaic function (which is curiously called `mosaic`)
# 
# KiliNP_Classic_RaoQ <- terra::mosaic(kili.rao.tiles)
# 
# # Export the final raster, and load it back in if necessary
# 
# writeRaster(
#   KiliNP_Classic_RaoQ,
#   file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ.tif"),
#   overwrite = TRUE
# )
# 
# KiliNP_Classic_RaoQ <- rast(file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ.tif")) # Load in the raster
# 
# plot(KiliNP_Classic_RaoQ) # Plot it! (Good for data exploration and checking that the raster loaded in as normal)

### Step 3: Demosaic the classical Rao's Q tiles
# For the gap-filled raster

kili.rao.files.sg <- list.files( # Gather up all the files
  file.path(KiliNP_Processed, "Kili_Tiles_NDVI_SG-Filtered", "Rao-utputs"), # Amend if outputs are in a subfolder
  pattern = "KiliNP_Classic-RaoQ_SG_Tile-",
  full.names = TRUE
)
KiliNP_Classic_RaoQ_SG <- terra::mosaic(sprc(lapply(kili.rao.files.sg, rast)))
writeRaster(KiliNP_Classic_RaoQ_SG, file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ_SG.tif"), overwrite = TRUE)
KiliNP_Classic_RaoQ_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ_SG.tif"))
plot(KiliNP_Classic_RaoQ_SG)

# For the masked raster

kili.rao.files.masked <- list.files( # Gather up all the files
  file.path(KiliNP_Processed, "Kili_Tiles_NDVI_Masked", "Rao-utputs"), # Amend if outputs are in a subfolder
  pattern = "KiliNP_Classic-RaoQ_Masked_Tile-",
  full.names = TRUE
)

KiliNP_Classic_RaoQ_Masked <- terra::mosaic(sprc(lapply(kili.rao.files.masked, rast))) # Load in and demosaic the computed tiles
KiliNP_Classic_RaoQ_Masked <- terra::extend(KiliNP_Classic_RaoQ_Masked, KiliNP_Classic_RaoQ_SG) # Extend the masked raster, filling the uncomputed boundary tiles with NAs
writeRaster(KiliNP_Classic_RaoQ_Masked, file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ_Masked.tif"), overwrite = TRUE) # Save it for later
KiliNP_Classic_RaoQ_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_Classic-RaoQ_Masked.tif")) # Load it back in
plot(KiliNP_Classic_RaoQ_Masked)

### 3. Rao's Q with TWDTW ####

message("Calculating Rao's Q with TWDTW distance for Kilimanjaro...")

## Step 1: I'll have to tile this as well because it is too large to compute as a single object
# This tiling script is copied from Step 1 of the classical Rao's Q analysis
# Some objects like "kili.tiling.grid" are assumed to be loaded

# # Create a directory to put the tiles in
# 
# kili.twdtw.rao.dir <- file.path(KiliNP_Processed,"Timeseries NDVI tiles")
# dir.create(kili.twdtw.tile.dir, recursive = TRUE, showWarnings = FALSE)
# 
# # Create the timeseries tiles
# 
# kili.twdtw.tiles <- makeTiles(
#   trim(KiliNP_Timeseries_Clean), # Trimmed for easier computation
#   y = kili.tiling.grid, # Make sure this is still loaded in from the previous step!
#   buffer = kili.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
#   filename = file.path(kili.tile.dir, "KiliNP_2017-2021_NDVI_Tile-.tif"),
#   overwrite = TRUE
# )

# For the masked raster

kili.twdtw.tile.dir.masked <- file.path(KiliNP_Processed, "Kili_TS_Tiles_NDVI_Masked")
dir.create(kili.twdtw.tile.dir.masked, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Timeseries_Clean), # Trim it for easier computation
  y = kili.tiling.grid, # Make sure this is still loaded in from the previous step!
  buffer = kili.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
  filename = file.path(kili.twdtw.tile.dir.masked, "KiliNP_TS_NDVI_Masked_Tile-.tif"),
  overwrite = TRUE
)

# For the gap-filled raster

kili.twdtw.tile.dir.sg <- file.path(KiliNP_Processed, "Kili_TS_Tiles_NDVI_SG-Filtered")
dir.create(kili.twdtw.tile.dir.sg, recursive = TRUE, showWarnings = FALSE)

makeTiles(
  trim(KiliNP_Timeseries_SG), 
  y = kili.tiling.grid, 
  buffer = kili.tile.overlap, 
  filename = file.path(kili.twdtw.tile.dir.sg, "KiliNP_TS_NDVI_SG_Tile-.tif"),
  overwrite = TRUE
)

## Step 2: Submit the tiles for processing on University of Marburg's supercomputer
# Please see script 03.2A_Kilimanjaro_TWDTW-RaoQ_MaRC3a.R for the actual job scripts used
# If you wish to attempt computation on your local machine, then please see the code below

######### Beginning of not actually used section
# Kili_Rao_TWDTW <- paRao(
#   x = KiliNP_Timeseries_Clean,
#   time_vector = Kili.dates,
#   window = 3,
#   alpha = 2,
#   na.tolerance = 0,
#   simplify = 2,
#   np = detectCores() -1,
#   progBar = TRUE,
#   method = "multidimension",
#   dist_m = "twdtw",
#   midpoint = 6,          # Midpoint of annual cycle (June)
#   stepness = -0.5,
#   cycle_length = "year",
#   time_scale = "month"   # Now explicitly monthly data
# )
# 
# writeRaster(
#   Kili_Rao_TWDTW$window.3$alpha.2,
#   filename = file.path(KiliNP_Results, "KiliNP_RaoQ_TWDTW.tif"),
#   overwrite = TRUE
# )
######### End of not actually used section

## Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster

# ## Gather up all the files
# 
# kili.twdtw.rao.files <- list.files(
#   file.path(kili.twdtw.rao.dir, "TWDTW Rao-utputs"),
#   pattern = "KiliNP_2017-2021_TWDTW-RaoQ_Tile-",
#   full.names = TRUE
# )
# 
# # Tell R to apply the `rast` function to them 
# 
# kili.twdtw.rao.files <- lapply(kili.twdtw.rao.files, rast)
# 
# # Convert them into a spatial raster collection:
# 
# kili.twdtw.rao.files <- sprc(kili.twdtw.rao.files)
# 
# # Run the demosaic function (which is curiously called `mosaic`)
# 
# KiliNP_TWDTW_RaoQ <- terra::mosaic(kili.twdtw.rao.files)
# 
# # Export the final raster, and load it back in if necessary
# 
# writeRaster(
#   KiliNP_TWDTW_RaoQ,
#   file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ.tif"),
#   overwrite = TRUE
# )
# 
# KiliNP_TWDTW_RaoQ <- rast(file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ.tif")) # Load in the raster
# 
# plot(KiliNP_TWDTW_RaoQ)

## Step 3: Demosaic the raster tiles to create a final TWDTW Rao's Q raster

# For the gap-filled raster

kili.twdtw.rao.files.sg <- list.files(
  file.path(KiliNP_Processed, "Kili_TS_Tiles_NDVI_SG-Filtered", "Rao-utputs"), 
  pattern = "KiliNP_TWDTW-RaoQ_SG_Tile-",
  full.names = TRUE
)
KiliNP_TWDTW_RaoQ_SG <- terra::mosaic(sprc(lapply(kili.twdtw.rao.files.sg, rast)))
writeRaster(KiliNP_TWDTW_RaoQ_SG, file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ_SG.tif"), overwrite = TRUE) # Save it for later
KiliNP_TWDTW_RaoQ_SG <- rast(file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ_SG.tif")) # Load it back in
plot(KiliNP_TWDTW_RaoQ_SG)

# For the masked raster

kili.twdtw.rao.files.masked <- list.files(
  file.path(KiliNP_Processed, "Kili_TS_Tiles_NDVI_Masked", "Rao-utputs"), 
  pattern = "KiliNP_TWDTW-RaoQ_Masked_Tile-",
  full.names = TRUE
)
KiliNP_TWDTW_RaoQ_Masked <- terra::mosaic(sprc(lapply(kili.twdtw.rao.files.masked, rast)))
KiliNP_TWDTW_RaoQ_Masked <- terra::extend(KiliNP_TWDTW_RaoQ_Masked, KiliNP_TWDTW_RaoQ_SG) # Extend the masked raster, filling the uncomputed boundary tiles with NAs
writeRaster(KiliNP_TWDTW_RaoQ_Masked, file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ_Masked.tif"), overwrite = TRUE) # Save it for later
KiliNP_TWDTW_RaoQ_Masked <- rast(file.path(KiliNP_Results, "Kilimanjaro_TWDTW-RaoQ_Masked.tif")) # Load it back in
plot(KiliNP_TWDTW_RaoQ_Masked)

### Export all rasters for comparison ####

KiliNP_Comparison_Rasters <- c(
  KiliNP_Mean_Masked, # Trimmed so that it matches the extent of the other rasters
  KiliNP_ShannonH_Masked,
  KiliNP_Classic_RaoQ_Masked,
  KiliNP_TWDTW_RaoQ_Masked,
  KiliNP_Mean_SG,
  KiliNP_ShannonH_SG,
  KiliNP_Classic_RaoQ_SG,
  KiliNP_TWDTW_RaoQ_SG
)

names(KiliNP_Comparison_Rasters) <- c(
  "Sentinel-2_MeanNDVI_Masked",
  "ShannonH_Masked",
  "RaosQ_Classic_Masked",
  "RaosQ_TWDTW_Masked",
  "Sentinel-2_MeanNDVI_SG_Gap-Filled",
  "ShannonH_SG",
  "RaosQ_Classic_SG",
  "RaosQ_TWDTW_SG"
)

writeRaster( # So I don't have to compute it every time
  KiliNP_Comparison_Rasters,
  filename = file.path(KiliNP_Results, "KiliNP_NDVI_Diversity_Comparison.tif"),
  overwrite = TRUE
)

KiliNP_Comparison_Rasters <- rast(file.path(KiliNP_Results, "KiliNP_NDVI_Diversity_Comparison.tif")) # Load it back in

png(file.path(KiliNP_Results, "KiliNP_NDVI-Based_Indices_Comparison.png"), # Exported for the paper
    width = 2560, height = 1440, res = 150)

plot(KiliNP_Comparison_Rasters) # Plot it in a file for export

dev.off()

plot(KiliNP_Comparison_Rasters) # Plot it to see how it looks

### Assess index performance using vegetation ground truth ####

message("Assessing diversity indices against vegetation ground truth...")

# Load KiliNP_LandCover_Vector boundary (this is our land cover ground truth data)

KiliNP_LandCover_Vector <- 
  vect(file.path(KiliNP_Input,
                 "/Kili Ground Truthing Land Cover Classifications/VegAug1_KILI_SES_withnewcof.shp")) # Load in the ground truth data again just in case I haven't earlier

# Ensure CRS matches

if (crs(KiliNP_Comparison_Rasters) != crs(KiliNP_LandCover_Vector)){
  KiliNP_LandCover_Vector <- project(KiliNP_LandCover_Vector, crs(KiliNP_Comparison_Rasters))
}

# Crop and mask the ENTIRE comparison stack simultaneously to match the ground truth extent

masked.KiliNP_Comparison_Rasters <- mask(
  crop(KiliNP_Comparison_Rasters, KiliNP_LandCover_Vector), 
  KiliNP_LandCover_Vector
)

## Rasterise vegetation class
# First I need to update the ground truth vector to use the proper category names
# I'll use a lookup table

kili.land.cover.lookup <- c(
  "0"  = "Snow/glacier",
  "1"  = "Agriculture (MAI)",
  "2"  = "Savannah (SAV)",
  "3"  = "Swamp",
  "4"  = "Overgrown clearing",
  "7"  = "Forest plantation",
  "9"  = "Riverine",
  "10" = "Upper montane Erica excelsa forest (FPO Podocarpus disturbed)",
  "11" = "Subalpine Erica trimera bushland (FED incl FER (Erica forest and bushland))",
  "12" = "Podocarpus forest (FPO)",
  "13" = "Subalpine tussock grassland",
  "14" = "Chagga homegardens (HOM)",
  "15" = "Alpine Helichrysum vegetation (HEL)",
  "16" = "Ocotea forest (FOC)",
  "17" = "Bare rock",
  "18" = "Sub/lower montane rainforest (FLM)",
  "19" = "Coffee plantations (COF)"
)

# Assign readable names to the grid code vector of the ground truth

KiliNP_LandCover_Vector$grid_code <- kili.land.cover.lookup[as.character(KiliNP_LandCover_Vector$grid_code)]

# Rasterise the land cover vector

KiliNP_LandCover_Raster <- rasterize(
  KiliNP_LandCover_Vector,
  KiliNP_Comparison_Rasters,
  field = "grid_code"
)

### Convert the index rasters to a dataframe for performance analysis ####
# Bind the cropped comparison stack with the newly rasterised ground truth

KiliNP_Indices_Comparison_Raster <- c(
  masked.KiliNP_Comparison_Rasters,
  KiliNP_LandCover_Raster
)

# Extract the existing 7 names and append the ground truth name

names(KiliNP_Indices_Comparison_Raster) <- c(
  names(masked.KiliNP_Comparison_Rasters),
  "Veg_GroundTruth"
)

KiliNP_Indices_Comparison_DF <- as.data.frame(
  KiliNP_Indices_Comparison_Raster,
  na.rm = TRUE
)

### PERMANOVA ####
## These datasets are too large to conduct a PERMANOVA upon (37121.9GB RAM required)
## Instead, I will use a random representative subset of the data
# 1. Strip out the NA and 0 value background pixels (cloud masks, dead space, boundaries)

KiliNP_NDVI_Filtered <- KiliNP_Indices_Comparison_DF %>%
  filter(!is.na(Veg_GroundTruth)) %>% # Keep rows where the ground truth is valid
  filter(if_all(-Veg_GroundTruth, ~ . != 0 & !is.na(.))) # Ensure none of the index columns (everything except Veg_GroundTruth) equal 0 or NA

# 2. Dynamically identify all valid vegetation classes on the mountain

target_classes <- unique(KiliNP_NDVI_Filtered$Veg_GroundTruth)

# 3. Calculate the equal split for the 10,000 pixel target

total_target_pixels <- 10000
n_classes <- length(target_classes)
pixels_per_class <- floor(total_target_pixels / n_classes) 

message("Sampling up to ", pixels_per_class, " pixels across ", n_classes, " vegetation classes.")

# 4. Execute the Stratified Subsample

subset.KiliNP_Indices_Comparison_DF <- KiliNP_NDVI_Filtered %>%
  filter(Veg_GroundTruth %in% target_classes) %>%
  group_by(Veg_GroundTruth) %>%
  slice_sample(n = pixels_per_class, replace = FALSE) %>%  # slice_sample safely takes all available pixels if a rare class has fewer than pixels_per_class
  ungroup() %>%
  as.data.frame() # Cast back to base R data.frame for the PERMANOVA

# 5. Save for later

saveRDS(subset.KiliNP_Indices_Comparison_DF,
        file = file.path(KiliNP_Processed, "Kili_NDVI_Comparison_DF_Filtered.rds"))

subset.KiliNP_Indices_Comparison_DF <- readRDS(file.path(KiliNP_Processed, "Kili_NDVI_Comparison_DF_Filtered.rds"))

### NOTE FOR REVIEWERS: ###
# The following PERMANOVAs can be run below on a powerful machine, 
# or by script 03.4_Kilimanjaro_NDVI_PERMANOVAs, which uses the exact same code but can be run externally on a server.

## Conduct a series of PERMANOVAs
# For the masked raster

PERMANOVA_ShannonsH_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_Classic_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_TWDTW_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

# For the gap-filled raster

PERMANOVA_ShannonsH_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_Classic_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_TWDTW_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

# Put the PERMANOVA results into a dataframe for effective presentation

KiliNP_PERMANOVA_Results <- data.frame(
  Pipeline = c("Masked", "Masked", "Masked", "SG Filtered", "SG Filtered", "SG Filtered"),
  Index = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q", "Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2 = c(
    PERMANOVA_ShannonsH_Masked$R2[1],
    PERMANOVA_RaosQ_Classic_Masked$R2[1],
    PERMANOVA_RaosQ_TWDTW_Masked$R2[1],
    PERMANOVA_ShannonsH_SG$R2[1],
    PERMANOVA_RaosQ_Classic_SG$R2[1],
    PERMANOVA_RaosQ_TWDTW_SG$R2[1]),
  F = c(
    PERMANOVA_ShannonsH_Masked$F[1],
    PERMANOVA_RaosQ_Classic_Masked$F[1],
    PERMANOVA_RaosQ_TWDTW_Masked$F[1],
    PERMANOVA_ShannonsH_SG$F[1],
    PERMANOVA_RaosQ_Classic_SG$F[1],
    PERMANOVA_RaosQ_TWDTW_SG$F[1]),   
  p_value = c(
    PERMANOVA_ShannonsH_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_Masked$`Pr(>F)`[1],
    PERMANOVA_ShannonsH_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_SG$`Pr(>F)`[1])
)

write.csv(KiliNP_PERMANOVA_Results, file = file.path(KiliNP_Results, "KiliNP_NDVI_PERMANOVA_Summary.csv"), row.names = FALSE)
message("Results successfully saved to: ", file.path(KiliNP_Results, "KiliNP_NDVI_PERMANOVA_Summary.csv"))

print(KiliNP_PERMANOVA_Results)

message("Kilimanjaro NDVI analyses complete!")