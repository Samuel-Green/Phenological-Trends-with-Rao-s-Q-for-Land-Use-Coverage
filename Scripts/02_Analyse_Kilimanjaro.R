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

# Export and reload the rounded raster so I don't have to calculate it every time

writeRaster(KiliNP_Mean_Raster2dec, file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Cropped-Masked-Rounded2DP.tif"), overwrite = TRUE)
KiliNP_Mean_Raster2dec <- rast(file.path(KiliNP_Processed, "KiliNP_MeanNDVI_Cropped-Masked-Rounded2DP.tif"))

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

# Export the raster or reload it if necessary

writeRaster(
  KiliNP_ShannonH_Raster,
  filename = file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_raster.tif"),
  overwrite = TRUE
)

KiliNP_ShannonH_Raster <- rast(file.path(KiliNP_Results, "Kilimanjaro_2017-2021_ShannonH_raster.tif"))

### 2. Classic Rao's Q  ####
## Due to the large size of the raster, I need to tile it so that it can be run
## The tiles will be stitched back together once they're computed

message("Calculating classical Rao's Q for Kilimanjaro...")

### Step 1: Create a grid to define zones for tiling

# Optional but recommended: trim outer NA borders

trimmed.KiliNP_Mean_Raster <- trim(KiliNP_Mean_Raster)

# We want approximately 72 tiles (not too many, not too few)

kili.total.tiles <- 72
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

# Window size (I think 3 is the default, but this can be changed as necessary)

RaoQ.window.size <- 3
kili.tile.overlap <- floor(RaoQ.window.size / 2)

# Create a directory to put the tiles in

kili.tile.dir <- "Data 🔢/Processed Data/Kilimanjaro/tiles"
dir.create(kili.tile.dir, recursive = TRUE, showWarnings = FALSE)

## Finally, create the tiles

kili.tiles <- makeTiles(
  trimmed.KiliNP_Mean_Raster, # Trimmed for easier computation
  y = kili.tiling.grid, # Specifies how the tiles should be allocated
  buffer = kili.tile.overlap, # Adds a little buffer so Rao's Q can compute without edge NAs
  filename = file.path(kili.tile.dir, "KiliNP_MeanNDVI_Tile-.tif"),
  overwrite = TRUE
)

### Step 2: Compute classic Rao's Q for each tile
## Firstly, I need to setup the environment for parallelisation
# Create a subfolder to store the classic Rao's Q output tiles

kili.rao.dir  <- file.path(kili.tile.dir, "rao-utputs") 
dir.create(kili.rao.dir, recursive = TRUE, showWarnings = FALSE)

## Create a computing cluster to parallelise the calculation at the tile level
# Set the number of cores to be used by the cluster

kili.cores <- max(1, detectCores() - 2)

# Initialise a log file so I can actually see what's going on

kili.log.file <- file.path(kili.rao.dir, "KiliNP_RaoQ_processing_log.txt")

# If the log file doesn't exist already, create one

if(!file.exists(kili.log.file)) file.create(kili.log.file)

# Create the cluster (alliterative and punny names are mandatory)

kili.cluster <- makeCluster(kili.cores)

clusterEvalQ(kili.cluster, {
  library(terra)
  library(rasterdiv)
})

clusterExport(kili.cluster, c(
  "kili.tiles",
  "kili.rao.dir",
  "RaoQ.window.size",
  "kili.log.file"
))

# Identify tiles still needing processing (so resources aren't wasted processing tiles already done)

tile.outputs <- file.path(
  kili.rao.dir,
  paste0("KiliNP_Classic-RaoQ_Tile-", seq_along(kili.tiles), ".tif")
)

tiles.to.process <- which(!file.exists(tile.outputs))

cat(length(tiles.to.process), "tiles remaining.\n")

## Now actually run the code

kili.classic.rao.results <- parLapply( # Function call
  kili.cluster,
  tiles.to.process,
  function(i){
    
    library(terra)
    library(rasterdiv)
    
    log_file <- kili.log.file
    
    log_msg <- function(msg){
      cat(
        paste0(Sys.time(), " | Worker ", Sys.getpid(), " | ", msg, "\n"),
        file = log_file,
        append = TRUE
      )
    }
    
    out.file <- file.path(
      kili.rao.dir,
      paste0("KiliNP_Classic-RaoQ_Tile-", i, ".tif")
    )
    
    if(file.exists(out.file)){
      log_msg(paste("Tile", i, "already exists — skipped"))
      return(NULL)
    }
    
    log_msg(paste("Tile", i, "STARTED"))
    
    tmp.tile <- rast(kili.tiles[i]) # Load in the raster for processing
    
    tmp.result <- paRao(
      tmp.tile,
      window = RaoQ.window.size,
      alpha = 2,
      simplify = 2, # This is necessary to maintain consistency with the Shannon's H test (keeps just 2 decimal places)
      method = "classic", # Because this is not looking at timeseries Rao's Q, just regular unidimensional Rao's Q
      np = 1 # Explicitly prevents nested parallelisation (or set above 1 if you want to melt your CPU)
    )
    
    tmp.rao_raster <- tmp.result[[1]][[1]] # Subsetting avoids hardcoding "$window.3$alpha.2"
    
    writeRaster(
      tmp.rao_raster,
      filename = out.file,
      overwrite = TRUE
    )
    
    rm(tmp.tile,tmp.result,tmp.rao_raster)
    gc()
    
    log_msg(paste("Tile №", i, "'s classic Rao's Q calculated successfully."))
    
    return(NULL) # So that each worker doesn't fill up R's memory with bloat upon completion
  }
)

### Step 3: Demosaic the classical Rao's Q tiles
## Gather up all the files

kili.rao.files <- list.files(
  kili.rao.dir,
  pattern = "Classic-RaoQ",
  full.names = TRUE
)

# Tell R to apply the `rast` function to them 

kili.rao.tiles <- lapply(kili.rao.files, rast)

# Convert them into a spatial raster collection:

kili.rao.tiles <- sprc(kili.rao.tiles)

# Run the demosaic function

KiliNP_Classic_RaoQ <- terra::mosaic(kili.rao.tiles)

# Export the final raster

writeRaster(
  KiliNP_Classic_RaoQ,
  "Results 📈📉/Kilimanjaro/Kilimanjaro_Classic-RaoQ.tif",
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