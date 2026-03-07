###############################################################################
# run_kili_twdtw_tile.R
# Processes ONE tile for TWDTW Rao's Q
###############################################################################

library(terra)
library(rasterdiv)
library(twdtw)

# Tile index passed by SLURM

tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Starting tile:", tile_id, "\n")

# directories

tile_dir <- "Data/Processed Data/Kilimanjaro/tiles"
out_dir  <- "Data/Processed Data/Kilimanjaro/tiles/twdtw_outputs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(
  tile_dir,
  paste0("KiliNP_2017-2021_NDVI_Tile-", tile_id, ".tif")
)

out_file <- file.path(
  out_dir,
  paste0("KiliNP_2017-2021_TWDTW-RaoQ_Tile-", tile_id, ".tif")
)

# Skip already processed tiles

if(file.exists(out_file)){
  cat("Tile", tile_id, "already processed. Skipping.\n")
  quit(save="no")
}

cat("Loading tile...\n")

tile_raster <- rast(tile_file)

# Construct time vector

n_layers <- nlyr(tile_raster)

years <- rep(2017:2021, each=12)
months <- rep(1:12, times=5)

dates <- as.Date(paste(years, months, "01", sep="-"))

cat("Running TWDTW Rao's Q...\n")

res <- paRao(
  x = tile_raster,
  time_vector = dates,
  window = 3,
  alpha = 2,
  simplify = 2,
  na.tolerance = 0,
  np = 1,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 6,
  stepness = -0.5,
  cycle_length = "year",
  time_scale = "month"
)

rao_raster <- res[[1]][[1]]

cat("Writing output...\n")

writeRaster(
  rao_raster,
  out_file,
  overwrite = TRUE
)

cat("Tile", tile_id, "completed.\n")