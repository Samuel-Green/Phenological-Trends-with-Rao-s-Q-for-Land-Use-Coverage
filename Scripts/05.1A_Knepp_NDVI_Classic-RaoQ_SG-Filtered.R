###############################################################################
# 05.1A_Knepp_NDVI_Classic-RaoQ_SG-Filtered.R
# Processes ONE tile for classical Rao's Q on the 2019-2024 SG-filtered NDVI stack
###############################################################################

library(terra)
library(rasterdiv)

# 1. Grab the tile index passed by the SLURM --array scheduler (1 to 750)
tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Starting Knepp SG Filtered NDVI Classic Rao tile:", tile_id, "\n")
cat("Start time:", format(Sys.time()), "\n")
flush.console()

# 2. Define directories (Update base_dir if your Knepp folder on MaRC3a is named differently)
base_dir <- "~/TWDTW_Paper"
tile_dir <- file.path(base_dir, "Data/Processed Data/Knepp Estate/Knepp_Tiles_NDVI_SG-Filtered")
out_dir  <- file.path(tile_dir, "Rao-utputs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Construct input and output tile paths matching Step 2 of 05_Analyse_Knepp_Estate_Sentinel_NDVI.R
tile_file <- file.path(
  tile_dir,
  paste0("Knepp_MeanNDVI_SG_Tile-", tile_id, ".tif")
)

# If makeTiles() skipped this tile because it fell entirely in the NA corner outside the 5km buffer,
# exit cleanly with status = 0 so SLURM records it as COMPLETED instead of FAILED.
if (!file.exists(tile_file)) {
  cat("Tile file missing (Empty Boundary Tile):", tile_file, "- Exiting cleanly.\n")
  quit(save = "no", status = 0)
}

out_file <- file.path(
  out_dir,
  paste0("Knepp_NDVI_Classic-RaoQ_SG_Tile-", tile_id, ".tif")
)

# Checkpoint check: skip tiles that have already been computed
if (file.exists(out_file)) {
  cat("Knepp SG NDVI Tile", tile_id, "already processed. Skipping.\n")
  flush.console()
  quit(save = "no", status = 0)
}

# 4. Load the pre-scaled float tile directly (No /10000 division needed!)
cat("Loading tile...\n")
flush.console()

tile_raster <- rast(tile_file)

# 5. Run Classical Rao's Q
cat("Running classical Rao's Q...\n")
cat("Rao computation start:", format(Sys.time()), "\n")
flush.console()

res <- paRao(
  x = tile_raster,
  window = 3,
  alpha = 2,
  na.tolerance = 0, # Skips NA border pixels cleanly
  simplify = 2,     # Rounds floats to 2 decimal places (consistent with Shannon's H)
  method = "classic",
  np = 1,           # 1 core per tile since SLURM parallelises across 750 tiles
  progBar = FALSE   # Prevents stdout buffer crashes on headless compute nodes
)

cat("Rao computation finished:", format(Sys.time()), "\n")
flush.console()

# Extract the computed SpatRaster layer from the nested list
rao_raster <- res[[1]][[1]]

# Free memory immediately
rm(res, tile_raster)
gc()

# 6. Export completed tile to Rao-utputs
cat("Writing output...\n")
flush.console()

writeRaster(
  rao_raster,
  out_file,
  overwrite = TRUE
)

cat("Knepp SG NDVI Classic Tile", tile_id, "completed.\n")
cat("End time:", format(Sys.time()), "\n")
flush.console()