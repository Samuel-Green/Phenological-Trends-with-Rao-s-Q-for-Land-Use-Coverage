###############################################################################
# 03.2B_Kilimanjaro_TWDTW-RaoQ_SG-Filtered
# Processes ONE tile for TWDTW Rao's Q on the Savitzky-Golay filtered raster
###############################################################################

library(terra)
library(rasterdiv)
library(twdtw)

# Tile index passed by SLURM

tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Starting SG Filtered TWDTW tile:", tile_id, "\n")
cat("Start time:", format(Sys.time()), "\n")
flush.console()

# Directories (Using Linux-safe paths)

tile_dir <- "~/TWDTW_Paper/Data/Processed_Data/Kilimanjaro/Kili_TS_Tiles_NDVI_SG-Filtered"
out_dir  <- file.path(tile_dir, "Rao-utputs")
results_dir <- "~/TWDTW_Paper/Results/Kilimanjaro"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(
  tile_dir,
  paste0("KiliNP_TS_NDVI_SG_Tile-", tile_id, ".tif")
)

if(!file.exists(tile_file)){
  cat("Tile file missing (Empty Boundary Tile):", tile_file, "- Exiting cleanly.\n")
  quit(save="no", status = 0)
}

out_file <- file.path(
  out_dir,
  paste0("KiliNP_TWDTW-RaoQ_SG_Tile-", tile_id, ".tif")
)

# Skip already processed tiles

if(file.exists(out_file)){
  cat("SG Tile", tile_id, "already processed. Skipping.\n")
  flush.console()
  quit(save="no")
}

cat("Loading optimal TWDTW \U03B1 (stepness)...\n")
optimal_alpha_path <- file.path(results_dir, "Kili_NDVI_Optimal_Alpha_NoBin.rds")

if(!file.exists(optimal_alpha_path)) {
  stop("Optimal alpha RDS file not found! Ensure grid search has completed.")
}
opt_alpha <- readRDS(optimal_alpha_path)
cat("Optimal \U03B1 loaded:", opt_alpha, "\n")

cat("Loading tile...\n")
flush.console()
tile_raster <- rast(tile_file)

# Construct time vector

years <- rep(2017:2021, each=12)
months <- rep(1:12, times=5)
dates <- as.Date(paste(years, months, "01", sep="-"))

cat("Running TWDTW Rao's Q...\n")
cat("Rao computation start:", format(Sys.time()), "\n")
flush.console()

res <- paRao(
  x = tile_raster,
  time_vector = dates,
  window = 3,
  alpha = 2,
  na.tolerance = 0,
  simplify = 2,
  np = 1,
  progBar = FALSE,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 6,
  stepness = opt_alpha, # Dynamically assigned
  cycle_length = "year",
  time_scale = "month"
)

cat("Rao computation finished:", format(Sys.time()), "\n")
flush.console()

rao_raster <- res[[1]][[1]]

# Free memory used by large objects

rm(res, tile_raster)
gc()

cat("Writing output...\n")
flush.console()

writeRaster(
  rao_raster,
  out_file,
  overwrite = TRUE
)

cat("SG Tile", tile_id, "completed.\n")
cat("End time:", format(Sys.time()), "\n")
flush.console()