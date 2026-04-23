###############################################################################
# 05.2A_Knepp_Estate_TWDTW-RaoQ_MaRC3a
# Processes ONE tile for TWDTW Rao's Q
###############################################################################

library(terra)
library(rasterdiv)
library(twdtw)

tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
year    <- Sys.getenv("YEAR")

cat("Year:", year, "| Tile:", tile_id, "\n")
cat("Start:", format(Sys.time()), "\n")
flush.console()

# Directories
tile_dir <- file.path("/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Data/Processed Data/Knepp Estate", paste0("Knepp_", year, "_NDVI_Timeseries_Tiles"))
out_dir  <- file.path(tile_dir, "Rao-utputs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(tile_dir, paste0("Knepp_", year, "_TS_Tile-", tile_id, ".tif"))

if (!file.exists(tile_file)) {
  stop(paste("Missing tile:", tile_file))
}

out_file <- file.path(out_dir, paste0("Knepp_", year, "_TWDTW_Rao_Tile-", tile_id, ".tif"))

if (file.exists(out_file)) {
  cat("Already done. Skipping.\n")
  quit(save = "no")
}

cat("Loading tile...\n")
flush.console()

r <- rast(tile_file)

# Reconstruct time vector
n_layers <- nlyr(r)

# 2. THE ARRAY SEVERANCE TRICK
# This breaks the ALTREP bindings and forces strict double-precision numeric memory
r <- terra::rast(terra::as.array(r * 1.0), crs = terra::crs(r), ext = terra::ext(r))

# -------------------------------
# LOAD TRUE TIME VECTOR
# -------------------------------

# 1. Protect against trailing blank lines causing NAs in C++
time_lines <- readLines(file.path(tile_dir, "time_vector.txt"))
time_lines <- time_lines[time_lines != ""] 
time_vector <- as.Date(time_lines)

cat("Running TWDTW Rao...\n")
flush.console()

# 2. Force the raster completely into RAM to prevent lazy-evaluation segfaults
r <- r + 0 

res <- paRao(
  x = r,
  time_vector = time_vector,
  window = 3,
  alpha = 2,
  na.tolerance = 0, # CRITICAL: Stops empty background pixels from crashing the C++ code
  simplify = 2,
  method = "multidimension",
  dist_m = "twdtw",
  midpoint = 6,
  stepness = -0.5,
  cycle_length = "year",
  time_scale = "month",
  np = 1,
  progBar = FALSE # 3. CRITICAL: Prevents stdout crashes in SLURM/RStudio!
)

rao <- res[[1]][[1]]

rm(res); gc()

cat("Writing output...\n")
flush.console()

writeRaster(rao, out_file, overwrite = TRUE)

cat("Done:", format(Sys.time()), "\n")
flush.console()