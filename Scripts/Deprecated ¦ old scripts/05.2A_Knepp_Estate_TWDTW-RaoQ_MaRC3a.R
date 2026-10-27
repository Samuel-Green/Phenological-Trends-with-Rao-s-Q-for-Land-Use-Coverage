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

# -------------------------------
# LOAD TRUE TIME VECTOR
# -------------------------------

# 1. Extract the layer names directly from the loaded raster tile 'r'
current_names <- names(r)

# 2. Strip the "_NDVI" suffix and append "-15" to create a YYYY-MM-DD string
# (This matches the naming convention you set up during your data prep)
dynamic_dates_string <- paste0(sub("_NDVI$", "", current_names), "-15")

# 3. Convert the strings into actual R Date objects
dynamic_time_vector <- as.Date(dynamic_dates_string, format = "%Y-%m-%d")

# Reconstruct time vector
n_layers <- nlyr(r)

# 2. THE ARRAY SEVERANCE TRICK
# Extract the array to break ALTREP bindings
r_array <- terra::as.array(r)

# Flag Sentinel-2 NoData (-9999 or -32768) as NA to prevent mathematical overflow
r_array[r_array <= -9000] <- NA

# Scale the integers back to a normal float range (-1 to 1)
# Multiplying by 1.0 forces strict double-precision numeric memory for the C++ code
r_array <- (r_array / 10000) * 1.0

# Rebuild the raster with the cleaned, scaled data
r <- terra::rast(r_array, crs = terra::crs(r), ext = terra::ext(r))

## Run the TWDTW Rao's Q using the newly generated, perfectly matched time vector
res <- paRao(
  x = r,
  time_vector = dynamic_time_vector,
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
  progBar = FALSE # CRITICAL: Prevents stdout crashes in SLURM/RStudio!
)

rao <- res[[1]][[1]]

rm(res); gc()

cat("Writing output...\n")
flush.console()

writeRaster(rao, out_file, overwrite = TRUE)

cat("Done:", format(Sys.time()), "\n")
flush.console()