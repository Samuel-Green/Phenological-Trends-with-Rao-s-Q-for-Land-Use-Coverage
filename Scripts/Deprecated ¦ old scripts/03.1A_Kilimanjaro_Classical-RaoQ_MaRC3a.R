###############################################################################
# 03.1A_Kilimanjaro_Classical-RaoQ_MaRC3a
# Processes ONE tile for classical Rao's Q
###############################################################################

library(terra)
library(rasterdiv)

# Tile index passed by SLURM

tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Starting tile:", tile_id, "\n")
cat("Start time:", format(Sys.time()), "\n")   # *** ADDED ***
flush.console()                                # *** ADDED ***

# directories

tile_dir <- "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Data/Processed Data/Kilimanjaro/Mean NDVI tiles"
out_dir  <- "/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Data/Processed Data/Kilimanjaro/Mean NDVI tiles/rao-utputs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(
  tile_dir,
  paste0("KiliNP_MeanNDVI_Tile-", tile_id, ".tif")
)

if(!file.exists(tile_file)){
  stop(paste("Tile file missing:", tile_file))
}

out_file <- file.path(
  out_dir,
  paste0("KiliNP_Classic-RaoQ_Tile-", tile_id, ".tif")
)

# Skip already processed tiles

if(file.exists(out_file)){
  cat("Tile", tile_id, "already processed. Skipping.\n")
  flush.console()                              # *** ADDED ***
  quit(save="no")
}

cat("Loading tile...\n")
flush.console()                                # *** ADDED ***

tile_raster <- rast(tile_file)

cat("Running classical Rao's Q...\n")
cat("Rao computation start:", format(Sys.time()), "\n")  # *** ADDED ***
flush.console()                                          # *** ADDED ***

res <- paRao(
  tile_raster,
  window = 3,
  alpha = 2,
  simplify = 2,
  method = "classic",
  np = 1
)

cat("Rao computation finished:", format(Sys.time()), "\n")  # *** ADDED ***
flush.console()                                             # *** ADDED ***

rao_raster <- res[[1]][[1]]

# *** ADDED: free memory used by large object ***
rm(res)
gc()

cat("Writing output...\n")
flush.console()                                # *** ADDED ***

writeRaster(
  rao_raster,
  out_file,
  overwrite = TRUE
)

cat("Tile", tile_id, "completed.\n")
cat("End time:", format(Sys.time()), "\n")     # *** ADDED ***
flush.console()                                # *** ADDED ***