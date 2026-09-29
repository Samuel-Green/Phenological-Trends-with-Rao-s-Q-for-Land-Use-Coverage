###############################################################################
# 03.1A_Kilimanjaro_Classical-RaoQ_Masked
# Processes ONE tile for classical Rao's Q on the masked dataset
###############################################################################

library(terra)
library(rasterdiv)

# Tile index passed by SLURM

tile_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Starting Masked tile:", tile_id, "\n")
cat("Start time:", format(Sys.time()), "\n")
flush.console()

# Directories (Using Linux-safe paths)

tile_dir <- "~/TWDTW_Paper/Data/Processed_Data/Kilimanjaro/Kili_Tiles_NDVI_Masked"
out_dir  <- file.path(tile_dir, "Rao-utputs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(
  tile_dir,
  paste0("KiliNP_MeanNDVI_Masked_Tile-", tile_id, ".tif")
)

if(!file.exists(tile_file)){
  cat("Tile file missing (Empty Boundary Tile):", tile_file, "- Exiting cleanly.\n")
  quit(save="no", status = 0)
}

out_file <- file.path(
  out_dir,
  paste0("KiliNP_Classic-RaoQ_Masked_Tile-", tile_id, ".tif")
)

# Skip already processed tiles

if(file.exists(out_file)){
  cat("Masked Tile", tile_id, "already processed. Skipping.\n")
  flush.console()
  quit(save="no")
}

cat("Loading tile...\n")
flush.console()

tile_raster <- rast(tile_file)

cat("Running classical Rao's Q...\n")
cat("Rao computation start:", format(Sys.time()), "\n")
flush.console()

res <- paRao(
  tile_raster,
  window = 3,
  alpha = 2,
  simplify = 2,
  method = "classic",
  np = 1
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

cat("Masked Tile", tile_id, "completed.\n")
cat("End time:", format(Sys.time()), "\n")
flush.console()