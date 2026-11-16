###############################################################################
# 06.1A_Knepp_Estate_Classical-RaoQ_MaRC3a
# Processes ONE tile for classical Rao's Q
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
tile_dir <- file.path("/home/shayle/TWDTW Paper (B3 Hackathon Postprint)/Data/Processed Data/Knepp Estate", paste0("Knepp_", year, "_MeanPPI_Tiles"))
out_dir  <- file.path(tile_dir, "Rao-utputs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tile_file <- file.path(tile_dir, paste0("Knepp_", year, "_Mean_Tile-", tile_id, ".tif"))

if (!file.exists(tile_file)) {
  stop(paste("Missing tile:", tile_file))
}

out_file <- file.path(out_dir, paste0("Knepp_", year, "_ClassicRao_Tile-", tile_id, ".tif"))

if (file.exists(out_file)) {
  cat("Already done. Skipping.\n")
  quit(save = "no")
}

cat("Loading tile...\n")
flush.console()

r <- rast(tile_file)

cat("Running Classical Rao...\n")
flush.console()

res <- paRao(
  r,
  window = 3,
  alpha = 2,
  simplify = 2,
  method = "classic",
  np = 1
)

rao <- res[[1]][[1]]

rm(res); gc()

cat("Writing output...\n")
flush.console()

writeRaster(rao, out_file, overwrite = TRUE)

cat("Done:", format(Sys.time()), "\n")
flush.console()