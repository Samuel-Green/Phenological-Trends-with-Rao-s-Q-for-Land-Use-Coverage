############################################################################ ###
### Elliot Samuel Shayle - University of Marburg - 24/11/2026                  #
### 06.3_Knepp_PPI_α_Grid-search.R                                             #
### TWDTW α (Steepness) Optimization for Knepp PPI via Consensus UKCEH (No Bin)#
############################################################################ ###

library(rasterdiv)
library(twdtw)
library(vegan)
library(terra)
library(dplyr)
library(here)     # Ensures dynamic, drive-agnostic pathing
library(parallel)

### 1. Drive-Agnostic Path Definitions ####
## Automatically resolves Keller (spaces/emojis) vs. MaRC3a (underscores) paths

if (dir.exists(here::here("Data", "Input Data"))) {
  InputData <- here::here("Data", "Input Data")
} else {
  InputData <- here::here("Data", "Input_Data")
}

if (dir.exists(here::here("Data", "Processed Data"))) {
  ProcessedData <- here::here("Data", "Processed Data")
} else {
  ProcessedData <- here::here("Data", "Processed_Data")
}

if (dir.exists(here::here("Results \U1F4C8\U1F4C9"))) {
  Results <- here::here("Results \U1F4C8\U1F4C9")
} else {
  Results <- here::here("Results")
}

# Resolve Knepp-specific subdirectories
Knepp_Input     <- if (dir.exists(file.path(InputData, "Knepp Estate"))) file.path(InputData, "Knepp Estate") else file.path(InputData, "Knepp_Estate")
Knepp_Processed <- if (dir.exists(file.path(ProcessedData, "Knepp Estate"))) file.path(ProcessedData, "Knepp Estate") else file.path(ProcessedData, "Knepp_Estate")
Knepp_Results   <- if (dir.exists(file.path(Results, "Knepp Estate")) || !dir.exists(file.path(Results, "Knepp_Estate"))) file.path(Results, "Knepp Estate") else file.path(Results, "Knepp_Estate")
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Allocate safe parallel core counts (respects SLURM allocation if running on MaRC3a)
slurm_cpus      <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
available_cores <- if (!is.na(slurm_cpus)) slurm_cpus else max(1, detectCores() - 2)
parao.cores     <- min(12, available_cores)
permanova.cores <- min(8, available_cores) # Hard cap at 8 cores to prevent adonis2 OOM crashes

message("Using ", parao.cores, " cores for paRao and an ", permanova.cores, "-core PSOCK cluster for PERMANOVA.")

### 2. Import Spatial Vectors & Cleaned 2019-2024 PPI Timeseries ####
message("Importing spatial vectors and 2019-2024 SG-filtered PPI timeseries...")

knepp.tiling.grid               <- vect(file.path(Knepp_Processed, "Knepp_Tiling_Grid_Polygons.geojson"))
KneppEstate_Boundaries          <- vect(file.path(Knepp_Input, "Schulte-to-Bühne_Knepp_Site_Data", "knepp_boundary_lcc_pretty.shp"))
KneppEstate_Buffered_Boundaries <- vect(file.path(Knepp_Processed, "Knepp_Estate_Buffered_Boundary.shp"))

Knepp_PPI_Timeseries_Clean <- rast(file.path(Knepp_Processed, "Knepp_Buffered_PPI_Cleaned_SG-method.tif"))

# Ensure all vectors share the exact CRS of the PPI raster stack
if (crs(knepp.tiling.grid) != crs(Knepp_PPI_Timeseries_Clean)) {
  knepp.tiling.grid <- project(knepp.tiling.grid, crs(Knepp_PPI_Timeseries_Clean))
}
if (crs(KneppEstate_Boundaries) != crs(Knepp_PPI_Timeseries_Clean)) {
  KneppEstate_Boundaries <- project(KneppEstate_Boundaries, crs(Knepp_PPI_Timeseries_Clean))
}
if (crs(KneppEstate_Buffered_Boundaries) != crs(Knepp_PPI_Timeseries_Clean)) {
  KneppEstate_Buffered_Boundaries <- project(KneppEstate_Buffered_Boundaries, crs(Knepp_PPI_Timeseries_Clean))
}

# Dynamically extract the 72-month Date vector from the layer names (e.g., "2019-01_PPI")
current.knepp.names <- names(Knepp_PPI_Timeseries_Clean)
knepp.dates.string  <- paste0(sub("_PPI$", "", current.knepp.names), "-15")
Knepp.dates         <- as.Date(knepp.dates.string, format = "%Y-%m-%d")

if (any(is.na(Knepp.dates)) || length(Knepp.dates) != nlyr(Knepp_PPI_Timeseries_Clean)) {
  message("Layer names stripped on disk; reconstructing 2019-2024 monthly calendar...")
  Knepp.dates <- seq(as.Date("2019-01-15"), by = "month", length.out = nlyr(Knepp_PPI_Timeseries_Clean))
}

message("Loaded ", length(Knepp.dates), " monthly PPI layers spanning ", min(Knepp.dates), " to ", max(Knepp.dates), ".")

### 3. Stratified Random Tile Subsampling (10 Inside Estate, 15 Outside Buffer) ####
message("Selecting stratified random subsample of 25 valid tiles (identical seed to NDVI grid search)...")

knepp.tiling.grid$Tile_ID <- seq_len(nrow(knepp.tiling.grid))

# Keep only tiles inside the 5 km buffered boundary (excludes empty NA corners)
valid_tiles <- knepp.tiling.grid[KneppEstate_Buffered_Boundaries, ]

# Separate core estate tiles from surrounding 5 km buffer tiles
inner_tiles <- valid_tiles[KneppEstate_Boundaries, ]
outer_tiles <- valid_tiles[!valid_tiles$Tile_ID %in% inner_tiles$Tile_ID, ]

set.seed(42)
n_inner <- min(10, nrow(inner_tiles))
n_outer <- min(15, nrow(outer_tiles))

sampled_inner <- inner_tiles[sample(seq_len(nrow(inner_tiles)), n_inner), ]
sampled_outer <- outer_tiles[sample(seq_len(nrow(outer_tiles)), n_outer), ]

sample_tiles <- rbind(sampled_inner, sampled_outer)
sample_poly  <- aggregate(sample_tiles)

message("Cropping, masking, and trimming PPI timeseries to the 25 sampled tiles...")
sample_ts <- crop(Knepp_PPI_Timeseries_Clean, sample_poly)
sample_ts <- mask(sample_ts, sample_poly)
sample_ts <- trim(sample_ts)

### 4. The Array Severance Trick (Full 10m Resolution - NO BINNING) ####
message("Severing lazy C++ pointers and mapping strict numeric memory (No Binning)...")

sample_array <- terra::as.array(sample_ts)
sample_array <- sample_array * 1.0 # Blocks ALTREP logical-NA corruption

sample_ts_RAM <- terra::rast(
  sample_array,
  crs = crs(sample_ts),
  ext = ext(sample_ts)
)
names(sample_ts_RAM) <- names(sample_ts)

# Bake numeric subset raster to disk for safe parallel worker transmission
tmp.tif_path <- file.path(Knepp_Processed, "Knepp_PPI_Grid_Search_Sample_FullRes.tif")
writeRaster(sample_ts_RAM, tmp.tif_path, datatype = "FLT8S", overwrite = TRUE)

rm(sample_array, sample_ts_RAM)
gc()

sample_ts   <- terra::rast(tmp.tif_path)
sample_mean <- app(sample_ts, fun = mean, na.rm = TRUE)

### 5. Pre-Loop Strict Consensus UKCEH Land Cover Preparation (2019-2024) ####

UKCEH_LC_Classes <- c(
  "1"  = "Broadleaved woodland",
  "2"  = "Coniferous woodland",
  "3"  = "Arable and horticulture",
  "4"  = "Improved grassland",
  "5"  = "Neutral grassland",
  "6"  = "Calcareous grassland",
  "7"  = "Acid grassland",
  "8"  = "Fen, marsh and swamp",
  "9"  = "Heather",
  "10" = "Heather grassland",
  "11" = "Bog",
  "12" = "Inland rock",
  "13" = "Saltwater",
  "14" = "Freshwater",
  "15" = "Supralittoral rock",
  "16" = "Supralittoral sediment",
  "17" = "Littoral rock",
  "18" = "Littoral sediment",
  "19" = "Saltmarsh",
  "20" = "Urban",
  "21" = "Suburban"
)

# All 7 UKCEH classes actually present in the Knepp consensus dataset (including Freshwater)
target_classes <- c(
  "Broadleaved woodland",
  "Coniferous woodland",
  "Arable and horticulture",
  "Improved grassland",
  "Freshwater",
  "Urban",
  "Suburban"
)

consensus_lc_path <- file.path(Knepp_Processed, "Knepp_2019-2024_UKCEH_Consensus_LandCover.tif")

if (file.exists(consensus_lc_path)) {
  message("Loading pre-computed 2019-2024 UKCEH Consensus Land Cover raster...")
  lc_full <- rast(consensus_lc_path)
  Knepp_LandCover_Consensus <- project(lc_full, sample_mean, method = "near")
  Knepp_LandCover_Consensus <- mask(Knepp_LandCover_Consensus, sample_mean)
  names(Knepp_LandCover_Consensus) <- "Veg_GroundTruth"
} else {
  message("Pre-computed consensus raster not found; building from raw 2019-2024 UKCEH files...")
  
  ceh_base_dir <- file.path(Knepp_Input, "CEH Land-cover data")
  if (!dir.exists(ceh_base_dir)) ceh_base_dir <- file.path(Knepp_Input, "CEH_Land-cover_data")
  
  ceh_tif_files <- list.files(ceh_base_dir, pattern = "^LCM.*\\.tif$", recursive = TRUE, full.names = TRUE)
  ceh_years     <- sub(".*(20\\d{2}).*", "\\1", ceh_tif_files)
  valid_ceh_idx <- which(ceh_years %in% as.character(2019:2024))
  
  ceh_tif_files <- ceh_tif_files[valid_ceh_idx]
  ceh_years     <- ceh_years[valid_ceh_idx]
  names(ceh_tif_files) <- ceh_years
  
  if (length(ceh_tif_files) == 0) stop("No UKCEH LCM .tif files found for 2019-2024 in: ", ceh_base_dir)
  
  Knepp_LC_Aligned_List <- lapply(ceh_years, function(y) {
    lc_raw     <- rast(ceh_tif_files[[y]])[[1]]
    lc_aligned <- project(lc_raw, sample_mean, method = "near")
    mask(lc_aligned, sample_mean)
  })
  
  Knepp_LC_Stack <- rast(Knepp_LC_Aligned_List)
  consensus_mask <- (min(Knepp_LC_Stack, na.rm = FALSE) == max(Knepp_LC_Stack, na.rm = FALSE))
  
  Knepp_LandCover_Consensus <- mask(Knepp_LC_Stack[[1]], consensus_mask, maskvalues = 0)
  names(Knepp_LandCover_Consensus) <- "Veg_GroundTruth"
}

### 6. Helper Function: Run paRao & Stratified PERMANOVA for a Single Alpha ####
run_knepp_ppi_alpha_step <- function(a_val, sample_ts, dates_vec, lc_consensus, ukceh_lookup, target_classes, n_cores_rao, perm_cluster) {
  
  message(paste0(Sys.time(), " | Testing PPI \U03B1 = ", a_val, " (midpoint = 2)..."))
  
  tmp.RaoQ <- paRao(
    x            = sample_ts,
    time_vector  = dates_vec,
    window       = 3,
    alpha        = 2,
    na.tolerance = 0,
    simplify     = 2,
    np           = n_cores_rao,
    progBar      = FALSE,
    method       = "multidimension",
    dist_m       = "twdtw",
    midpoint     = 2,     # 2-month temperate UK seasonal inflection point
    stepness     = a_val,
    cycle_length = "year",
    time_scale   = "month"
  )
  
  message(paste0(Sys.time(), " | paRao completed for \U03B1 = ", a_val))
  
  tmp.raster       <- tmp.RaoQ$window.3$alpha.2
  tmp.stack        <- c(tmp.raster, lc_consensus)
  names(tmp.stack) <- c("RaosQ", "Veg_GroundTruth")
  
  tmp.df <- as.data.frame(tmp.stack, na.rm = TRUE)
  
  # Map raw integer UKCEH codes (1:21) to descriptive strings in base R
  tmp.df$Veg_GroundTruth <- ukceh_lookup[as.character(tmp.df$Veg_GroundTruth)]
  
  # Filter out 0/NA artifacts and subset to the 7 target classes
  tmp.df_filtered <- tmp.df %>%
    filter(!is.na(Veg_GroundTruth) & Veg_GroundTruth %in% target_classes) %>%
    filter(!is.na(RaosQ) & RaosQ != 0)
  
  present_classes  <- unique(tmp.df_filtered$Veg_GroundTruth)
  pixels_per_class <- floor(10000 / length(present_classes))
  
  set.seed(42)
  tmp.df_sub <- tmp.df_filtered %>%
    group_by(Veg_GroundTruth) %>%
    slice_sample(n = pixels_per_class, replace = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  
  message(paste0(Sys.time(), " | Running PERMANOVA on ", nrow(tmp.df_sub),
                 " stratified consensus pixels across ", length(present_classes), " classes..."))
  
  tmp.permanova <- adonis2(
    tmp.df_sub$RaosQ ~ tmp.df_sub$Veg_GroundTruth,
    permutations = 999,
    parallel     = perm_cluster
  )
  
  step_result <- data.frame(
    Alpha   = a_val,
    R2      = tmp.permanova$R2[1],
    F_stat  = tmp.permanova$F[1],
    p_value = tmp.permanova$`Pr(>F)`[1]
  )
  
  rm(tmp.RaoQ, tmp.raster, tmp.stack, tmp.df, tmp.df_filtered, tmp.df_sub, tmp.permanova)
  gc()
  
  return(step_result)
}

### 7. Two-Stage Grid Search for Optimal TWDTW α ####
message("Starting two-stage PPI grid search for optimal TWDTW \U03B1 (midpoint = 2)...")

log_csv <- file.path(Knepp_Results, "Knepp_PPI_Alpha_GridSearch_NoBin_Log.csv")

alpha_coarse_grid <- c(-0.1, -0.3, -0.5, -0.7, -0.9)
alpha_results     <- data.frame(Alpha = numeric(), R2 = numeric(), F_stat = numeric(), p_value = numeric())

# Initialise one reusable PSOCK cluster for all adonis2 calls (avoids mclapply SIGPIPE crashes)
perm_cluster <- makeCluster(permanova.cores)

## STAGE 1: Coarse Grid
for (a in alpha_coarse_grid) {
  step_result <- run_knepp_ppi_alpha_step(
    a_val          = a,
    sample_ts      = sample_ts,
    dates_vec      = Knepp.dates,
    lc_consensus   = Knepp_LandCover_Consensus,
    ukceh_lookup   = UKCEH_LC_Classes,
    target_classes = target_classes,
    n_cores_rao    = parao.cores,
    perm_cluster   = perm_cluster
  )
  
  alpha_results <- rbind(alpha_results, step_result)
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  message(paste0("Logged \U03B1 = ", a, " | R2 = ", round(step_result$R2, 5)))
}

best.coarse.alpha <- alpha_results$Alpha[which.max(alpha_results$R2)]
message(paste("Best coarse \U03B1 is:", best.coarse.alpha))

## STAGE 2: Fine Grid
alpha_fine_grid <- seq(best.coarse.alpha + 0.1, best.coarse.alpha - 0.1, by = -0.05)
alpha_fine_grid <- round(alpha_fine_grid[alpha_fine_grid < 0], 2)

for (a in alpha_fine_grid) {
  if (a %in% round(alpha_coarse_grid, 2)) next
  
  step_result <- run_knepp_ppi_alpha_step(
    a_val          = a,
    sample_ts      = sample_ts,
    dates_vec      = Knepp.dates,
    lc_consensus   = Knepp_LandCover_Consensus,
    ukceh_lookup   = UKCEH_LC_Classes,
    target_classes = target_classes,
    n_cores_rao    = parao.cores,
    perm_cluster   = perm_cluster
  )
  
  alpha_results <- rbind(alpha_results, step_result)
  write.csv(alpha_results, file = log_csv, row.names = FALSE)
  message(paste0("Logged \U03B1 = ", a, " | R2 = ", round(step_result$R2, 5)))
}

# Shut down the PSOCK cluster cleanly
stopCluster(perm_cluster)
gc()

### 8. Sort and Export Final Results ####
alpha_results <- alpha_results[order(-alpha_results$R2), ]
print("Full Knepp PPI Grid Search Complete. Results:")
print(alpha_results)

Knepp_PPI_Optimal_Alpha <- alpha_results$Alpha[1]
message(paste("The absolute optimal PPI \U03B1 value is:", Knepp_PPI_Optimal_Alpha))

# Save optimal alpha to the exact RDS filename expected by 06.2A_Knepp_PPI_TWDTW-RaoQ_SG-Filtered.R
saveRDS(Knepp_PPI_Optimal_Alpha, file.path(Knepp_Results, "Knepp_PPI_Optimal_Alpha.rds"))
saveRDS(Knepp_PPI_Optimal_Alpha, file.path(Knepp_Results, "Knepp_PPI_Optimal_Alpha_NoBin.rds"))

write.csv(alpha_results, file = log_csv, row.names = FALSE)
message("All Knepp PPI grid search results exported successfully!")