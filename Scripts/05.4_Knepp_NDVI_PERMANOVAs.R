############################################################################ ###
# Elliot Samuel Shayle - University of Marburg - 28/10/2026                    #
# 05.4_Knepp_NDVI_PERMANOVAs.R                                                 #
# Standalone script to run 3-way PERMANOVA on 2019-2024 Knepp NDVI indices     #
############################################################################ ###

library(vegan)
library(parallel)
library(here) # Ensures dynamic, drive-agnostic pathing

message(paste0(Sys.time(), ": Initialising Knepp Estate NDVI PERMANOVA environment..."))

### 1. Drive-Agnostic Path Definitions ####
## Automatically detects whether the current machine uses spaces/emojis (Keller)
## or underscores (MaRC3a) relative to the .Rproj root directory.

# Resolve ProcessedData directory

if (dir.exists(here::here("Data", "Processed Data"))) {
  ProcessedData <- here::here("Data", "Processed Data")
} else {
  ProcessedData <- here::here("Data", "Processed_Data")
}

# Resolve Knepp_Processed directory

if (dir.exists(file.path(ProcessedData, "Knepp Estate"))) {
  Knepp_Processed <- file.path(ProcessedData, "Knepp Estate")
} else {
  Knepp_Processed <- file.path(ProcessedData, "Knepp_Estate")
}

# Resolve Results directory

if (dir.exists(here::here("Results \U1F4C8\U1F4C9"))) {
  Results <- here::here("Results \U1F4C8\U1F4C9")
} else {
  Results <- here::here("Results")
}

# Resolve Knepp_Results directory

if (dir.exists(file.path(Results, "Knepp Estate")) || !dir.exists(file.path(Results, "Knepp_Estate"))) {
  Knepp_Results <- file.path(Results, "Knepp Estate")
} else {
  Knepp_Results <- file.path(Results, "Knepp_Estate")
}
dir.create(Knepp_Results, showWarnings = FALSE, recursive = TRUE)

# Allocate cores

permanova.cores <- max(1, detectCores() /2.5)
message("Allocated ", permanova.cores, " CPU cores for adonis2 permutation processing.")

### 2. Load the Filtered, Stratified Dataframe ####

subset.Knepp_Indices_Comparison_DF <- readRDS(
  file.path(Knepp_Processed, "Knepp_NDVI_Comparison_DF_Filtered.rds")
)

### 3. Conduct a series of PERMANOVAs ####
## Create the PERMANOVAs subfolder inside the Knepp results directory

Knepp_PERMANOVA_Dir <- file.path(Knepp_Results, "PERMANOVAs")
dir.create(Knepp_PERMANOVA_Dir, showWarnings = FALSE, recursive = TRUE)

message(paste0(Sys.time(), ": Starting SG-Filtered Knepp NDVI raster PERMANOVAs..."))

# 1. SG-Filtered Shannon's H (NDVI) 

message(paste0(Sys.time(), ": [1/3] Running PERMANOVA_Knepp_ShannonsH_SG..."))
PERMANOVA_Knepp_ShannonsH_SG <- adonis2(
  subset.Knepp_Indices_Comparison_DF$`ShannonH SG (NDVI)` ~ subset.Knepp_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = permanova.cores
)
saveRDS(PERMANOVA_Knepp_ShannonsH_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_NDVI_ShannonsH_SG.rds"))
message(paste0(Sys.time(), ": [1/3] Completed and saved PERMANOVA_Knepp_NDVI_ShannonsH_SG.rds"))
gc()

# 2. SG-Filtered Classic Rao's Q (NDVI) 

message(paste0(Sys.time(), ": [2/3] Running PERMANOVA_Knepp_RaosQ_Classic_SG..."))
PERMANOVA_Knepp_RaosQ_Classic_SG <- adonis2(
  subset.Knepp_Indices_Comparison_DF$`RaosQ Classic SG (NDVI)` ~ subset.Knepp_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = permanova.cores
)
saveRDS(PERMANOVA_Knepp_RaosQ_Classic_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_NDVI_RaosQ_Classic_SG.rds"))
message(paste0(Sys.time(), ": [2/3] Completed and saved PERMANOVA_Knepp_NDVI_RaosQ_Classic_SG.rds"))
gc()

# 3. SG-Filtered TWDTW Rao's Q (NDVI)

message(paste0(Sys.time(), ": [3/3] Running PERMANOVA_Knepp_RaosQ_TWDTW_SG..."))
PERMANOVA_Knepp_RaosQ_TWDTW_SG <- adonis2(
  subset.Knepp_Indices_Comparison_DF$`RaosQ TWDTW SG (NDVI)` ~ subset.Knepp_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999,
  parallel     = permanova.cores
)
saveRDS(PERMANOVA_Knepp_RaosQ_TWDTW_SG, file.path(Knepp_PERMANOVA_Dir, "PERMANOVA_Knepp_NDVI_RaosQ_TWDTW_SG.rds"))
message(paste0(Sys.time(), ": [3/3] Completed and saved PERMANOVA_Knepp_NDVI_RaosQ_TWDTW_SG.rds"))
gc()

### 4. Compile and export final PERMANOVA Summary Table ####

message(paste0(Sys.time(), ": Compiling Knepp NDVI PERMANOVA results..."))

Knepp_NDVI_PERMANOVA_Results <- data.frame(
  Pipeline = c("SG Filtered", "SG Filtered", "SG Filtered"),
  Index    = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2       = c(
    PERMANOVA_Knepp_ShannonsH_SG$R2[1],
    PERMANOVA_Knepp_RaosQ_Classic_SG$R2[1],
    PERMANOVA_Knepp_RaosQ_TWDTW_SG$R2[1]
  ),
  F        = c(
    PERMANOVA_Knepp_ShannonsH_SG$F[1],
    PERMANOVA_Knepp_RaosQ_Classic_SG$F[1],
    PERMANOVA_Knepp_RaosQ_TWDTW_SG$F[1]
  ),
  p_value  = c(
    PERMANOVA_Knepp_ShannonsH_SG$`Pr(>F)`[1],
    PERMANOVA_Knepp_RaosQ_Classic_SG$`Pr(>F)`[1],
    PERMANOVA_Knepp_RaosQ_TWDTW_SG$`Pr(>F)`[1]
  )
)

write.csv(
  Knepp_NDVI_PERMANOVA_Results,
  file      = file.path(Knepp_Results, "Knepp_NDVI_PERMANOVA_Summary.csv"),
  row.names = FALSE
)
saveRDS(
  Knepp_NDVI_PERMANOVA_Results,
  file = file.path(Knepp_Results, "Knepp_PERMANOVA_Summary_NDVI.rds")
)

print(Knepp_NDVI_PERMANOVA_Results)
message("Knepp Estate 2019-2024 NDVI analysis complete!") # End of script ####