###############################################################################
# 03.4_Kilimanjaro_NDVI_PERMANOVAs
# Executes 6-way PERMANOVA across Masked and SG-Filtered NDVI datasets
###############################################################################

library(vegan)
library(parallel)
library(here)

message("Initialising Kilimanjaro NDVI PERMANOVA environment...")

# Define Linux-safe paths for the Keller Computer

InputData <- here::here("Data/Input Data")
ProcessedData <- here::here("Data/Processed Data")
Results <- here::here("Results 📈📉")

KiliNP_Input <- file.path(InputData, "Kilimanjaro")
KiliNP_Processed <- file.path(ProcessedData, "Kilimanjaro")
KiliNP_Results <- file.path(Results, "Kilimanjaro")

dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

# Dynamically detect available cores assigned by SLURM

kili.cores <- max(1, detectCores())
message("Allocated ", kili.cores, " cores for permutation processing.")

# Load the filtered, stratified dataset

message("Loading stratified dataset...")
subset.KiliNP_Indices_Comparison_DF <- readRDS(file.path(KiliNP_Processed, "Kili_NDVI_Comparison_DF_Filtered.rds"))

# Create the PERMANOVAs subfolder inside the Kilimanjaro results directory
KiliNP_PERMANOVA_Dir <- file.path(KiliNP_Results, "PERMANOVAs")
dir.create(KiliNP_PERMANOVA_Dir, showWarnings = FALSE, recursive = TRUE)

message(paste0(Sys.time(), ": Starting Masked NDVI raster PERMANOVAs..."))

# --- 1. Masked Shannon's H (NDVI) ---
message(paste0(Sys.time(), ": [1/6] Running PERMANOVA_ShannonsH_Masked..."))
PERMANOVA_ShannonsH_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`ShannonH Masked (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_ShannonsH_Masked, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_ShannonsH_Masked.rds"))
message(paste0(Sys.time(), ": [1/6] Completed and saved PERMANOVA_NDVI_ShannonsH_Masked.rds"))
gc()

# --- 2. Masked Classic Rao's Q (NDVI) ---
message(paste0(Sys.time(), ": [2/6] Running PERMANOVA_RaosQ_Classic_Masked..."))
PERMANOVA_RaosQ_Classic_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`RaosQ Classic Masked (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_RaosQ_Classic_Masked, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_RaosQ_Classic_Masked.rds"))
message(paste0(Sys.time(), ": [2/6] Completed and saved PERMANOVA_NDVI_RaosQ_Classic_Masked.rds"))
gc()

# --- 3. Masked TWDTW Rao's Q (NDVI) ---
message(paste0(Sys.time(), ": [3/6] Running PERMANOVA_RaosQ_TWDTW_Masked..."))
PERMANOVA_RaosQ_TWDTW_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`RaosQ TWDTW Masked (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_RaosQ_TWDTW_Masked, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_RaosQ_TWDTW_Masked.rds"))
message(paste0(Sys.time(), ": [3/6] Completed and saved PERMANOVA_NDVI_RaosQ_TWDTW_Masked.rds"))
gc()

message(paste0(Sys.time(), ": Starting SG-Filtered NDVI raster PERMANOVAs..."))

# --- 4. SG-Filtered Shannon's H (NDVI) ---
message(paste0(Sys.time(), ": [4/6] Running PERMANOVA_ShannonsH_SG..."))
PERMANOVA_ShannonsH_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`ShannonH SG (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_ShannonsH_SG, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_ShannonsH_SG.rds"))
message(paste0(Sys.time(), ": [4/6] Completed and saved PERMANOVA_NDVI_ShannonsH_SG.rds"))
gc()

# --- 5. SG-Filtered Classic Rao's Q (NDVI) ---
message(paste0(Sys.time(), ": [5/6] Running PERMANOVA_RaosQ_Classic_SG..."))
PERMANOVA_RaosQ_Classic_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`RaosQ Classic SG (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_RaosQ_Classic_SG, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_RaosQ_Classic_SG.rds"))
message(paste0(Sys.time(), ": [5/6] Completed and saved PERMANOVA_NDVI_RaosQ_Classic_SG.rds"))
gc()

# --- 6. SG-Filtered TWDTW Rao's Q (NDVI) ---
message(paste0(Sys.time(), ": [6/6] Running PERMANOVA_RaosQ_TWDTW_SG..."))
PERMANOVA_RaosQ_TWDTW_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$`RaosQ TWDTW SG (NDVI)` ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)
saveRDS(PERMANOVA_RaosQ_TWDTW_SG, file.path(KiliNP_PERMANOVA_Dir, "PERMANOVA_NDVI_RaosQ_TWDTW_SG.rds"))
message(paste0(Sys.time(), ": [6/6] Completed and saved PERMANOVA_NDVI_RaosQ_TWDTW_SG.rds"))
gc()

message(paste0(Sys.time(), ": Compiling results..."))

# Put the PERMANOVA results into a dataframe for effective presentation

KiliNP_PERMANOVA_Results <- data.frame(
  Pipeline = c("Masked", "Masked", "Masked", "SG Filtered", "SG Filtered", "SG Filtered"),
  Index = c("Shannon H", "Classic Rao Q", "TWDTW Rao Q", "Shannon H", "Classic Rao Q", "TWDTW Rao Q"),
  R2 = c(
    PERMANOVA_ShannonsH_Masked$R2[1],
    PERMANOVA_RaosQ_Classic_Masked$R2[1],
    PERMANOVA_RaosQ_TWDTW_Masked$R2[1],
    PERMANOVA_ShannonsH_SG$R2[1],
    PERMANOVA_RaosQ_Classic_SG$R2[1],
    PERMANOVA_RaosQ_TWDTW_SG$R2[1]),
  F = c(
    PERMANOVA_ShannonsH_Masked$F[1],
    PERMANOVA_RaosQ_Classic_Masked$F[1],
    PERMANOVA_RaosQ_TWDTW_Masked$F[1],
    PERMANOVA_ShannonsH_SG$F[1],
    PERMANOVA_RaosQ_Classic_SG$F[1],
    PERMANOVA_RaosQ_TWDTW_SG$F[1]),
  p_value = c(
    PERMANOVA_ShannonsH_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_Masked$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_Masked$`Pr(>F)`[1],
    PERMANOVA_ShannonsH_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_Classic_SG$`Pr(>F)`[1],
    PERMANOVA_RaosQ_TWDTW_SG$`Pr(>F)`[1])
)

print(KiliNP_PERMANOVA_Results)

# Final export for safekeeping

export_path <- file.path(KiliNP_Results, "KiliNP_NDVI_PERMANOVA_Summary.csv")
write.csv(KiliNP_PERMANOVA_Results, file = export_path, row.names = FALSE)
message("Results successfully saved to: ", export_path)

message(paste0(Sys.time(), ": Kilimanjaro NDVI analyses complete!"))