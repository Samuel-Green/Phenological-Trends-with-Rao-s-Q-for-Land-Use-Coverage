###############################################################################
# 03.4_Kilimanjaro_NDVI_PERMANOVAs
# Executes 6-way PERMANOVA across Masked and SG-Filtered NDVI datasets
###############################################################################

library(vegan)
library(parallel)

message("Initialising Kilimanjaro NDVI PERMANOVA environment...")

# Define Linux-safe paths for the MaRC3a server

KiliNP_Processed <- "~/TWDTW_Paper/Data/Processed_Data/Kilimanjaro"
KiliNP_Results <- "~/TWDTW_Paper/Results/Kilimanjaro"
dir.create(KiliNP_Results, showWarnings = FALSE, recursive = TRUE)

# Dynamically detect available cores assigned by SLURM

kili.cores <- max(1, detectCores())
message("Allocated ", kili.cores, " cores for permutation processing.")

# Load the filtered, stratified dataset

message("Loading stratified dataset...")
subset.KiliNP_Indices_Comparison_DF <- readRDS(file.path(KiliNP_Processed, "Kili_NDVI_Comparison_DF_Filtered.rds"))

message(paste0(Sys.time(), ": Starting Masked raster PERMANOVAs..."))

PERMANOVA_ShannonsH_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_Classic_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_TWDTW_Masked <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_Masked ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

message(paste0(Sys.time(), ": Starting SG-Filtered raster PERMANOVAs..."))

PERMANOVA_ShannonsH_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$ShannonH_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_Classic_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_Classic_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

PERMANOVA_RaosQ_TWDTW_SG <- adonis2(
  subset.KiliNP_Indices_Comparison_DF$RaosQ_TWDTW_SG ~ subset.KiliNP_Indices_Comparison_DF$Veg_GroundTruth,
  permutations = 999, parallel = kili.cores
)

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