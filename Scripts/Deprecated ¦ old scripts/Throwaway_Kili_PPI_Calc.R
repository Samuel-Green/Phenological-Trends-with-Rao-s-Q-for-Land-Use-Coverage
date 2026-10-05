# Throwaway_PPI_Calc.R
library(terra)
library(here)
library(parallel)

message("Loading paths and rasters...")
ProcessedData <- here::here("Data/Processed_Data")
KiliNP_Processed <- file.path(ProcessedData, "Kilimanjaro")
KiliNP_PPI_Processed <- "~/Kilimanjaro Geodata/Processed" 

# Load only the Red and NIR bands required for DVI
KiliNP_B04_Files <- list.files(KiliNP_PPI_Processed, pattern = "^KiliNP_\\d{4}_B04_Cropped\\.tif$", full.names = TRUE)
KiliNP_Red_Timeseries <- rast(KiliNP_B04_Files)

KiliNP_B08_Files <- list.files(KiliNP_PPI_Processed, pattern = "^KiliNP_\\d{4}_B08_Cropped\\.tif$", full.names = TRUE)
KiliNP_NIR_Timeseries <- rast(KiliNP_B08_Files)

# Extract and assign dates for the timeseries[cite: 16]
Kili.years <- sub(".*_(\\d{4})_B\\d{2}_Cropped\\.tif", "\\1", basename(sources(KiliNP_Red_Timeseries)))
Kili.layer.names <- unlist(lapply(Kili.years, function(y) paste0(y, " - ", month.name)))
names(KiliNP_Red_Timeseries) <- Kili.layer.names
names(KiliNP_NIR_Timeseries) <- Kili.layer.names

message("Calculating SZA and DVI...")
# Calculate Astronomical Solar Zenith Angles[cite: 16]
kili_ext <- ext(KiliNP_Red_Timeseries)
centroid_geom <- vect(matrix(c(mean(c(kili_ext$xmin, kili_ext$xmax)), mean(c(kili_ext$ymin, kili_ext$ymax))), ncol=2), crs = crs(KiliNP_Red_Timeseries))
centroid_latlon <- project(centroid_geom, "EPSG:4326")
kili_lat_rad <- geom(centroid_latlon)[,"y"] * (pi / 180)

dates <- seq(as.Date("2017-01-15"), as.Date("2021-12-15"), by = "1 month")
day_of_year <- as.numeric(strftime(dates, format = "%j"))
declination_rad <- (23.45 * sin((2 * pi / 365) * (day_of_year - 81))) * (pi / 180)
hour_angle_rad <- (-1.5 * 15) * (pi / 180)

cos_theta_z <- sin(kili_lat_rad) * sin(declination_rad) + cos(kili_lat_rad) * cos(declination_rad) * cos(hour_angle_rad)
sza_rad_vector <- acos(cos_theta_z)

# Calculate DVI[cite: 16]
KiliNP_DVI <- KiliNP_NIR_Timeseries - KiliNP_Red_Timeseries

message("Applying PPI formula...")
# Define NA-tolerant PPI math[cite: 16]
calc_ppi <- function(dvi_vals, sza_vector) {
  if(all(is.na(dvi_vals))) return(rep(NA_real_, length(dvi_vals)))
  M <- max(dvi_vals, na.rm = TRUE) + 0.005
  ppi_out <- rep(NA_real_, length(dvi_vals))
  valid <- !is.na(dvi_vals)
  if(any(valid)) {
    v_dvi <- dvi_vals[valid]
    v_sza <- sza_vector[valid]
    d_c <- 0.0336 + 0.0477 / cos(v_sza)
    Q_E <- d_c + (1 - d_c) * 0.5 / cos(v_sza)
    K <- 1 / (4 * Q_E) * (1 + M) / (1 - M)
    log_arg <- (M - v_dvi) / (M - 0.09)
    log_arg[log_arg <= 0] <- NA_real_
    ppi_out[valid] <- -K * log(log_arg)
  }
  return(ppi_out)
}

# Apply across Z-dimension using all available SLURM cores
kili.cores <- detectCores()
KiliNP_PPI_Timeseries <- app(KiliNP_DVI, fun = calc_ppi, sza_vector = sza_rad_vector, cores = kili.cores)
names(KiliNP_PPI_Timeseries) <- names(KiliNP_DVI)

message("Exporting final raster...")
writeRaster(
  KiliNP_PPI_Timeseries, 
  filename = file.path(KiliNP_Processed, "KiliNP_PPI_2017-2021_Timeseries.tif"), 
  overwrite = TRUE
)
message("PPI Calculation Complete.")