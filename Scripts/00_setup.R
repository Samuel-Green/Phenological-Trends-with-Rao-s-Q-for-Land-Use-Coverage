############################################################## ###
# Elliot Samuel Shayle - University of Marburg - 24/11/2025   ####
# 00_setup.R - Setup TWDTW Project                            ####
############################################################## ###

### Libraries ####
## Install packages
install.packages("terra")
install.packages("sf")
install.packages("here")
install.packages("dplyr")
install.packages("stringr")
#install.packages("rasterdiv")
install.packages("remotes")
remotes::install_github("mattmar/rasterdiv") # [20/02/2026]: Package author recommends using the pre-release version as it has a bug fix
install.packages("twdtw")
install.packages("vegan")
install.packages("snow")
install.packages("pROC")

## Load libraries
library(terra)
library(sf)
library(here)
library(dplyr)
library(stringr)
library(rasterdiv)
library(twdtw)
library(vegan)
library(snow)
library(pROC)

### Reproducibility settings ####
options(scipen = 999)
set.seed(123)

### Path definitions ####
InputData <- here::here("Data 🔢/Input Data")
ProcessedData <- here::here("Data 🔢/Processed Data")
Results <- here::here("Results 📈📉")

### Site-specific paths ####
Macchia_Input <- file.path(InputData, "Macchia Sacra")
Knepp_Input <- file.path(InputData, "Knepp Estate")
KiliNP_Input <- file.path(InputData, "Kilimanjaro")

Macchia_Processed <- file.path(ProcessedData, "Macchia Sacra")
Knepp_Processed <- file.path(ProcessedData, "Knepp Estate")
KiliNP_Processed <- file.path(ProcessedData, "Kilimanjaro")

message("Setup complete. Paths and environment ready.")
