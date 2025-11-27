# Time-weighted dynamic time warping of Rao's Q:

## Incorporating Phenology Into Mapping Land-use From Space

### Overview:



This project preprocesses Landsat 7 and Sentinel 2 satellite data surface reflectance data for the three study sites used in this analysis:

* Macchia Sacra (Italy)
* Knepp Estate (UK)
* Kilimanjaro National Park (Tanzania)



The workflow follows the British Ecological Society’s Guide to Reproducible Code and uses a renv environment to ensure package reproducibility.



#### Directory structure:

* Data 🔢/Raw Data			# Unmodified downloaded imagery
* Data 🔢/Processed Data		# Cropped, masked, cleaned rasters created by code
* Scripts/			# All R scripts used in preprocessing
* Results 📈📉/				# Figures, summaries, metadata



Each site has its own folder under Raw Data and Processed Data.



### Workflow:

**00\_setup.R**

* Loads packages, defines paths, ensures reproducibility, and sets up the session.



**01\_Preprocess\_Sites.R**

* Loads Landsat 7 scenes
* Loads Knepp AOI (MultiPolygon shapefile)
* Reprojects AOI if necessary
* Crops + masks all scenes
* Saves output to Data 🔢/Processed Data/Knepp Estate/



### Reproducibility:

All packages are version-controlled via the renv



No manual data manipulation (everything done via scripts)



All output filenames are constructed programmatically



No setwd() used, paths are handled with the here package



### Requirements:

R ≥ 4.2



QGIS installed if performing any optional spatial checks



Sufficient RAM (16 GB recommended)



### Contact:

Maintainer: Elliot Samuel Shayle

University of Marburg

