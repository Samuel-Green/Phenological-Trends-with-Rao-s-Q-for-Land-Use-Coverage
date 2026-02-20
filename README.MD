# Time-weighted dynamic time warping of Rao's Q:

## Incorporating Phenology Into Mapping Land-use From Space

### Overview:



This project preprocesses Landsat 7 and Sentinel 2 satellite data surface reflectance data for the three study sites used in this analysis:

* Macchia Sacra (Italy)
* Knepp Estate (UK)
* Kilimanjaro National Park (Tanzania)



The workflow follows the British Ecological Society’s Guide to Reproducible Code and uses a renv environment to ensure package reproducibility.



#### Directory structure:

* Data 🔢/Input Data			# These are the data which I've received (either from satellite imagery providers or collaborators) and remain unaltered
* Data 🔢/Processed Data		# Data which were in some way modified, such as cleaned datasets, buffered shapefiles, processed rasters
* Scripts/			# All R scripts used in preprocessing
* Results 📈📉/				# Figures, summaries, metadata, and finalised raster files



Each site has its own folder under Input Data and Processed Data.



### Workflow:

**00\_setup.R**

* Loads packages, defines paths, ensures reproducibility, and sets up the session.



**01\_Analyse\_Macchia\_Sacra.R**

* Loads Macchia Sacra PPI and NDVI data.
* Conducts Shannon-Wiener Index analysis.
* Conducts classical Rao's Q Index analysis.
* Conducts multidimensional  Rao's Q Index analysis.
* Conducts PERMANOVA performance analysis of indices' performance.
* Tests comparative ability of the different approaches to remove artefacts in the satellite data.



### Reproducibility:

All packages are version-controlled via the renv



All input data is processed reproducibly using R scripts



All output filenames are constructed programmatically



No setwd() used, paths are handled with the here package so the project can be easily rerun on any other machine



### Requirements:

R ≥ 4.2



QGIS installed if performing any optional spatial checks



Sufficient RAM (16 GB recommended)



### Contact:

Maintainer: Elliot Samuel Shayle

University of Marburg

Email 📧: [elliot.shayle.18@alumni.ucl.ac.uk](mailto:elliot.shayle.18@alumni.ucl.ac.uk?subject=Rao's%20Q%20and%20TWDTW%20Project%20Inquiry) || [elliot.shayle@staff.uni-marburg.de](mailto:elliot.shayle@staff.uni-marburg.de?subject=Rao's%20Q%20and%20TWDTW%20Project%20Inquiry)

