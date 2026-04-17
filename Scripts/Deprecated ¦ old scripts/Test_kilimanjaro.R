#
# library(devtools)
# PD <- "/home/matteo/git/rasterdiv"
# setwd(PD)

#
# load_all(PD)
#
remotes::install_github("mattmar/rasterdiv")
library(rasterdiv)
library(terra)

# kili <- rast("/home/matteo/own_data/PoD/hackathon/Elliot/Data/kili/KM_EPSG32737_SEN2_L3_BOA_NDV_2017_FBM.tif")
kili <- rast("F:/Elliot Shayle/Kilimanjaro Geodata/KM_EPSG32737_SEN2_L3_BOA_NDV_2017_FBM.tif")
kili.100 <- aggregate(kili, fact=100, fun="mean")
kili.100Mean <- mean(kili.100) 

#P: No error; user time elapsed 2.333 secs with 10 [8 for Elliot's laptop] cores
system.time(tmp.resultCP <- paRao(
	kili.100Mean,
	window = 9,
	alpha = 2,
	simplify = 0, 
	method = "classic",
	np = 8
	))

#S: No error; user time elapsed 2.333 secs with 10 cores
system.time(tmp.resultCS <- paRao(
	kili.100Mean,
	window = 9,
	alpha = 2,
	simplify = 0, 
	method = "classic",
	np = 8
	))

par(mfrow=c(3,1))
plot(kili.100Mean)
plot(tmp.resultCP[[1]][[1]])
plot(tmp.resultCS[[1]][[1]])

#S: No error; user time elapsed 395.614 secs with 20 cores
system.time(tmp.resultCP <- paRao(
	kili.100Mean,
	window = 9,
	alpha = 2,
	simplify = 0, 
	method = "classic",
	np = 8
	))

time_vector <- 1:12
#P: No error; user time elapsed 395.614 secs with 20 cores
system.time(tmp.resultTD <- paRao(
	kili.100,
	window = 9,
	alpha = 2,
	simplify = 0, 
	method = "multidimension",
	dist_m="twdtw",
	time_vector = time_vector,
	cycle_length = "year",
	time_scale = "month",
	np = 8,
	progBar=TRUE
	))

#S: No error; user time elapsed 3194.216 secs
system.time(tmp.result <- paRao(
	kili.100,
	window = 9,
	alpha = 2,
	simplify = 0, 
	method = "multidimension",
	dist_m="twdtw",
	time_vector = time_vector,
	cycle_length = "year",
	time_scale = "month",
	np = 1,
	progBar=TRUE
	))

par(mfrow=c(2,2))
plot(tmp.resultCS[[1]][[1]])
plot(tmp.resultCP[[1]][[1]])
plot(tmp.result[[1]][[1]])
plot(tmp.resultTD[[1]][[1]])
