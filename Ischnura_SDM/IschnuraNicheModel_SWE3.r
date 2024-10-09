#install.packages(c("maptools", "dismo", "sp", "rgeos", "rgdal", "ncdf", "reshape"))

library(dismo)
library(maptools)
library(sp)
library(rgeos)
library(rgdal)
library(ncdf)
library(reshape)
setwd("~/Desktop/Ischnura niche model")

ig <- read.csv("IschnuraPointsClean.csv")

#make shapefile of Norway, Finland, and Sweden, called "Fenno"
countries<-getData("countries")
SWE <- getData('GADM', country='SWE', level=0) #level 0 means just the country outline
#plot(SWE)
FIN <- getData('GADM', country='Finland', level=0) 
#plot(FIN)
NOR <- getData('GADM', country='Norway', level=0) 
#plot(NOR)
Fenno1 <- gUnion(FIN, NOR)
Fenno <- gUnion(Fenno1, SWE)
plot(Fenno)

#reimport bioclim variables from file:
for (i in 1:19){
	name <- paste("bio", i, sep="")
	#name2 <- paste("bio", i, ".mask", sep="")
	file <- raster(paste("bioclim raster files Europe/bio", i, ".grd", sep=""))
	file2 <- crop(file, Fenno) #crops file to rectangle that fits Fenno
	#file3 <- mask(x=file2, mask=Fenno.r) #crops to exact shape of Fenno
	assign(name, file2)
	#assign(name2, file3)
	}

#make a raster file that represents extent of Fenno
crop <- setValues(bio1, NA)
Fenno.r <- rasterize(Fenno, crop)

#clip bio data to mask:
bio1.mask <- mask(x=bio1, mask=Fenno.r) #crop bioclim map to shape of Fenno
bio2.mask <- mask(x=bio2, mask=Fenno.r)
bio3.mask <- mask(x=bio3, mask=Fenno.r)
bio4.mask <- mask(x=bio4, mask=Fenno.r)
bio5.mask <- mask(x=bio5, mask=Fenno.r)
bio6.mask <- mask(x=bio6, mask=Fenno.r)
bio7.mask <- mask(x=bio7, mask=Fenno.r)
bio8.mask <- mask(x=bio8, mask=Fenno.r)
bio9.mask <- mask(x=bio9, mask=Fenno.r)
bio10.mask <- mask(x=bio10, mask=Fenno.r)
bio11.mask <- mask(x=bio11, mask=Fenno.r)
bio12.mask <- mask(x=bio12, mask=Fenno.r)
bio13.mask <- mask(x=bio13, mask=Fenno.r)
bio14.mask <- mask(x=bio14, mask=Fenno.r)
bio15.mask <- mask(x=bio15, mask=Fenno.r)
bio16.mask <- mask(x=bio16, mask=Fenno.r)
bio17.mask <- mask(x=bio17, mask=Fenno.r)
bio18.mask <- mask(x=bio18, mask=Fenno.r)
bio19.mask <- mask(x=bio19, mask=Fenno.r)

####workspace saved!!!!

# occurrence points:
ischpts <- ig[c(10,9)]
ischpts$species <- c("Ischnura elegans")
ischpts <- ischpts[,c(3, 1:2)]
dups <- duplicated(ischpts[,1:2])
ischpts2 <- ischpts[!dups,]
colnames(ischpts2) <- c("species", "lon", "lat")
ipts <- SpatialPointsDataFrame(data=ischpts2, coords = ischpts2[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

#data(wrld_simpl)
# plot(wrld_simpl, xlim=c(5, 35), ylim=c(54,72), axes=TRUE)

#select occurrence points that are in Fenno only:
ov <- overlay(ipts, Fenno)
i <- which(!is.na(ov))
ipts.fen <- ipts[i,]
dim(ipts.fen) 
#[1] 1767    3

plot(wrld_simpl, xlim=c(-15, 45), ylim=c(30,75), axes=TRUE)
points(ipts.fen, col="red", cex=0.75)
points(bg.fen, col="green", cex=0.75)

#select target background points (= Zygoptera), Targeted background sampling using sightings of a wider group of similar organisms, to control for sampling effort (Phillips et al. 2009)
cal.rec <- gbif("calopteryx", "*", ext=Fenno)
#epa.rec <- gbif("epallage", "*", ext=Fenno) # no records
les.rec <- gbif("lestes", "*", ext=Fenno)
sym.rec <- gbif("sympecma", "*", ext=Fenno)
#cer.rec <- gbif("cercion", "*", ext=Fenno) # no records
coe.rec <- gbif("coenagrion", "*", ext=Fenno)
ery.rec <- gbif("erythromma", "*", ext=Fenno)
#nah.rec <- gbif("nahalennia", "*", ext=Fenno) #no records
pyr.rec <- gbif("pyrrhosoma", "*", ext=Fenno)
ena.rec <- gbif("enallagma", "*", ext=Fenno)
isch.rec <- gbif("ischnura", "*", ext=Fenno)
#cer.rec <- gbif("ceragrion", "*", ext=Fenno) #no records
pla.rec <- gbif("platycnemis", "*", ext=Fenno)

cal <- cal.rec
cal <- cal[,c(8,7)]
cal$species <- c("Calopteryx")
cal <- cal[,c(3,1:2)]
cal <- cal[!duplicated(cal),]
#cal <- SpatialPointsDataFrame(data=cal, coords = cal[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

les <- les.rec
les <- les[,c(8,7)]
les$species <- c("Lestes")
les <- les[,c(3,1:2)]
les <- les[!duplicated(les),]
#les <- SpatialPointsDataFrame(data=les, coords = les[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

#bg <- rbind(cal, les)

# #workspace saved

sym <- sym.rec
sym <- sym[,c(8,7)]
sym$species <- c("Sympecma")
sym <- sym[,c(3,1:2)]
sym <- sym[!duplicated(sym),]
head(sym)

coe <- coe.rec
coe <- coe[,c(8,7)]
coe$species <- c("Coenagrion")
coe <- coe[,c(3,1:2)]
coe <- coe[!duplicated(coe),]
head(coe)

ery <- ery.rec
ery <- ery[,c(8,7)]
ery$species <- c("Erythromma")
ery <- ery[,c(3,1:2)]
ery <- ery[!duplicated(ery),]
head(ery)

pyr <- pyr.rec
pyr <- pyr[,c(8,7)]
pyr$species <- c("Pyrrhosoma")
pyr <- pyr[,c(3,1:2)]
pyr <- pyr[!duplicated(pyr),]
head(pyr)

ena <- ena.rec
ena <- ena[,c(8,7)]
ena$species <- c("Enallagma")
ena <- ena[,c(3,1:2)]
ena <- ena[!duplicated(ena),]
head(ena)

isc <- isch.rec
isc <- isc[,c(8,7)]
isc$species <- c("Ischnura")
isc <- isc[,c(3,1:2)]
isc <- isc[!duplicated(isc),]
head(isc)


pla <- pla.rec
pla <- pla[,c(8,7)]
pla$species <- c("Platycnemis")
pla <- pla[,c(3,1:2)]
pla <- pla[!duplicated(pla),]
head(pla)

bg <- rbind(cal, les, sym, coe, ery, pyr, ena, isc, pla)
bg <- SpatialPointsDataFrame(data=bg, coords = bg[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

ov2 <- overlay(bg, Fenno)
i <- which(!is.na(ov2))
bg.fen <- bg[i,]
dim(bg.fen)
#[1] 23799     3


#bg.target <- rbind(cal.les, coordinates(sym), coordinates(coe), coordinates(ery), coordinates(pyr), coordinates(ena), coordinates (isc), coordinates(pla))
#dim(bg.target)
# # [1] 42883     2 (42883 background points)



################################
# data layers (besides bioclim)

# distance from coast. NASA, 0.01 degree resolution (~30 arcsec) GeoTIFF (interpolated from 0.04 degree resolution data set). <0 is land, >0 is ocean.
# http://oceancolor.gsfc.nasa.gov/DOCS/DistFromCoast/
dcoast <- raster("GMT_intermediate_coast_distance_01d/GMT_intermediate_coast_distance_01d.tif")
dcoast.crop <- crop(dcoast, Fenno)
crop2 <- setValues(dcoast.crop, NA)
Fenno.r2 <- rasterize(Fenno, crop2)
dcoast.mask <- mask(dcoast.crop, Fenno.r2)
dcoast.mask2 <- resample(dcoast.mask, bio1.mask)

predictors <- stack(c(elev.mask, bio1.mask))
predictors <- stack(c(predictors, dcoast.mask2))

#Harmonized World Soil Database available water storage capacity (mm/m of soil unit) 30 arcsec, lat-long.
H20cap <- raster("~/Desktop/Spatial data from NCEAS computer/spatial datasets/HWSD_/HWSD_AWC_CLASS.tif")
H2Ocap.crop<- crop(H20cap, Fenno)
crop3 <- setValues(H2Ocap.crop, NA)
Fenno.r3 <- rasterize(Fenno, crop3)
H2Ocap.mask <- mask(H2Ocap.crop, Fenno.r3)
##clearly NOT in 30 arcsec resolution!
#HWSD classes 1-7 need to be converted to mm/m values.
fun <- function(x) {x[x==0] <- NA;
	x[x==1] <- 150;
	x[x==2] <- 125;
	x[x==3] <- 100;
	x[x==4] <- 75;
	x[x==5] <- 50;
	x[x==6] <- 15;
	x[x==7] <- 0;
	return(x)} 
H2OcapVal.mask <- calc(H2Ocap.mask, fun)
H2OcapVal.mask2 <- projectRaster(H2OcapVal.mask, bio1.mask )

predictors <- stack(c(predictors, H2OcapVal.mask2))

#workspace saved
#GMTED2010 (30 arcsecond, lat-long), elev in meters(?)
elev <- raster("mn30_grd/mn30_grd")
elev.crop <- crop(elev, Fenno)
elev.crop2 <- projectRaster(elev.crop, bio1.mask)
#crop4 <- setValues(elev.crop, NA)
#projection(elev.crop)
#Fenno.r4 <- rasterize(Fenno, crop4) 
#projection(Fenno.r4) <- projection(elev.crop) #this fixed it!
elev.mask <- mask(elev.crop2, Fenno.r) 

#globcover 2009: 150m accuracy. Plate-carree projection. categorical values (11-220).
globcover <- raster("Globcover2009_V2/GLOBCOVER_L4_200901_200912_V2.3.tif")
##globcover.longlat <- projectRaster(globcover, crop4) #crop and reproject #already in longlat!! 
globcover2 <- resample(globcover, bio1)
#globcov.crop <- crop(globcover2, Fenno)
#crop5 <- setValues(globcover2, NA)
#Fenno.r5 <- rasterize(Fenno, crop5)
globcov.mask <- mask(globcover2, Fenno.r)
writeRaster(globcov.mask, filename="globcovMask.grd")

##workspace saved

#WWF terrestrial ecoregions (long-lat, 0.017 degree(?) resolution)
#http://worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
ogrInfo("wwf_terr_bioreg", "wwf_terr_ecos") #directory location, file name (w/out extension) #this command gives info about shapefile
ecoreg <- readOGR("wwf_terr_bioreg", "wwf_terr_ecos")
ecoreg.r <- rasterize(ecoreg, Fenno.r)
#crop6 <- setValues(ecoreg.r, NA)
#Fenno.r6 <- rasterize(Fenno, crop6)
ecoreg.mask <- mask(ecoreg.r, Fenno.r)

#NDVI vegetation index average for July. Averaged over 1982-2000. Values rescaled from (-1 -> +1) to (1 -> 255).
#lat-long proj., 0.1 degree resolution
#compiled by European distributed institute of technology
#processed at Clark Labs (www.clarklabs.org)
#http://edit.csic.es/Soil-Vegetation-LandCover.html
meanNDVI <- raster("mean_esriascii/JULAV18.asc")
meanNDVI.crop <- crop(meanNDVI, Fenno)
crop7 <- setValues(meanNDVI.crop, NA)
Fenno.r7 <- rasterize(Fenno, crop7)
meanNDVI.mask <- mask(meanNDVI.crop, Fenno.r7)
#meanNDVI.mask2 <- meanNDVI.mask
#projection(meanNDVI.mask2) <- projection(bio1.mask)
meanNDVI.mask2 <- resample(meanNDVI.mask, bio1.mask)

#NDVI vegetation index standard deviation for July, over 18 year period (as above)
sdNDVI <- raster("stand_dev_esriascii/JULSD18.asc")
sdNDVI.crop <- crop(sdNDVI, Fenno)
sdNDVI.mask <- mask(sdNDVI.crop, Fenno.r7)
sdNDVI.mask2 <- resample(sdNDVI.mask, bio1.mask)

predictors <- stack(c(predictors, globcov.mask, ecoreg.mask, meanNDVI.mask2, sdNDVI.mask2))

#Net primary productivity: WHAT ARE THE UNITS?
#Imhoff, M.L., L. Bounoua, T. Ricketts, C. Loucks, R. Harriss, and W.T. Lawrence. 2004. HANPP Collection: Global Patterns in Net Primary Productivity (NPP). Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). http://sedac.ciesin.columbia.edu/data/set/hanpp-net-primary-productivity. Accessed DAY MONTH YEAR. (April 3, 2013)
#http://sedac.ciesin.columbia.edu/data/set/hanpp-net-primary-productivity
#0.25 degree resolution, latlong proj
npp <- raster("npp_geotiff/npp_geotiff.tif")
npp.crop <- crop(npp, Fenno)
crop8 <- setValues(npp.crop, NA)
Fenno.r8 <- rasterize(Fenno, crop8)
npp.mask <- mask(npp.crop, Fenno.r8)
npp.mask2 <- resample(npp.mask, bio1.mask)

#European soils database raster 1kmx1km --CATEGORICAL--
#http://eusoils.jrc.ec.europa.eu/ESDB_Archive/ESDB_data_1k_raster_intro/ESDB_1k_raster_data_intro.html
# ETRS89 Lambert Azimuthal Equal Area (ETRS_LAEA)  projection
###dominant parent material (European soils) (PAR-MAT-DOM) 
#Many categories coded between 1000-9300. Category legends stored as text file in Niche Model folder.
parent <- raster("~/Desktop/Spatial data from NCEAS computer/spatial datasets/PARMADO_directory/parmado")
#http://www.nceas.ucsb.edu/scicomp/recipes/projections
parent2 <- parent
projection(parent2) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
# parent.geo <- projectRaster(parent2, crop5)
# Fenno.r5 <- rasterize(Fenno, crop5)
# parent.mask <- mask(parent.geo, Fenno.r5) #cannot open connection Fenno.r5!?!?!
# parent.geo1 <- parent.geo
# projection(parent.geo1) <- projection(Fenno.r)
# extent(parent.geo1) <-extent(Fenno.r)
# resolution(parent.geo1) <- resolution(Fenno.r)
# parent.mask <- mask(parent.geo1, Fenno.r) #nope didn't work
# parent.geo2 <- projectRaster(parent.geo, bio1.mask) #fuck, now can't open parent.geo
parent.geo3 <- projectRaster(parent2, bio1.mask)
parent.mask <- mask(parent.geo3, Fenno.r)
#need to convert 0 values to NA:
fun <- function(x) {x[x==0] <- NA;
	return(x)} 
parent2.mask <- calc(parent.mask, fun)

#AVHRR continuous fields tree cover
#http://glcf.umd.edu/data/treecover/index.shtml
#values: 10-80 = percent tree cover. 254= non-vegetated. 255= tree cover less than 10%. Need to convert 254 and 255 to 0.
treecov <- raster("ea-latlong-treecover.grd")
projection(treecov) <- projection(bio1.mask)
treecov.crop <- crop(treecov, Fenno)
fun <- function(x) {x[x>200] <- 0; return(x)} #set values >200 to 0.
treecov.crop.clean <- calc(treecov.crop, fun)
crop9 <- setValues(treecov.crop.clean, NA)
Fenno.r9 <- rasterize(Fenno, crop9)
treecov.mask <- mask(treecov.crop.clean, Fenno.r9)

##Added later:
#topsoil water content at field capacity; EFSA_THETA_FC_TOP.asc; http://eusoils.jrc.ec.europa.eu/library/data/efsa/download/index.cfm
#projection = ETRS 89 LAEA, resolution = 1km. units = m^3*m^-3
topsoilH2O <- raster("EFSA_V1_1_20121130/EFSA_THETA_FC_TOP.asc")
projection(topsoilH2O) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
topsoilH2O <- projectRaster(topsoilH2O, bio1.mask)
topsoilH2O.mask <- mask(topsoilH2O, Fenno.r)
#need to convert -9000 values to NA:
fun <- function(x) {x[x==-9000] <- NA;
	return(x)} 
topsoilH2O.mask2 <- calc(topsoilH2O.mask, fun)


predictors <- stack(c(predictors, npp.mask2, parent2.mask, treecov.mask, bio2.mask, bio3.mask, bio4.mask, bio5.mask, bio6.mask, bio7.mask, bio8.mask, bio9.mask, bio10.mask, bio11.mask, bio12.mask, bio13.mask, bio14.mask, bio15.mask, bio16.mask, bio17.mask, bio18.mask, bio19.mask))

names(predictors) <- c("elev", "bio1", "dcoast", "H2Ocap", "globcov", "ecoreg", "meanJulNDVI", "sdJulNDVI", "npp", "soilParent", "treecov", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")



###############################################
#############################################


# extract raster layers for both occurrence and background points:
ischvals<-extract(predictors, ipts.fen)
absencevals <- extract(predictors, bg.fen)
pa <- c(rep(1, nrow(ischvals)), rep(0, nrow(absencevals)))
padata <- data.frame(cbind(pa, rbind(ischvals, absencevals)))
padata <- na.omit(padata)

#look for correlations
c <- cor(padata[, c(2:30)], padata[, c(2:30)], method = "pearson", use = "complete.obs")
m <- melt(c) # create dataframe from correlation output
m <- m[order(- abs(m$value)), ] # sort by descending absolute correlation
m

#at 80% cutoff for correlations across background and isch observations, the following bioclim variables were selected: 1,2,3,4,5,8,10,12,15.
predictors <- stack(c(elev.mask, dcoast.mask2, H2OcapVal.mask2, globcov.mask, ecoreg.mask, meanNDVI.mask2, sdNDVI.mask2, npp.mask2, parent2.mask, treecov.mask, bio1.mask, bio2.mask, bio3.mask, bio4.mask, bio5.mask, bio8.mask, bio10.mask, bio12.mask, bio15.mask))

names(predictors) <- c("elev", "dcoast", "H2Ocap", "globcov", "ecoreg", "meanJulNDVI", "sdJulNDVI", "npp", "soilParent", "treecov", "bio1", "bio2", "bio3", "bio4", "bio5", "bio8", "bio10", "bio12", "bio15")

################
############# run maxent

#move maxent.jar to dismo/java folder:
#Type: system.file("java", package="dismo")  to see path to this folder.
#Then go to terminal and type: sudo mv drag-and-drop-maxent.jar-file paste-path-to-dismo-java-folder

#check that it is there:
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep="")
file.exists(jar)
#[1] TRUE

##allocate 1GB of memory to Java before loading dismo:
#uncheck dismo in package manager.
options(java.parameters = "-Xmx1g")
library(dismo)

#withhold 20% sample for testing:
fold <- kfold(ipts.fen, k=5)
#see group assignments:
#fold[1:10]
#unique(fold)

#single partition for preliminary tests:
itest <- ipts.fen[fold == 1,]
itrain <- ipts.fen[fold != 1,]

#maxent args
#http://code.google.com/p/api-maxent-programming/source/browse/trunk/src/edu/berkeley/mvz/amp/MaxEntRun.java?spec=svn15&r=15
# "-J" = jackknife, "-P" = turn on response curves, "-K" = turn on picture making. (-K was 'not understood')
# see help menu in Maxent for more.

#run maxent with training data:
me <- maxent(predictors, itrain, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P"))
##Warning message:
#In .local(x, p, ...) :
#  240 (19.62%) of the presence points have NA predictor values

#typing 'me' opens up html to view results:
me
#results:
#1. analysis of omission-- omission rate should be similar to the predicted omission rate. fraction of background predicted should decline steeply.
#2. ROC curve with area under the curve (AUC). AUC is 0.931 here.
#3. response curves: how each variable affects maxent prediction. If the curve is high at low values, this means that low values of that variable were particularly influential, whereas larger values of the variable did not have much impact. Marginal response curve holds all other predictors costnat. Maxent resonse curves disregard all other variables.
#4. precent contribution and permuation imporance. Percent contribution depends on the path that the algorithm took to its solution. Permutation importance is independent of the idiosyncracies of the path, and is more reliable.
#5. jackknife estimates of variable importance. If any of the light blue bars are longer than the red bar, this means that the model improves when the variable *isn't* used. 
###Files saved under Maxent training results May6_2013 in the Ischnura Niche Model folder.

#plot percent contribution (same information as table in html results)
plot(me)

#plot marginal response curves (not as good-looking as in html version)
response(me)

#predict to entire dataset: (see ?maxent for additional arguments, including saving as .grd)
r <- predict(me, predictors)

#view:
plot(r)
points(ipts.fen)

#evaluate:
e1 <- evaluate(me, p=itest, a=bg.fen, x=predictors)
e1 # AUC = 0.64. Not too good... TPR= true positive rate, TNR = true negative rate.

#not sure what this does:
threshold(e1)

#some 'args' from a posting on message board (no confirmation that these are 'good'):
#me<-maxent(rasters,training,removeDuplicates=TRUE, args=c("jacknife=TRUE","outputdirectory=MaxEntOut","replicates=5", "warnings=FALSE","askoverwrite=FALSE","writeplotdata=TRUE","outputgrids=FALSE","plots=TRUE","appendtoresultsfile=TRUE","applythresholdrule=Fixed cumulative value 5"))
 
# #use model to predict to study region
# r <- predict(me, rasters, progress='window', args=c("doclamp=TRUE", "writeclampgrid=TRUE", "writemess=TRUE")

#try testing data with isch removed from bg.fen:
bg2 <- rbind(cal, les, sym, coe, ery, pyr, ena, pla)
bg2 <- SpatialPointsDataFrame(data=bg2, coords = bg2[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
ov3 <- overlay(bg2, Fenno) #use 'over' instead
i2 <- which(!is.na(ov3))
bg.fen2 <- bg[i2,]
dim(bg.fen2)
#[1] 21891     3

e2 <- evaluate(me, p=itest, a=bg.fen2, x=predictors)
e2 #AUC = 0.64 still. didn't make a difference.

#try a smaller number of background pts:
bg.fen3 <- bg.fen[sample(nrow(bg.fen), 2000),]
e3 <- evaluate(me, p=itest, a=bg.fen3, x=predictors)
e3 # AUC = 0.63. Didn't matter.

#try including background data in maxent model:
me2 <- maxent(predictors, p=itrain, a=bg.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P"))
e4 <- evaluate(me2, p=itest, a=bg.fen, x=predictors)
e4 # AUC = 0.57
r2 <- predict(me2, predictors)
plot(r2) #looks horrible! definitely don't include bg pts in maxent()

#try evaluating me with random background points:
bg.rand <- randomPoints(predictors, 1000)
e5 <- evaluate(me, p=itest, a=bg.rand, x=predictors)
#AUC = 0.95. OK! Need an alternate background point system. *****

#try fewer predictors:
predictors2 <- stack(c(elev.mask, dcoast.mask2, globcov.mask, ecoreg.mask, treecov.mask, bio1.mask, bio2.mask, bio3.mask, bio4.mask, bio5.mask, bio8.mask, bio10.mask, bio12.mask, bio15.mask))
names(predictors2) <- c("elev", "dcoast", "globcov", "ecoreg", "treecov", "bio1", "bio2", "bio3", "bio4", "bio5", "bio8", "bio10", "bio12", "bio15")
me3 <- maxent(predictors2, itrain, factors = c("globcov", "ecoreg"), args=c("-J", "-P")) #AUC ~ 91 or so. Looks fine.
e6 <- evaluate(me3, p=itest, a=bg.rand, x=predictors)
# AUC = 0.95, using the random background pts.

###########graph range
install.packages("alphahull")
library(alphahull)

ipts.fenNS <- as.data.frame(ipts.fen)
ahull <- ahull(ipts.fenNS$lon, ipts.fenNS$lat, a=1)
ahull2 <- ahull(ig$lon, ig$lat, a=1)

plot(wrld_simpl, xlim=c(-5, 35), ylim=c(35, 75), axes=TRUE)
plot(ahull2, add=TRUE, col="blue", lty=4, wpoints=FALSE)

pdf("IschnuraRange.pdf")
plot(countries, xlim=c(10, 20), ylim=c(55,66), axes=TRUE)
plot(ahull2, add=TRUE, col="blue", lty=4, wpoints=FALSE)
points(x=15.466666, y=59.73333, col="red", pch=19)
points(x=13.1931, y=55.7028, col="red", pch=19)
dev.off()

######################

me <- maxent(predictors, ipts.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))
#Random background points from throughout Fenno (not sure how many?). AUC ~ 0.93.

bg.rand2 <- randomPoints(predictors, 10000)
me2 <- maxent(predictors, p=ipts.fen, a=bg.rand2, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))
#similar to above. Fever bg points (~1000) resulted in lower AUC.

#test for spatial bias in sampling:
sb <- ssb(itest, bg.rand, itrain)
sb[,1]/sb[,2]
#0.01992249
#ssb should be close to 1 if there is no bias. This value is very close to 0, indicating high spatial bias. (= bad).

#sample within a .7 degree radius of occurrence points
x <- circles(ipts.fenNS[,2:3], d=0.7, lonlat=FALSE) #couldon't get lonlat=TRUE to work.

#fuse circles into single spatial polygon:
pol <- gUnaryUnion(x@polygons)

#crop circles polygon to coastline of Fenno:
#proj4string(pol) <- CRS(proj4string(Fenno))
# #proj4string(pol) <- CRS(" +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#pol2 <- spTransform(pol, CRS(" +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
Fenno.temp <- Fenno
Fenno.temp2 <- spTransform(Fenno.temp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
pol2 <- gIntersection(pol, Fenno.temp2, byid=TRUE) #doesn't work!!!
## finally got the proj4 strings to match by reprojecting Fenno without the space in front of the first plus. It did not work to try to ADD the space to pol's string!

plot(pol2, axes=TRUE)

##Trying to extract one point per circle, but couldn't figure out how to extract individual circles:
# y@polygons <- SpatialPolygonsDataFrame(x@polygons)
# projection(x@polygons) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# bg.circ <- data.frame()
# for i in 1:10 {
	# point <- spsample(x@polygons[i, extract=TRUE], 1, type="random", iter=10)
	# bg.circ <- rbind(bg.circ, point)
# }

#sample 10000 background points from polygon describing 0.7 degree radius around observations:
samp1 <- spsample(pol2, 10000, type="random", iter=25)
cells <- cellFromXY(Fenno.r, samp1)
length(cells)
cells <- unique(cells)
bg.circ <- xyFromCell(Fenno.r, cells)

plot(bg.circ, axes=TRUE)

#test for spatial bias in sampling:
sb <- ssb(itest, bg.circ, itrain)
sb[,1]/sb[,2]
# 0.1571349 seems ok?

me3 <- maxent(predictors, p=ipts.fen, a=bg.circ, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))
#AUC = 0.869. Pretty good, but not as high as when sampling from throughout Fenno.

#Now try pairwise sampling from within bg.circ:
pwd <- data.frame()
for (i in 1:5) {
itest <- ipts.fen[fold == i,]
itrain <- ipts.fen[fold != i,]
i <- pwdSample(itest, bg.circ, itrain, n=4, tr=0.2)
pwd <- rbind(pwd, i)
}

pwd1 <- as.data.frame(pwd[,1], colnames=FALSE)
1 <- col
pwd2 <- as.data.frame(pwd[,2], colnames=FALSE)
pwd3 <- as.data.frame(pwd[,3], colnames=FALSE)
pwd4 <- as.data.frame(pwd[,4], colnames=FALSE)
pwdn <- rbind(pwd1[1:1767,], pwd2[1:1767,], pwd3[1:1767,], pwd4[1:1767,])

me4 <- maxent(predictors, p=ipts.fen, a=bg.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))
#AUC = 0.778, but predictors don't look good in jackknife
r.me4 <- predict(me4, predictors)
pdf(file="me4 Fenno prediction map.pdf")
plot(r.me4, axes=TRUE)
dev.off()

#which predictor has that big hole in the middle of Sweden??
pdf(file="predictors_maps.pdf")
plot(predictors, axes=TRUE)
dev.off()
# H2Ocap

###################################
#### OK, now start model selection for real. I am choosing to use random background points, since I want to test all of Fenno as potentially suitable habitat-- I am not interested in comparing habitt suitability of I. elegans compared to other damselflies, or in comparion to only the 50km region immediately surrounding their known occurrence.

#workspace saved

me <- maxent(predictors, ipts.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5")). 
#AUC= 0.92, no variable looks like a clear candidate for removal. H2Ocap only explains 2%.

#getting layer for soil moisture content to use instead of H2Ocap (topsoil water content at field capacity; EFSA_THETA_FC_TOP.asc; eurosoils.jrc.ec.europa.eu). In the meantime, just try without H2Ocap:
predictors3 <- stack(c(elev.mask, dcoast.mask2, globcov.mask, ecoreg.mask, meanNDVI.mask2, sdNDVI.mask2, npp.mask2, parent2.mask, treecov.mask, bio1.mask, bio2.mask, bio3.mask, bio4.mask, bio5.mask, bio8.mask, bio10.mask, bio12.mask, bio15.mask))
names(predictors3) <- c("elev", "dcoast", "globcov", "ecoreg", "meanJulNDVI", "sdJulNDVI", "npp", "soilParent", "treecov", "bio1", "bio2", "bio3", "bio4", "bio5", "bio8", "bio10", "bio12", "bio15")
me5 <- maxent(predictors3, ipts.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))
#AUC = 0.92. Looks similar.

#try with only predictors that have <1% permutation importance in me5:
predictors4 <- stack(c(elev.mask, dcoast.mask2, ecoreg.mask, bio1.mask, bio2.mask, bio8.mask, bio10.mask))
names(predictors4) <- c("elev", "dcoast", "ecoreg", "bio1", "bio2", "bio8", "bio10")
me6 <- maxent(predictors4, ipts.fen, factors = c("ecoreg"), args=c("-J", "-P", "replicates=5"))

#add topsoilH2O.mask (instead of H2OcapVal.mask2, which had that big hole)
predictors5 <- stack(c(elev.mask, dcoast.mask2, topsoilH2O.mask2, globcov.mask, ecoreg.mask, meanNDVI.mask2, sdNDVI.mask2, npp.mask2, parent2.mask, treecov.mask, bio1.mask, bio2.mask, bio3.mask, bio4.mask, bio5.mask, bio8.mask, bio10.mask, bio12.mask, bio15.mask))
names(predictors5) <- c("elev", "dcoast", "topsoilH2O", "globcov", "ecoreg", "meanJulNDVI", "sdJulNDVI", "npp", "soilParent", "treecov", "bio1", "bio2", "bio3", "bio4", "bio5", "bio8", "bio10", "bio12", "bio15")
me7 <- maxent(predictors5, ipts.fen, factors = c("globcov", "ecoreg", "soilParent"), args=c("-J", "-P", "replicates=5"))

#it's between me6 and me7. Look at them:
r.me6 <- predict(me6, predictors4)
r.me7 <- predict(me7, predictors5)
r.me6.mean <- calc(r.me6, mean)
r.me7.mean <- calc(r.me7, mean)

plot(r.me6.mean, axes=TRUE) #looks great!
writeRaster(r.me6.mean, file="me6_sdm_raster.grd")

pdf(file="me6.pdf")
plot(me6, axes=TRUE, col=colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(40))
dev.off()

plot(r.me7.mean, axes=TRUE) #looks crappy! 

#soil moisture layer is all stipply in Norway, does not give a smooth image there. Try full niche model without soil moisture at all:
r.me5 <- predict(me5, predictors3)
r.mr5.mean <- calc(r.me5, mean)





