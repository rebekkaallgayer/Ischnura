library(data.table)
library(geodata)
library(terra)
library(sdm)
library(dismo)
library(maps)
library(CoordinateCleaner)
library(rgbif)
library(corrplot)
library(mecofun)
library(usdm)
library(mgcv)
library(predicts)


setwd("C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/")

#read in env data
terr_cls<- terrain.colors(100, rev=T)
clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")
clim_proj<- terra::project(clim_fenno, "EPSG:3006") #roughly 500x500km
plot(clim_proj, col=terr_cls)

#cleaned and combined ischnura points, no duplicates, with extracted env data, spatially thinned!
#and it has already been split into the 70:30 train:test 
#isch_env_dup<- read.table("data/isch_presabs_thinned10.txt", sep="\t", header=T)
isch_env_dup<- read.table("data/isch_train.txt", sep="\t", header=T)

#from sdm package tutorial
v <- vifstep(isch_env_dup[,-c(1:4,24:25)]) #i think threshold is by default 10
#ensemble models from sdm package
#Phisically exclude the collinear variables which are identified 
#using vifcor or vifstep from a set of variables.
clim_fenno_ex <- exclude(clim_proj, v)

#read in other layers
fenno_elev<- rast("data/elev_Fenno.tif")
names(fenno_elev)<- "elevation"
fenno_land<- rast("data/land_cover_Fenno.tif")
names(fenno_land)<- "land_cover"
fenno_water<- rast("data/perc_cover_freshwater_Fenno.tif")
names(fenno_water)<- "water_cover"

#make extents match
ext(clim_fenno_ex)<- ext(fenno_water)
ext(fenno_land)<- ext(fenno_water)
ext(fenno_elev)<- ext(fenno_water)
#make resolutions match
clim_layers<- terra::resample(clim_fenno_ex, fenno_water)
fenno_land<- terra::resample(fenno_land, fenno_water)
fenno_elev<- terra::resample(fenno_elev, fenno_water)
plot(clim_layers)

fenno_layers<- c(clim_layers,fenno_elev, fenno_land, fenno_water)

ipts_train<- vect(isch_env_dup[isch_env_dup$occ==1,1:2], crs=crs(clim_fenno))
ipts_train$species=1
ipts_train<- terra::project(ipts_train,"EPSG:3006" )
#isch_sp<- vect(isch_env_dup[,1:3], crs=crs(clim_fenno))
#isch_sp_proj<- terra::project(isch_sp,"EPSG:3006")

#~. means it takes species column and considers the rest as predictors, don't need to individually specify

d <- sdmData(species~., ipts_train, predictors= fenno_layers, 
             bg = list(method='gRandom',n=1000))
d
m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m

p1 <- predict(m, fenno_layers,filename='isch_predict_bio_land.img', overwrite=T)
#p1<- rast("pr.img")
p1<- rast("isch_predict_bio_land.img")
plot(p1, col=terr_cls)
p1
names(p1)