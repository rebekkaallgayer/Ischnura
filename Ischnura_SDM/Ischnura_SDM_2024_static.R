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
library(sp)


setwd("C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/")

#read in env data
terr_cls<- terrain.colors(100, rev=T)
# chelsa_mean_lesley<- rast("data/chelsa_yr_mean_crop_lesley_ll.tif")
# #this is CHELSA data, averaged over 2000-2013 for Norway, Finland and Sweden
# #lat-lon, only bio vars 1 (anual mean temp*100), 2(Mean Diurnal Range (Mean of monthly (max temp - min temp))), 
# #8 (Mean Temperature of Wettest Quarter) and 10(Mean Temperature of Warmest Quarter)
# names(chelsa_mean_lesley)<- c("temp", "diurnal_range", "temp_wet", "temmp_warm")
# #plot(chelsa_mean_lesley)
# 
# #read in other layers
# #elevation
# fenno_elev<- rast("data/elev_Fenno_ll.tif")
# names(fenno_elev)<- "elevation"
# #land cover types
# fenno_land<- rast("data/land_cover_Fenno_ll.tif")
# names(fenno_land)<- "land_cover"
# #percentage freshwater cover per cell
# fenno_water<- rast("data/perc_cover_freshwater_Fenno_ll.tif")
# names(fenno_water)<- "water_cover"
# #distance to coast per cell
# fenno_coast<- rast( "data/distance_to_coast_ll.tif")
# names(fenno_coast)<- "distance_to_coast"
# 
# #make resolutions match
# crs(chelsa_mean_lesley)<- crs(fenno_coast)
# chelsa_layers<- terra::resample(chelsa_mean_lesley, fenno_coast)
# writeRaster(chelsa_layers, "data/final_fenno/chelsa_layers_resampled.tif")
# chelsa_layers<- rast("data/final_fenno/chelsa_layers_resampled.tif")
# crs(fenno_water)<- crs(fenno_coast)
# fenno_water<- terra::resample(fenno_water, fenno_coast)
# 
# #make extents match
# ext(chelsa_layers)<- ext(fenno_coast)
# 
# fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)
# writeRaster(fenno_layers, "data/final_fenno/fenno_layers.tif", overwrite=T)


fenno_layers<- rast("data/final_fenno/static/fenno_layers.tif")
plot(fenno_layers)

#data has already been split into the 70:30 train:test 
ipts_train<- read.table("data/sdm_isch_train.txt", sep="\t", header=T)
#this is all the presence points from 2000-2013
ipts_train<- vect(ipts_train, crs=crs(fenno_coast))
ipts_train$species=1

#plot together
plot(chelsa_mean_lesley,1)
points(ipts_train, cex=.1)

#~. means it takes species column and considers the rest as predictors, don't need to individually specify

d <- sdmData(species~., ipts_train, predictors= fenno_layers,
             bg = list(method='gRandom',n=1000))
d

m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m
write.sdm(m, filename='data/final_fenno/mean_fenno_sdm')
m<-  read.sdm('data/final_fenno/mean_fenno_sdm.sdm')

p1 <- predict(m, fenno_layers,filename='data/final_fenno/isch_predict_chelsa_mean.img', overwrite=T)
p1<- rast("data/final_fenno/isch_predict_chelsa_mean.img")
plot(p1, col=terr_cls)
# p1
# names(p1)

#but if you have already run the predict() as above, you can use the outputs 
en1 <- ensemble(m, p1, filename='data/final_fenno/ensemble_chelsa_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en1<- rast("data/ensemble_chelsa_lesley.img")
names(en1)<- "ensemble_weighted"

plot(en1, col=terr_cls)
points(ipts_train, cex=0.1)

##evaluate with model df
#need to extract some background points
df1 <- as.data.frame(d)
df <- data.frame(species=df1$species,coords(d))
colnames(df)<-c("species", "lon", "lat")
bg<- df[df$species==0,]
ipts_test<- read.table("data/sdm_isch_test.txt", sep="\t", header=T)
ipts_test$species<-1
ipts_test<- ipts_test[,c('species', 'lon', 'lat')]
ipts_eval<- rbind(ipts_test, bg)
xy <- as.matrix(ipts_eval[,c('lon','lat')])
head(xy)
p <- terra::extract(en1,xy)$ensemble_weighted
head(p)
nrow(ipts_eval)
length(p)
ev <- evaluates(ipts_eval$species,p)
ev@statistics

##need to apply the large model to Sweden
#env layers
swe_layers<- rast("data/final_sweden/sweden_layers_WGS84.tif" )
#points
swe_pts<- read.table("data/final_sweden/swe_isch.txt", header=T)

#take the big model
m<- read.sdm('data/final_fenno/mean_fenno_sdm.sdm')
#use it to predict onto smaller scale...i think
p_swe <- predict(m, swe_layers,filename='data/final_sweden/isch_predict_swe.img', overwrite=T)
p_swe<- rast("data/final_sweden/isch_predict_swe.img")
#but if you have already run the predict() as above, you can use the outputs 
en_swe <- ensemble(m, p_swe, filename='data/final_sweden/ensemble_sweden.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en_swe<- rast("data/final_sweden/ensemble_sweden.img")
names(en_swe)<- "ensemble_weighted"
en_swe

plot(en_swe, col=terr_cls)
points(swe_pts, cex=0.1)

#so for the evaluates function, i am extracting for the test points
#from the ensemble model, which gives me the likelihood of presence based on the model
#what i need is a matrix of "species", "lon", "lat"

#create presence data
df_swe<- terra::extract(x=swe_layers,y=swe_pts, cells=T, xy=T)
#remove NAs
df_swe_na<- df_swe[complete.cases(df_swe),] #to make sure we're only using points we have data for
df_swe_pres<- cbind(species=rep(1,nrow(df_swe_na)), df_swe_na[,11:12])

#create background data
swe_bg<- background(swe_layers, 1000, "gRandom")
df_swe_abs<- cbind(species=rep(0,nrow(swe_bg)), swe_bg[,1:2])

#combine
df_swe1<- rbind(df_swe_pres, df_swe_abs)
colnames(df_swe1)<- c("species", "lon", "lat")
swe_xy<- as.matrix(df_swe1[,c('lon','lat')])
head(swe_xy)
p <- terra::extract(en_swe,swe_xy)$ensemble_weighted
head(p)
nrow(df_swe1)
length(p)
ev <- evaluates(df_swe1$species,p)
ev@statistics

##need to project this new layer to be in m!
en_swe_proj<- terra::project(en_swe, "EPSG:3006")
plot(en_swe_proj)
