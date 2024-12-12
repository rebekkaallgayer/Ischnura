####this is to compare sdm versions#######

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

#things that are common between all
terr_cls<- terrain.colors(100, rev=T)
#read in other layers
fenno_elev<- rast("data/elev_Fenno_ll.tif")
names(fenno_elev)<- "elevation"
fenno_land<- rast("data/land_cover_Fenno_ll.tif")
names(fenno_land)<- "land_cover"
fenno_water<- rast("data/perc_cover_freshwater_Fenno_ll.tif")
names(fenno_water)<- "water_cover"
fenno_coast<- rast( "data/distance_to_coast_ll.tif")
names(fenno_coast)<- "distance_to_coast"
#make resolutions match
crs(fenno_water)<- crs(fenno_coast)
fenno_water<- terra::resample(fenno_water, fenno_coast)

##################################################
##mean
chelsa_mean_lesley<- rast("data/chelsa_yr_mean_crop_lesley_ll.tif")
ipts_train<- read.table("data/sdm_isch_train.txt", sep="\t", header=T)
#make resolutions match
crs(chelsa_mean_lesley)<- crs(fenno_coast)
chelsa_mean_layers<- terra::resample(chelsa_mean_lesley, fenno_coast)
#make extents match
ext(chelsa_mean_layers)<- ext(fenno_coast)
fenno_mean_layers<- c(chelsa_mean_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)
#read in presence
ipts_train<- vect(ipts_train, crs=crs(fenno_coast))
ipts_train$species=1
#run the model
d_mean <- sdmData(species~., ipts_train, predictors= fenno_mean_layers,
                  bg = list(method='gRandom',n=1000))
d_mean
write.sdm(d_mean, filename='data/isch_dmean_chelsa_lesley')
d_mean<- read.sdm('data/isch_dmean_chelsa_lesley.sdd')
m_mean <- sdm(species~., d_mean, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
              test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_mean
write.sdm(m_mean, filename='data/isch_mmean_chelsa_lesley')
m_mean<- read.sdm('data/isch_mmean_chelsa_lesley.sdm')
p_mean <- predict(m_mean, fenno_mean_layers,filename='data/isch_predict_chelsa_mean_lesley.img', overwrite=T)
p_mean<- rast("data/isch_predict_chelsa_mean_lesley.img")
en_mean <- ensemble(m_mean, p_mean, filename='data/ensemble_chelsa_mean_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en_mean<- rast("data/ensemble_chelsa_mean_lesley.img")
names(en_mean)<- "ensemble_weighted"
plot(en_mean, col=terr_cls)
points(ipts_train, cex=0.1)
##evaluate with model df
#need to extract some background points
df_mean <- as.data.frame(d_mean)
df_mean <- data.frame(species=df_mean$species,coords(d_mean))
colnames(df_mean)<-c("species", "lon", "lat")
bg_mean<- df_mean[df_mean$species==0,]
ipts_test<- read.table("data/sdm_isch_test.txt", sep="\t", header=T)
ipts_test$species<-1
ipts_test<- ipts_test[,c('species', 'lon', 'lat')]
ipts_mean_eval<- rbind(ipts_test, bg_mean)
xy_mean <- as.matrix(ipts_mean_eval[,c('lon','lat')])
head(xy_mean)
ex_mean <- terra::extract(en_mean,xy_mean)$ensemble_weighted
#head(p_mean)
#nrow(ipts_mean_eval)
#length(p_mean)
ev_mean <- evaluates(ipts_mean_eval$species,ex_mean)
ev_mean@statistics

##################################################
##2013
chelsa_lesley<- rast("data/chelsa_lesley_ll.tif")
ipts_2013<- read.table("data/sdm_isch_2013.txt", sep="\t", header=T)
#make resolutions match
crs(chelsa_lesley)<- crs(fenno_coast)
chelsa_layers<- terra::resample(chelsa_lesley, fenno_coast)
#make extents match
ext(chelsa_layers)<- ext(fenno_coast)
fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)
#read in presence
ipts_2013<- vect(ipts_2013, crs=crs(fenno_coast))
ipts_2013$species=1
#run the model
d_2013 <- sdmData(species~., ipts_2013, predictors= fenno_layers,
                  bg = list(method='gRandom',n=1000))
d_2013
write.sdm(d_2013, filename='data/isch_d2013_chelsa_lesley')
d_2013<- read.sdm('data/isch_d2013_chelsa_lesley.sdd')
m_2013 <- sdm(species~., d_2013, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
              test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_2013
write.sdm(m_2013, filename='data/isch_m2013_chelsa_lesley')
m_2013<- read.sdm('data/isch_m2013_chelsa_lesley.sdm')
p_2013 <- predict(m_2013, fenno_layers,filename='data/isch_predict_chelsa_2013_lesley.img', overwrite=T)
p_2013<- rast("data/isch_predict_chelsa_2013_lesley.img")
en_2013 <- ensemble(m_2013, p_2013, filename='data/ensemble_chelsa_2013_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en_2013<- rast("data/ensemble_chelsa_2013_lesley.img")
names(en_2013)<- "ensemble_weighted"
plot(en_2013, col=terr_cls)
points(ipts_2013, cex=0.1)
##evaluate with model df
#need to extract some background points
df_2013 <- as.data.frame(d_2013)
df_2013 <- data.frame(species=df_2013$species,coords(d_2013))
colnames(df_2013)<-c("species", "lon", "lat")
xy_2013 <- as.matrix(df_2013[,c('lon','lat')])
head(xy_2013)
ex_2013 <- terra::extract(en_2013,xy_2013)$ensemble_weighted
# head(p)
# nrow(ipts_eval)
# length(p)
ev_2013 <- evaluates(df_2013$species,ex_2013)
ev_2013@statistics

##################################################
##hybrid (2013 env but all presence points)
chelsa_lesley<- rast("data/chelsa_lesley_ll.tif")
ipts_train<- read.table("data/sdm_isch_train.txt", sep="\t", header=T)
#make resolutions match
crs(chelsa_lesley)<- crs(fenno_coast)
chelsa_layers<- terra::resample(chelsa_lesley, fenno_coast)
#make extents match
ext(chelsa_layers)<- ext(fenno_coast)
fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)
#read in presence
ipts_train<- vect(ipts_train, crs=crs(fenno_coast))
ipts_train$species=1
#run the model
d <- sdmData(species~., ipts_train, predictors= fenno_layers,
             bg = list(method='gRandom',n=1000))
d
write.sdm(d, filename='data/isch_d_chelsa_lesley', overwrite=T)
d<- read.sdm('data/isch_d_chelsa_lesley.sdd')
m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m
write.sdm(m, filename='data/isch_m_chelsa_lesley', overwrite=T)
m<- read.sdm('data/isch_m_chelsa_lesley.sdm')
p1 <- predict(m, fenno_layers,filename='data/isch_predict_chelsa_lesley.img', overwrite=T)
p1<- rast("data/isch_predict_chelsa_lesley.img")
#but if you have already run the predict() as above, you can use the outputs 
en1 <- ensemble(m, p1, filename='data/ensemble_chelsa_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
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
#head(xy)
ex <- terra::extract(en1,xy)$ensemble_weighted
#head(p)
#nrow(ipts_eval)
#length(p)
ev <- evaluates(ipts_eval$species,ex)
ev@statistics


##i think the mean model is the best so that's what we'll use going forward
ipts_train<- read.table("data/sdm_isch_train.txt", sep="\t", header=T)
ipts_train<- vect(ipts_train, crs=crs(fenno_coast))
ipts_train$species=1
d_mean<- read.sdm('data/isch_dmean_chelsa_lesley.sdd')
m_mean<- read.sdm('data/isch_mmean_chelsa_lesley.sdm')
p_mean<- rast("data/isch_predict_chelsa_mean_lesley.img")
en_mean<- rast("data/ensemble_chelsa_mean_lesley.img")
names(en_mean)<- "ensemble_weighted"
plot(en_mean, col=terr_cls)
points(ipts_train, cex=0.1)

#read in sweden-specific layers
#the %water layer is SWEREF99 TM projection so I'm going to project them all?
swe_water<- rast("data/perc_cover_freshwater.tif")
swe_land<- rast("data/SWE_msk_cov_tif/SWE_msk_cov.tif")
swe_elev<- rast("data/SWE_msk_alt_tif/SWE_msk_alt.tif")

fenno_coast<- rast( "data/distance_to_coast_ll.tif")
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_ex<- extract(fenno_coast, SWE, xy=T)
swe_coast<- rast(swe_ex[,c(3,4,2)])
swe_coast<- resample(swe_coast, swe_land)
#ext(swe_coast)<- ext(swe_land)
#swe_coast<- extend(swe_coast, swe_land)
names(swe_coast)<- "distance_to_coast"


chelsa_mean_lesley<- rast("data/chelsa_yr_mean_crop_lesley_ll.tif")
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
chelsa_mean_lesley<- terra::resample(chelsa_mean_lesley, swe_coast)

chelsa_ex<- extract(chelsa_mean_lesley, SWE, xy=T)
chelsa_ras<- rast(chelsa_ex[,c(6,7,2)])
for(l in 2:nlyr(chelsa_mean_lesley)){
  chelsa_lyr<- rast(chelsa_ex[,c(6,7,l+1)])
  chelsa_ras<- c(chelsa_ras, chelsa_lyr)
}
chelsa_ras<- terra::resample(chelsa_ras, swe_coast)

#make resolutions match
#crs(chelsa_mean_lesley)<- crs(swe_coast)
#chelsa_mean_lesley<- terra::resample(chelsa_mean_lesley, swe_coast)
#make extents match
#ext(chelsa_mean_lesley)<- ext(swe_coast)
swe_mean_layers<- c(chelsa_ras,swe_elev, swe_land, swe_coast)
swe_proj_layers<- terra::project(swe_mean_layers, "EPSG:3006")
plot(swe_proj_layers)

swe_wat_resample<- resample(swe_water,swe_proj_layers)
swe_proj_layers<- c(swe_proj_layers, swe_wat_resample)
