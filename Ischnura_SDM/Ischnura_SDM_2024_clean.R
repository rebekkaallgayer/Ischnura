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
#clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")
# clim_proj<- terra::project(clim_fenno, "EPSG:3006") #roughly 500x500km
# plot(clim_proj, col=terr_cls)

chelsa_lesley<- rast("data/chelsa_lesley_ll.tif")
chelsa_mean_lesley<- rast("data/chelsa_yr_mean_crop_lesley_ll.tif")

#cleaned and combined ischnura points, no duplicates, with extracted env data, spatially thinned!
#and it has already been split into the 70:30 train:test 
#isch_env_dup<- read.table("data/isch_presabs_thinned10.txt", sep="\t", header=T)
#isch_env<- read.table("data/isch_train.txt", sep="\t", header=T)
ipts_2013<- read.table("data/sdm_isch_2013.txt", sep="\t", header=T)
ipts_train<- read.table("data/sdm_isch_train.txt", sep="\t", header=T)
#from sdm package tutorial
#v <- vifstep(isch_env[,-c(1:4,24:25)]) #i think threshold is by default 10
#ensemble models from sdm package
#Phisically exclude the collinear variables which are identified 
#using vifcor or vifstep from a set of variables.
# clim_fenno_ex <- exclude(clim_proj, v)
#clim_fenno_ex <- exclude(clim_fenno, v)

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
# crs(clim_fenno_ex)<- crs(fenno_coast)
crs(chelsa_lesley)<- crs(fenno_coast)
crs(chelsa_mean_lesley)<- crs(fenno_coast)
chelsa_layers<- terra::resample(chelsa_lesley, fenno_coast)
chelsa_layers<- terra::resample(chelsa_mean_lesley, fenno_coast)
crs(fenno_water)<- crs(fenno_coast)
fenno_water<- terra::resample(fenno_water, fenno_coast)

#make extents match
ext(chelsa_layers)<- ext(fenno_coast)
#plot(clim_layers)

fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)

# ipts_train<- vect(isch_env_dup[isch_env_dup$occ==1,1:2], crs=crs(clim_fenno))
ipts_2013<- vect(ipts_2013, crs=crs(fenno_coast))
ipts_2013$species=1
ipts_train<- vect(ipts_train, crs=crs(fenno_coast))
ipts_train$species=1
#ipts_train<- terra::project(ipts_train,"EPSG:3006" )
#isch_sp<- vect(isch_env_dup[,1:3], crs=crs(clim_fenno))
#isch_sp_proj<- terra::project(isch_sp,"EPSG:3006")

#~. means it takes species column and considers the rest as predictors, don't need to individually specify

d <- sdmData(species~., ipts_train, predictors= fenno_layers,
             bg = list(method='gRandom',n=1000))
d
d_mean <- sdmData(species~., ipts_train, predictors= fenno_layers,
             bg = list(method='gRandom',n=1000))
d_mean
d_2013 <- sdmData(species~., ipts_2013, predictors= fenno_layers,
             bg = list(method='gRandom',n=1000))
d_2013
m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m
m_mean <- sdm(species~., d_mean, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_mean
m_2013 <- sdm(species~., d_2013, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
              test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_2013

p1 <- predict(m, fenno_layers,filename='data/isch_predict_chelsa_lesley.img', overwrite=T)
p_mean <- predict(m_mean, fenno_layers,filename='data/isch_predict_chelsa_mean_lesley.img', overwrite=T)
p_2013 <- predict(m_2013, fenno_layers,filename='data/isch_predict_chelsa_2013_lesley.img', overwrite=T)
#p1<- rast("pr.img")
p1<- rast("data/isch_predict_chelsa_lesley.img")
p1_mean<- rast("data/isch_predict_chelsa_mean_lesley.img")
p1_mean<- rast("data/isch_predict_chelsa_2013_lesley.img")
plot(p1, col=terr_cls)
p1
names(p1)

#but if you have already run the predict() as above, you can use the outputs 
en1 <- ensemble(m, p1, filename='data/ensemble_chelsa_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en1_mean <- ensemble(m_mean, p_mean, filename='data/ensemble_chelsa_mean_lesley.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en1<- rast("data/ensemble_chelsa_lesley.img")
en1_mean<- rast("data/ensemble_chelsa_mean_lesley.img")
names(en1)<- "ensemble_weighted"
names(en1_mean)<- "ensemble_weighted"

plot(en1, col=terr_cls)
points(ipts_train)
plot(en1_mean, col=terr_cls)
points(ipts_train)

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
