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
#elevation
fenno_elev<- rast("data/elev_Fenno_ll.tif")
names(fenno_elev)<- "elevation"

#land cover types
fenno_land<- rast("data/land_cover_Fenno_ll.tif")
names(fenno_land)<- "land_cover"

#distance to coast per cell
fenno_coast<- rast( "data/distance_to_coast_ll.tif")
names(fenno_coast)<- "distance_to_coast"

#percentage freshwater cover per cell
fenno_water<- rast("data/final_fenno/perc_cover_freshwater_Fenno.tif")
names(fenno_water)<- "water_cover"

env_ext<- read.table("data/final_fenno/dynamic_dataset.txt", header=T)
#data has already been split into the 70:30 train:test
isch_train<- read.table("data/final_fenno/dynamic_train.txt", header=T)
isch_test<- read.table("data/final_fenno/dynamic_test.txt", header=T)

#this is all the presence points from the entire dataset, 2000-2018
ipts_train<- vect(isch_train, geom=c("x", "y"), crs=crs(fenno_coast))
ipts_train$species=1

#~. means it takes species column and considers the rest as predictors, don't need to individually specify
train_df<- isch_train[,-c(2,3,12)]
d<- sdmData(occ~., train=train_df)

d
write.sdm(d, filename='data/final_fenno/dynamic_fenno_d')
d<- read.sdm('data/final_fenno/dynamic_fenno_d')

m <- sdm(occ~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m
write.sdm(m, filename='data/final_fenno/dynamic_fenno_sdm')
m<-  read.sdm('data/final_fenno/dynamic_fenno_sdm.sdm')

#this is the raster stack of the mean of 2000-2013!
mean_layers<- rast("data/final_fenno/static/fenno_layers.tif")
#this is the data from 2013 only
layers_2013<- rast("data/final_fenno/2013_layers.tif")

plot(mean_layers)
p1 <- predict(m, mean_layers,filename='data/final_fenno/dynamic_predict_mean.img', overwrite=T)
p1<- rast("data/final_fenno/dynamic_predict_mean.img")
plot(p1, col=terr_cls)
p13 <- predict(m, layers_2013,filename='data/final_fenno/dynamic_predict_2013.img', overwrite=T)
p13<- rast("data/final_fenno/dynamic_predict_2013.img")
plot(p13, col=terr_cls)
# p1
# names(p1)

#but if you have already run the predict() as above, you can use the outputs 
en1 <- ensemble(m, p1, filename='data/final_fenno/ensemble_dynamic.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en1<- rast("data/ensemble_dynamic.img")
names(en1)<- "ensemble_weighted"
en13 <- ensemble(m, p13, filename='data/final_fenno/ensemble_dynamic_2013.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en13<- rast("data/ensemble_dynamic_2013.img")
names(en13)<- "ensemble_weighted"

plot(en1, col=terr_cls)
plot(en13, col=terr_cls)
points(ipts_train, cex=0.1)

isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_2013<- isch_dat[(isch_dat$year==2013)&(!is.na(isch_dat$year)),]
vect_2013<- vect(isch_2013[,c(5,4)], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
points(vect_2013, cex=.4)

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
#2013
swe_2013<- rast("data/final_sweden/swe_2013_layers.tif")
plot(swe_2013)
#points
#swe_pts<- read.table("data/final_sweden/swe_isch.txt", header=T)

#take the big model
m<- read.sdm('data/final_fenno/dynamic_fenno_sdm.sdm')
#use it to predict onto smaller scale...i think
p_swe <- predict(m, swe_2013,filename='data/final_sweden/2013_predict_swe.img', overwrite=T)
p_swe<- rast("data/final_sweden/2013_predict_swe.img")
#but if you have already run the predict() as above, you can use the outputs 
en_swe <- ensemble(m, p_swe, filename='data/final_sweden/ensemble_sweden_2013.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en_swe<- rast("data/final_sweden/ensemble_sweden_2013.img")
names(en_swe)<- "ensemble_weighted"
en_swe

plot(en_swe, col=terr_cls)

#########################


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
