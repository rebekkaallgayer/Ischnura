##
#this script is for running sdms on damselfly data for Sweden, Norway and Finland
##Dr. Rebekka Allgayer, August 2025

#necessary libraries
library(terra)
library(sdm)


setwd("C:/Users/Rey/Documents/Ischnura/Handover")

#colour pallette for plotting maps
terr_cls<- terrain.colors(100, rev=T)

##FENNO
land_layers<- rast("data/land_layers_fenno.tif")

#reading in the presence data with extracted environmental variables and duplicates removed
#data has already been split into the 70:30 train:test
isch_train<- read.table("data/sdm_fenno/dynamic_train.txt", header=T)
isch_test<- read.table("data/sdm_fenno/dynamic_test.txt", header=T)

#make data a vect
ipts_train<- vect(isch_train, geom=c("x", "y"), crs=crs(land_layers))
ipts_train$species=1

#using the sdm package, build the ensemble model
train_df<- isch_train[,-c(2,3)]
test_df<- isch_test[,-c(2,3)]

#prepare the data in the right format
#~. means it takes species column and considers the rest as predictors, don't need to individually specify
d<- sdmData(occ~.+ coords(x+y), train=train_df, test=test_df)

d
write.sdm(d, filename='data/sdm_fenno/dynamic_fenno_d', overwrite=T)
d<- read.sdm('data/sdm_fenno/dynamic_fenno_d.sdd')

#prepare the model: THIS TAKES A WHILE
m <- sdm(occ~.+ coords(x+y), d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))

m
write.sdm(m, filename='data/sdm_fenno/dynamic_fenno_sdm', overwrite=T)
m<-  read.sdm('data/final_fenno/dynamic_fenno_sdm.sdm')

#predict across environmental layers for Fenno (this is for illustrative purposes to compare 2013 and 2022)
#NB: THIS TAKES A WHILE, so don't run it if you don't need the maps
layers_2013<- rast("data/sdm_fenno/2013_layers.tif")
layers_2022<- rast("data/sdm_fenno/2022_layers.tif")

p13 <- predict(m, layers_2013,filename='data/sdm_fenno/dynamic_predict_2013.img', overwrite=T)
p13<- rast("data/sdm_fenno/dynamic_predict_2013.img")
#plot(p13, col=terr_cls)

p22 <- predict(m, layers_2022,filename='data/sdm_fenno/dynamic_predict_2022.img', overwrite=T)
p22<- rast("data/sdm_fenno/dynamic_predict_2022.img")

#if you have already run the predict() as above, you can use the outputs directly
en13 <- ensemble(m, p13, filename='data/sdm_fenno/ensemble_dynamic_2013.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en13<- rast("data/total_timeseries/ensemble_dynamic_2013.img")
names(en13)<- "ensemble_weighted"
en22 <- ensemble(m, p22, filename='data/sdm_fenno/ensemble_dynamic_2022.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en22<- rast("data/total_timeseries/ensemble_dynamic_2022.img")
names(en22)<- "ensemble_weighted"

par(mfrow=c(1,2))
plot(en13, col=terr_cls, main="2013")
plot(en22, col=terr_cls, main="2022")

#calculate differences between two years
ch <- en22 - en13
cl2 <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
par(mfrow=c(1,1))
plot(ch,col=cl2(200))

#############################################
###running sdms for Sweden + Finland
##############
##do the same for sweden+finland
#take the big model
m<- read.sdm('data/sdm_fenno/dynamic_fenno_sdm.sdm')
years<- seq(2013,2023,1)
for(y in 2:length(years)){
  print(years[y])
  yr_layer<- rast(paste("data/sdm_swefin/swefin_", years[y], "_layers.tif", sep=""))
  p_yr <- predict(m, yr_layer,filename=paste('data/sdm_swefin/', years[y], "_predict_swefin.img", sep=""), overwrite=T)
  en_yr <- ensemble(m, p_yr, filename=paste("data/sdm_swefin/ensemble_swefin_", years[y], 
                                            ".img", sep=""),setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
  rm(yr_layer)
  rm(p_yr)
  rm(en_yr)
  gc()
}

#############################################
###running sdms for Sweden only
##############
m<- read.sdm('data/sdm_fenno/dynamic_fenno_sdm.sdm')
years<- seq(2013,2023,1)
for(y in 1:length(years)){
  print(years[y])
  yr_layer<- rast(paste("data/sdm_swe/swe_", years[y], "_layers.tif", sep=""))
  p_yr <- predict(m, yr_layer,filename=paste('data/sdm_swe/', years[y], "_predict_swe.img", sep=""), overwrite=T)
  en_yr <- ensemble(m, p_yr, filename=paste("data/sdm_swe/ensemble_swe_", years[y], 
                                            ".img", sep=""),setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
  rm(yr_layer)
  rm(p_yr)
  rm(en_yr)
  gc()
}

swe_2013<- rast("data/sdm_swe/ensemble_swe_2013.img")
swefin_2013<- rast("data/sdm_swefin/ensemble_swefin_2013.img")
par(mfrow=c(1,2))
plot(swe_2013)
plot(swefin_2013)
