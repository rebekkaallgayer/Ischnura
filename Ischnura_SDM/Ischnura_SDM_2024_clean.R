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
fenno_elev<- rast("data/final_fenno/elev_Fenno_ll.tif")
names(fenno_elev)<- "elevation"

#land cover types
fenno_land<- rast("data/final_fenno/fenno_land_wupdate.tif")
names(fenno_land)<- "land_cover"

#distance to coast per cell
fenno_coast<- rast( "data/final_fenno/distance_to_coast_ll.tif")
names(fenno_coast)<- "distance_to_coast"

#percentage freshwater cover per cell
fenno_water<- rast("data/final_fenno/perc_cover_freshwater_Fenno.tif")
names(fenno_water)<- "water_cover"

#reading in the presence data with extracted environmental variables
env_ext_dup<- read.table("data/total_timeseries/dynamic_dataset.txt", header=T)
#data has already been split into the 70:30 train:test
isch_train<- read.table("data/total_timeseries/dynamic_train.txt", header=T)
isch_test<- read.table("data/total_timeseries/dynamic_test.txt", header=T)

#this is all the presence points from the entire dataset, 2000-2018
ipts_train<- vect(isch_train, geom=c("x", "y"), crs=crs(fenno_coast))
ipts_train$species=1

#~. means it takes species column and considers the rest as predictors, don't need to individually specify
train_df<- isch_train[,-c(2,3)]
test_df<- isch_test[,-c(2,3)]
d<- sdmData(occ~.+ coords(x+y), train=train_df, test=test_df)

d
write.sdm(d, filename='data/total_timeseries/dynamic_fenno_d', overwrite=T)
d<- read.sdm('data/total_timeseries/dynamic_fenno_d.sdd')

m <- sdm(occ~.+ coords(x+y), d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))

m
write.sdm(m, filename='data/total_timeseries/dynamic_fenno_sdm', overwrite=T)
m<-  read.sdm('data/final_fenno/dynamic_fenno_sdm.sdm')

#this is the raster stack of the mean of 2000-2013!
#mean_layers<- rast("data/final_fenno/static/fenno_layers.tif")
#this is the data from 2013 only
layers_2013<- rast("data/total_timeseries/2013_layers.tif")
layers_2022<- rast("data/total_timeseries/2022_layers.tif")

plot(layers_2013)
# p1 <- predict(m, mean_layers,filename='data/final_fenno/dynamic_predict_mean.img', overwrite=T)
# p1<- rast("data/final_fenno/dynamic_predict_mean.img")
# plot(p1, col=terr_cls)
p13 <- predict(m, layers_2013,filename='data/total_timeseries/dynamic_predict_2013.img', overwrite=T)
p13<- rast("data/total_timeseries/dynamic_predict_2013.img")
plot(p13, col=terr_cls)

p22 <- predict(m, layers_2022,filename='data/total_timeseries/dynamic_predict_2022.img', overwrite=T)
p22<- rast("data/total_timeseries/dynamic_predict_2022.img")
plot(p22, col=terr_cls)
# p1
# names(p1)

#but if you have already run the predict() as above, you can use the outputs 
# en1 <- ensemble(m, p1, filename='data/final_fenno/ensemble_dynamic.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
# en1<- rast("data/ensemble_dynamic.img")
# names(en1)<- "ensemble_weighted"
en13 <- ensemble(m, p13, filename='data/total_timeseries/ensemble_dynamic_2013.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en13<- rast("data/total_timeseries/ensemble_dynamic_2013.img")
names(en13)<- "ensemble_weighted"
en22 <- ensemble(m, p22, filename='data/total_timeseries/ensemble_dynamic_2022.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en22<- rast("data/total_timeseries/ensemble_dynamic_2022.img")
names(en22)<- "ensemble_weighted"

# plot(en1, col=terr_cls)
par(mfrow=c(1,2))
plot(en13, col=terr_cls, main="2013")
plot(en22, col=terr_cls, main="2022")

library(mapview)
cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))
mapview(c(en13, en22), col.regions=cl(200))

#-----------
ch <- en22 - en13
cl2 <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
par(mfrow=c(1,1))
plot(ch,col=cl2(200))
library(rasterVis)
levelplot(ch,par.settings=RdBuTheme(), margin=F)
levelplot(ch,par.settings=viridisTheme(), margin=F)

# isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
# isch_2013<- isch_dat[(isch_dat$year==2013)&(!is.na(isch_dat$year)),]
# vect_2013<- vect(isch_2013[,c(5,4)], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
# points(vect_2013, cex=.4)

##evaluate with model df
#need to extract some background points
df1 <- as.data.frame(d)
#df1$species<- "Ischnura_elegans"
df <- data.frame(species=df1$occ,coords(d))
colnames(df)<-c("species", "lon", "lat")
bg<- df[df$species==0,]
#ipts_test<- read.table("data/sdm_isch_test.txt", sep="\t", header=T)
ipts_test<- read.table("data/final_fenno/dynamic_test.txt", header=T)
ipts_test<-ipts_test[ipts_test$occ==1,]
#ipts_test$species<-1
ipts_test<- ipts_test[,c('occ', 'x', 'y')]
colnames(ipts_test)<- c("species", "lon", "lat")
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
swe_2013<- rast("data/total_timeseries/sweden/swe_2013_layers.tif")
plot(swe_2013)
#points
#swe_pts<- read.table("data/final_sweden/swe_isch.txt", header=T)

#take the big model
m<- read.sdm('data/total_timeseries/dynamic_fenno_sdm.sdm')
#use it to predict onto smaller scale...i think
p_swe <- predict(m, swe_2013,filename='data/total_timeseries/sweden/2013_predict_swe.img', overwrite=T)
p_swe<- rast("data/final_sweden/2013_predict_swe.img")
#but if you have already run the predict() as above, you can use the outputs 
en_swe <- ensemble(m, p_swe, filename='data/total_timeseries/sweden/ensemble_sweden_2013.img',setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
en_swe<- rast("data/total_timeseries/sweden/ensemble_sweden_2013.img")
names(en_swe)<- "ensemble_weighted"
en_swe

plot(en_swe, col=terr_cls)

###create layers for each year
years<- seq(2013,2022,1)
for(y in 1:length(years)){
  print(years[y])
  yr_layer<- rast(paste("data/total_timeseries/sweden/swe_", years[y], "_layers.tif", sep=""))
  p_yr <- predict(m, yr_layer,filename=paste('data/total_timeseries/sweden/', years[y], "_predict_swe.img", sep=""), overwrite=T)
  en_yr <- ensemble(m, p_yr, filename=paste("data/total_timeseries/sweden/ensemble_sweden_", years[y], 
                                            ".img", sep=""),setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
}


##############
##do the same for sweden+finland
years<- seq(2013,2022,1)
for(y in 1:length(years)){
  print(years[y])
  yr_layer<- rast(paste("data/total_timeseries/sweden_finland/swefin_", years[y], "_layers.tif", sep=""))
  p_yr <- predict(m, yr_layer,filename=paste('data/total_timeseries/sweden_finland/', years[y], "_predict_swefin.img", sep=""), overwrite=T)
  en_yr <- ensemble(m, p_yr, filename=paste("data/total_timeseries/sweden_finland/ensemble_swefin_", years[y], 
                                            ".img", sep=""),setting=list(method='weighted',stat='tss',opt=2), overwrite = T)
}


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
