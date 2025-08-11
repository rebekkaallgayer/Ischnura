###
#This script is for preparing presence/absence for use in SDMs
#Sweden, Finland and Norway

##Dr. Rebekka Allgayer, August 2025

#necessary libraries
library(terra)
library(sdm)
library(spThin)
library(dynamicSDM)
library(caret)

# library(geodata)
# 
# library(maps)
# library(CoordinateCleaner)
# library(rgbif)
# library(spThin)
# library(igraph)
# library(sf)

# library(usdm)

#set your working directory
setwd("C:/Users/Rey/Documents/Ischnura/Handover/")

############################################
##creating presence points for a DYNAMIC sdm (Regional: Fenno)
##########################################

#get presence records
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
#get points from 2000 onwards, lat lon, country
xy_2000s<- isch_dat[(isch_dat$year>=2000) &(!is.na(isch_dat$year)),c(5,4, 10,12)]


############################
#specific year for env variables
#use land_layers as examples
land_layers<- rast("data/land_layers_fenno.tif")
Fenno<- vect("data/Fenno.gpkg")
crs(Fenno)<- crs(land_layers)
years<- sort(unique(xy_2000s$year))

#now i need to match the presence points (per year) to the environment of that year
#all variables need to be in a single stack so that you can pass it to extract()
env_ext<-c()
bg_ext<- c()
for(y in 1:length(years)){
  print(years[y])
  yr_dat<- xy_2000s[xy_2000s$year==years[y],] #get the position data
  #turn it into a terra::vect
  yr_vect<- vect(yr_dat[,1:2], geom=c("longitude", "latitude"), crs=crs(land_layers))
  #read in that year's environmental/habitat data
  yr_env<- rast(paste("data/sdm_fenno/", years[y], "_layers.tif", sep=""))
  #match coordinates with env data
  yr_ext<- terra::extract(yr_env, yr_vect, xy=T, cells=T)[,c(10:12, 2:9)]
  #yr_ext<- yr_ext[complete.cases(yr_ext), ]
  yr_ext<- cbind(occ=rep(1,nrow(yr_ext)), year=rep(years[y],nrow(yr_ext)), yr_ext)
  #get pseudo-absence points 
  yr_bg<- sdm::background(yr_env, 100, "gRandom" )#[,-c(1:2)]
  yr_bg<- cbind(occ=rep(0,nrow(yr_bg)),year=rep(years[y],nrow(yr_bg)),
                cell=cellFromXY(yr_env, as.matrix(cbind(yr_bg$x,yr_bg$y))),yr_bg )
  #add to overall dataset
  if(y==1){env_ext<- yr_ext; bg_ext<- yr_bg}
  else{
    env_ext<- rbind(env_ext, yr_ext)
    bg_ext<- rbind(bg_ext, yr_bg)
  }
  rm(yr_env)
  gc() #garbage collector
}
write.table(env_ext, "data/sdm_fenno/dynamic_pres_dataset.txt", quote=F, row.names=F)
write.table(bg_ext, "data/sdm_fenno/dynamic_bg_dataset.txt", quote=F, row.names=F)

env_ext<- read.table("data/sdm_fenno/dynamic_pres_dataset.txt", header=T)
bg_ext<- read.table("data/sdm_fenno/dynamic_bg_dataset.txt", header=T)

################################
##filling NAs near the coast by imputation
##################################

#most NAs are located near the coast, so want to fill in the NAs

NA_impute <- caret::preProcess(env_ext[,-c(1:5)],method="bagImpute")
new_pnt <- terra::predict(NA_impute,env_ext[,-c(1:5)])
env_no_NA <- cbind(env_ext[,1:5],new_pnt)
env_no_NA<- env_no_NA[complete.cases(env_no_NA),]
write.table(env_no_NA, "data/sdm_fenno/dynamic_pres_impute_dataset.txt", quote=F, row.names=F)

#in the case of Land cover, where values are categorical, this has taken means and created new values
#could round up or down, but that doesn't mean much so I have just left it 
#there are only 66 of them so won't skew the model too much

######################################
#spatial autocorrelation
######################################

#get rid of duplicates
env_ext_dup<- env_no_NA[!duplicated(env_no_NA[c('cell')]), ]
write.table(env_ext_dup, "data/sdm_fenno/dynamic_pres_dup_dataset.txt", quote=F, row.names=F)

# test spatial autocorrelation
occ<- env_ext_dup[,c(4,5,2, 6:13)]
colnames(occ)<- c("x","y", "year", "Temp","Diurnal_range","Temp_wet","Temp_warm","Elevation",
                  "Land_cover","Perc_water","Dist_coast")

variablenames<- colnames(occ)[-(1:3)]
#NB: this takes a few minutes:
autocorrelation <- spatiotemp_autocorr(occ,
                                       varname = variablenames,
                                       plot = TRUE,
                                       temporal.level = c("year"))
autocorrelation
# $Statistical_tests$Temp$Spatial_autocorrelation
# observed      expected          sd p.value
# 1 -0.1594322 -0.0002555584 0.000170611       0

# 
# $Statistical_tests$Diurnal_range$Spatial_autocorrelation
#     observed      expected           sd p.value
# 1 -0.0293214 -0.0002555584 0.0001706231       0
# 
# $Statistical_tests$Temp_wet$Spatial_autocorrelation
# observed      expected          sd p.value
# 1 -0.01009252 -0.0002555584 0.000170629       0
# 
# $Statistical_tests$Temp_warm$Spatial_autocorrelation
#      observed      expected           sd p.value
# 1 -0.02989776 -0.0002555584 0.0001706002       0
# 
# $Statistical_tests$Elevation$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.01180911 -0.0002555584 0.0001705965       0
# 
# $Statistical_tests$Land_cover$Spatial_autocorrelation
#       observed      expected           sd      p.value
# 1 -0.003217194 -0.0002555584 0.0001706406 1.777229e-67
# 
# $Statistical_tests$Perc_water$Spatial_autocorrelation
# observed      expected           sd       p.value
# 1 -0.004101675 -0.0002555584 0.0001706123 1.572157e-112
# 
# $Statistical_tests$Dist_coast$Spatial_autocorrelation
#      observed      expected           sd p.value
# 1 -0.04831741 -0.0002555584 0.0001705428       0

#combine presence and pseudo-absence
all_data<- rbind(env_ext_dup, bg_ext)
##split for train and test
train_i <- sample(seq_len(nrow(all_data)), size=round(0.7*nrow(all_data)))

# Then, we can subset the training and testing data
isch_train <- all_data[train_i,]
isch_test <- all_data[-train_i,]
write.table(isch_train, "data/sdm_fenno/dynamic_train.txt", quote=F, row.names=F)
write.table(isch_test, "data/sdm_fenno/dynamic_test.txt", quote=F, row.names=F)
isch_train<- read.table("data/sdm_fenno/dynamic_train.txt", header=T)
isch_test<- read.table("data/sdm_fenno/dynamic_test.txt", header=T)

