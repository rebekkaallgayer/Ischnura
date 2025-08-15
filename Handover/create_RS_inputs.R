##This script prepares input files for RangeShiftR

##Dr. Rebekka Allgayer, August 2025

#required libraries
library(terra)

#set working directory
setwd("C:/Users/Rey/Documents/Ischnura/Handover")

################################
##adjusting resolution
################################

#the sdm outputs are all in WGS84 but RS needs inputs to be in metres
##so we need to project habitat suitability layers to metres

#create a blank for wanted resolution
RS_res<- 1000 #change as needed
#read in EPSG3006 blank for Sweden
swe_blank<- rast("data/sdm_swe/swe_EPSG3006_blank.tif")
#and for Sweden_Finland
swefin_blank<- rast("data/sdm_swefin/swefin_EPSG3006_blank.tif")

#create a new blank for sweden with given resolution
res_blank<- rast(ext(swe_blank), resolution=c(RS_res,RS_res))
swe_res_blank<- resample(swe_blank, res_blank)
plot(swe_res_blank)
  
#and for sweden_finland
res_blank<- rast(ext(swefin_blank), resolution=c(RS_res,RS_res))
swefin_res_blank<- resample(swefin_blank, res_blank)
plot(swefin_res_blank)

#use these to project all years of sdm projections 
years<- seq(2013,2023,1)
for(y in 1:length(years)){
  print(years[y])
  #for sweden
  ens_yr<- rast(paste("data/sdm_swe/ensemble_swe_", years[y], 
                      ".img", sep=""))
  proj_yr<- terra::project(ens_yr,"EPSG:3006")
  res_yr<- resample(proj_yr, swe_res_blank)
  writeRaster(res_yr, paste("data/RS_swe/swe_en_EPSG3006_",
                            years[y], "_1km.tif", sep=""), overwrite=T)
  
  #for sweden_finland
  ens_yr<- rast(paste("data/sdm_swefin/ensemble_swefin_", years[y], 
                      ".img", sep=""))
  proj_yr<- terra::project(ens_yr,"EPSG:3006")
  res_yr<- resample(proj_yr, swefin_res_blank)
  writeRaster(res_yr, paste("data/RS_swefin/swefin_en_EPSG3006_",
                            years[y], "_1km.tif", sep=""), overwrite=T)
  
}


##################################
##turning sdm predictions into habitat suitability/densities
##################################

##calculating catch per unit effort (CPUE) estimates from empirical data
#read in density data

dens_data<- read.csv("data/Ischnura_densities_Scandinavia.csv", header=T)

#using data from sweden and finland to maximise chance of finding a pattern
swefin_dens<- dens_data[(dens_data$Transect=="Sweden") | (dens_data$Transect=="Finland"),c(4:5,7:8)]
years<- unique(swefin_dens$Year)

#to use as a reference for creating vector data with WGS84 
#since survey data was in latlon
swefin_2013<- rast("data/sdm_swefin/ensemble_swefin_2013.img")

dens_res<- c()
for(y in 1:length(years)){
  #read in year's data
  dens_yr<- swefin_dens[swefin_dens$Year==years[y],]
  #turn it into a vector
  dens_vect<- vect(dens_yr[,2:1], geom=c("Longitude.adj", "Latitude.adj"), crs=crs(swefin_2013))
  
  #read in predicted probability of presence
  en_data<- rast(paste("data/sdm_swefin/ensemble_swefin_", years[y],".img", sep=""))
  
  #extract values
  en_dens<- extract(en_data,dens_vect)
  dens_yr<- cbind(dens_yr, hab_suit=en_dens[,c(2)])
  if(y==1){
    dens_res=dens_yr
  }
  else{
    dens_res<- rbind(dens_res, dens_yr)
  }
  
}

#run a linear model on the raw data
lm_dens_suit<- lm(Captures_per_hour ~ hab_suit, data=dens_res)
summary(lm_dens_suit)
#try it when surveys of 0 captures are excluded
dens_res_0<- dens_res[dens_res$Captures_per_hour>0,]
lm_dens_suit_0<- lm(Captures_per_hour ~ hab_suit, data=dens_res_0)
summary(lm_dens_suit_0)

#plot the results
plot(x=dens_res$hab_suit,y=dens_res$Captures_per_hour, xlab="Habitat Suitability",
     ylab="Captures per hour")
abline(lm_dens_suit, col = "blue", lwd = 2)
abline(lm_dens_suit_0, col = "red", lwd = 2)
text(x=0.4, y=50,paste("R2 =", round(summary(lm_dens_suit)$adj.r.squared,2), sep=""), col="blue")
text(x=0.4, y=40,paste("R2 =", round(summary(lm_dens_suit_0)$adj.r.squared,2), sep=""), col="red")


##using the maximum and mean captures per hour in each hab suit BIN
hab_suits<- seq(0,1,0.05) #create bins of habitat suitability
hab_suits_res<- as.data.frame(matrix(NA, nrow=length(hab_suits), ncol=3))
colnames(hab_suits_res)<- c("Hab_Suit", "Max_CPUE", "Mean_CPUE")
hab_suits_res[,1]<- hab_suits

#fill in the data
for(h in 1:(length(hab_suits)-1)){
  hab_r<- subset(dens_res, hab_suit>=hab_suits[h] 
                 & hab_suit<hab_suits[h+1])$Captures_per_hour
  
  if(length(hab_r)>0){
    hab_suits_res[h, 2]<- max(hab_r)
    hab_suits_res[h,3]<- mean(hab_r)
  }
}
#run a linear model for maximum and mean CPUE 
lm_max_dens<- lm(Max_CPUE ~ Hab_Suit, data=hab_suits_res)
summary(lm_max_dens)
lm_mean_dens<- lm(Mean_CPUE ~ Hab_Suit, data=hab_suits_res)
summary(lm_mean_dens)

#plot the results
plot(x=hab_suits_res$Hab_Suit,y=hab_suits_res$Max_CPUE, xlab="Habitat Suitability",
     ylab="Captures per hour", col="blue")
points(x=hab_suits_res$Hab_Suit,y=hab_suits_res$Mean_CPUE, col="red")
abline(lm_max_dens, col = "blue", lwd = 2)
abline(lm_mean_dens, col = "red", lwd = 2)
text(x=0.4, y=50,paste("R2 =", round(summary(lm_max_dens)$adj.r.squared,2), sep=""), col="blue")
text(x=0.4, y=40,paste("R2 =", round(summary(lm_mean_dens)$adj.r.squared,2), sep=""), col="red")

##from these relationships, I've pulled the following assumptions:
#1) anything with habitat suitability <0.15 is unsuitable
#2) anything with habitat suitability 0.15-0.5 has low carrying capacity
#3) anything with habitat suitability 0.5-0.8 has medium carrying capacity
#4) anything with habitat suitability >0.8 has high carrying capacity

#assugning K to habitat type
hab_K<- as.data.frame(matrix(NA, nrow=4, ncol=3))
colnames(hab_K)<- c("Min_suit", "Hab_type", "K")
hab_K[,1]<- c(0, 0.15, 0.5, 0.8)
hab_K[,2]<- c(1,2,3,4)
hab_K[,3]<- c(0,100,250,450)

#creating habitat type layers for RS input
#remember that this is now in EPSG3006 because RS needs it in metres
years<- seq(2013,2023,1)
for(y in 1:length(years)){
  print(years[y])
  #for sweden
  ens_yr<- rast(paste("data/RS_swe/swe_en_EPSG3006_", years[y], 
                      "_1km.tif", sep=""))

  ens_yr[(ens_yr>=0.8)]<- 4
  ens_yr[(ens_yr<0.15)]<- 1
  ens_yr[] = ifelse(ens_yr[]>=0.15 & ens_yr[]<0.5, 2, ens_yr[])
  ens_yr[] = ifelse(ens_yr[]>=0.5 & ens_yr[]<0.8, 3, ens_yr[])
  
  #for some reason there are a few NaNs in here
  ens_yr<- classify(ens_yr, cbind(NaN, NA))

  #plot(ens_yr)
  
  writeRaster(ens_yr, paste("data/RS_swe/hab_type_",
                            years[y], "_1km.tif", sep=""), overwrite=T)
  writeRaster(ens_yr, paste("data/RS_swe/hab_type_", years[y], "_1km.asc", 
                            sep=""), NAflag=-99, overwrite=TRUE)
  
  
  #for sweden_finland
  ens_yr<- rast(paste("data/RS_swefin/swefin_en_EPSG3006_", years[y], 
                      "_1km.tif", sep=""))

  ens_yr[(ens_yr>=0.8)]<- 4
  ens_yr[(ens_yr<0.15)]<- 1
  ens_yr[] = ifelse(ens_yr[]>=0.15 & ens_yr[]<0.5, 2, ens_yr[])
  ens_yr[] = ifelse(ens_yr[]>=0.5 & ens_yr[]<0.8, 3, ens_yr[])
  
  #for some reason there are a few NaNs in here
  ens_yr<- classify(ens_yr, cbind(NaN, NA))
  
  writeRaster(ens_yr, paste("data/RS_swefin/hab_type_",
                            years[y], "_1km.tif", sep=""), overwrite=T)
  writeRaster(ens_yr, paste("data/RS_swefin/hab_type_", years[y], "_1km.asc", 
                            sep=""), NAflag=-99, overwrite=TRUE)
  
}


swefin_hab_2013<-rast("data/RS_swefin/hab_type_2013_1km.tif")
swefin_hab_2023<-rast("data/RS_swefin/hab_type_2023_1km.tif")
par(mfrow=c(1,2))
plot(swefin_hab_2013)
plot(swefin_hab_2023)

####################################################
##create initial species distribution 
##############################################

#read in the presence points
pres_points<- read.table("data/sdm_fenno/dynamic_pres_impute_dataset.txt", header=T)
#take points up to 2013
pres_2013<- pres_points[pres_points$year<=2013,] #3968 observation

#create a raster to hold species presence data
swefin_en_proj<- rast("data/RS_swefin/swefin_en_EPSG3006_2013_1km.tif")
presence_ras<- swefin_en_proj

#to use as a reference for creating vector data with WGS84 
#since pres data was in latlon
swefin_2013<- rast("data/sdm_swefin/ensemble_swefin_2013.img")

#i have presence in lat lon
pres_2013_vect<- vect(pres_2013[,4:5], geom=c("x", "y"), crs=crs(swefin_2013))
pres_2013_proj<- terra::project(pres_2013_vect,"EPSG:3006" )
#extract cell number
extract_2013<- terra::extract(presence_ras,pres_2013_proj, cells=T)
extract_2013<- extract_2013[complete.cases(extract_2013),] #3184
extract_2013_dup<- extract_2013[!duplicated(extract_2013[c('cell')]), ] #1589

#assign cells 1 and 0 for presence
presence_ras[!(is.na(presence_ras))]<- 0
presence_ras[extract_2013$cell]=1
plot(presence_ras)

SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
swfin<- terra::union(SWE, FIN)
swfin_proj<- terra::project(swfin, "EPSG:3006")
plot(presence_ras, col=c("white", "red"))
plot(swfin_proj, add=T)

#to make sure they match the habitat files exactly
swefin_hab_2013<-rast("data/RS_swefin/hab_type_2013_1km.tif")
presence_ras_crop<- crop(presence_ras, swefin_hab_2013)
#try it with 0s instead of NAs
presence_ras_crop[is.na(presence_ras_crop)]<- 0
plot(presence_ras_crop)
writeRaster(presence_ras_crop, "data/RS_swe/initial_dist_1km_0.asc", NAflag=-99, overwrite=TRUE)

#save files
writeRaster(presence_ras, "data/RS_swefin/initial_dist_1km.tif", overwrite=T)
writeRaster(presence_ras, "data/RS_swefin/initial_dist_1km.asc", NAflag=-99, overwrite=TRUE)


##########################
#get it for Sweden only
swe_en_proj<- rast("data/RS_swe/swe_en_EPSG3006_2013_1km.tif")
swe_presence_ras<- swe_en_proj

#extract cell number
extract_2013<- terra::extract(swe_presence_ras,pres_2013_proj, cells=T)
extract_2013<- extract_2013[complete.cases(extract_2013),] #2745
extract_2013_dup<- extract_2013[!duplicated(extract_2013[c('cell')]), ] #1328

#assign cells 1 and 0 for presence
swe_presence_ras[!(is.na(swe_presence_ras))]<- 0
swe_presence_ras[extract_2013$cell]=1
plot(swe_presence_ras)

#plot
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_out_proj<- terra::project(SWE, "EPSG:3006")
plot(swe_presence_ras, col=c("white", "red"))
plot(swe_out_proj, add=T)

#to make sure they match the habitat files exactly
swe_hab_2013<-rast("data/RS_swe/hab_type_2013_1km.tif")
swe_presence_ras_crop<- crop(swe_presence_ras, swe_hab_2013)

#try it with 0s instead of NAs
swe_presence_ras_crop[is.na(swe_presence_ras_crop)]<- 0
plot(swe_presence_ras_crop)

#save files
writeRaster(swe_presence_ras_crop, "data/RS_swe/initial_dist_1km.tif", overwrite=T)
writeRaster(swe_presence_ras_crop, "data/RS_swe/initial_dist_1km.asc", NAflag=-99, overwrite=TRUE)



#####################
##cutting off the south of sweden
#####################
years<- seq(2013,2023,1)
for(y in 1:length(years)){
  print(years[y])
  hab_yr<- rast(paste("data/RS_swe/hab_type_",
                      years[y], "_1km.tif", sep=""))

  
  # Define extent to exclude rows (e.g., remove ~300 rows)
  new_extent <- ext(xmin(hab_yr), xmax(hab_yr), 6450000, ymax(hab_yr))  # Adjust ymin/ymax as needed
  cropped_raster <- crop(hab_yr, new_extent)
  
  writeRaster(cropped_raster, paste("data/RS_swe/hab_type_",
                            years[y], "_1km_crop.tif", sep=""), overwrite=T)
  writeRaster(cropped_raster, paste("data/RS_swe/hab_type_", 
                                    years[y], "_1km_crop.asc", 
                            sep=""), NAflag=-99, overwrite=TRUE)
}

######################
##creating initial distribution for the range front
#####################

##Sweden
spdist_base<- rast("data/RS_swe/initial_dist_1km.tif")
hab_base<- rast("data/RS_swe/hab_type_2013_1km_crop.tif")
spdist_crop<- crop(spdist_base, hab_base)
plot(spdist_crop)
writeRaster(spdist_crop, "data/RS_swe/initial_dist_1km_crop.tif", overwrite=T)
writeRaster(spdist_crop, "data/RS_swe/initial_dist_1km_crop.asc", NAflag=-99, overwrite=TRUE)


###############################################
##edit files for running RSR with dynamic landscapes
##############################################

#it's important that all the NAs in the dynamic landscapes have to match
#I know that yeas 2013-2018 have no issues, but 2019 has new NAs 
#after checking the rest of the years, 2019-2023 have the same NAs
#this is most likely due to using era5 data for these years!

years<- seq(2013,2023,1)
hab_13<- rast("data/RS_swe/hab_type_2013_1km_crop.tif")
hab_19<- rast("data/RS_swe/hab_type_2019_1km_crop.tif")
for(y in 1:length(years)){
  print(years[y])
  hab_yr<- rast(paste("data/RS_swe/hab_type_",
                      years[y], "_1km_crop.tif", sep=""))
  hab_matched<- hab_yr
  if(years[y]<2019){
    hab_matched <- mask(hab_yr, hab_19)
  }
  else{
    hab_matched <- mask(hab_yr, hab_13)
  }
  na_cells_matched <- which(is.na(values(hab_matched)))
  print(length(na_cells_matched))
  
  writeRaster(hab_matched, paste("data/RS_swe/hab_type_",
                                    years[y], "_1km_cropNA.tif", sep=""), 
              datatype="INT4S", overwrite=T)
  writeRaster(hab_matched, paste("data/RS_swe/hab_type_", 
                                    years[y], "_1km_cropNA.asc", 
                                    sep=""), NAflag=-99, 
              datatype="INT4S",overwrite=TRUE)
}
