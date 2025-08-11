###
#This script is for preparing environmental layers and maps 
#for use in SDMs for damselfly across Sweden, Finland and Norway

##Dr. Rebekka Allgayer, July 2025

#necessary libraries
library(geodata)
library(terra)
library(sdm)
library(maps)
library(CoordinateCleaner)
library(rgbif)
library(spThin)
library(igraph)
library(sf)
library(dynamicSDM)
library(usdm)

#set your working directory
setwd("C:/Users/Rey/Documents/Ischnura/Handover/")

#make shapefile of Norway, Finland, and Sweden, called "Fenno"
#SWE <- gadm(country='SWE', level=0, path="data") #level 0 means just the country outline
#if you've already downloaded them, simply:
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
#plot(SWE)

# FIN <- gadm(country='Finland', level=0, path="data") 
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
# #plot(FIN)

# NOR <- gadm(country='Norway', level=0, path="data") 
NOR<- readRDS("data/gadm/gadm41_NOR_0_pk.rds")
# #plot(NOR)

#join the countries together
#Fenno1<- terra::union(SWE, FIN)
#Fenno<- union(Fenno1, NOR)
#plot(Fenno)
#writeVector(Fenno, "data/Fenno.gpkg")
Fenno<- vect("data/Fenno.gpkg")

##exploring land cover
nor_land<- rast("data/NOR_msk_cov_tif/NOR_msk_cov.tif")
#plot(nor_land)
swe_land<- rast("data/SWE_msk_cov_tif/SWE_msk_cov.tif")
#plot(swe_land)
fin_land<- rast("data/FIN_msk_cov_tif/FIN_msk_cov.tif")
#plot(fin_land)
fenno_land<- merge(swe_land, fin_land)
fenno_land<- merge(fenno_land, nor_land)
plot(fenno_land)
#terra::writeRaster(fenno_land, "data/land_cover_Fenno_ll.tif", overwrite=T)
#fenno_land<- rast("data/land_cover_Fenno_ll.tif")

#combining some landcover types
fenno_land[(fenno_land<=6) | (fenno_land==10)]<- 1 #Tree cover
fenno_land[(fenno_land==7) | (fenno_land==8)]<- 2 #tree cover flooded
fenno_land[(fenno_land==9) | (fenno_land==17)| (fenno_land==18)]<- 3 #Mosaic
fenno_land[(fenno_land>=11) & (fenno_land<=14)]<- 4 #shrubland
fenno_land[(fenno_land==15)]<- 5 #shrubland flooded
fenno_land[(fenno_land==16)]<- 6 #cultivated/managed
fenno_land[(fenno_land==19)]<- 7 #Bare areas
fenno_land[(fenno_land==20)]<- 8 #water bodies
fenno_land[(fenno_land==21)]<- 9 #snow and ice
fenno_land[(fenno_land==22)]<- 10 #artificial
terra::writeRaster(fenno_land, "data/land_cover_Fenno_ll_combined.tif", overwrite=T)

terra::plot(fenno_land, col=c("darkgreen", "lightgreen",
                              "darkviolet","darkolivegreen3",
                              "darkseagreen", "darkgoldenrod",
                              "burlywood", "deepskyblue3",
                              "darkslategray1", "brown1"))
terra::plot(fenno_land)
fenno_land<- rast("data/land_cover_Fenno_ll_combined.tif")

############################################
#elevation
nor_elev<- rast("data/NOR_msk_alt_tif/NOR_msk_alt.tif")
plot(nor_elev)
swe_elev<- rast("data/SWE_msk_alt_tif/SWE_msk_alt.tif")
fin_elev<- rast("data/FIN_msk_alt_tif/FIN_msk_alt.tif")
plot(fin_elev)

fenno_elev<- merge(swe_elev, fin_elev)
fenno_elev<- merge(fenno_elev, nor_elev)
plot(fenno_elev)
terra::writeRaster(fenno_elev, "data/elev_Fenno_ll.tif", overwrite=T)
fenno_elev<- rast("data/elev_Fenno_ll.tif")

########################################################
##distance to coast
#NB: THIS TAKES A LONG TIME TO RUN!(>1 day)

#need Finland's neighbour to really calculate distance from coast
RUS <- gadm(country='RUS', level=0, path="data") #level 0 means just the country outline
#if you've already downloaded them, simply:
RUS<- readRDS("data/gadm/gadm41_RUS_0_pk.rds")
#join the countries together
Fenno_Rus<- terra::union(Fenno, RUS)
Fenno_Rus_crop<- crop(Fenno_Rus, ext(c(xmin(Fenno),40, 
                                       ymin(Fenno), ymax(Fenno))))
plot(Fenno_Rus_crop)

#create a blank raster to hold values
r <- rast(ext(Fenno_Rus_crop), resolution = res(fenno_elev))  # Adjust resolution as needed
# Rasterize the country outline (fill it)
filled_raster <- rasterize(Fenno_Rus_crop, r, field = 1)  # 'field' can be adjusted for attributes
filled_raster[is.na(filled_raster)]<-2 #where there's NA, make it =2, water
plot(filled_raster)
#rather counter-intuitively, turn the land back to NA
filled_raster[filled_raster==1]<-NA
plot(filled_raster)
#calculate distance of every NA cell to closest neighbouring valid cell
#this can take HOURS!
fenno_dist<- terra::distance(filled_raster)
fenno_dist[fenno_dist==0]<- NA
plot(fenno_dist)
writeRaster(fenno_dist, "data/distance_to_coast_RUS_ll.tif", overwrite=T)
fenno_dist<- rast("data/distance_to_coast_RUS_ll.tif")
#now extract the area we actually want
crs(Fenno)<- crs(fenno_elev)
crs(fenno_dist)<- crs(fenno_elev)
fenno_dist_crop<- mask(fenno_dist, Fenno)
plot(fenno_dist_crop)
writeRaster(fenno_dist_crop, "data/distance_to_coast_ll.tif", overwrite=T)
fenno_coast<- rast("data/distance_to_coast_ll.tif")

#############################################################
###exploring Norway and FInland water shapefiles
#downloaded data from https://diva-gis.org/data.html
fin_water<- vect("data/FIN_wat/FIN_water_areas_dcw.shp")
fin_line<- vect("data/FIN_wat/FIN_water_lines_dcw.shp")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")

r <- rast(ext(fenno_coast), crs=crs(fenno_coast),resolution = res(fenno_coast))  # Adjust resolution as needed
fin_water_ras<- terra::rasterize(fin_water, r, cover=T)
#plot(fin_water_ras)
fin_lines_ras<- terra::rasterize(fin_line, r, cover=T)
#plot(fin_lines_ras)
fin_comb<- merge(fin_lines_ras, fin_water_ras)
plot(fin_comb)

nor_water<- vect("data/NOR_wat/NOR_water_areas_dcw.shp")
nor_line<- vect("data/NOR_wat/NOR_water_lines_dcw.shp")
NOR<- readRDS("data/gadm/gadm41_NOR_0_pk.rds")
nor_water_ras<- terra::rasterize(nor_water, r, cover=T)
#plot(nor_water_ras)
nor_lines_ras<- terra::rasterize(nor_line, r, cover=T)
#plot(nor_lines_ras)
nor_comb<- merge(nor_lines_ras, nor_water_ras)
plot(nor_comb)

#getting equivalent sweden data
swe_wat<- vect("data/SWE_wat/SWE_water_areas_dcw.shp")
swe_line<- vect("data/SWE_wat/SWE_water_lines_dcw.shp")
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")

swe_water_ras<- terra::rasterize(swe_wat, r, cover=T)
#plot(swe_water_ras)
swe_lines_ras<- terra::rasterize(swe_line, r, cover=T)
#plot(swe_lines_ras)
swe_comb<- merge(swe_lines_ras, swe_water_ras)
plot(swe_comb)

#for some reason, can't merge multiple rasters so do this twice
fenno_water<- merge(swe_comb, fin_comb)
fenno_water<- merge(fenno_water, nor_comb)
plot(fenno_water)
fenno_filled<- fenno_coast
fenno_filled[!is.na(fenno_filled)]<- 0
fenno_filled<- merge(fenno_water, fenno_filled)
plot(fenno_filled)

terra::writeRaster(fenno_filled, "data/perc_cover_freshwater_Fenno_ll.tif")
fenno_water<- rast("data/perc_cover_freshwater_Fenno_ll.tif")

#update fenno_land with this waterbody info
#create a copy that has water being habitat type 8
fenno_w<- fenno_water
fenno_w[fenno_w>0]<- 8
fenno_w[fenno_w==0]<- NA
fenno_w<- resample(fenno_w, fenno_land)
fenno_land_w<- merge(fenno_w, fenno_land, first=T)


#this provides another 200k cells of water

terra::plot(fenno_land_w, col=c("darkgreen", "lightgreen",
                                "darkviolet","darkolivegreen3",
                                "darkseagreen", "darkgoldenrod",
                                "burlywood", "deepskyblue3",
                                "darkslategray1", "brown1"))
terra::plot(fenno_land, col=c("darkgreen", "lightgreen",
                              "darkviolet","darkolivegreen3",
                              "darkseagreen", "darkgoldenrod",
                              "burlywood", "deepskyblue3",
                              "darkslategray1", "brown1"))

terra::writeRaster(fenno_land_w, "data/fenno_land_wupdate.tif", overwrite=T)
fenno_land<- rast("data/fenno_land_wupdate.tif")

#concatenate land layers
#NB some extents vary slightly so need to standardise them
#make extents and dimensions match smallest, which is fenno_elev
fenno_land<- resample(fenno_land, fenno_elev)
fenno_coast<- resample(fenno_coast, fenno_elev)
fenno_water<- resample(fenno_water, fenno_elev)

land_layers<- c(fenno_elev, fenno_land, fenno_water, fenno_coast)
names(land_layers)<- c("Elevation", "Land_cover","Perc_water", "Dist_coast")
terra::writeRaster(land_layers, "data/land_layers_fenno.tif", overwrite=T)
land_layers<- rast("data/land_layers_fenno.tif")
  
####################################################
###getting ERA5 land monthly mean data for years 2019-2022
##CHELSA dataset only has until 2018, so we need new data and
#to calculate the biovars for these years
#the three things we need are precipitation, tmin and tmax

#downloaded data from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means?tab=download
#monthly averaged!
era5_land<- rast("data/era5/data_stream-moda.nc")
#the era5 data is just stacked and has no identifying names, so
#for ease of use, give the layers names
era5_names<- c()
months<- c("January", "February","March", "April", "May", "June", "July",
           "August", "September", "October", "November", "December")
years<- c( "2019", "2020", "2021", "2022", "2023")
for(y in 1:length(years)){
  for(m in 1:length(months)){
    era5_names<- c(era5_names, paste(months[m], "_", years[y],sep=""))
  }
}
era5_names
names(era5_land)<- rep(era5_names,2)

#this dataset has both temperature and precipitation
#BUT we need a different temporal resolution for temperature so we mostly
#use precipitation from this file

#era5_temp<- subset(era5_land, 1:60)
era5_precip<- subset(era5_land, 61:120)
counter=1
for(y in 1:length(years)){
  year_sub<- subset(era5_precip,counter:(counter+11))
  writeRaster(year_sub, paste("data/era5/precip_", years[y],".tif",sep=""))
  counter=counter+12
}

#Monthly averaged reanalysis by hour of day 
#for temperature only
era5_temp<- rast("data/era5/data_stream-mnth.nc")
hours<- seq(0,23,1)
temp_names<- c()
for(y in 1:length(years)){
  for(m in 1:length(months)){
    for(h in 1:length(hours)){
      temp_names<- c(temp_names, paste(months[m], "_", years[y],"_", hours[h],sep=""))
    }
  }
}
names(era5_temp)<- temp_names
counter=1

for(y in 1:length(years)){
  print(years[y])
  tmin<- subset(era5_temp,1) #create a blank
  tmax<- subset(era5_temp,1) #create a blank
  for(m in 1:length(months)){
    print(months[m])
    month_sub<- subset(era5_temp, counter:(counter+23)) #subset the right hours for the right month
    month_range<- range(month_sub) #calculate min and max
    if(m==1){ #if it's january
      tmin<- subset(month_range,1) #overwrite the blanks
      tmax<- subset(month_range,2)
    }
    else{
      tmin<- c(tmin, subset(month_range,1)) #otherwise, add a layer
      tmax<- c(tmax, subset(month_range,2))
    }
    counter=counter+24 #increase counter for next month
  }
  names(tmin)<- months #overwrite the names
  names(tmax)<- months
  writeRaster(tmin, paste("data/era5/tmin_", years[y],".tif", sep="")) #write the files
  writeRaster(tmax, paste("data/era5/tmax_", years[y], ".tif", sep=""))
}

#calculate biovars
for(y in 1:length(years)){ #for each year
  #read in the files
  #NB this is using the raster package because dismo needs it to be a raster::raster
  y_tmin<- as(rast(paste("data/era5/tmin_", years[y], ".tif", sep="")), 'Raster')
  y_tmax<- as(rast(paste("data/era5/tmax_", years[y], ".tif", sep="")), 'Raster')
  y_precip<- as(rast(paste("data/era5/precip_", years[y], ".tif", sep="")), 'Raster')
  #calculate biovars
  y_biovars<- dismo::biovars(y_precip, y_tmin, y_tmax)
  #output the files
  writeRaster(y_biovars, paste("data/era5/biovars_",years[y], ".tif", sep=""))
}

#############################################
##formatting CHELSA data
#NB: pre-processing and calculating biovars was done by Camila Rocobado-Penanco

####variable selection
#NB: some variable selection has been carried over from previous
#SDM work done by Lesley but this script still does collinearity analysis
#I use 2013 as a baseline here to do this analysis

chelsa_yr<- rast(paste("D:/Post_doc/Ischnura_SDM/CHELSA/chelsa_Biovars_", 2013, "_extv3.tif", sep=""))
chelsa_yr_crop<- mask(chelsa_yr, Fenno)

#testing correlation and selecting variables from CHELSA
#Calculates variance inflation factor (VIF) for a set of variables and exclude 
#the highly correlated variables from the set through a stepwise procedure. 
#if a variable has a strong linear relationship with at least one other variables, 
#the correlation coefficient would be close to 1, and VIF for that variable would be large.
#A VIF greater than 10 is a signal that the model has a collinearity problem

##check collinearity of certain variables
#Bio1: annual temp 
#Bio2: diurnal range (temp) 
#Bio8: temp of wettest quarter 
#Bio10: temp of warmest quarter

chelsa_sub<- subset(chelsa_yr_crop,c(1,2,8,10))
chelsa_sub_vif<- vif(chelsa_sub)
# Variables      VIF
# 1   layer.1 5.069003
# 2   layer.2 1.123501
# 3   layer.8 1.429671
# 4  layer.10 5.649928
#so there are no collinearity problems with the variables we've chosen

#creating the environmental raster stack for each year
years<- seq(2000,2018,1)
#NB: I suggest doing this in blocks of 5, ie for (y in 1:5) then for(y in 6:10) etc
#because this builds up a lot of memory in R for some reason, even with gc() and removing the big data files

for(y in 1:19){
  print(years[y])
  chelsa_yr<- rast(paste("D:/Post_doc/Ischnura_SDM/CHELSA/chelsa_Biovars_", years[y], "_extv3.tif", sep=""))
  #chelsa_fenno<- subset(crop(chelsa_yr, Fenno),c(1,2,8,10))
  crs(chelsa_yr)<- crs(fenno_coast)
  chelsa_yr<- terra::resample(chelsa_yr, fenno_coast)
  #make extents match
  ext(chelsa_yr)<- ext(fenno_coast)
  chelsa_yr<- mask(chelsa_yr, fenno_coast)

  #write raster
  writeRaster(chelsa_yr, paste("data/CHELSA/",years[y], "_layers.tif", sep=""), overwrite=T)
  rm(chelsa_yr)
  gc()
}

#make era5 match CHELSA and create layers for Fenno
#need Fenno outline
crs(Fenno)<- crs(fenno_coast)
years<- seq(2000,2023,1)
land_layers<- rast("data/land_layers_fenno.tif")

for(y in 1:length(years)){
  print(years[y])
  if(years[y]<=2018){ #for the first 19 years, use CHELSA
    env_y<- rast(paste("data/CHELSA/",years[y], "_layers.tif", sep=""))
    env_sub<- subset(env_y, c(1,2,8,10))/10-273.15
    # original values from CHELSA were K/10
    #these values are stored *10 for memory reasons
    env_sub[[2]]<- env_sub[[2]]+273.15 #this is a range, so didn't need to be converted from K
  }
  else{ #otherwise use era5
    env_y<- rast(paste("data/era5/biovars_",years[y], ".tif", sep=""))
    env_sub<- subset(env_y,c(1,2,8,10))-273.15 #convert from K to C
    env_sub[[2]]<- env_sub[[2]]+273.15 #this is a range, so didn't need to be converted from K
  }
  
  #make sure env layers match land layers
  env_resam<- terra::resample(env_sub, fenno_coast)
  crs(env_resam)<- crs(fenno_coast)
  env_crop<- mask(env_resam, Fenno)

  #create stack
  fenno_layers<- c(env_crop, land_layers)
  names(fenno_layers)<- c("Temp", "Diurnal_range", "Temp_wet", "Temp_warm", "Elevation", "Land_cover",
                          "Perc_water", "Dist_coast")
  #write file
  writeRaster(fenno_layers, paste("data/sdm_fenno/", years[y], "_layers.tif", sep=""), overwrite=T)
  #free up memory
  rm(env_y)
  rm(env_sub)
  rm(env_resam)
  rm(env_crop)
  gc()
}

##############################################
##take a subset for Sweden + Finland
##############################################

##create the outline
##create for sweden+finland
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
#combine
swe_fin<- terra::union(SWE, FIN)
land_layers<- rast("data/land_layers_fenno.tif")
crs(swe_fin)<- crs(land_layers)

#create a blank outline of Sweden+Finland
swefin_blank<- mask(land_layers[[1]], swe_fin)
swefin_blank[!is.na(swefin_blank)]<- 0
swefin_blank<- crop(swefin_blank, swe_fin)
#save it
writeRaster(swefin_blank, "data/sdm_swefin/swefin_WGS84_blank.tif", overwrite=T)
swefin_blank<- rast("data/sdm_swefin/swefin_WGS84_blank.tif")
plot(swefin_blank)

#create a blank outline of just Sweden
crs(SWE)<- crs(land_layers)
swe_blank<- mask(land_layers[[1]], SWE)
swe_blank[!is.na(swe_blank)]<- 0
swe_blank<- crop(swe_blank, SWE)
#save it
writeRaster(swe_blank, "data/sdm_swe/swe_WGS84_blank.tif", overwrite=T)
swe_blank<- rast("data/sdm_swe/swe_WGS84_blank.tif")
plot(swe_blank)

#################
#we have Sweden-specific data for freshwater cover, very detailed!
###need to create layers for just sweden to correspond to large layers

##we have more detailed information for Sweden, specifically
#THE CRS IS DIFFERENT! it's a Swedish one
swe_blank<- rast("data/sdm_swe/swe_WGS84_blank.tif")
swe_proj<- terra::project(swe_blank,"EPSG:3006")
plot(swe_proj)
writeRaster(swe_proj, "data/sdm_swe/swe_EPSG3006_blank.tif", overwrite=T)
swe_proj<- rast("data/sdm_swe/swe_EPSG3006_blank.tif")

##and since we'll incorporate it into the sweden+finland model runs
swefin_blank<- rast("data/sdm_swefin/swefin_WGS84_blank.tif")
swefin_proj<- terra::project(swefin_blank,"EPSG:3006")
plot(swefin_proj)
writeRaster(swefin_proj, "data/sdm_swefin/swefin_EPSG3006_blank.tif", overwrite=T)
swefin_proj<- rast("data/sdm_swefin/swefin_EPSG3006_blank.tif")

##working with .gml files for freshwater layer
#where to find the freshwater files (131 of them)
gml_dir<- "data/lakes_gml_inspire/"
gml_files<- list.files(gml_dir, pattern=".gml")

#combine them all together into one file
#NB: this take A WHILE! don't run if you don't need to
swe_water<- swe_proj
for(f in 1:length(gml_files)){
  print(f)
  sp_dat<- st_read(paste(gml_dir,gml_files[f], sep=""))
  sp_vect<- vect(sp_dat$geometry)
  sp_ras<- terra::rasterize(sp_vect, swe_proj, update=T,cover=T)
  if(f==1){
    swe_water<- sp_ras
  }
  else{
    swe_water<- merge(swe_water, sp_ras)
  }
  rm(sp_dat)
  rm(sp_vect)
  rm(sp_ras)
  gc()
}
swe_water
plot(swe_water)
terra::writeRaster(swe_water, "data/sdm_swe/swe_perc_freshwater_EPSG3006.tif", overwrite=T)
swe_water<- rast("data/sdm_swe/swe_perc_freshwater_EPSG3006.tif")
#fill in all the NAs with 0 if they're on land
swe_filled<- cover(swe_water, swe_proj)
writeRaster(swe_filled, "data/final_sweden/sweden_perc_water_EPSG3006.tif", overwrite=T)
plot(swe_filled)

#put it back into WGS84 to merge with layers used for the ensemble sdm
swe_water_wgs84<- terra::project(swe_filled, land_layers[[1]])
plot(swe_water_wgs84)
writeRaster(swe_water_wgs84, "data/sdm_swe/swe_perc_freshwater_WGS84.tif", overwrite=T)
swe_water_wgs84<- rast("data/sdm_swe/swe_perc_freshwater_WGS84.tif")
#swe_water_wgs84<- resample(swe_water_wgs84, swe_blank)

# SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
# crs(SWE)<- crs(fenno_coast)
# swe_blank<- mask(fenno_coast, SWE)
# swe_blank[!is.na(swe_blank)]<- 0
# #swe_blank <- resample(swe_blank, fenno_coast)
# writeRaster(swe_blank, "data/final_sweden/sweden_WGS84_blank.tif", overwrite=T)
# swe_blank<- rast("data/final_sweden/sweden_WGS84_blank.tif")

##now we need to subset out Sweden+Finland and Sweden on its own
#as well as add in Sweden's water info
years<- seq(2000,2023,1)
#need to make extents of the blanks match the land layers or further steps won't work
swe_blank_large<- resample(swe_blank, land_layers[[1]])
swefin_blank_large<- resample(swefin_blank, land_layers[[1]])
for(y in 1:length(years)){
  print(years[y])
  #read in year's env data
  yr_env<- rast(paste("data/sdm_fenno/", years[y], "_layers.tif", sep=""))
  #subset sweden
  swe_msk<- mask(yr_env, swe_blank_large)
  #create a mosaic that takes the maximum water percentage and overwrite
  swe_msk[["Perc_water"]]<- mosaic(swe_water_wgs84, swe_msk[["Perc_water"]], fun="max")
  writeRaster(swe_msk, paste("data/sdm_swe/swe_", years[y], "_layers.tif", sep=""), overwrite=T)
  
  #do the same for sweden_finland
  swefin_msk<-mask(yr_env, swefin_blank_large) 
  swefin_msk[["Perc_water"]]<- mosaic(swe_water_wgs84, swefin_msk[["Perc_water"]], fun="max")
  writeRaster(swefin_msk, paste("data/sdm_swefin/swefin_", years[y], "_layers.tif", sep=""), overwrite=T)
  
  rm(yr_env)
  gc()
}








# swe_2013<- rast("data/final_sweden/swe_2013_layers.tif")
swe_2013<- rast("data/total_timeseries/sweden/swe_2013_layers.tif")
plot(swe_2013)







swefin_blank_proj<- project(swefin_blank, "EPSG:3006")
writeRaster(swefin_blank_proj, "data/total_timeseries/sweden_finland/swefin_EPSG3006_blank.tif", overwrite=T)
res_blank<- rast(ext(swefin_blank_proj), resolution=c(1000,1000))
crs(res_blank)<- crs(swefin_blank_proj)

years<- seq(2013,2022,1)
for(y in 1:length(years)){
  print(years[y])
  swefin_yr<- rast(paste("data/total_timeseries/sweden_finland/ensemble_swefin_", years[y], 
                         ".img", sep=""))
  proj_yr<- terra::project(swefin_yr,"EPSG:3006")
  res_yr<- resample(proj_yr, res_blank)
  writeRaster(res_yr, paste("data/total_timeseries/sweden_finland/swefin_en_proj_",
                            years[y], "_1km.tif", sep=""), overwrite=T)
}
swefin_en_2013_proj<- rast("data/total_timeseries/sweden_finland/swefin_en_proj_2013_1km.tif")










###########################################
##Take a subset for Sweden only 
##########################################

