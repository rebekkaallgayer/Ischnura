###
#This script is for preparing environmental layers and maps 
#for use in SDMs for damselfly presence/absence across
#Sweden, Finland and Norway

##Dr. Rebekka Allgayer, June 2025

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
fenno_land<- rast("data/land_cover_Fenno_ll.tif")

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
#NB: THIS TAKES A LONG TIME TO RUN!

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
fenno_w<- merge(swe_comb, fin_comb)
fenno_w<- merge(fenno_w, nor_comb)
fenno_w[fenno_w>0]<-8
fenno_land_w<- merge(fenno_w, fenno_land)
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

terra::writeRaster(fenno_land_w, "data/fenno_land_wupdate.tif")
fenno_land<- rast("data/fenno_land_wupdate.tif")



#################
#getting env data from CHELSA
#NB the files were pre-processed by Camila so the raw calculations
#are not included here

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
for(y in 14:19){
  print(years[y])
  chelsa_yr<- rast(paste("D:/Post_doc/Ischnura_SDM/CHELSA/chelsa_Biovars_", years[y], "_extv3.tif", sep=""))
  chelsa_fenno<- subset(crop(chelsa_yr, Fenno),c(1,2,8,10))
  crs(chelsa_fenno)<- crs(fenno_coast)
  chelsa_layers<- terra::resample(chelsa_fenno, fenno_coast)
  chelsa_layers<- mask(chelsa_layers, fenno_coast)
  #make extents match
  ext(chelsa_layers)<- ext(fenno_coast)
  fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_water, fenno_coast)
  
  #create stack
  names(fenno_layers)<- c("Temp", "Diurnal_range", "Temp_wet", "Temp_warm", "Elevation", "Land_cover",
                          "Perc_water", "Dist_coast")
  #write raster
  writeRaster(fenno_layers, paste("data/CHELSA/",years[y], "_layers.tif", sep=""), overwrite=T)
  rm(chelsa_yr)
  gc()
}

###getting ERA5 land monthly mean data 
# #API token 8ea23d42-450b-41ac-abef-c7e6b47d6a69

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
    month_range<- range(month_sub) #calcualte min and max
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

#make it match CHELSA
CHELSA_2013<- rast("data/CHELSA/2013_layers.tif")
land_layers<- subset(CHELSA_2013, c(5,6,7,8))
Fenno<- vect("data/Fenno.gpkg")
crs(Fenno)<- crs(CHELSA_2013)

for(y in 1:length(years)){
  y_biovars<- rast(paste("data/era5/biovars_",years[y], ".tif", sep=""))
  y_sub<- subset(y_biovars,c(1,2,8,10))*10 #these values are stored 
  y_sub_resam<- terra::resample(y_sub, subset(CHELSA_2013,1))
  crs(y_sub_resam)<- crs(CHELSA_2013)
  y_sub_crop<- mask(y_sub_resam, Fenno)
  y_layers<- c(y_sub_crop, land_layers)
  names(y_layers)<- names(CHELSA_2013)
  writeRaster(y_layers, paste("data/total_timeseries/", years[y], "_layers.tif", sep=""), overwrite=T)
}


##creating a DYNAMIC sdm (Regional)
#need to link presence points with the environment in their year!
#and need to take all the information from 2000-2022

#make a stack of all env layers for the years instead of taking averages
#all variables need to be in a single stack so that you can pass it to extract()

#get presence records
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
#lon, lat, country, year of records of 2013 and not NA
xy_2013<- isch_dat[(isch_dat$year==2013)&(!is.na(isch_dat$year)),c(5,4, 10,12)]
#just lon, lat
xy_crop<- xy_2013[,1:2]



############################
#specific year for env variables
years<- sort(unique(xy_2000s$year))
years<- seq(2000,2018,1) 

#now i need to match the presence points (per year) to the environment of that year
#year specific data
yr_dat<- xy_2000s[xy_2000s$year==years[1],]
yr_vect<- vect(yr_dat[,1:2], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
env_layers<- rast("data/final_fenno/2000_layers.tif")
env_ext<- terra::extract(env_layers, yr_vect, xy=T, cells=T)[,c(10:12, 2:9)]
#env_ext<- terra::extract(env_layers, yr_vect)[,-1]
env_ext<- env_ext[complete.cases(env_ext), ]
env_ext<- cbind(occ=rep(1,nrow(env_ext)), year=rep(2000,nrow(env_ext)),env_ext)
bg_ext<- sdm::background(env_layers, 100, "gRandom" )#[,-c(1:2)]
#yr_bg$cell<- cellFromXY(env_layers, as.matrix(cbind(yr_bg$x,yr_bg$y)))
bg_ext<- cbind(occ=rep(0,nrow(bg_ext)),year=rep(2000,nrow(bg_ext)),
              cell=cellFromXY(env_layers, as.matrix(cbind(bg_ext$x,bg_ext$y))), bg_ext )

#env_ext<- rbind(env_ext, yr_bg)

for(y in 20:length(years)){
  print(years[y])
  yr_dat<- xy_2000s[xy_2000s$year==years[y],] #get the position data
  yr_vect<- vect(yr_dat[,1:2], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
  yr_env<- rast(paste("data/total_timeseries/", years[y], "_layers.tif", sep=""))
  yr_ext<- terra::extract(yr_env, yr_vect, xy=T, cells=T)[,c(10:12, 2:9)]
  yr_ext<- yr_ext[complete.cases(yr_ext), ]
  yr_ext<- cbind(occ=rep(1,nrow(yr_ext)), year=rep(years[y],nrow(yr_ext)), yr_ext)
  yr_bg<- sdm::background(yr_env, 100, "gRandom" )#[,-c(1:2)]
  yr_bg<- cbind(occ=rep(0,nrow(yr_bg)),year=rep(years[y],nrow(yr_bg)),
                cell=cellFromXY(yr_env, as.matrix(cbind(yr_bg$x,yr_bg$y))),yr_bg )
  env_ext<- rbind(env_ext, yr_ext)
  bg_ext<- rbind(bg_ext, yr_bg)
  rm(yr_env)
  gc() #garbage collector
}
write.table(env_ext, "data/total_timeseries/dynamic_pres_dataset.txt", quote=F, row.names=F)
write.table(bg_ext, "data/total_timeseries/dynamic_bg_dataset.txt", quote=F, row.names=F)
env_ext<- read.table("data/total_timeseries/dynamic_pres_dataset.txt", header=T)


#spatial autocorrelation : get rid of duplicates
env_ext_dup<- env_ext[!duplicated(env_ext[c('cell')]), ]
pres<- vect(env_ext[,c(4,5)], geom=c("x", "y"),crs=crs(fenno_coast))
dup<- vect(env_ext_dup[,c(4,5)], geom=c("x", "y"),crs=crs(fenno_coast))
#env_crop<- crop(env_layers, y=ext(4.6, 31.7, 55.2, 60))
plot(env_layers,1)
points(pres, cex=.2)
points(dup, cex=.2, "red")


# Remove adjacent cells of presence/background data:
library(spThin)
occ<- env_ext_dup[,c(4,5,2, 6:13)]
colnames(occ)<- c("x","y", "year", "Temp","Diurnal_range","Temp_wet","Temp_warm","Elevation",
                  "Land_cover","Perc_water","Dist_coast")

variablenames<- colnames(occ)[-(1:3)]
autocorrelation <- spatiotemp_autocorr(occ,
                                       varname = variablenames,
                                       plot = TRUE,
                                       temporal.level = c("year"))
autocorrelation
# $Statistical_tests$Temp$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.2164336 -0.0003652301 0.0002524802       0

# $Statistical_tests$Diurnal_range$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.03443414 -0.0003652301 0.0002524907       0

# $Statistical_tests$Temp_wet$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.01074117 -0.0003652301 0.0002525123       0

# $Statistical_tests$Temp_warm$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.03308299 -0.0003652301 0.0002524563       0

# $Statistical_tests$Elevation$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.01279481 -0.0003652301 0.0002524707       0

# $Statistical_tests$Land_cover$Spatial_autocorrelation
# observed      expected           sd      p.value
# 1 -0.007119536 -0.0003652301 0.0002525302 1.35504e-157

# $Statistical_tests$Perc_water$Spatial_autocorrelation
# observed      expected           sd      p.value
# 1 -0.00291247 -0.0003652301 0.0002524747 6.173895e-24

# $Statistical_tests$Dist_coast$Spatial_autocorrelation
# observed      expected           sd p.value
# 1 -0.04630107 -0.0003652301 0.0002523877       0


#combine presence and pseudo-absence
all_data<- rbind(env_ext_dup, bg_ext)
##split for train and test
train_i <- sample(seq_len(nrow(all_data)), size=round(0.7*nrow(all_data)))

# Then, we can subset the training and testing data
isch_train <- all_data[train_i,]
isch_test <- all_data[-train_i,]
write.table(isch_train, "data/total_timeseries/dynamic_train.txt", quote=F, row.names=F)
write.table(isch_test, "data/total_timeseries/dynamic_test.txt", quote=F, row.names=F)
isch_train<- read.table("data/total_timeseries/dynamic_train.txt", header=T)
isch_test<- read.table("data/total_timeseries/dynamic_test.txt", header=T)

##pass to sdm 
# train_df<- isch_train[,-c(2,3,12)]
# train_df<- env_ext_dup[,-c(2,3,14)]
# d<- sdmData(occ~., train=sp_thinned)
# d<- sdmData(occ~., train=train_df)

#checking memory usage
# sort( sapply(ls(),function(x){object.size(get(x))})) 



##SWEDEN####
###need to create layers for just sweden to correspond to large layers


##working with .gml files for freshwater layer
#library(igraph)
library(sf)
# Fenno<- vect("data/Fenno.gpkg")
# fenno_proj<- terra::project(Fenno,"EPSG:3006")
# plot(Fenno)
# plot(fenno_proj)
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_proj<- terra::project(SWE,"EPSG:3006")
plot(swe_proj)

gml_dir<- "../IschnuraRangeShifter/lakes_gml_inspire/"
gml_files<- list.files(gml_dir, pattern=".gml")

##read in a single file 
sp_dat<- st_read(paste(gml_dir,gml_files[1], sep=""))
plot(sp_dat$geometry)
plot(fenno_proj, add=T)
plot(swe_proj, add=T)


#how to create a new layer, if .gml file overlaps, give cell value of 1
#swe_ras<- rast(swe_proj)

##get environmental data
terr_cls<- terrain.colors(100, rev=T)
fenno_elev<- rast("data/final_fenno/elev_Fenno_ll.tif")
#plot(clim_fenno, col=terr_cls)
clim_proj<- terra::project(fenno_elev, "EPSG:3006") #roughly 500x500km

# #get raster just for sweden? 
# swe_clim<- crop(clim_proj, swe_proj)
# swe_clim[!is.na(swe_clim)]<- 1
# swe_clim
# plot(swe_clim)
# plot(sp_dat$geometry, add=T)
# 
# #rasterise water data
# sp_vect<- vect(sp_dat$geometry)
# plot(sp_vect)
# sp_ras<- terra::rasterize(sp_vect, swe_clim)
# plot(sp_ras)
# swe_water<- swe_clim + sp_ras
# plot(swe_water, col=terr_cls)
# 


swe_blank<- crop(clim_proj, swe_proj)
swe_blank[!is.na(swe_blank)]<-0
swe_water<- swe_blank
for(f in 1:length(gml_files)){
  print(f)
  sp_dat<- st_read(paste(gml_dir,gml_files[f], sep=""))
  sp_vect<- vect(sp_dat$geometry)
  sp_ras<- terra::rasterize(sp_vect, swe_blank, update=T,cover=T)
  if(f==1){
    swe_water<- sp_ras
  }
  else{
    swe_water<- merge(swe_water, sp_ras)
  }
}
swe_water
plot(swe_water)
terra::writeRaster(swe_water, "data/perc_cover_freshwater.tif")
swe_water<- rast("data/perc_cover_freshwater.tif")
plot(swe_water)

swe_blank<- rast("data/final_sweden/static/sweden_blank.tif")
swe_filled<- cover(swe_water, swe_blank)
writeRaster(swe_filled, "data/final_sweden/sweden_perc_water_EPSG3006.tif", overwrite=T)
plot(swe_filled)
swe_elev<- rast("data/SWE_msk_alt_tif/SWE_msk_alt.tif")
swe_water_proj<- terra::project(swe_filled, swe_elev)
plot(swe_water_proj)
writeRaster(swe_water_proj, "data/final_sweden/static/sweden_perc_water_WGS84.tif", overwrite=T)
swe_water_proj<- rast("data/final_sweden/static/sweden_perc_water_WGS84.tif")

SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
crs(SWE)<- crs(fenno_coast)
swe_blank<- mask(fenno_coast, SWE)
swe_blank[!is.na(swe_blank)]<- 0
#swe_blank <- resample(swe_blank, fenno_coast)
writeRaster(swe_blank, "data/final_sweden/sweden_WGS84_blank.tif", overwrite=T)
swe_blank<- rast("data/final_sweden/sweden_WGS84_blank.tif")
years<- seq(2000,2022,1)
swe_water<- rast("data/final_sweden/static/sweden_perc_water_WGS84.tif")
swe_water<- resample(swe_water, swe_blank)

for(y in 1:length(years)){
  print(years[y])
  # yr_env<- rast(paste("data/final_fenno/", years[y], "_layers.tif", sep=""))
  yr_env<- rast(paste("data/total_timeseries/", years[y], "_layers.tif", sep=""))
  swe_msk<- mask(yr_env, swe_blank)
  swe_msk[["Perc_water"]]<- swe_water
  # writeRaster(swe_msk, paste("data/final_sweden/swe_", years[y], "_layers.tif", sep=""))
  writeRaster(swe_msk, paste("data/total_timeseries/sweden/swe_", years[y], "_layers.tif", sep=""), overwrite=T)
  rm(yr_env)
  gc()
}

# swe_2013<- rast("data/final_sweden/swe_2013_layers.tif")
swe_2013<- rast("data/total_timeseries/sweden/swe_2013_layers.tif")
plot(swe_2013)

##create for sweden+finland
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
#combine
swe_fin<- terra::union(SWE, FIN)
fenno_coast<- rast( "data/final_fenno/distance_to_coast_ll.tif")
crs(swe_fin)<- crs(fenno_coast)

swefin_blank<- mask(fenno_coast, swe_fin)
swefin_blank[!is.na(swefin_blank)]<- 0
#swe_blank <- resample(swe_blank, fenno_coast)
writeRaster(swefin_blank, "data/total_timeseries/swefin_WGS84_blank.tif", overwrite=T)
swefin_blank<- rast("data/total_timeseries/swefin_WGS84_blank.tif")
plot(swefin_blank)

#need to update sweden water
swe_water<- rast("data/final_sweden/static/sweden_perc_water_WGS84.tif")
#make extents match
water_test<- resample(swe_water, swefin_blank)
test1<- cover(water_test, fenno_coast)


years<- seq(2000,2022,1)
for(y in 1:length(years)){
  print(years[y])
  # yr_env<- rast(paste("data/final_fenno/", years[y], "_layers.tif", sep=""))
  yr_env<- rast(paste("data/total_timeseries/", years[y], "_layers.tif", sep=""))
  swefin_msk<- mask(yr_env, swefin_blank)
  swefin_msk[["Perc_water"]]<- cover(water_test, swefin_msk$Perc_water)
  # writeRaster(swe_msk, paste("data/final_sweden/swe_", years[y], "_layers.tif", sep=""))
  writeRaster(swefin_msk, paste("data/total_timeseries/sweden_finland/swefin_", years[y], "_layers.tif", sep=""), overwrite=T)
  rm(yr_env)
  gc()
}


##project habitat suitability layers to metres

#create a blank for 1000m
swe_filled<- rast("data/final_sweden/static/sweden_perc_water_EPSG3006.tif")
res_blank<- rast(ext(swe_filled), resolution=c(1000,1000))
crs(res_blank)<- crs(swe_filled)

years<- seq(2013,2022,1)
for(y in 1:length(years)){
  print(years[y])
  swe_yr<- rast(paste("data/total_timeseries/sweden/ensemble_sweden_", years[y], 
                       ".img", sep=""))
  proj_yr<- terra::project(swe_yr,"EPSG:3006")
  res_yr<- resample(proj_yr, res_blank)
  writeRaster(res_yr, paste("data/total_timeseries/sweden/swe_en_proj_",
                             years[y], "_1km.tif", sep=""), overwrite=T)
}
swe_en_2013_proj<- rast("data/total_timeseries/sweden/swe_en_proj_2013_1km.tif")

##for sweden+finland
swefin_blank<- rast("data/total_timeseries/sweden_finland/swefin_WGS84_blank.tif")
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


####
#convert habitat suitability to discrete habitat types for carrying capacity

#read in density data

dens_data<- read.csv("data/Ischnura_densities_Scandinavia.csv", header=T)

swe_dens<- dens_data[dens_data$Transect=="Sweden",c(4:5,7:8)]
years<- unique(swe_dens$Year)

dens_res<- c()

for(y in 1:length(years)){
  #read in year's data
  dens_yr<- swe_dens[swe_dens$Year==years[y],]
  #turn it into a vector
  dens_vect<- vect(dens_yr[,2:1], geom=c("Longitude.adj", "Latitude.adj"), crs=crs(swe_2013))
  #project it
  #dens_proj<- terra::project(dens_vect,"EPSG:3006" )
  
  #read in habitat suitability
  if(y==1){
    en_data<- rast(paste("data/total_timeseries/sweden/ensemble_sweden_", years[y], 
                         ".img", sep=""))
  }
  else{
    en_data<- rast(paste("data/total_timeseries/sweden/ensemble_sweden_", years[y]-1, 
                         ".img", sep=""))
  }

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

lm_dens_suit<- lm(Captures_per_hour ~ hab_suit, data=dens_res)
summary(lm_dens_suit)
dens_res_0<- dens_res[dens_res$Captures_per_hour>0,]
lm_dens_suit_0<- lm(Captures_per_hour ~ hab_suit, data=dens_res_0)
summary(lm_dens_suit_0)

plot(x=dens_res$hab_suit,y=dens_res$Captures_per_hour, xlab="Habitat Suitability",
     ylab="Captures per hour")
abline(lm_dens_suit, col = "blue", lwd = 2)
abline(lm_dens_suit_0, col = "red", lwd = 2)
text(x=0.4, y=50,paste("R2 =", round(summary(lm_dens_suit)$adj.r.squared,2), sep=""), col="blue")
text(x=0.4, y=40,paste("R2 =", round(summary(lm_dens_suit_0)$adj.r.squared,2), sep=""), col="red")

########################
#deciding how to discretise the suitability
hab_K<- as.data.frame(matrix(NA, nrow=10, ncol=4))
colnames(hab_K)<- c("Min_suit", "Cap_hr", "Cap_hr_0", "Hab_type")
hab_K[,1]<- seq(0,90,10)
hab_K[,4]=0

for(r in 1:nrow(hab_K)){
  hab_K[r,2]<- lm_dens_suit$coefficients[2]*(hab_K[r,1]/100)+lm_dens_suit$coefficients[1]
  hab_K[r,3]<- lm_dens_suit_0$coefficients[2]*(hab_K[r,1]/100)+lm_dens_suit_0$coefficients[1]
  if(hab_K[r,2]>0 || hab_K[r,3] >0){
    hab_K[r,4]=hab_K[r-1,4]+1
  }  
}


# library(ggplot2)
# ggplot(dens_yr, aes(x=hab_suit, y=Captures_per_hour)) +
#   geom_point() +
#   geom_smooth(method=lm) #add linear trend line


#deciding how to discretise the suitability
hab_K<- as.data.frame(matrix(NA, nrow=10, ncol=3))
colnames(hab_K)<- c("Min_suit", "Hab_type", "K")
hab_K[,1]<- seq(0,90,10)
hab_K[,2]<- seq(0,9,1)
hab_K[,3]<- c(0,0,10,100,500, 1500,3000,3500,4000, 4200)

plot(x=hab_K[,1], y=hab_K[,3], xlab="Suitability %", ylab="K (ind/ha)")

swe_hab<- rast("data/total_timeseries/sweden/swe_en_proj_2013_1km.tif")
swe_hab_2022<- rast("data/total_timeseries/sweden/swe_en_proj_2022_1km.tif")

swe_hab[(swe_hab<0.6)]<- 0
swe_hab[(swe_hab>0)]<- floor(swe_hab[(swe_hab>0)]*10-5)
plot(swe_hab)

swe_hab_2022[(swe_hab_2022<0.6)]<- 0
swe_hab_2022[(swe_hab_2022>0)]<- floor(swe_hab_2022[(swe_hab_2022>0)]*10-5)
plot(swe_hab_2022)


par(mfrow=c(1,2))
plot(swe_en_2013)
plot(swe_en_2013_hab)


##create initial species distribution 
env_ext_dup<- read.table("data/total_timeseries/dynamic_pres_dataset.txt", header=T)
pres_2013<- env_ext_dup[env_ext_dup$year<=2013,] #1912 observation

swe_en_proj<- rast("data/total_timeseries/sweden/swe_en_proj_2013_1km.tif")
presence_ras<- swe_en_proj
presence_ras[!(is.na(presence_ras))]<- 0

#i have presence in lat lon
pres_2013_vect<- vect(pres_2013[,4:5], geom=c("x", "y"), crs=crs(swe_2013))
pres_2013_proj<- terra::project(pres_2013_vect,"EPSG:3006" )
#extract cell number

extract_2013<- terra::extract(swe_en_proj,pres_2013_proj, cells=T)
extract_2013<- extract_2013[complete.cases(extract_2013),]

#assign cells 1 and 0 for presence
presence_ras[extract_2013$cell]=1
plot(presence_ras)
writeRaster(presence_ras, "data/total_timeseries/sweden/initial_dist.tif", overwrite=T)

SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_out_proj<- terra::project(SWE, "EPSG:3006")
plot(swe_out_proj, add=T)
plot(presence_ras, col=c("white", "red"))

######
##for sweden+finland
swefin_en_2013<- rast("data/total_timeseries/sweden_finland/ensemble_swefin_2013.img")
swefin_en_2013_proj<- rast("data/total_timeseries/sweden_finland/swefin_en_proj_2013_1km.tif")
#get presence points for sweden and finland
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_2000s<- isch_dat[(isch_dat$year>=2000)&(!is.na(isch_dat$year)),]
xy_2000s<- isch_2000s[,c(5,4, 10,12)]
xy_swefin<- xy_2000s[(xy_2000s$country!="Norway" & xy_2000s$year<=2013),]
#project them to EPSG:3006
swefin_vect<- vect(xy_swefin[,1:2], geom=c("longitude", "latitude"), 
                   crs=crs(swefin_en_2013))
swefin_pres_proj<- terra::project(swefin_vect,"EPSG:3006" )
plot(swefin_en_2013_proj)
points(swefin_vect)

#create an empty raster to be distribution
presence_ras<- swefin_en_2013_proj
presence_ras[!(is.na(presence_ras))]<- 0
swefin_extract<- terra::extract(swefin_en_2013_proj,swefin_pres_proj, cells=T)
swefin_ex_2013<- swefin_extract[complete.cases(swefin_extract),]
swefin_ex_dup<- swefin_ex_2013[!duplicated(swefin_ex_2013[c('cell')]), ]

#assign cells 1 and 0 for presence
presence_ras[swefin_ex_2013$cell]=1
plot(presence_ras)
writeRaster(presence_ras, "data/total_timeseries/sweden_finland/initial_dist.tif", overwrite=T)

##plot
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
#combine
swe_fin<- terra::union(SWE, FIN)
swefin_out_proj<- terra::project(swe_fin, "EPSG:3006")
plot(swefin_out_proj, add=T)
plot(presence_ras, col=c("white", "red"))




