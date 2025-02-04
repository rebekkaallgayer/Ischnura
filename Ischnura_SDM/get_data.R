library(geodata)
library(terra)
library(sdm)
#library(dismo)
library(maps)
library(CoordinateCleaner)
library(rgbif)
library(spThin)
library(igraph)
library(sf)

setwd("C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/")
##working with env layers
#make shapefile of Norway, Finland, and Sweden, called "Fenno"
#countries<-getData("countries")
SWE <- gadm(country='SWE', level=0, path="data") #level 0 means just the country outline
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
plot(SWE)
# FIN <- gadm(country='Finland', level=0, path="data") 
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
# #plot(FIN)
# NOR <- gadm(country='Norway', level=0, path="data") 
NOR<- readRDS("data/gadm/gadm41_NOR_0_pk.rds")
# #plot(NOR)
Fenno1<- terra::union(SWE, FIN)
Fenno<- union(Fenno1, NOR)
plot(Fenno)
writeVector(Fenno, "data/Fenno.gpkg")
Fenno<- vect("data/Fenno.gpkg")

##exploring land cover
nor_land<- rast("data/NOR_msk_cov_tif/NOR_msk_cov.tif")
plot(nor_land)
swe_land<- rast("data/SWE_msk_cov_tif/SWE_msk_cov.tif")
plot(swe_land)
swe_land_proj<- terra::project(swe_land,"EPSG:3006")
terra::writeRaster(swe_land_proj, "data/land_cover_swe.tif", overwrite=T)

fin_land<- rast("data/FIN_msk_cov_tif/FIN_msk_cov.tif")
plot(fin_land)

fenno_land<- merge(swe_land, fin_land)
fenno_land<- merge(fenno_land, nor_land)
plot(fenno_land)
fenno_land_proj<- terra::project(fenno_land,"EPSG:3006")
terra::writeRaster(fenno_land, "data/land_cover_Fenno_ll.tif", overwrite=T)

#elevation
nor_elev<- rast("data/NOR_msk_alt_tif/NOR_msk_alt.tif")
plot(nor_elev)
swe_elev<- rast("data/SWE_msk_alt_tif/SWE_msk_alt.tif")
swe_elev_proj<- terra::project(swe_elev,"EPSG:3006")
terra::writeRaster(swe_elev_proj, "data/elev_swe.tif", overwrite=T)
fin_elev<- rast("data/FIN_msk_alt_tif/FIN_msk_alt.tif")
plot(fin_elev)

fenno_elev<- merge(swe_elev, fin_elev)
fenno_elev<- merge(fenno_elev, nor_elev)
fenno_elev_proj<- terra::project(fenno_elev,"EPSG:3006")
plot(fenno_elev)
terra::writeRaster(fenno_elev, "data/elev_Fenno_ll.tif", overwrite=T)

##distance to coast
fenno_elev<- rast("data/elev_Fenno_ll.tif")
fenno_1<- fenno_elev
fenno_1[!is.na(fenno_1)]<- 1
#plot(fenno_1)
fenno_2<-fenno_1
fenno_2[is.na(fenno_2)]<- 2
#plot(fenno_2)
fenno_NA<- fenno_2
fenno_NA[fenno_NA==1]<-NA
plot(fenno_NA)
fenno_dist<- terra::distance(fenno_NA)
plot(fenno_dist)
fenno_dist_NA<- fenno_dist
fenno_dist_NA[fenno_dist_NA==0]<- NA
plot(fenno_dist_NA)
writeRaster(fenno_dist_NA, "data/distance_to_coast_ll.tif", overwrite=T)

###exploring Norway and FInland water shapefiles
#downloaded data from https://diva-gis.org/data.html
fin_water<- vect("data/FIN_wat/FIN_water_areas_dcw.shp")
#fin_wat_proj<- terra::project(fin_water,"EPSG:3006")
#plot(fin_wat_proj)
fin_line<- vect("data/FIN_wat/FIN_water_lines_dcw.shp")
#fin_line_proj<- terra::project(fin_line,"EPSG:3006")
#plot(fin_line_proj)
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
#fin_proj<- terra::project(FIN,"EPSG:3006")
#plot(fin_proj)
fin_blank<- crop(clim_fenno, FIN)
fin_blank[!is.na(fin_blank)]<- 0
fin_water_ras<- terra::rasterize(fin_water, fin_blank, cover=T)
plot(fin_water_ras)
fin_lines_ras<- terra::rasterize(fin_line, fin_blank, cover=T)
plot(fin_lines_ras)
fin_comb<- merge(fin_lines_ras, fin_water_ras)
plot(fin_comb)
plot(fin_comb, add=T)

nor_water<- vect("data/NOR_wat/NOR_water_areas_dcw.shp")
#nor_wat_proj<- terra::project(nor_water,"EPSG:3006")
#plot(nor_wat_proj)
nor_line<- vect("data/NOR_wat/NOR_water_lines_dcw.shp")
#nor_line_proj<- terra::project(nor_line,"EPSG:3006")
#plot(nor_line_proj)
NOR<- readRDS("data/gadm/gadm41_NOR_0_pk.rds")
#nor_proj<- terra::project(NOR,"EPSG:3006")
#plot(nor_proj)
nor_blank<- crop(clim_fenno, nor_water)
nor_blank[!is.na(nor_blank)]<- 0
nor_water_ras<- terra::rasterize(nor_water, nor_blank, cover=T)
plot(nor_water_ras)
nor_lines_ras<- terra::rasterize(nor_line, nor_blank, cover=T)
plot(nor_lines_ras)
nor_comb<- merge(nor_lines_ras, nor_water_ras)
plot(nor_comb)
plot(nor_comb, add=T)

#getting equivalent sweden data
swe_wat<- vect("data/SWE_wat/SWE_water_areas_dcw.shp")
#swe_wat_proj<- terra::project(swe_wat,"EPSG:3006")
#plot(swe_wat_proj)
swe_line<- vect("data/SWE_wat/SWE_water_lines_dcw.shp")
#swe_line_proj<- terra::project(swe_line,"EPSG:3006")
#plot(swe_line_proj)
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
#swe_proj<- terra::project(SWE,"EPSG:3006")
#plot(swe_proj)
swe_blank<- crop(clim_fenno, swe_wat)
swe_blank[!is.na(swe_blank)]<- 0
swe_water_ras<- terra::rasterize(swe_wat, swe_blank, cover=T)
plot(swe_water_ras)
swe_lines_ras<- terra::rasterize(swe_line, swe_blank, cover=T)
plot(swe_lines_ras)
swe_comb<- merge(swe_lines_ras, swe_water_ras)
plot(swe_comb)
plot(swe_comb, add=T)

fenno_water<- merge(swe_comb, fin_comb)
plot(fenno_water, add=T)
fenno_water<- merge(fenno_water, nor_comb)
plot(fenno_water)

terra::writeRaster(fenno_water, "data/perc_cover_freshwater_Fenno_ll.tif")


##creating a DYNAMIC sdm (Regional)
#i need to link presence points with the environment in their year!
#and need to take all the information from 2000-2022

#make a stack of all env layers for the years instead of taking averages
#all variables need to be in a single stack so that i can pass it to extract()
#so either one mega stack of everything or an array of stacks? does that work?

#get the points for 2000-2022
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_2000s<- isch_dat[(isch_dat$year>=2000)&(!is.na(isch_dat$year)),]

xy_2000s<- isch_2000s[,c(5,4, 10,12)]

#split per year
#extract() function taking year-specific env data for each set of presence points
#combine into one big dataframe to pass to sdmData()
Fenno<- vect("data/Fenno.gpkg")
plot(Fenno)

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
fenno_water<- rast("data/perc_cover_freshwater_Fenno_ll.tif")
names(fenno_water)<- "water_cover"
crs(fenno_water)<- crs(fenno_coast)
fenno_water<- terra::resample(fenno_water, fenno_coast)
#fill in NAs in water cover
fenno_blank<- crop(fenno_coast, Fenno)
fenno_blank[!is.na(fenno_blank)]<- 0
fenno_blank <- resample(fenno_blank, fenno_coast)
plot(fenno_blank)
fenno_filled<- cover(fenno_water, fenno_blank)
plot(fenno_filled)
writeRaster(fenno_filled, "data/final_fenno/perc_cover_freshwater_Fenno.tif")

#specific year for env variables
years<- sort(unique(xy_2000s$year))
years<- seq(2000,2018,1) #i only have Chelsa data until 2018!
for(y in 14:length(years)){
  print(years[y])
  chelsa_yr<- rast(paste("D:/Post_doc/Ischnura_SDM/CHELSA/chelsa_Biovars_", years[y], "_extv3.tif", sep=""))
  chelsa_fenno<- subset(crop(chelsa_yr, Fenno),c(1,2,8,10))
  crs(chelsa_fenno)<- crs(fenno_coast)
  chelsa_layers<- terra::resample(chelsa_fenno, fenno_coast)
  chelsa_layers<- mask(chelsa_layers, fenno_coast)
  #make extents match
  ext(chelsa_layers)<- ext(fenno_coast)
  fenno_layers<- c(chelsa_layers,fenno_elev, fenno_land, fenno_filled, fenno_coast)
  
  #create stack
  names(fenno_layers)<- c("Temp", "Diurnal_range", "Temp_wet", "Temp_warm", "Elevation", "Land_cover",
                          "Perc_water", "Dist_coast")
  #write raster
  writeRaster(fenno_layers, paste("data/final_fenno/",years[y], "_layers.tif", sep=""), overwrite=T)
  rm(chelsa_yr)
}


#year specific data
yr_dat<- xy_2000s[xy_2000s$year==years[1],]
yr_vect<- vect(yr_dat[,1:2], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
env_layers<- rast("data/final_fenno/2000_layers.tif")
env_ext<- terra::extract(env_layers, yr_vect, xy=T)[,c(10:11, 2:9)]
#env_ext<- terra::extract(env_layers, yr_vect)[,-1]
env_ext<- env_ext[complete.cases(env_ext), ]
env_ext<- cbind(occ=rep(1,nrow(env_ext)), env_ext)
yr_bg<- sdm::background(env_layers, 100, "gRandom" )#[,-c(1:2)]
yr_bg<- cbind(occ=rep(0,nrow(yr_bg)),yr_bg )

env_ext<- rbind(env_ext, yr_bg)

for(y in 2:length(years)){
  print(years[y])
  yr_dat<- xy_2000s[xy_2000s$year==years[y],] #get the position data
  yr_vect<- vect(yr_dat[,1:2], geom=c("longitude", "latitude"), crs=crs(fenno_coast))
  yr_env<- rast(paste("data/final_fenno/", years[y], "_layers.tif", sep=""))
  yr_ext<- terra::extract(yr_env, yr_vect, xy=T)[,c(10:11, 2:9)]
  yr_ext<- yr_ext[complete.cases(yr_ext), ]
  yr_ext<- cbind(occ=rep(1,nrow(yr_ext)), yr_ext)
  yr_bg<- sdm::background(yr_env, 100, "gRandom" )#[,-c(1:2)]
  yr_bg<- cbind(occ=rep(0,nrow(yr_bg)),yr_bg )
  env_ext<- rbind(env_ext, yr_ext, yr_bg)
  rm(yr_env)
  gc() #garbage collector
}
write.table(env_ext, "data/final_fenno/dynamic_dataset.txt", quote=F, row.names=F)
env_ext<- read.table("data/final_fenno/dynamic_dataset.txt", header=T)

##spatial thinning
# thin() expects that the data.frame contains a column with the species name
env_ext$sp <- 'Ischnura_elegans'

# Remove adjacent cells of presence/background data:
#library(spThin)
# xy <- thin(env_ext, lat.col='y',long.col='x',spec.col='sp', 
#            thin.par=1,reps=1, write.files=F,locs.thinned.list.return=T)
# 
# # Keep the coordinates with the most presence records
# xy_keep <- xy[[1]]
# 
# # Thin the dataset - here, we first extract the cell numbers for the thinned coordinates and then use these to subset our data frame.
# cells_thinned <- terra::cellFromXY(clim_fenno, xy_keep)
# sp_thinned <- sp_env[sp_env$cell %in% cells_thinned,]

##split for train and test
train_i <- sample(seq_len(nrow(env_ext)), size=round(0.7*nrow(env_ext)))

# Then, we can subset the training and testing data
isch_train <- env_ext[train_i,]
isch_test <- env_ext[-train_i,]
write.table(isch_train, "data/final_fenno/dynamic_train.txt", quote=F, row.names=F)
write.table(isch_test, "data/final_fenno/dynamic_test.txt", quote=F, row.names=F)
isch_train<- read.table("data/final_fenno/dynamic_train.txt", header=T)
isch_test<- read.table("data/final_fenno/dynamic_test.txt", header=T)

##pass to sdm 
train_df<- isch_train[,-c(2,3,12)]
d<- sdmData(occ~., train=train_df)

#checking memory usage
sort( sapply(ls(),function(x){object.size(get(x))})) 



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
plot(SWE)

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
clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")
plot(clim_fenno, col=terr_cls)
clim_proj<- terra::project(clim_fenno$bio1, "EPSG:3006") #roughly 500x500km

#get raster just for sweden? 
swe_clim<- crop(clim_proj, swe_proj)
swe_clim[!is.na(swe_clim)]<- 1
swe_clim
plot(swe_clim)
plot(sp_dat$geometry, add=T)

#rasterise water data
sp_vect<- vect(sp_dat$geometry)
plot(sp_vect)
sp_ras<- terra::rasterize(sp_vect, swe_clim)
plot(sp_ras)
swe_water<- swe_clim + sp_ras
plot(swe_water, col=terr_cls)

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
swe_filled<- swe_water+swe_blank

#perc water
##more detailed water data, this is in EPSG:3006
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_proj<- terra::project(SWE,"EPSG:3006")
swe_water<- rast("data/perc_cover_freshwater.tif")
names(swe_water)<- "water_cover"
plot(swe_water)
swe_blank<- rast("data/final_sweden/sweden_blank.tif")
swe_filled<- cover(swe_water, swe_blank)
writeRaster(swe_filled, "data/final_sweden/sweden_perc_water_EPSG3006.tif", overwrite=T)
plot(swe_filled)
swe_water_proj<- terra::project(swe_filled, swe_elev)
plot(swe_water_proj)
writeRaster(swe_water_proj, "data/final_sweden/sweden_perc_water_WGS84.tif", overwrite=T)


SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
crs(SWE)<- crs(fenno_coast)
swe_blank<- mask(fenno_coast, SWE)
swe_blank[!is.na(swe_blank)]<- 0
#swe_blank <- resample(swe_blank, fenno_coast)
writeRaster(swe_blank, "data/final_sweden/sweden_WGS84_blank.tif")
swe_blank<- rast("data/final_sweden/sweden_WGS84_blank.tif")
years<- seq(2000,2018,1)
swe_water<- rast("data/final_sweden/static/sweden_perc_water_WGS84.tif")
swe_water<- resample(swe_water, swe_blank)

for(y in 1:length(years)){
  print(years[y])
  yr_env<- rast(paste("data/final_fenno/", years[y], "_layers.tif", sep=""))
  swe_msk<- mask(yr_env, swe_blank)
  swe_msk[["Perc_water"]]<- swe_water
  writeRaster(swe_msk, paste("data/final_sweden/swe_", years[y], "_layers.tif", sep=""))
  rm(yr_env)
  gc()
}

swe_2013<- rast("data/final_sweden/swe_2013_layers.tif")
plot(swe_2013)