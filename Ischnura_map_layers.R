library(terra)
library(sp)
#working on reprojecting the occurence points to match the CRS in GIS
setwd("Ischnura/IschnuraRangeShifter/IschnuraRangeShifter/")
#read in the occurence data
occ<- read.csv("Occurrence_data/Ischnura_elegans_2022_10_26.csv",header=T)
sw_occ<- occ[occ$country=="Sweden",]
sites<- read.csv("Occurrence_data/allsites_Ischnura.csv", header=T)
sw_sites<- sites[sites$transect=="Sweden",]
sw_sites$corr_lat<- sw_sites$latitude
sw_sites$corr_lon<- sw_sites$longitude

sw_sites$corr_lat[57:66]<- sw_sites$longitude[57:66]
sw_sites$corr_lon[57:66]<- sw_sites$latitude[57:66]

sw_2013<- sw_sites[sw_sites$year==2013,]
sw_2023<- sw_sites[sw_sites$year==2023,]

sites_latlon<-sw_2023[,c(7,6)]
#colnames(sites_latlon)<- c("Longitude","Latitude")
colnames(sites_latlon)<- c("lon","lat")
# df.sp<- SpatialPointsDataFrame(coords=sites_latlon,data=sw_sites, 
#                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
df.sp<- terra::vect(sites_latlon, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# SWEREF99_Transverse_Mercator
# EPSG:3006

o <- new("Spatial")
proj4string(o) <- CRS("epsg:3006")
sp_reproj<- project(df.sp,"+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" )
#crds(sp_reproj)
writeVector(sp_reproj, "occ_2023.gpkg")

##read in sweden shapefile
sweden<- terra::vect("Shapefiles/se_shp/se.shp")
plot(sweden)
plot(df.sp, col="blue", add=T)

