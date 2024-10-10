library(geodata)
library(terra)
library(sdm)
#library(dismo)
library(maps)
library(CoordinateCleaner)
library(rgbif)

setwd("C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/")

##read in occurrence data
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')

summary(isch_dat)

#filter for country, if needed
isch_sw<- isch_dat[isch_dat$country=="Sweden",]
isch_yr<- isch_dat[isch_dat$year>=2000,]

# occurrence points:
ischpts <- isch_yr[c(5,4)]
#ischpts$species <- c("Ischnura_elegans")
#ischpts <- ischpts[,c(3, 1:2)]
dups <- duplicated(ischpts[,1:2])
ischpts2 <- ischpts[!dups,]
# colnames(ischpts2) <- c("species", "lon", "lat")
colnames(ischpts2) <- c("lon", "lat")
write.table(ischpts2, "data/isch_points_2000.txt", quote=F,row.names=F,sep="\t" )
ischpts2<- read.table("data/isch_points_2000.txt", header=T, sep="\t")
#ipts <- SpatialPointsDataFrame(data=ischpts2, coords = ischpts2[,c(2:3)], proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
ipts <- terra::vect(ischpts2, geom=c("lon","lat"),crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
writeVector(ipts, filename="data/isch_points_2000.gpkg")
ipts<- vect("data/isch_points_2000.gpkg")
plot(Fenno)
points(ischpts2, col="red", cex=0.1)
#points(ipts, col="red", cex=0.1)

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

#select target background points (= Zygoptera), 
#Targeted background sampling using sightings of a wider group of 
#similar organisms, to control for sampling effort (Phillips et al. 2009)
e<- ext(Fenno)[1:4]

# Check for synonyms
#name_suggest(q="sympecma", rank='genus')

#see how many records we have
occ_count(scientificName= "sympecma", basisOfRecord='HUMAN_OBSERVATION', 
          hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","), 
          decimalLongitude=paste(e[[1]], e[[2]], sep=","))
# Check number of records - here filtered to those with coordinate information
isch_rec<- occ_search(scientificName= "Ischnura", basisOfRecord='HUMAN_OBSERVATION',
          hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","),
          decimalLongitude=paste(e[[1]], e[[2]], sep=","))
isch_gbif<- isch_rec$data[,c(1:4,40:44, 64)]
isch_gbif<- isch_gbif[isch_gbif$scientificName %like% "elegans", ]
write.table(isch_gbif, "data/background_ischnura.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
isch_gbif<- read.table("data/background_ischnura.txt", header=T, sep="\t")
isch_vect<- vect(isch_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(isch_vect, Fenno, "intersects")
isch_fenno<- isch_gbif[i,]
write.table(isch_fenno, "data/fenno_ischnura.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
isch_fenno<- vect(isch_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(isch_fenno, "data/fenno_ischnura.rds")
isch_fenno<- readRDS("data/fenno_ischnura.rds")

sym.rec<- occ_search(scientificName= "sympecma", basisOfRecord='HUMAN_OBSERVATION',
                     hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","),
                     decimalLongitude=paste(e[[1]], e[[2]], sep=","))
sym_gbif<- sym.rec$data[colnames(isch_gbif)]
write.table(sym_gbif, "data/background_sympecma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
sym_gbif<- read.table("data/background_sympecma.txt", header=T, sep="\t")
sym_vect<- vect(sym_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(sym_vect, Fenno, "intersects")
sym_fenno<- sym_gbif[i,]
write.table(sym_fenno, "data/fenno_sympecma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
sym_fenno<- vect(sym_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(sym_fenno, "data/fenno_sympecma.rds")
sym_fenno<- readRDS("data/fenno_sympecma.rds")


coe.rec <-occ_search(scientificName= "coenagrion", basisOfRecord='HUMAN_OBSERVATION',
                     hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","),
                     decimalLongitude=paste(e[[1]], e[[2]], sep=","))
coe_gbif<- coe.rec$data[colnames(isch_gbif)]
write.table(coe_gbif, "data/background_coenagrion.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
coe_gbif<- read.table("data/background_coenagrion.txt", header=T, sep="\t")
coe_vect<- vect(coe_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(coe_vect, Fenno, "intersects")
coe_fenno<- coe_gbif[i,]
write.table(coe_fenno, "data/fenno_coenagrion.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
coe_fenno<- vect(coe_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(coe_fenno, "data/fenno_coenagrion.rds")
coe_fenno<- readRDS("data/fenno_coenagrion.rds")

ery.rec <-occ_search(scientificName= "erythromma", basisOfRecord='HUMAN_OBSERVATION',
                     hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","),
                     decimalLongitude=paste(e[[1]], e[[2]], sep=","))
ery_gbif<- ery.rec$data[colnames(isch_gbif)]
write.table(ery_gbif, "data/background_erythromma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
ery_gbif<- read.table("data/background_erythromma.txt", header=T, sep="\t")
ery_vect<- vect(ery_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(ery_vect, Fenno, "intersects")
ery_fenno<- ery_gbif[i,]
write.table(ery_fenno, "data/fenno_erythromma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
ery_fenno<- vect(ery_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(ery_fenno, "data/fenno_erythromma.rds")
ery_fenno<- readRDS("data/fenno_erythromma.rds")


pyr.rec <-occ_search(scientificName= "pyrrhosoma", basisOfRecord='HUMAN_OBSERVATION',
                     hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","),
                     decimalLongitude=paste(e[[1]], e[[2]], sep=","))
pyr_gbif<- pyr.rec$data[colnames(isch_gbif)]
write.table(pyr_gbif, "data/background_pyrrhosoma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
pyr_gbif<- read.table("data/background_pyrrhosoma.txt", header=T, sep="\t")
pyr_vect<- vect(pyr_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(pyr_vect, Fenno, "intersects")
pyr_fenno<- pyr_gbif[i,]
write.table(pyr_fenno, "data/fenno_pyrrhosoma.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
pyr_fenno<- vect(pyr_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(pyr_fenno, "data/fenno_pyrrhosoma.rds")
pyr_fenno<- readRDS("data/fenno_pyrrhosoma.rds")

pla.rec <-occ_search(scientificName= "platycnemis", basisOfRecord='HUMAN_OBSERVATION', 
                     hasCoordinate=T, decimalLatitude=paste(e[[3]], e[[4]], sep=","), 
                     decimalLongitude=paste(e[[1]], e[[2]], sep=","))

pla_gbif<- pla.rec$data[colnames(isch_gbif)]
write.table(pla_gbif, "data/background_platycnemis.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
pla_gbif<- read.table("data/background_platycnemis.txt", header=T, sep="\t")
pla_vect<- vect(pla_gbif[,c(4,3)], geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
i <- is.related(pla_vect, Fenno, "intersects")
pla_fenno<- pla_gbif[i,]
write.table(pla_fenno, "data/fenno_platycnemis.txt", quote=F, sep="\t", row.names=F, col.names=colnames(isch_gbif))
pla_fenno<- vect(pla_gbif[i,],geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno))
saveRDS(pla_fenno, "data/fenno_platycnemis.rds")
pla_fenno<- readRDS("data/fenno_platycnemis.rds")


#plot
plot(Fenno)
points(isch_gbif$decimalLongitude, isch_gbif$decimalLatitude, col="red", cex=0.1)
points(sym_gbif$decimalLongitude, sym_gbif$decimalLatitude, col="blue", cex=0.1)
points(coe_gbif$decimalLongitude, sym_gbif$decimalLatitude, col="green", cex=0.1)
points(ery_gbif$decimalLongitude, sym_gbif$decimalLatitude, col="orange", cex=0.1)
points(pyr_gbif$decimalLongitude, sym_gbif$decimalLatitude, col="pink", cex=0.1)
points(pla_gbif$decimalLongitude, sym_gbif$decimalLatitude, col="purple", cex=0.1)

#create all background points
species<- c("ischnura", "sympecma","coenagrion","erythromma","pyrrhosoma","platycnemis")
bg_points<- c()
for(s in 1:length(species)){
  dat<- read.table(paste("data/fenno_", species[s], ".txt", sep=""),header=T, sep="\t")
  dat$species<- rep(species[s], nrow(dat))
  if(s==1){bg_points<- dat}
  else{bg_points<- rbind(bg_points, dat)}
}
write.table(bg_points,"data/background_fenno.txt",quote=F, sep="\t", row.names=F)
bg_points<- read.table("data/background_fenno.txt", header=T, sep="\t")
bg_vect<- vect(bg_points,geom=c("decimalLongitude", "decimalLatitude"), crs=crs(Fenno) )
saveRDS(bg_vect, "data/background_points.rds")
plot(Fenno)
points(bg_vect, "blue")


##environmental data 1970-2000
# Download global bioclimatic data from worldclim (you may have to set argument 'download=T' for first download, if 'download=F' it will attempt to read from file):
#clim <- geodata::worldclim_tile(var = 'bio', res = 0.5, lon=e[[1]], lat=e[[3]],download = T, path = 'data')
clim_swe <- geodata::worldclim_country(var = 'bio', res = 0.5, country=c("SWE"),download = T, path = 'data')
clim_nor <- geodata::worldclim_country(var = 'bio', res = 0.5, country=c("Norway"),download = T, path = 'data')
clim_fin <- geodata::worldclim_country(var = 'bio', res = 0.5, country=c("Finland"),download = T, path = 'data')

# Now, let's look at the data:
terr_cls<- terrain.colors(100, rev=T)
clim_fenno<- merge(clim_swe, clim_nor, clim_fin)
names(clim_fenno)<- c("bio1", "bio2", "bio3","bio4", "bio5", "bio6",
                      "bio7", "bio8", "bio9","bio10", "bio11", "bio12",
                      "bio13", "bio14", "bio15","bio16", "bio17", "bio18", "bio19")
plot(clim_fenno, col=terr_cls)
clim_na<- focal(clim_fenno,w=3,fun=mean)
terra::writeRaster(clim_na, "data/climate/wc2.1_country/clim_fenno.tif", overwrite=TRUE)
clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")

# We can extract the environmental data for the GBIF coordinates.
# Coordinates are always provided as x/y format, in our case lon/lat.
# We also extract the cellnumbers as this allows checking for duplicates later.
fenno_extract<- terra::extract(x=clim_fenno,
                               y=ischpts2, cells=T)

fenno_nona<- fenno_extract[!is.na(fenno_extract[,2]),]
ischpts_nona<- ischpts2[!is.na(fenno_extract[,2]),]
isch_env<- cbind(ischpts_nona, fenno_nona)
#sum(duplicated(isch_env$cell))
isch_env$Ischnura_elegans<- rep(1, nrow(isch_env))

#do the same for the background points
bg_extract<- terra::extract(x=clim_fenno,
                               y=bg_points[,c(4,3)], cells=T)

bg_nona<- bg_extract[!is.na(bg_extract[,2]),]
bgpts_nona<- bg_points[!is.na(bg_extract[,2]),]
bg_env<- cbind(bgpts_nona[,c(4,3)], bg_nona)
bg_env$Ischnura_elegans<- rep(0, nrow(bg_env))
colnames(bg_env)<- c("lon", "lat", colnames(bg_env[,-c(1:2)]))
write.table(bg_env, "data/bg_points_extract.txt", sep="\t", quote=F, row.names=F)
isch_env<- rbind(isch_env, bg_env)


# Only retain non-duplicated cells (will not work in this example as we don't have duplicates):
isch_env_dup <- isch_env[!duplicated(isch_env$cell),]
write.table(isch_env_dup, "data/isch_env.txt", quote=F,row.names=F,sep="\t" )
isch_env_dup<- read.table("data/isch_env.txt", sep="\t", header=T)


##create cpue landscape for background points
null_rast<- rast(ext(clim_fenno_ex), resolution=res(clim_fenno_ex))
crs(null_rast) <- crs(Fenno)
values(null_rast)<-0
pres_freq<- as.data.frame(table(bg_env[,23]))
pres_freq[,1] <- as.numeric(as.character(pres_freq[,1]))
for(i in 1:nrow(pres_freq)){
  null_rast[pres_freq[i,1]]<- pres_freq[i,2]
}
plot(null_rast, col=terr_cls)
null_coarse<- aggregate(null_rast, fact=10, fun="sum")
plot(null_coarse, col=terr_cls)
null_coarse50<- aggregate(null_rast, fact=50, fun="sum")
plot(null_coarse50, col=terr_cls)
writeRaster(null_rast, "data/cpue_background.tif")
writeRaster(null_coarse50, "data/cpue_background_coarse50.tif")

#using CHELSA

#using https://github.com/MirzaCengic/climatedata/?tab=readme-ov-file

# Get models with all 4 RCP scenarios
# all_models<- check_models()
# models_all_pmip<- all_models[all_models$scenario=="pmip3",]$model
# #"CCSM4"        "CNRM-CM5"     "FGOALS-g2"    "IPSL-CM5A-LR" "MIROC-ESM"    "MPI-ESM-P"   "MRI-CGCM3"
# my_output_directory <- "C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/data/CHELSA"
# 
# 
# chelsa_bioclim <- get_chelsa(output_dir = my_output_directory, model_string=models_all_pmip,
#                              scenario_string = "pmip3",
#                              period = "past")
# 
# chelsa_bioclim <- get_chelsa(output_dir = my_output_directory, model_string="IPSL-CM5A-LR",
#                              scenario_string = "pmip3",
#                              period = "past")


#using https://insileco.r-universe.dev/rchelsa/doc/manual.html#chl_global_mod
#  my_output_directory <- "C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/data/CHELSA"
# 
# chl_global_mod(
#   horizon = 2040,
#   var = "bio",
#   val = 10,
#   model = "ipsl-cm6a",
#   scenario = "ssp126",
#   path = my_output_directory
# )
# 
# # chl_global_obsv_daily(var, month = 1, day = 1, path = ".")
# # 
# # chl_global_obsv_monthly(var, month = 1, year = 1980, path = ".")
# 
# #exploring data
# chl_global_ras<- terra::rast("data/CHELSA/CHELSA_bio10_2011-2040_ipsl-cm6a-lr_ssp126_V.2.1.tif")
