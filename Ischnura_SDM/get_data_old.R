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
isch_env<- isch_env[!duplicated(isch_env$cell),]
write.table(isch_env, "data/isch_points_extract.txt", sep="\t", quote=F, row.names=F)

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

###generating background points using different pseudo-absence data (from Damaris' tutorials)
#load presence data
ischpts2<- read.table("data/isch_points_2000.txt", header=T, sep="\t")
isch_env<- read.table("data/isch_points_extract.txt", sep="\t", header=T)

plot(Fenno)
points(ischpts2, col="red",pch=19, cex=0.1)
clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")
isch_coords<- isch_env[,c("lon", "lat")]

# Make SpatVector:
presences <- terra::vect( as.matrix(isch_coords), crs=crs(clim_fenno))

# Then, place a buffer of 100 km radius around our presence points
v_buf <- terra::buffer(presences, width=100000)

# Create a background mask with target resolution and extent from climate layers
# Set all raster cells outside the buffer to NA.
bg <- clim_fenno[[1]]
values(bg)[!is.na(values(bg))] <- 1
region_buf <- terra::mask(bg, v_buf)
plot(bg, col='grey90', legend=F)
plot(region_buf, add=T, col='grey60', legend=F)

# Exclude presence locations:
sp_cells <- terra::extract(region_buf, presences, cells=T)$cell
region_buf_exclp <- region_buf
values(region_buf_exclp)[sp_cells] <- NA

# Randomly select background data within the buffer, excluding presence locations. We sample 10 times as many background data as we have presences. To ensure that we find enough samples, we use the argument exhaustive=T
bg_rand_buf <- terra::spatSample(region_buf_exclp, length(presences), "random", na.rm=T, as.points=TRUE, exhaustive=T)
points(bg_rand_buf, pch=19, cex=0.1)
points(presences, pch=19, cex=0.1, col='red')

# First, we prepare the presences data to contain a column indicating 1 for presence.
sp_env <- data.frame(isch_coords, occ=1)

# Second, we make sure the background data have the same columns, and indicate 0 for absence.
bg_rand_buf_df <- data.frame(terra::geom(bg_rand_buf)[,c('x','y')])
summary(bg_rand_buf_df)
names(bg_rand_buf_df) <- c('lon','lat')
bg_rand_buf_df$occ <- 0
summary(bg_rand_buf_df)
# Third, we bind these two data sets
sp_env <- rbind(sp_env, bg_rand_buf_df)
summary(sp_env)
# Last, we join this combined data set with the climate data.
sp_env <- cbind(sp_env, terra::extract(x = clim_fenno, y = sp_env[,c('lon','lat')], cells=T) )
summary(sp_env)

##spatial thinning

# The spThin package requires longitude/latitude coordinates, which we already have.

# thin() expects that the data.frame contains a column with the species name
sp_env$sp <- 'Ischnura_elegans'

# Remove adjacent cells of presence/background data:
xy <- thin(sp_env, lat.col='lat',long.col='lon',spec.col='sp', 
           thin.par=10,reps=1, write.files=F,locs.thinned.list.return=T)

# Keep the coordinates with the most presence records
xy_keep <- xy[[1]]

# Thin the dataset - here, we first extract the cell numbers for the thinned coordinates and then use these to subset our data frame.
cells_thinned <- terra::cellFromXY(clim_fenno, xy_keep)
sp_thinned <- sp_env[sp_env$cell %in% cells_thinned,]
write.table(sp_thinned, "data/isch_presabs_thinned10.txt", sep="\t", row.names=F, quote=F)
# Plot the map and data
plot(bg, col='grey90', legend=F)
points(sp_thinned[,1:2],pch=19,col=c('black','red')[as.factor(sp_thinned$occ)], cex=0.3)

# First, we randomly select 70% of the rows that will be used as training data
train_i <- sample(seq_len(nrow(sp_thinned)), size=round(0.7*nrow(sp_thinned)))

# Then, we can subset the training and testing data
isch_train <- sp_thinned[train_i,]
isch_test <- sp_thinned[-train_i,]

# We store the split information for later:
write(train_i, file='data/indices_traindata.txt')
write.table(isch_train, "data/isch_train.txt", sep="\t", quote=F, row.names=F)
write.table(isch_test, "data/isch_test.txt", sep="\t", quote=F, row.names=F)


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

##creating distance from coast
sea<-rast("data/swe_sea.tif")
plot(sea)
#make all cells that are not 8 NA
sea_8<- subst(sea, 1:7, NA)
plot(sea_8)
sea_dist<- distance(sea_8)
plot(sea_dist)

#but that's just sweden!
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




###need to remake the presence points for sdm package rather than Damaris' tutorials
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_start<- isch_dat[(isch_dat$year>=2000) &(isch_dat$year<=2013) &(!is.na(isch_dat$year)),]
isch_occ<- isch_start[,c(5,4)]
colnames(isch_occ)<- c("lon", "lat")
train_i <- sample(seq_len(nrow(isch_occ)), size=round(0.7*nrow(isch_occ)))

# Then, we can subset the training and testing data
isch_train <- isch_occ[train_i,]
isch_test <- isch_occ[-train_i,]
# We store the split information for later:
write(train_i, file='data/sdm_indices_traindata.txt')
write.table(isch_train, "data/sdm_isch_train.txt", sep="\t", quote=F, row.names=F)
write.table(isch_test, "data/sdm_isch_test.txt", sep="\t", quote=F, row.names=F)

isch_sp<- vect(isch_occ, crs=crs(clim_fenno))
isch_sp$species=1
ipts_sp<- terra::project(isch_sp,"EPSG:3006" )
points(ipts_sp)

##get presence points just for 2013
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_start<- isch_dat[(isch_dat$year==2013) &(!is.na(isch_dat$year)),]
isch_occ<- isch_start[,c(5,4)]
colnames(isch_occ)<- c("lon", "lat")
write.table(isch_occ, "data/sdm_isch_2013.txt", sep="\t", quote=F, row.names=F)


#exploring CHELSA data
chelsa_2013<- rast("D:/Post_doc/Ischnura_SDM/CHELSA/chelsa_Biovars_2013_extv3.tif")

plot(chelsa_ex)
chelsa_fenno<- crop(chelsa_2013, fenno_coast)
v_chelsa<- vifstep(chelsa_fenno)
chelsa_fenno_ex<- exclude(chelsa_fenno, v_chelsa)
writeRaster(chelsa_fenno_ex, "data/chelsa_fenno_ex_ll.tif", overwrite=T)
plot(chelsa_fenno_ex)


#using variables that Lesley used
#stack(c(elev.mask, dcoast.mask2, ecoreg.mask, bio1.mask, bio2.mask, bio8.mask, bio10.mask))

chelsa_2013<- rast("D:/Post_doc/Ischnura_SDM/CHELSA/upto_2013/chelsa_Biovars_2013_extv3.tif")
chelsa_fenno<- crop(chelsa_2013, fenno_coast)
chelsa_lesley<- subset(chelsa_fenno, c(1,2,8,10))
writeRaster(chelsa_lesley, "data/chelsa_lesley_ll.tif", overwrite=T)


#take the mean of 2000-2013
chelsa_ls<- list.files("D:/Post_doc/Ischnura_SDM/CHELSA/upto_2013/")
yr_stack<- rast("D:/Post_doc/Ischnura_SDM/CHELSA/upto_2013/chelsa_Biovars_2000_extv3.tif")
for(f in 2:length(chelsa_ls)){
  yr_ras<- rast(paste("D:/Post_doc/Ischnura_SDM/CHELSA/upto_2013/", chelsa_ls[f],sep=""))
  yr_stack<- yr_stack + yr_ras
  rm(yr_ras)
}
yr_ave<- yr_stack/length(chelsa_ls)
chelsa_mean_fenno<- crop(yr_ave, fenno_coast)
writeRaster(chelsa_mean_fenno, "data/chelsa_yrmean_crop_ll.tif", overwrite=T)
chelsa_mean_lesley<- subset(chelsa_mean_fenno, c(1,2,8,10))
writeRaster(chelsa_mean_lesley, "data/chelsa_yr_mean_crop_lesley_ll.tif", overwrite=T)


camila<- rast("D:/Post_doc/Ischnura_SDM/CHELSA/averaged_biovars.tif")

##SWEDEN####
###need to create layers for just sweden to correspond to large layers

##elevation and land cover
swe_elev<- rast("data/SWE_msk_alt_tif/SWE_msk_alt.tif")
names(swe_elev)<- "elevation"
swe_land<- rast("data/SWE_msk_cov_tif/SWE_msk_cov.tif")
names(swe_land)<- "land_cover"
#sweden outline
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_proj<- terra::project(SWE,crs(swe_elev))
writeVector(swe_proj, "data/final_sweden/sweden_outline.gpkg")
swe_proj<- vect("data/final_sweden/sweden_outline.gpkg")
# plot(swe_elev)
# plot(swe_proj, add=T)

##crop env layers to sweden
chelsa_layers<- rast("data/final_fenno/chelsa_layers_resampled.tif")
chelsa_mean_swe<- terra::crop(chelsa_layers, swe_proj, mask=T)
# plot(chelsa_mean_swe,1)
ext(chelsa_mean_swe)<- ext(swe_elev)
chelsa_mean_swe<- terra::resample(chelsa_mean_swe, swe_elev)
writeRaster(chelsa_mean_swe, "data/final_sweden/chelsa_mean_layers.tif", overwrite=T)

#distance to coast
fenno_coast<- rast("data/distance_to_coast_ll.tif")
swe_coast<- terra::crop(fenno_coast, swe_proj, mask=T)
plot(swe_coast)
ext(swe_coast)<- ext(swe_elev)
swe_coast<- terra::resample(swe_coast, swe_elev)
writeRaster(swe_coast, "data/final_sweden/swe_coast.tif", overwrite=T)

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

#make layers
swe_layers<- c(chelsa_mean_swe,swe_elev, swe_land, swe_water_proj, swe_coast)
plot(swe_layers)
writeRaster(swe_layers,"data/final_sweden/sweden_layers_WGS84.tif", overwrite=T )


##points for just sweden
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_start<- isch_dat[(isch_dat$year>=2000) &(isch_dat$year<=2013) &(!is.na(isch_dat$year)),]
isch_swe<- isch_start[isch_start$country=="Sweden",]
isch_occ<- isch_swe[,c(5,4)]
colnames(isch_occ)<- c("lon", "lat")
write.table(isch_occ, "data/final_sweden/swe_isch.txt", sep="\t", quote=F, row.names=F)
train_i <- sample(seq_len(nrow(isch_occ)), size=round(0.7*nrow(isch_occ)))

# Then, we can subset the training and testing data
isch_train <- isch_occ[train_i,]
isch_test <- isch_occ[-train_i,]
# We store the split information for later:
write(train_i, file='data/final_sweden/sdm_indices_traindata.txt')
write.table(isch_train, "data/final_sweden/swe_isch_train.txt", sep="\t", quote=F, row.names=F)
write.table(isch_test, "data/final_sweden/swe_isch_test.txt", sep="\t", quote=F, row.names=F)
isch_train<- vect(isch_train, crs=crs(swe_coast))
isch_train$species=1

plot(swe_layers,1)
points(isch_train, cex=0.1)


#ok now need to look at new years
#should i just add the locations from the new years to the whole dataset?

#first explore what the points data look like
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')
isch_swe<- isch_dat[isch_dat$country=="Sweden",]
isch_2013<- isch_swe[(isch_swe$year>=2000) &(isch_swe$year<=2013) 
                     &(!is.na(isch_swe$year)),][,c(5,4)]
colnames(isch_2013)<- c("lon", "lat")

isch_plus<- isch_swe[(isch_swe$year>2013) &(!is.na(isch_swe$year)),][,c(5,4)]
colnames(isch_plus)<- c("lon", "lat")

swe_proj<- vect("data/final_sweden/sweden_outline.gpkg")
plot(swe_proj)
pts_2013<- vect(isch_2013, crs=crs(swe_layers$temp))
pts_2013$species=1
pts_plus<- vect(isch_plus, crs=crs(swe_layers$temp))
pts_plus$species=1

plot(swe_proj)
points(pts_2013, cex=0.5, col="blue")
points(pts_plus, cex=0.5, col="red")


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
