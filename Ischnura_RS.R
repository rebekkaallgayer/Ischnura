##################setting up RangeShiftR#################
setwd("Ischnura/")
# install.packages("Rcpp")
# install.packages("devtools")
# install.packages("Rdpack")
# install.packages("exactextractr")

library(Rcpp)
library(devtools)
library(Rdpack)
#Rcpp::evalCpp("2+2") #should give 4, if not, toolchain not complete
# devtools::install_github("https://github.com/RangeShifter/RangeShiftR-pkg", 
#                          ref = "main", subdir="RangeShiftR")
library(RangeShiftR)
library(terra)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(exactextractr)

##running a tutorial

#Create a parameter master object with all the default settings and store it:
s <- RSsim()
#store path to working directory
dirpath = "RSR_map/"
#Create the RS folder structure, if it doesn’t yet exist:
# dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE)
# dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
# dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)

#RunRS(s,dirpath)
#s

#set simulation parameters
# sim <- Simulation(Simulation = 0,
#                   Years = 20,
#                   Replicates = 2,
#                   OutIntPop = 2)
# 
# 
# 
# ######run a single-cell population model######
# land <- ArtificialLandscape(Resolution = 1000,  # in meters, 1kmx1km gives 100 hectares
#                             K_or_DensDep = 150,  # ~ 15000 inds/cell
#                             propSuit = 1,
#                             dimX = 1, dimY = 1)
# 
# #set demographic parameters with non-overlapping generations
# 
# demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.5)
# 
# #dispersal parameters
# disp <-  Dispersal(Emigration = Emigration(EmigProb = 0), 
#                    Transfer   = DispersalKernel(Distances=1500),
#                    Settlement = Settlement() )
# 
# #set initialisation
# init<- Initialise(InitType= 0, #free initialisation in suitable
#                   FreeType=1, # all cells 
#                   InitDens=1 #at half K
# )
# 
# #put it all together
# s <- RSsim(simul = sim, land = land, 
#            demog = demo, dispersal = disp, 
#            init = init)
# 
# #validate 
# validateRSparams(s)
# 
# #run the model
# RunRS(s, dirpath)
# 
# #plot results
# range_df <- readRange(s, dirpath)
# 
# plotAbundance(range_df)
# plotAbundance(range_df, sd=T, replicates=F)

######for real landscape#######
#set land parameters with imported landscape
Swedenmap <- terra::rast(paste0(dirpath, "Inputs/hab_ras.txt"))
SpDist <- terra::rast(paste0(dirpath, "Inputs/spd_ras_10.1.txt"))
values(SpDist)[values(SpDist) < 1] <- NA

#check that spdist points are in suitable habitat
sp_df <- terra::as.data.frame(SpDist, xy = TRUE, na.rm = TRUE) 
e <- extract(Swedenmap, sp_df[,1:2])

# plot land cover map and highlight cells with initial species distribution - option 1:
plot(Swedenmap, col=brewer.pal(n = 6, name = "Spectral"), axes=F)
plot(as.polygons(SpDist, dissolve=F), add=T)

# plot land cover map and highlight cells with initial species distribution - option 2 with categorical legend:
Swedenmap.f <- as.factor(Swedenmap)
# add the land cover classes to the raster attribute table (RAT)
rat <- levels(Swedenmap.f)[[1]][-2]
rat[["landcover"]] <- c("Forest, not wetland", "Wetland", "Arable land", "Forest, wetland",
                        "Open land", "Artificial", "Inland water", "Marine")
levels(Swedenmap.f) <- rat

custom.pal <- c("#1A9850","#7AC5CD", "#D9EF8B","#91CF60",  "#FEE08B", "#777777", "#98F5FF", "#00008B")
plot(Swedenmap.f, col=custom.pal, axes=F)  
plot(as.polygons(SpDist, dissolve=F), border='red', col=NA, add=T)

# carrying capacitíes and landscape parameter object
carrycap <- c(0, 100, 0, 50, 0, 0, 150, 150)
land <- ImportedLandscape(LandscapeFile = "hab_ras.txt", 
                          Resolution = 1000, 
                          Nhabitats = 8, 
                          K_or_DensDep = carrycap, 
                          SpDistFile = "spd_ras_10.1.txt", 
                          SpDistResolution = 10000)

#demography
demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.5)

#set dispersal for emigration

#exploring a dispersal kernel with a long tail
# maybe want a double kernel
plotProbs(DispersalKernel(Distances = matrix(c(1500, 10000,.9),nrow=1), DoubleKernel=T))

disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.8), #, UseFullKern=T, only if densdep=F!
                   Transfer   = DispersalKernel(Distances = matrix(c(1500, 10000,.9),nrow=1), DoubleKernel=T),
                   Settlement = Settlement() )

#set initialisation
init<- Initialise(InitType= 1, #species distribution
                  SpType =0, # all cells in species distribution
                  InitDens=1 #at half K
                  )

#set up sim
sim_1 <- Simulation(Simulation = 1, 
                    Replicates = 2, 
                    Years = 20,
                    OutIntPop = 2,
                    OutIntOcc = 2,
                    OutIntRange = 2)

#run sim
#s <- RSsim(land = land, demog = demo, dispersal = disp, simul = sim_0, init = init)
s <- RSsim(land = land, demog = demo, dispersal = disp, simul = sim_1, init = init)
s
validateRSparams(s)
RunRS(s, dirpath)

#plot results
# read 'range' output into a data frame
range_df <- readRange(s, dirpath)

# plot trajectories of all individual runs and overlay with mean:
par(mfrow=c(1,2))
plotAbundance(range_df)
plotOccupancy(range_df)


# read population output file into a data frame
pop_df <- readPop(s, dirpath)

# Not all years have the same number of cells, since only cells that had ever established a population are recorded. 
# For later stacking, we need a common extent. This is a quick & dirty solution:
#ext <- terra::ext(c(min(pop_df$x)-500,max(pop_df$x)+500,min(pop_df$y)-500,max(pop_df$y)+500))

# Make stack of different raster layers for each year and for only one repetition (Rep==0):
pop_wide_rep0 <- reshape(subset(pop_df,Rep==0)[,c('Year','x','y','NInd')], timevar='Year', v.names=c('NInd'), idvar=c('x','y'), direction='wide')
##because for some reason RangeShiftR sets xllcorner and yllcorner coordinates to 0,0
pop_wide_rep0$x <- pop_wide_rep0$x+ext(Swedenmap)[1]
pop_wide_rep0$y <- pop_wide_rep0$y+ext(Swedenmap)[3]

r_years_rep0 <- terra::rast(pop_wide_rep0, extent=ext(Swedenmap), type="xyz")
plot(r_years_rep0[["NInd.20"]])

# Overlay with Sweden mask
#r_years_rep0 <- terra::extend(r_years_rep0, Swedenmap, "out")
values(r_years_rep0)[is.na(values(r_years_rep0))] <- 0
r_years_rep0 <- terra::mask(r_years_rep0, Swedenmap)

# Map abundance
par(mfrow=c(1,2))
plot(r_years_rep0[['NInd.0']], col=c('grey',rev(inferno(150))), axes=F)
plot(as.polygons(SpDist, dissolve=F), border='red', col=NA, add=T)
plot(r_years_rep0[['NInd.20']], col=c('grey',rev(inferno(150))), axes=F)
plot(as.polygons(SpDist, dissolve=F), border='red', col=NA, add=T)





#creating initial distribution
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
presence_res_5<- aggregate(presence_ras, 5)

#to use as a reference for creating vector data with WGS84 
#since pres data was in latlon
swefin_2013<- rast("data/sdm_swefin/ensemble_swefin_2013.img")

#i have presence in lat lon
pres_2013_vect<- vect(pres_2013[,4:5], geom=c("x", "y"), crs=crs(swefin_2013))
pres_2013_proj<- terra::project(pres_2013_vect,"EPSG:3006" )
#extract cell number
extract_2013<- terra::extract(presence_res_5,pres_2013_proj, cells=T)
extract_2013<- extract_2013[complete.cases(extract_2013),] #2498
extract_2013_dup<- extract_2013[!duplicated(extract_2013[c('cell')]), ] #781

#assign cells 1 and 0 for presence
presence_res_5[!(is.na(presence_res_5))]<- 0
presence_res_5[extract_2013$cell]=1
plot(presence_res_5)
#for some reason there are a few NaNs in here
presence_res_5<- classify(presence_res_5, cbind(NaN, NA))
#save files
writeRaster(presence_res_5, "data/RS_swefin/initial_dist_5km.tif", overwrite=T)
writeRaster(presence_res_5, "data/RS_swefin/initial_dist_5km.asc", NAflag=-99, overwrite=TRUE)


SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
swfin<- terra::union(SWE, FIN)
swfin_proj<- terra::project(swfin, "EPSG:3006")
plot(presence_res_5, col=c("white", "red"))
plot(swfin_proj, add=T)

#try adding original res with aggregated
#disaggregate again to get back to 1km so we can fill in the gaps
presence_res_1<- disagg(presence_res_5, 5)
presence_ras[!(is.na(presence_ras))]<- 0
presence_ras<- mosaic(presence_ras, presence_res_1, fun="max")
plot(presence_ras)
#to make sure they match the habitat files exactly
swefin_hab_2013<-rast("data/RS_swefin/hab_type_2013_1km.tif")
presence_ras_crop<- crop(presence_ras, swefin_hab_2013)
#for some reason there are a few NaNs in here
presence_ras_crop<- classify(presence_ras_crop, cbind(NaN, NA))
writeRaster(presence_ras_crop, "data/RS_swefin/initial_dist_1km.tif", overwrite=T)
writeRaster(presence_ras_crop, "data/RS_swefin/initial_dist_1km.asc", NAflag=-99, overwrite=TRUE)


##########################
#get it for Sweden only
swe_en_proj<- rast("data/RS_swe/swe_en_EPSG3006_2013_1km.tif")
swe_presence_ras<- swe_en_proj
#swe_presence_res_5<- aggregate(swe_presence_ras, 5)

#extract cell number
# extract_2013<- terra::extract(swe_presence_res_5,pres_2013_proj, cells=T)
extract_2013<- terra::extract(swe_presence_ras,pres_2013_proj, cells=T)
extract_2013<- extract_2013[complete.cases(extract_2013),] #2145
extract_2013_dup<- extract_2013[!duplicated(extract_2013[c('cell')]), ] #621

#assign cells 1 and 0 for presence
swe_presence_res_5[!(is.na(swe_presence_res_5))]<- 0
swe_presence_res_5[extract_2013$cell]=1
plot(swe_presence_res_5)
#for some reason there are a few NaNs in here
swe_presence_res_5<- classify(swe_presence_res_5, cbind(NaN, NA))
#save files
writeRaster(swe_presence_res_5, "data/RS_swe/initial_dist_5km.tif", overwrite=T)
writeRaster(swe_presence_res_5, "data/RS_swe/initial_dist_5km.asc", NAflag=-99, overwrite=TRUE)

#plot
SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
swe_out_proj<- terra::project(SWE, "EPSG:3006")
plot(swe_presence_res_5, col=c("white", "red"))
plot(swe_out_proj, add=T)

#try adding original res with aggregated
#disaggregate again to get back to 1km so we can fill in the gaps
swe_presence_res_1<- disagg(swe_presence_res_5, 5)
swe_presence_ras[!(is.na(swe_presence_ras))]<- 0
swe_presence_ras<- mosaic(swe_presence_ras, swe_presence_res_1, fun="max") #15525
plot(swe_presence_ras)

#to make sure they match the habitat files exactly
swe_hab_2013<-rast("data/RS_swe/hab_type_2013_1km.tif")
swe_presence_ras<- crop(swe_presence_ras, swe_hab_2013)
#for some reason there are a few NaNs in here
swe_presence_ras<- classify(swe_presence_ras, cbind(NaN, NA))
writeRaster(swe_presence_ras, "data/RS_swe/initial_dist_1km.tif", overwrite=T)
#write to txt file
writeRaster(swe_presence_ras, "data/RS_swe/initial_dist_1km.asc", NAflag=-99, overwrite=TRUE)

#try it with 0s instead of NAs
swe_presence_ras[is.na(swe_presence_ras)]<- 0
plot(swe_presence_ras)
writeRaster(swe_presence_ras, "data/RS_swe/initial_dist_1km_0.asc", NAflag=-99, overwrite=TRUE)
