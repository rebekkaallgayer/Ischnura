##This script is for running RSR simulations for Ischnura elegans

#Dr Rebekka Allgayer August 2025

##################setting up RangeShiftR#################
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

#set working directory
setwd("C:/Users/Rey/Documents/Ischnura/Handover")

#store path to working directory
dirpath = "Simulation1/"

######for real landscape#######
#set land parameters with imported landscape, Sweden 2013
Swedenmap <- terra::rast(paste0(dirpath, "Inputs/hab_type_2013_1km.asc"))
SpDist <- terra::rast(paste0(dirpath, "Inputs/initial_dist_1km.asc"))
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
rat[["landcover"]]<- c("Unsuitable", "Low quality", "Medium quality", "High quality")
levels(Swedenmap.f) <- rat

custom.pal <- c("#777777",  "#FEE08B","#D9EF8B","#91CF60")
custom.pal <- c("#777777",  "#FEE08B","#7AC5CD","#98F5FF")
plot(Swedenmap.f, col=custom.pal, axes=F)  
plot(as.polygons(SpDist, dissolve=F), border='red', col=NA, add=T)

# carrying capacitíes and landscape parameter object
carrycap <- c(0, 1000, 2500, 4500)
land <- ImportedLandscape(LandscapeFile = "hab_type_2013_1km.asc", 
                          Resolution = 1000, 
                          Nhabitats = 4, 
                          K_or_DensDep = carrycap, 
                          SpDistFile = "initial_dist_1km.asc", 
                          SpDistResolution = 1000)

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
