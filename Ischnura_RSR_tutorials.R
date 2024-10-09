##################setting up RangeShiftR#################
setwd("Ischnura/")
# install.packages("Rcpp")
# install.packages("devtools")
# install.packages("Rdpack")

library(Rcpp)
library(devtools)
library(Rdpack)
Rcpp::evalCpp("2+2") #should give 4, if not, toolchain not complete
devtools::install_github("https://github.com/RangeShifter/RangeShiftR-pkg", 
                         ref = "main", subdir="RangeShiftR")
library(RangeShiftR)

##running a tutorial

#Create a parameter master object with all the default settings and store it:
s <- RSsim()
#store path to working directory
dirpath = "RSR_Tutorial/"
#Create the RS folder structure, if it doesn’t yet exist:
dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)

RunRS(s,dirpath)

s

#This module is used to set general simulation
#parameters (e.g. simulation ID, number of 
#replicates, and number of years to simulate) 
#and to control output types (plus some more specific settings). 
#For this overview, we will stick to the defaults:
sim <- Simulation(Simulation = 2,
                  Years = 50,
                  Replicates = 2,
                  OutIntPop = 50)

#RangeShiftR can either import a map from an 
#ASCII raster file in the ‘Inputs’ folder or generate
#a random map to use in the simulation.
#For each option, there is a corresponding 
#function to create a Landscape parameter object
land <- ImportedLandscape() 
land <- ArtificialLandscape()

land <- ArtificialLandscape(Resolution = 10,  # in meters
                            K_or_DensDep = 1500,  # ~ 15 inds/cell
                            propSuit = 0.2,
                            dimX = 129, dimY = 257, 
                            fractal = T, hurst = 0.3,
                            continuous = F)

##set demographic parameters
#unstructured model/non-overlapping gens
demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.45)

#dispersal
disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.2), 
                   Transfer   = DispersalKernel(Distances = 50),
                   Settlement = Settlement() )

plotProbs(DispersalKernel(Distances = matrix(c(0,1,2,70,50,30),nrow = 3), StageDep = T))

#settlement
disp <-  disp + Settlement(SexDep = T, 
                           Settle = matrix(c(0,1,1,0), nrow = 2))

disp2 <-  disp + Emigration(IndVar = TRUE,
                            EmigProb = matrix(c(.7,.1), nrow = 1),
                            TraitScaleFactor = c(.1))

#initialisation
init <- Initialise(FreeType = 0, 
                   NrCells = 2250,
                   InitDens = 2, 
                   IndsHaCell = 3, 
                   PropStages = c(0,0.7,0.3))
init


##tutorial 1
carrycap <- c(0, 0, 0, 5, 0, 0)
land <- ImportedLandscape(LandscapeFile = "UKmap_1km.txt", 
                          Resolution = 1000, 
                          Nhabitats = 6, 
                          K_or_DensDep = carrycap, 
                          SpDistFile = "Species_Distribution_10km.txt", 
                          SpDistResolution = 10000)
demo <- Demography(Rmax = 1.5)
disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.1), 
                   Transfer = DispersalKernel(Distances = 2000), 
                   Settlement = Settlement() )
init <- Initialise(InitType = 1, # = initialisation from a loaded species distribution map
                   SpType = 0,   # = all suitable cells within all distribution presence cells
                   InitDens = 0) # = at carrying capacity
sim_0 <- Simulation(Simulation = 0, 
                    Replicates = 2, 
                    Years = 100,
                    OutIntPop = 5,
                    OutIntOcc = 5,
                    OutIntRange = 5)
s <- RSsim(land = land, demog = demo, dispersal = disp, simul = sim_0, init = init)
RunRS(s, dirpath)

# read population output file into a data frame
pop_df <- readPop(s, dirpath)

# Not all years have the same number of cells, since only cells that had ever established a population are recorded. 
# For later stacking, we need a common extent. This is a quick & dirty solution:
ext <- terra::ext(c(min(pop_df$x)-500,max(pop_df$x)+500,min(pop_df$y)-500,max(pop_df$y)+500))

# Make stack of different raster layers for each year and for only one repetition (Rep==0):
pop_wide_rep0 <- reshape(subset(pop_df,Rep==0)[,c('Year','x','y','NInd')], timevar='Year', v.names=c('NInd'), idvar=c('x','y'), direction='wide')
r_years_rep0 <- terra::rast(pop_wide_rep0, type="xyz")

# Overlay with UK mask
UKmap <- terra::rast(paste0(dirpath, "Inputs/UKmap_1km.txt"))
SpDist <- terra::rast(paste0(dirpath, "Inputs/Species_Distribution_10km.txt"))
values(SpDist)[values(SpDist) < 1] <- NA

r_years_rep0 <- terra::extend(r_years_rep0, UKmap)
values(r_years_rep0)[is.na(values(r_years_rep0))] <- 0
r_years_rep0 <- terra::mask(r_years_rep0, UKmap)

# Map abundance
plot(r_years_rep0[['NInd.90']], col=c('grey',rev(inferno(150))), axes=F)
plot(as.polygons(SpDist, dissolve=F), border='red', col=NA, add=T)

########################################
##tutorial2
