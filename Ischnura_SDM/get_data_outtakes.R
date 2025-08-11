swe_c<- swe_comb
swe_c[swe_c>0]<-8

#create habitat type copy
nor_c<- nor_comb
nor_c[nor_c>0]<- 8
plot(nor_c)

#to update land layer, create just habitat type 8
fin_c<- fin_comb
fin_c[fin_c>0]<- 8
plot(fin_c)


fenno_w[fenno_w!=8]<- NA
plot(fenno_w)

Fenno<- vect("data/Fenno.gpkg")
fenno_blank<- crop(fenno_land, Fenno)
fenno_blank[!is.na(fenno_blank)]<- 0
fenno_blank <- resample(fenno_blank, fenno_land)

fenno_filled<- cover(fenno_w, fenno_blank)
plot(fenno_filled)


Fenno<- vect("data/Fenno.gpkg")
plot(Fenno)

#elevation
fenno_elev<- rast("data/elev_Fenno_ll.tif")
names(fenno_elev)<- "elevation"

#land cover types
# fenno_land<- rast("data/land_cover_Fenno_ll.tif")
fenno_land<- rast("data/final_fenno/fenno_land_wupdate.tif")
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
fenno_water<- rast("data/final_fenno/perc_cover_freshwater_Fenno.tif")

#specific year for env variables
years<- sort(unique(isch_2000s$year))
#years<- seq(2000,2018,1) #i only have Chelsa data until 2018!
#split per year
#extract() function taking year-specific env data for each set of presence points
#combine into one big dataframe to pass to sdmData()

isch_2000s<- isch_dat[(isch_dat$year>=2000)&(!is.na(isch_dat$year)),c(5,4, 10,12)]
#just lon, lat
xy_crop<- isch_2000s[,1:2]

#need to read in the presence data
#xy_2013$species<-1
xy_2013_vect<- vect(xy_2013[,1:2], geom=c("longitude", "latitude"))
extract_2013<- terra::extract(chelsa_yr_crop,xy_2013_vect)

#calculates VIF for all variables, excludes the one with the highest VIF 
#(if it is greater than the threshold), repeat the procedure until no variables 
#with a VIF greater than th remains.
vif_2013_extract<- vif(extract_2013[,-(1)])
vif_2013_extract

v_2013 <- vifstep(extract_2013[,-1]) #i think threshold is by default 10
v_2013

#This would have excluded : layer.5 layer.4 layer.2 layer.11 layer.12 layer.13 layer.7 layer.1 layer.16 
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1    layer.3 2.107992
# 2    layer.6 2.354404
# 3    layer.8 1.942563
# 4    layer.9 2.754825
# 5   layer.10 2.943945
# 6   layer.14 3.402239
# 7   layer.15 3.751145
# 8   layer.17 7.968019
# 9   layer.18 5.262451
# 10  layer.19 8.079809

#chelsa_sub_extract<- extract(chelsa_sub, xy_2013_vect)
#chelsa_sub_extract_vif<- vif(chelsa_sub_extract)


####old spatial thinning code
#playing with dynamicSDM spatiotemporal autocorrelation
# library(ape)
# load("data/sample_explan_data.rda")
# variablenames<-c("eight_sum_prec","year_sum_prec","grass_crop_percentage")
# 
# autocorrelation <- spatiotemp_autocorr(sample_explan_data,
#                                        varname = variablenames,
#                                        plot = TRUE,
#                                        temporal.level = c("year")) # can choose month or day too
# 
# #> `geom_smooth()` using formula = 'y ~ x'
# autocorrelation
# 
# variablenames<- colnames(env_ext)[-(1:5)]
# autocorrelation <- spatiotemp_autocorr(env_ext_dup,
#                                        varname = variablenames,
#                                        plot = F,
#                                        temporal.level = c("year"))
# 

##spatial thinning
# thin() expects that the data.frame contains a column with the species name
# env_ext$sp <- 'Ischnura_elegans'
# env_ext_dup$sp <- 'Ischnura_elegans'

# xy <- thin(occ, lat.col='y',long.col='x',spec.col='sp',
#            thin.par=10,reps=1, write.files=F,locs.thinned.list.return=T)
# 
# # Keep the coordinates with the most presence records
# xy_keep <- xy[[1]]
# xy_vect<- vect(xy_keep, geom=c("Longitude", "Latitude"), crs=crs(fenno_coast))
# plot(fenno_coast)
# points(xy_vect, cex=0.2, col="red")
# 
# sp_thinned<- occ[rownames(xy_keep),]
# 
# autocorrelation <- spatiotemp_autocorr(sp_thinned,
#                                        varname = variablenames,
#                                        plot = F,
#                                        temporal.level = c("year"))
# 
# 
# ##alternative way of testing spatial autocorrelation
# #https://github.com/jorgeassis/spatialAutocorrelation/blob/master/mainScript.R
# library(raster)
# library(ggplot2)
# library(sdmpredictors)
# library(ecodist)
# library(sp)
# library(spThin)
# 
# source("functions.R")
# 
# 
# source("https://raw.githubusercontent.com/jorgeassis/marineforestsDB/master/sourceMe.R")
# 
# # Define the distance class 
# autocorrelationClassDistance <- 5
# 
# # Define the maximum distance at which autocorrelation will be determined
# autocorrelationMaxDistance <- 500
# 
# # Define the significance level of the test
# autocorrelationSignif <- 0.05
# 
# sa_layers<- stack(env_layers)
# distanceUncorr <- data.frame(Predictor=variablenames,Distance=NA)
# for( i in 1:length(variablenames)) {
#   distanceUncorr[i,2] <- spatialAutocorrelation(occurrenceRecords=env_ext_dup[,c(4,5)],subset(sa_layers,i),
#                                                 autocorrelationClassDistance,autocorrelationMaxDistance,
#                                                 autocorrelationSignif)
# }
# 
# distanceUncorrPlot <- ggplot(distanceUncorr[sort(distanceUncorr[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
#   geom_bar( aes(x= reorder(Predictor, Distance) , y=Distance), stat="identity", fill="black", alpha=0.5) +
#   coord_flip() + theme(
#     axis.text=element_text(size=12),
#     axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
#     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
#     panel.border = element_blank(),
#     panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
#     panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
#     panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
#   ) + labs(x = "Predictor") + 
#   labs(y = "Spatial correlation (km)") + geom_hline(aes(yintercept=round(mean(distanceUncorr[,2]),digits=2)   ),color="Black", linetype="dashed", size=0.3) +
#   annotate("text", y = round(mean(distanceUncorr[,2]),digits=2) + 2 , x = 1 , label = paste0( round(mean(distanceUncorr[,2]),digits=2)," Km") , hjust = 0)
# 
# distanceUncorrPlot
# 
# meanCorrDistance <- mean(distanceUncorr[,2])
# 
# # Prune occurrence recors 
# #occurrenceRecords <- spatialThinning(occ[,c(4,5)],sa_layers, meanCorrDistance )
# env_ext_dup$sp <- 'Ischnura_elegans'
# xy <- thin(env_ext_dup, lat.col='y',long.col='x',spec.col='sp',
#            thin.par=meanCorrDistance,reps=1, write.files=F,locs.thinned.list.return=T)
# xy <- thin(env_ext_dup, lat.col='y',long.col='x',spec.col='sp',
#            thin.par=3,reps=1, write.files=F,locs.thinned.list.return=T)
# xy_keep <- xy[[1]]
# 
# # Thin the dataset - here, we first extract the cell numbers for the thinned coordinates and then use these to subset our data frame.
# cells_thinned <- terra::cellFromXY(env_layers, xy_keep)
# sp_thinned <- env_ext_dup[env_ext_dup$cell %in% cells_thinned,]

# 
# yr_vect<- vect(occurrenceRecords, geom=c("Lon", "Lat"), crs=crs(fenno_coast))
# yr_ext<- terra::extract(env_layers, yr_vect, xy=T)[,2:11]
# 
# autocorrelation <- spatiotemp_autocorr(yr_ext,
#                                        varname = variablenames,
#                                        plot = F,
#                                        temporal.level = c("year"))
# 
# 



fenno_coast1<- resample(fenno_coast, CHELSA_2013)
years=seq(2000,2018,1)
for(y in 1:length(years)){
  chelsa_y<- rast(paste("../Ischnura_SDM/data/final_fenno/",years[y], "_layers.tif", sep=""))
  
  chelsa_y[["Dist_coast"]]<- fenno_coast1
  writeRaster(chelsa_y, paste("data/CHELSA/",years[y], "_layers.tif", sep=""), overwrite=T)
}