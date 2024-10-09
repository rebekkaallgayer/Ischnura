#code copied from Damaris Zurell's tutorial 
#https://damariszurell.github.io/EEC-MGC/b3_SDM_intro.html#1_Introduction

library(terra)
library(corrplot)
library(mecofun)


setwd("C:/Users/Rey/Documents/Ischnura/SDM_tutorials/Damaris_new")
bg <- terra::rast('Prac3_data/Prac3_UK_mask.tif')
sp_dat <- read.table('Prac3_data/Prac3_RingOuzel.txt',header=T)
summary(sp_dat)

# Plot GB land mass
terra::plot(bg,col='grey',axes=F,legend=F)

# Plot presences in red and absences in black
plot(extend(terra::rast(sp_dat[,1:3], crs=crs(bg), type='xyz'), bg), col=c('black','red'), legend=F,add=T)

# We first fit a GLM for the bio1 variable assuming a linear relationship:
m1 <- glm(Turdus_torquatus ~ bio1, family="binomial", data= sp_dat)

# We can get a summary of the model:
summary(m1)

# Fit a quadratic relationship with bio1:
m1_q <- glm(Turdus_torquatus ~ bio1 + I(bio1^2), family="binomial", data= sp_dat)
summary(m1_q)

# Or use the poly() function:
summary( glm(Turdus_torquatus ~ poly(bio1,2) , family="binomial", data= sp_dat) )

# Fit two linear variables:
summary( glm(Turdus_torquatus ~ bio1 + bio8, family="binomial", data= sp_dat) )

# Fit three linear variables:
summary( glm(Turdus_torquatus ~ bio1 + bio8 + bio17, family="binomial", data= sp_dat) )

# Fit three linear variables with up to three-way interactions
summary( glm(Turdus_torquatus ~ bio1 * bio8 * bio17, family="binomial", data= sp_dat) )

# Fit three linear variables with up to two-way interactions
summary( glm(Turdus_torquatus ~ bio1 + bio8 + bio17 + 
               bio1:bio8 + bio1:bio17 + bio8:bio17, 
             family="binomial", data= sp_dat) )

# We first estimate a correlation matrix from the predictors. 
# We use Spearman rank correlation coefficient, as we do not know 
# whether all variables are normally distributed.
cor_mat <- cor(sp_dat[,-c(1:3)], method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Run select07()
var_sel <- select07(X=sp_dat[,-c(1:3)], 
                    y=sp_dat$Turdus_torquatus, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel

# How many presence points do we have?
sum(sp_dat$Turdus_torquatus)
#According to the rule of thumb, we can only include 4 parameters and want 
#to include quadratic terms as well. Thus, for this example, we should only 
#use two predictors. In our case, the best two predictors in terms of 
#univariate AIC were bio1 and bio8.

# Fit the full model:
m_full <- glm( Turdus_torquatus ~ bio1 + I(bio1^2) + bio8 + I(bio8^2), 
               family='binomial', data=sp_dat)

# Inspect the model:
summary(m_full)

# Explained deviance:
expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m_full$fitted)

m_step <- step(m_full) 
# Inspect the model:
summary(m_step)

# Explained deviance:
expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m_step$fitted)
#The explained deviance of the final model is a tiny bit 
#lower than for the full model, but the final model is more parsimonious.

# Exercise:
#   
#   Load the second species data set for the Yellowhammer.
#bg <- terra::rast('Prac3_data/Prac3_UK_mask.tif')
yh_dat <- read.table('Prac3_data/Prac3_YellowHammer.txt',header=T)
summary(yh_dat)

# Plot the presences and absences of the Yellowhammer.
# Plot GB land mass
terra::plot(bg,col='grey',axes=F,legend=F)

# Plot presences in red and absences in black
plot(extend(terra::rast(yh_dat[,1:3], crs=crs(bg), type='xyz'), bg), col=c('black','red'), legend=F,add=T)

# Check for collinearity and reduce the data set to only weakly correlated variables.
cor_mat <- cor(yh_dat[,-c(1:3)], method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)
# Run select07()
var_sel <- select07(X=yh_dat[,-c(1:3)], 
                    y=yh_dat$Emberiza_citrinella, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel
sum(yh_dat$Emberiza_citrinella)

# Define a full model with all weakly correlated variables including their linear and quadratic terms.
# Fit the full model:
m_full_yh <- glm( Emberiza_citrinella ~ bio12 + I(bio12^2) + bio9 + I(bio9^2)
                  + bio1 + I(bio1^2)+ bio8 + I(bio8^2)+ bio3 + I(bio3^2), 
               family='binomial', data=yh_dat)

# Inspect the model:
summary(m_full_yh)
# Run the full model.
# Simplify the model using stepwise variable selection step()
# Explained deviance:
expl_deviance(obs = yh_dat$Emberiza_citrinella,
              pred = m_full_yh$fitted)

m_step_yh <- step(m_full_yh) 
# Inspect the model:
summary(m_step_yh)
# Compare the full model and the reduced model in terms of AIC and explained deviance.
# Explained deviance:
expl_deviance(obs = yh_dat$Emberiza_citrinella,
              pred = m_step_yh$fitted)


#getting species data
library(rgbif)

# Check out the number of occurrences found in GBIF:
occ_count()
# number of observations:
occ_count(basisOfRecord='OBSERVATION')

# number of observations reported for Germany:
occ_count(
  country = rgbif::rgb_country_codes("Germany")$iso2, 
  basisOfRecord = 'OBSERVATION'
)
# Check for synonyms
name_suggest(q='Sorex alpinus', rank='species')
# Check number of records - here filtered to those with coordinate information
occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, limit = 10)
occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, basisOfRecord='HUMAN_OBSERVATION', limit = 10)

gbif_shrew <- occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, basisOfRecord='HUMAN_OBSERVATION', limit = 600) 

# We are just interested in the data frame containing the records
gbif_shrew <- gbif_shrew$data

library(maps)
maps::map('world',xlim=c(5,30), ylim=c(40,55))
points(gbif_shrew$decimalLongitude, gbif_shrew$decimalLatitude, col='red',  pch=19)

library(CoordinateCleaner)

# We use only those data entries with coordinate information - Note that you don't need this if you have used the hasCoordinate=T in the occ_search() function:
#gbif_shrew <- subset(gbif_shrew, !is.na(decimalLatitude))

# We now clean the coordinates and check for outliers - see ?clean_coordinates for more options
gbif_shrew_cleaned_coord <- clean_coordinates(gbif_shrew, lon="decimalLongitude", lat="decimalLatitude", countries="countryCode", tests=c("centroids", "outliers", "duplicates", "institutions"), inst_rad = 1000)

# Plot world map
maps::map('world',xlim=c(5,30), ylim=c(40,55))
# Plot all gbif points downloaded
points(gbif_shrew$decimalLongitude, gbif_shrew$decimalLatitude, col='red',  pch=19)
# Plot all remaining points after cleaning
points(gbif_shrew$decimalLongitude[gbif_shrew_cleaned_coord$.summary], gbif_shrew$decimalLatitude[gbif_shrew_cleaned_coord$.summary], col='blue',  pch=18)

# Store the cleaned point locations in a new object
gbif_shrew_cleaned <- gbif_shrew[gbif_shrew_cleaned_coord$.summary,]
save(gbif_shrew_cleaned,file='data/gbif_shrew_cleaned.RData')

###################
##environmental data
##################

# Download global bioclimatic data from worldclim (you may have to set argument 'download=T' for first download, if 'download=F' it will attempt to read from file):
clim <- geodata::worldclim_global(var = 'bio', res = 10, download = T, path = 'data')

# Now, let's look at the data:
clim
terr_cls<- terrain.colors(100, rev=T)
plot(clim,col=terr_cls )

clim_ag<- terra::aggregate(clim[[1]], fact=6, fun="mean")
plot(clim_ag, col=terr_cls)

# Download future climate scenario from 'ACCESS-ESM1-5' climate model.
# Please note that you have to set download=T if you haven't downloaded the data before:
clim_fut <- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2041-2060', var='bioc', download=F, res=10, path='data')

# Inspect the SpatRaster object:
clim_fut
plot(clim_fut, col=terr_cls)

# In this case, let's keep the names of the future climate layers
names(clim) <- names(clim_fut)

# Download fractional tree cover at 30-sec resolution:
# Please note that you have to set download=T if you haven't downloaded the data before:
trees_30sec <- geodata::landcover(var='trees', path='data', download=F)

# map the tree cover
plot(trees_30sec)

# Aggregate tree cover to 10-min spatial resolution
trees_10min <- terra::aggregate(trees_30sec, fact=20, fun='mean')

# Map the 10-min tree cover
plot(trees_10min)

# This produces an error that spatial extents do not match:
env_cur <- c(clim, trees_10min)
## Error: [rast] extents do not match

# Which SpatRaster object has the larger extent?
terra::ext(clim)
## SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)

terra::ext(trees_10min)
## SpatExtent : -180, 179.99999999999, -59.999999999996, 84 (xmin, xmax, ymin, ymax)

# As the climate data have the larger extent, we now have to "extend" our land cover extent
terra::extend(trees_10min, clim)
# Produce the multi-layer environmental data object with matching extents:
env_cur <- c(clim, terra::extend(trees_10min, clim))

# I suggest to keep the following columns for now:
gbif_shrew2 <- gbif_shrew_cleaned[,c("key", "scientificName", "decimalLatitude", 
                                     "decimalLongitude", "basisOfRecord", "speciesKey", 
                                     "species", "year")]

# We can extract the environmental data for the GBIF coordinates.
# Coordinates are always provided as x/y format, in our case lon/lat.
# We also extract the cellnumbers as this allows checking for duplicates later.
head(terra::extract(x = env_cur, 
                    y = data.frame(gbif_shrew2[,c('decimalLongitude','decimalLatitude')]), cells=T ))

# Finally, we put species and environmental data into the same data frame:
gbif_shrew2 <- cbind(gbif_shrew2, terra::extract(x = env_cur, y = data.frame(gbif_shrew2[,c('decimalLongitude','decimalLatitude')]), cells=T ))
sum(duplicated(gbif_shrew2$cell))

# Only retain non-duplicated cells (will not work in this example as we don't have duplicates):
gbif_shrew_env <- gbif_shrew2[!duplicated(gbif_shrew2$cell),]

save(gbif_shrew2, gbif_shrew_cleaned,file='data/gbif_shrew_cleaned.RData')