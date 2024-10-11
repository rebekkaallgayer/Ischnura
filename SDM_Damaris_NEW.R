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

##practical 4: assessment and prediction
m1 <- step(glm( Turdus_torquatus ~ bio1 + I(bio1^2) + bio8 + I(bio8^2), family='binomial', data=sp_dat))
summary(m1)

#visualising response curves
# If we do not provide "newdata", then predict() should simply return the fitted values: 
head(predict(m1, type='response'))
head(m1$fitted)
# We want to make predictions for all combinations of the two predictor variables
# and along their entire environmental gradients:
xyz <- expand.grid(
  # We produce a sequence of environmental values within the predictor ranges:
  bio1 = seq(min(sp_dat$bio1),max(sp_dat$bio1),length=50),
  bio8 = seq(min(sp_dat$bio8),max(sp_dat$bio8),length=50)
)

# Now we can make predictions to this new data frame
xyz$z <- predict(m1, newdata=xyz, type='response')
summary(xyz)
# As result, we have a 3D data structure and want to visualise this.
# Here, I first set a color palette
library(RColorBrewer)
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# Finally, we plot the response surface using the wireframe function from the lattice package
library(lattice)
wireframe(z ~ bio1 + bio8, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1))
# We can also rotate the axes to better see the surface
wireframe(z ~ bio1 + bio8, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          screen=list(z = -160, x = -70, y = 3))

library(mecofun)

# Names of our variables:
my_preds <- c('bio1', 'bio8')

# We want two panels next to each other:
par(mfrow=c(1,2))

# Plot the partial responses
partial_response(m1, predictors = sp_dat[,my_preds], ylab='Occurrence probability')

expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m1$fitted)
preds_cv <- crossvalSDM(m1, traindat = sp_dat, colname_species = 'Turdus_torquatus', colname_pred = my_preds)
plot(m1$fitted.values, preds_cv, xlab='Fitted values', ylab='Predicted values from CV')
abline(0,1,col='red',lwd=2)

library(PresenceAbsence)
#Threshold-dependent performance measures
# We first prepare our data:
# Prepare cross-validated predictions:
thresh_dat <- data.frame(
  ID = seq_len(nrow(sp_dat)), 
  obs = sp_dat$Turdus_torquatus,
  pred = preds_cv)

# Then, we find the optimal thresholds:     
(thresh_cv <- PresenceAbsence::optimal.thresholds(DATA= thresh_dat))
(cmx_maxSSS <- PresenceAbsence::cmx(DATA= thresh_dat, threshold=thresh_cv[3,2]))

# From such a confusion matrix, we can calculate different evaluation criteria. For example,
# - the proportion of correctly classified test observations pcc
# - the proportion of correctly classified presences, also called sensitivity or true positive rate
# - the proportion of correctly classified absences, also called specificity or true negative rate

# Proportion of correctly classified observations
PresenceAbsence::pcc(cmx_maxSSS, st.dev=F)
# Sensitivity = true positive rate
PresenceAbsence::sensitivity(cmx_maxSSS, st.dev=F)
# Specificity = true negative rate
PresenceAbsence::specificity(cmx_maxSSS, st.dev=F)
# Kappa
PresenceAbsence::Kappa(cmx_maxSSS, st.dev=F)
# True skill statistic
TSS(cmx_maxSSS) 

#Threshold-independent performance measures
library(AUC)

# Let's have a look a the ROC curve:
roc_cv <- roc(preds_cv, as.factor(sp_dat$Turdus_torquatus))
plot(roc_cv, col = "grey70", lwd = 2)
# Compute the AUC:
AUC::auc(roc_cv)
evalSDM(sp_dat$Turdus_torquatus, preds_cv)

#Spatio-temporal predictions
library(geodata)

# Please note that you have to set download=T if you haven't downloaded the data before:
bio_curr <- geodata::worldclim_country('GB', var = 'bio', res = 5, download = F, path = 'data')

# Please note that you have to set download=T if you haven't downloaded the data before:
bio_fut <- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2041-2060', var='bioc', download=F, res=5, path='data')
# the approx. spatial extent of UK in lon/lat coordinates
extent_uk <- c(-12, 3, 48, 62)

# Crop and reproject current climate
bio_curr <- terra::crop(bio_curr, extent_uk)
bio_curr <- terra::project(bio_curr, bg)
bio_curr <- terra::mask(bio_curr, bg)

# Crop and reproject future climate
bio_fut <- terra::crop(bio_fut, extent_uk)
bio_fut <- terra::project(bio_fut, bg)
bio_fut <- terra::mask(bio_fut, bg)
# Change names of climate layers
names(bio_curr) <- names(bio_fut) <- names(sp_dat)[-c(1:3)]

#Make predictions to the environmental layers
# Prepare data frames
bio_curr_df <- data.frame(crds(bio_curr),as.points(bio_curr))
bio_fut_df <- data.frame(crds(bio_fut),as.points(bio_fut))

# Make continuous predictions:
bio_curr_df$pred_glm <- predict(m1, newdata= bio_curr_df, type="response")
bio_fut_df$pred_glm <- predict(m1, newdata= bio_fut_df, type="response")

par(mfrow=c(1,2))

# Make raster of predictions to current environment:
r_pred_curr <- terra::rast(bio_curr_df[,c('x','y','pred_glm')], type='xyz', crs=crs(bg))
plot(r_pred_curr, axes=F, main='Occ. prob. - today', col=terr_cls)

# Make raster stack of predictions to future environment:
r_pred_fut <- terra::rast(bio_fut_df[,c('x','y','pred_glm')], type='xyz', crs=crs(bg))
plot(r_pred_fut, axes=F, main='Occ. prob. - 2050', col=terr_cls)

# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm >= thresh_cv[3,2], 1, 0)
bio_fut_df$bin_glm <- ifelse(bio_fut_df$pred_glm >= thresh_cv[3,2], 1, 0)

# Make raster stack of predictions to current environment:
r_pred_curr <- terra::rast(bio_curr_df[,c('x','y','pred_glm','bin_glm')], type='xyz', crs=crs(bg))
plot(r_pred_curr, axes=F, col=terr_cls)
# Make raster stack of predictions to future environment:
r_pred_fut <- terra::rast(bio_fut_df[,c('x','y','pred_glm','bin_glm')], type='xyz', crs=crs(bg))
plot(r_pred_fut, axes=F, col=terr_cls)

#Assessing novel environments
library(predicts)
#Novel environments can be assessed in different ways. 
#MESS (Multivariate environmental similarity surface) maps assess for each 
#environmental variables separately whether the projection data contain novel 
#conditions beyond the sampled range.

# MESS maps from the predicts package:
r_mess <- mess(bio_fut[[my_preds]], sp_dat[,my_preds])
plot(r_mess, axes=F, col=terr_cls)
# Negative values indicate dissimilar=novel environments:
r_mess_mask <- r_mess < 0
plot(r_mess_mask, axes=F)
#we can already see that novel environments should not be 
#any issue for the ring ouzel as novel environments could mainly arise in the
#South while the Ring Ouzel is a northern distributed species
# Predictions to analogous climates:
r_analog_fut <- terra::extend(r_pred_fut, bg)
values(r_analog_fut)[values(r_mess)<0] <- NA
plot(r_analog_fut, axes=F, col=terr_cls)

# Predictions to novel climates:
r_novel_fut <- terra::extend(r_pred_fut, bg)
values(r_novel_fut)[values(r_mess)>=0] <- NA
plot(r_novel_fut, axes=F, col=terr_cls)

###practical 5: Pseudo-absence and background data
# Our study region and species data
region <- terra::rast('Prac5_data/Prac5_Europe5min.grd') # mask of study region
sp <- read.table('Prac5_data/Prac5_presences.txt', header=T) # species presences

# Plot the map and data
plot(region,col='grey',legend=F)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(sp[sp$sp=='Populus_spp',1:2],pch='+',cex=0.3,col='grey20')

#Random selection of points within study area but excluding the presence location
# Randomly select background points from study region
# The argument na.rm=T ensures that we only sample points within the masked regions (not in the ocean)
# The argument as.points=T indicates that a SpatVector of coordinates should be returned
bg_rand <- terra::spatSample(region, 500, "random", na.rm=T, as.points=TRUE)

# Plot the map and data
plot(region,col='grey',legend=F)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(bg_rand,pch=19,cex=0.3)
# Make a new regional mask that contains NAs in presence locations:
sp_cells <- terra::extract(region, sp[sp$sp=='Populus_imagines',1:2], cells=T)$cell
region_exclp <- region
values(region_exclp)[sp_cells] <- NA

# Randomly select background data but excluding presence locations
bg_rand_exclp <- terra::spatSample(region_exclp, 500, "random", na.rm=T, as.points=TRUE)

# Plot the map and data
plot(region,col='grey',legend=F)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(bg_rand_exclp,pch=19,cex=0.3)

#Also, we can define an extent from where random points should be drawn.

# Define extent object:
e <- ext(8,24,46,57)

# Randomly select background data within a restricted extent and excluding presence locations:
bg_rand_exclp_e <- terra::spatSample(region_exclp, 500, "random", na.rm=T, as.points=TRUE, ext=e)

# Plot the map and data
plot(region,col='grey',legend=F)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(bg_rand_exclp_e,pch=19,cex=0.3)
lines(e, col='red')

#Last, we could also restrict the random samples to within a certain buffer distance.
# Create SpatVector object of known occurrences:
pop_imag <- terra::vect( as.matrix( sp[sp$sp=='Populus_imagines',1:2]) , crs=crs(region))

# Then, place a buffer of 200 km radius around our presence points
v_buf <- terra::buffer(pop_imag, width=200000)

# Set all raster cells outside the buffer to NA
region_buf <- terra::mask(region, v_buf)

# Randomly select background data within the buffer
bg_rand_buf <- terra::spatSample(region_buf, 500, "random", na.rm=T, as.points=TRUE)

# Plot the map and data
plot(region,col='grey',legend=F)
plot(region_buf, legend=F, add=T)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(bg_rand_buf,pch=19,cex=0.3)

#Random selection of points outside of study area
#Barbet-Massin et al. (2012) also suggested a method to sample pseudo-absences 
#only beyond a minimum distance from the presence points. This is more of a 
#macroecological approach, suitable for characterising the climatic limits of species. 
#We can also use the buffering approach from above to achieve this.


# Place a buffer of 200 km radius around our presence points
v_buf <- terra::buffer(pop_imag, width=200000)

# Set all raster cells inside the buffer to NA using the argument inverse=T
region_outbuf <- terra::mask(region, v_buf, inverse=T)
region_buf <- terra::mask(region, v_buf)  # we use this only for illustrative purpose below

# Randomly select background data outside the buffer
bg_rand_outbuf <- terra::spatSample(region_outbuf, 500, "random", na.rm=T, as.points=TRUE)

# Plot the map and data
plot(region,col='grey',legend=F)
plot(region_buf,legend=F, add=T)
points(sp[sp$sp=='Populus_imagines',1:2],pch='+',col='red')
points(bg_rand_outbuf,pch=19,cex=0.3)