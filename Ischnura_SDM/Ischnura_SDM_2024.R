library(data.table)
# library(randomForest)
# library(lattice)
# library(RColorBrewer)
# library(PresenceAbsence)
library(geodata)
library(terra)
library(sdm)
#library(dismo)
library(maps)
library(CoordinateCleaner)
library(rgbif)
library(corrplot)
library(mecofun)
library(usdm)

setwd("C:/Users/Rey/Documents/Ischnura/Ischnura_SDM/")

##read in occurrence data
isch_dat <- read.table('data/Ischnura_elegans_2022_10_26.csv', header=T, sep=',')

summary(isch_dat)

#filter for country, if needed
isch_sw<- isch_dat[isch_dat$country=="Sweden",]


#make shapefile of Norway, Finland, and Sweden, called "Fenno"
# SWE<- readRDS("data/gadm/gadm41_SWE_0_pk.rds")
# plot(SWE)
#  FIN<- readRDS("data/gadm/gadm41_FIN_0_pk.rds")
# # #plot(FIN)
#  NOR<- readRDS("data/gadm/gadm41_NOR_0_pk.rds")
#  Fenno1<- terra::union(SWE, FIN)
#  Fenno<- union(Fenno1, NOR)
Fenno<- vect("data/Fenno.gpkg")
plot(Fenno)

# occurrence points:
ischpts2<- read.table("data/isch_points_2000.txt", header=T, sep="\t")
ipts<- vect("data/isch_points_2000.gpkg")
plot(Fenno)
points(ischpts2, col="red", cex=0.1)
points(ipts, col="blue", cex=0.1)


#select target background points (= Zygoptera), 
#Targeted background sampling using sightings of a wider group of 
#similar organisms, to control for sampling effort (Phillips et al. 2009)
e<- ext(Fenno)[1:4]

# Check for synonyms
#name_suggest(q="sympecma", rank='genus')

#occurance background points of similar species
isch_fenno<- readRDS("data/fenno_ischnura.rds")
sym_fenno<- readRDS("data/fenno_sympecma.rds")
coe_fenno<- readRDS("data/fenno_coenagrion.rds")
ery_fenno<- readRDS("data/fenno_erythromma.rds")
pyr_fenno<- readRDS("data/fenno_pyrrhosoma.rds")
pla_fenno<- readRDS("data/fenno_platycnemis.rds")
bg_vect<- readRDS("data/background_points.rds")
bg_points<- read.table("data/background_fenno.txt", header=T, sep="\t")


maps::map('world',xlim=c(1,40), ylim=c(50,71))
plot(Fenno)
#points(ischpts2, "red", cex=0.1) #ischnura
points(ipts, "red", cex=0.1) #ischnura
#points(isch_fenno, "blue", cex=0.1) #ischnura
points(sym_fenno, col="blue", cex=0.1)
points(coe_fenno, col="green", cex=0.1)
points(ery_fenno, col="orange", cex=0.1)
points(pyr_fenno, col="pink", cex=0.1)
points(pla_fenno, col="purple", cex=0.1)

plot(Fenno)
points(ipts, "red", cex=0.1)
points(bg_vect, "blue", cex=0.1)



##get environmental data
terr_cls<- terrain.colors(100, rev=T)
clim_fenno<- rast("data/climate/wc2.1_country/clim_fenno.tif")
plot(clim_fenno, col=terr_cls)

#cleaned and combined ischnura points, no duplicates, with extracted env data
isch_env_dup<- read.table("data/isch_env.txt", sep="\t", header=T)

# Check for collinearity and reduce the data set to only weakly correlated variables.
cor_mat <- cor(isch_env_dup[,-c(1:3,23:24)], method='spearman')
# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

#from sdm package tutorial
vif(clim_fenno)
v <- vifstep(isch_env_dup[,-c(1:3,23:24)]) #i think threshold is by default 10

# Run select07()
var_sel <- select07(X=isch_env_dup[,-c(1:3,23:24)], 
                    y=isch_env_dup$Ischnura_elegans, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel
sum(isch_env_dup$Ischnura_elegans)

#the variables in pred_sel vs v are slightly different!
#the selection from v i think is better so using that going forward

# Define a full model with all weakly correlated variables including their linear and quadratic terms.
# Fit the full model:
m_full_isch <- glm( Ischnura_elegans ~ bio3 + I(bio3^2) + bio8 + I(bio8^2)
                  + bio9 + I(bio9^2)+ bio10 + I(bio10^2)+ bio11 + I(bio11^2)
                  + bio15 + I(bio15^2)+ bio18 + I(bio18^2),
                  family='binomial', data=isch_env_dup)

# Inspect the model:
summary(m_full_isch)
# Run the full model.
# Simplify the model using stepwise variable selection step()
# Explained deviance:
expl_deviance(obs = isch_env_dup$Ischnura_elegans,
              pred = m_full_isch$fitted)

m_step_isch <- step(m_full_isch) 
# Inspect the model:
summary(m_step_isch)
# Compare the full model and the reduced model in terms of AIC and explained deviance.
# Explained deviance:
expl_deviance(obs = isch_env_dup$Ischnura_elegans,
              pred = m_step_isch$fitted)

