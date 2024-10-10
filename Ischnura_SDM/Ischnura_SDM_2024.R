library(data.table)
# library(randomForest)
# library(lattice)
# library(RColorBrewer)
# library(PresenceAbsence)
library(geodata)
library(terra)
library(sdm)
library(dismo)
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
#these models are terrible! 

#ensemble models
#Phisically exclude the collinear variables which are identified 
#using vifcor or vifstep from a set of variables.
clim_fenno_ex <- exclude(clim_fenno, v)

ipts_k<- kfold(ischpts2, k=3)
ipts_test<- ischpts2[ipts_k==3,]
ipts_train<- ischpts2[!ipts_k==3,]
ipts_train$species=1
ipts_train<- vect(ipts_train)

#~. means it takes species column and considers the rest as predictors, don't need to individually specify

d <- sdmData(species~., ipts_train, predictors= clim_fenno_ex, 
             bg = list(method='gRandom',n=1000))
d
m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m

p1 <- predict(m, clim_fenno_ex,filename='isch_predict_worldclim.img')
#p1<- rast("pr.img")
p1<- rast("isch_predict_worldclim.img")
plot(p1, col=terr_cls)
p1
names(p1)

#try it with set background points
bg_env<- read.table("data/bg_points_extract.txt", header=T, sep="\t")
bg_env1<- bg_env[,c(1:2, 6,11:14, 18,21)]
colnames(bg_env1)<- c("x", "y", colnames(bg_env[,-c(1:2)]))
d_bg <- sdmData(species~., ipts_train, predictors= clim_fenno_ex, bg = bg_env1)

m_bg <- sdm(species~., d_bg, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_bg
## the AUC and especially the TSS are WAY lower using the species-specific background points!

#maybe i can weight the background points by creating a CPUE landscape
#from the species specific data
d_wbg <- sdmData(species~., ipts_train, predictors= clim_fenno_ex, 
             bg = list(method='gRandom',n=1000,bias=null_rast))
d_wbg
m_wbg <- sdm(species~., d_wbg, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
            test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_wbg
#weighted by species specific background is exactly the same as using the points themselves
null_coarse50<- rast("data/cpue_background_coarse50.tif")
d_wcbg <- sdmData(species~., ipts_train, predictors= clim_fenno_ex, 
                 bg = list(method='gRandom',n=1000,bias=null_coarse50))
d_wcbg
m_wcbg <- sdm(species~., d_wcbg, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
             test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))
m_wcbg
#weighted by species specific background is only very slightly better at aggregation factor 10
#it does get better with a x50 aggregation factor but the TSS is still really low
