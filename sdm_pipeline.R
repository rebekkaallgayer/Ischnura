setwd("Ischnura/")

library(lubridate)
library(sdm)
library(usdm)
library(terra)
library(sdmpredictors)
library(parallel)
library(dplyr)
library(viridis)
library(elsa)
library(psych) # principal() function
library(dismo)


# Study area
indochina <- terra::vect("Data/Sabah/01_Spatial object_Indochina inset.shp")
study_area <- crop(indochina, c(114, 121, 3, 9))

Swedenmap <- terra::rast(paste0(dirpath, "Inputs/hab_ras.txt"))


# Environmental variables (EV)
ev.dry <- rast("Data/Environment/4_Predictor_stacks/6_Raster_9_dry.season.tif")

# CPUE map
cpue.dry <- rast("Data/CPUE/01_total/cpue_bss.grd")



# Test for colinnearity of EVs
vif(ev.dry)
vifstep(ev.dry) # None have collinearity problems


#############################





# Read in the CPUE file for the wet season
cpue.dry <- rast("Data/CPUE/01_total/cpue_bss.grd")

# Extract non-zero and non-NA cells
cell.dry <- which(!is.na(values(cpue.dry)) & values(cpue.dry) != 0)

# Set a seed for reproducibility
set.seed(1234)

# Initialize lists to store model objects and prediction rasters
model_data_list <- list()
model_list <- list()


# Loop to sample presence points, input data using sdmData())= and run models using sdm()
for (i in 1:30) {
    
    # Feedback text
    cat("Running model", i, "out of 30\n")
    
    tryCatch({
        
        # 1) Subsample the map to get presence points
        presence.dry <- sample(cell.dry,
                               size = 301,
                               prob = values(cpue.dry)[cell.dry],
                               replace = TRUE)
        
        # Get the coordinates of the selected cells
        pres_coord.dry <- xyFromCell(cpue.dry, presence.dry)
        
        # Create a data frame with coordinates and species presence
        bycatch.dry <- data.frame(pres_coord.dry, species = 1)
        
        # Convert the presence data to spatial vector format
        presence_vect <- vect(bycatch.dry, c("x", "y"))
        
        
        # ### The following commands are for keeping only a single presence point per cell
        
        # Get the cell numbers corresponding to these coordinates
        presence_cells <- cellFromXY(cpue.dry, as.matrix(crds(presence_vect)))
        
        # Create a data frame with cell indices and coordinates
        bycatch.dry$cell <- presence_cells
        
        # Keep only the first occurrence of each unique cell
        unique_bycatch_dry <- bycatch.dry[!duplicated(bycatch.dry$cell), ]
        
        unique_bycatch_dry_vect <- vect(unique_bycatch_dry[ , ],
                                        geom = c("x", "y"),
                                        crs=crs(cpue.dry))
        
        # Print the number of unique cells with presence points
        cat("Number of unique cells with presence points: ", nrow(unique_bycatch_dry_vect), "\n")
        
        # Prepare the data for SDM
        model_data.d <- sdmData(formula = species ~ .,
                                train = unique_bycatch_dry_vect,
                                predictors = ev.dry,
                                impute = TRUE,
                                bg = list(n = nrow(unique_bycatch_dry),
                                          method = "gRandom",
                                          bias = bias_map_2)
        )
        
        # Store the model in the list
        model_data_list[[i]] <- model_data.d
        
        # Save sdmData as .sdd
        file_name <- paste0("Data/SDM_results/bss/total/model_data/sdmData_", i, "_total.sdd")
        write.sdm(model_data.d, file_name, overwrite = TRUE)
        
        # Confirm sdmData creation
        cat("Model_data_", i, "saved\n")
        
        # Run the SDM
        model.d <- sdm(formula = species ~ .,
                       data = model_data.d,
                       methods = c("gam", "brt", "rf", "maxent", "maxnet"),
                       replication = c("cv", "boot"),
                       n = 5,
                       test.percent = 30,
                       parallelSetting = list(ncore = 11, method = "parallel"))
        
        # Store the model in the list
        model_list[[i]] <- model.d
        
        # Save the model using write.sdm
        file_name <- paste0("Data/SDM_results/bss/total/model_iterations/model_", i, "_total.sdm")
        write.sdm(model.d, file_name, overwrite = TRUE)
        
        cat("Model", i, "saved as", file_name, "\n")
        
    }, error = function(e) {
        
        cat("Error encountered in iteration", i, ": ", e$message, "\n")
        
    })
}


cat("All iterations completed\n")



# Add up all generated models
model_list

model_list2 <- model_list

null_rm<- c() # Create vector of NULL models (models that failed to run)

for(t in 1:length(model_list2)){
    
    if(is.null(model_list2[[t]])) {
        
        null_rm<- c(null_rm, t)
    }
}

null_rm

model_list2 <- model_list2[-null_rm]

# Combine all successfully created models (of the 30 iterations)
# i.e., model_combined <- model[[ 1 ]] + model[[ 2 ]] + .... + model[[ i ]]

model_combined.d <- Reduce(`+`, model_list2)
model_combined.d
write.sdm(model_combined.d, "Data/SDM_results/bss/total/bss_sdm_total.sdm", overwrite = TRUE)

gui(model_combined.d)

eval <- getEvaluation(model_combined.d,
                      stat = c("AUC", "TSS", "COR"),
                      wtest = "test.dep")
head(eval)
nrow(eval)

mean(eval$AUC)
mean(eval$TSS)
mean(eval$COR)


plot(getVarImp(model_combined.d))





# After all models are trained and predictions are generated, combine them into an ensemble model
ensemble_model.d <- ensemble(model_combined.d,
                             newdata = ev.dry,
                             setting = list(method = "weighted",
                                            stat = "auc",
                                            expr = "auc > 0.75 & tss > 0.5"),
                             filename = "Data/SDM_results/bss/total/maps/bss_ensemble_total.tif",
                             overwrite = TRUE
)


# Plot the ensemble model
plot(ensemble_model.d)
plot(study_area, col = "grey90", add = TRUE)



# Example niche model 
niche(x = ev.dry, # predictors
      h = ensemble_model.d, # ensemble model
      n = c('salinity','reef') # EVs to plot
)
