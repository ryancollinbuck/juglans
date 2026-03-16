qrsh -l h_rt=24:00:00,h_vmem=25G -pe shared 10
module load gcc/10.2.0
module load gdal/3.1.3
module load geos/3.11.1
module load proj/7.1.1
module load sqlite/3.33.0
module load curl/8.4.0
module load R/4.3.0
module load pandoc/2.17.1.1
module load cmake/3.30.0
cd project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY/Jcal2Walnuts




SUBsnp <- read.table("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.forR", header = T, row.names = 1)
library(parallel)
cl<-makeCluster(detectCores())
SUBsnp.imp <- parApply(cl,SUBsnp, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #remove NAs and replace with most common genotype

library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)


library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

sample.coord.sp <- vect(Coordinates, geom = c("x", "y"), crs = crs(ras_8110))
clim.points_8110 <- extract(ras_8110, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
clim.points_8110<-data.frame(clim.points_8110)


library(gradientForest)
env.gf <- cbind(clim.points_8110[,-1])
maxLevel <- log2(0.368*nrow(env.gf)/2)
gf <- gradientForest(cbind(env.gf, SUBsnp.imp), predictor.vars=colnames(env.gf), response.vars=colnames(SUBsnp.imp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)







#mask current clim by donor range, then GF predict across it
library(adehabitatLT)
library(extendedForest) 
library(fields)
library(gradientForest)
library(reshape)
library(parallel)
library(doParallel)
library(readr)
library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp45_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C45_2040 <- rast(tif_files)
layer_names <- gsub("CNRMrcp45_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C45_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C45_2040 <- project(ras_C45_2040, new_crs, method = "bilinear")
ras_C45_2040[ras_C45_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp45_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C45_2070 <- rast(tif_files)
layer_names <- gsub("CNRMrcp45_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C45_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C45_2070 <- project(ras_C45_2070, new_crs, method = "bilinear")
ras_C45_2070[ras_C45_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp85_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C85_2040 <- rast(tif_files)
layer_names <- gsub("CNRMrcp85_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C85_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C85_2040 <- project(ras_C85_2040, new_crs, method = "bilinear")
ras_C85_2040[ras_C85_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp85_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C85_2070 <- rast(tif_files)
layer_names <- gsub("CNRMrcp85_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C85_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C85_2070 <- project(ras_C85_2070, new_crs, method = "bilinear")
ras_C85_2070[ras_C85_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp45_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H45_2040 <- rast(tif_files)
layer_names <- gsub("HADES_rcp45_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H45_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H45_2040 <- project(ras_H45_2040, new_crs, method = "bilinear")
ras_H45_2040[ras_H45_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp45_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H45_2070 <- rast(tif_files)
layer_names <- gsub("HADES_rcp45_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H45_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H45_2070 <- project(ras_H45_2070, new_crs, method = "bilinear")
ras_H45_2070[ras_H45_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp85_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H85_2040 <- rast(tif_files)
layer_names <- gsub("HADES_rcp85_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H85_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H85_2040 <- project(ras_H85_2040, new_crs, method = "bilinear")
ras_H85_2040[ras_H85_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp85_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H85_2070 <- rast(tif_files)
layer_names <- gsub("HADES_rcp85_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H85_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H85_2070 <- project(ras_H85_2070, new_crs, method = "bilinear")
ras_H85_2070[ras_H85_2070 == -9999] <- NA

remove.NAs.stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_H45_2040 <- remove.NAs.stack(ras_H45_2040)
ras_H45_2070 <- remove.NAs.stack(ras_H45_2070)
ras_H85_2040 <- remove.NAs.stack(ras_H85_2040)
ras_H85_2070 <- remove.NAs.stack(ras_H85_2070)
ras_C45_2040 <- remove.NAs.stack(ras_C45_2040)
ras_C45_2070 <- remove.NAs.stack(ras_C45_2070)
ras_C85_2040 <- remove.NAs.stack(ras_C85_2040)
ras_C85_2070 <- remove.NAs.stack(ras_C85_2070)


###################################################################################
#################################################################################
###############################################################################
Jcal_2p5 <- vect("Jcal_2p5km_ALL_Clip.shp")
Jcal_2p5 <- project(Jcal_2p5, crs(ras_8110))
Jcalmask <- mask(crop(ras_8110, Jcal_2p5), Jcal_2p5)
currClim <- as.data.frame(Jcalmask,cells = TRUE)
currClim <- na.omit(currClim)
transClim<-predict(gf,currClim[,-(1)])

# SPATIO-TEMPORAL OFFSETS -------------------------------------------------
currClim1<- as.data.frame(Jcalmask,xy=TRUE,cells = TRUE)
currClim.coords <- na.omit(currClim1)
donOff_xy<-data.frame(currClim.coords[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there

## GF transform climate data
# Current transformed climate (coordinates required for mapping - not part of this script)
currClim_trans_xy <- data.frame(donOff_xy, transClim)
currClim_trans_xy <- currClim_trans_xy[order(currClim_trans_xy$cell),]
# Remove xy coordinates and cell index
currClim_trans<-currClim_trans_xy[,-(1:3)]


CDL_Walnuts <- vect("CDL_Cali_Walnuts_Polygon.shp")
CDL_Walnuts <- project(CDL_Walnuts, crs(ras_8110))


#to C852070
Walnutsmask_C852070 <- mask(crop(ras_C85_2070, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_C852070_points <- as.data.frame(Walnutsmask_C852070,xy=TRUE,cells = TRUE)
Walnutsmask_C852070_points <- na.omit(Walnutsmask_C852070_points)
pred_C85_2070 <- predict(gf, Walnutsmask_C852070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_C85_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_C85_2070
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C85_2070$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_C85_2070$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_C852070.tif", filetype="GTiff", overwrite=TRUE)



#to C852040
Walnutsmask_C852040 <- mask(crop(ras_C85_2040, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_C852040_points <- as.data.frame(Walnutsmask_C852040,xy=TRUE,cells = TRUE)
Walnutsmask_C852040_points <- na.omit(Walnutsmask_C852040_points)
pred_C85_2040 <- predict(gf, Walnutsmask_C852040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_C85_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_C85_2040
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C85_2040$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_C85_2040$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_C852040.tif", filetype="GTiff", overwrite=TRUE)



#to C452070
Walnutsmask_C452070 <- mask(crop(ras_C45_2070, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_C452070_points <- as.data.frame(Walnutsmask_C452070,xy=TRUE,cells = TRUE)
Walnutsmask_C452070_points <- na.omit(Walnutsmask_C452070_points)
pred_C45_2070 <- predict(gf, Walnutsmask_C452070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_C45_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_C45_2070
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C45_2070$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_C45_2070$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_C452070.tif", filetype="GTiff", overwrite=TRUE)



#to C452040
Walnutsmask_C452040 <- mask(crop(ras_C45_2040, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_C452040_points <- as.data.frame(Walnutsmask_C452040,xy=TRUE,cells = TRUE)
Walnutsmask_C452040_points <- na.omit(Walnutsmask_C452040_points)
pred_C45_2040 <- predict(gf, Walnutsmask_C452040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_C45_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_C45_2040
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C45_2040$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_C45_2040$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_C452040.tif", filetype="GTiff", overwrite=TRUE)



#to H852070
Walnutsmask_H852070 <- mask(crop(ras_H85_2070, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_H852070_points <- as.data.frame(Walnutsmask_H852070,xy=TRUE,cells = TRUE)
Walnutsmask_H852070_points <- na.omit(Walnutsmask_H852070_points)
pred_H85_2070 <- predict(gf, Walnutsmask_H852070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_H85_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_H85_2070
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H85_2070$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_H85_2070$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_H852070.tif", filetype="GTiff", overwrite=TRUE)



#to H852040
Walnutsmask_H852040 <- mask(crop(ras_H85_2040, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_H852040_points <- as.data.frame(Walnutsmask_H852040,xy=TRUE,cells = TRUE)
Walnutsmask_H852040_points <- na.omit(Walnutsmask_H852040_points)
pred_H85_2040 <- predict(gf, Walnutsmask_H852040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_H85_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_H85_2040
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H85_2040$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_H85_2040$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_H852040.tif", filetype="GTiff", overwrite=TRUE)



#to H452070
Walnutsmask_H452070 <- mask(crop(ras_H45_2070, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_H452070_points <- as.data.frame(Walnutsmask_H452070,xy=TRUE,cells = TRUE)
Walnutsmask_H452070_points <- na.omit(Walnutsmask_H452070_points)
pred_H45_2070 <- predict(gf, Walnutsmask_H452070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_H45_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_H45_2070
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H45_2070$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_H45_2070$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_H452070.tif", filetype="GTiff", overwrite=TRUE)



#to H452040
Walnutsmask_H452040 <- mask(crop(ras_H45_2040, CDL_Walnuts), CDL_Walnuts)
Walnutsmask_H452040_points <- as.data.frame(Walnutsmask_H452040,xy=TRUE,cells = TRUE)
Walnutsmask_H452040_points <- na.omit(Walnutsmask_H452040_points)
pred_H45_2040 <- predict(gf, Walnutsmask_H452040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(currClim_trans), cut(1:nrow(currClim_trans), 100, labels = FALSE))
Offset_H45_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    currClim_trans[breakItAll[[i]], ],
    pred_H45_2040
  )
  # Return only what you need
  data.frame(
    cell = currClim$cell[breakItAll[[i]]],
    avgDonOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- Jcalmask
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H45_2040$avgDonOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[donOff_xy$cell]<-Offset_H45_2040$avgDonOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGdonorOffset_toWalnuts_H452040.tif", filetype="GTiff", overwrite=TRUE)
