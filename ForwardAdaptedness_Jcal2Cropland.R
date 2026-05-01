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
cd project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY/Jcal2Cropland




# =============================================================================
# ForwardAdaptedness_Jcal2Cropland.R
# "Forward" genomic offset: how well-adapted is J. californica to future
# conditions at current walnut cropland locations in California?
#
# Same workflow as ForwardAdaptedness_Jhin2Cropland.R but uses J. californica
# adaptive outlier SNPs and sample coordinates.
#
# Approach:
#   1. Load RDA-identified adaptive outlier SNPs for J. californica.
#   2. Impute missing genotypes (modal allele, parallel apply).
#   3. Extract current climate at sample coordinates via terra.
#   4. Train Gradient Forest (GF) on genotype × climate.
#   5. For each future scenario (CNRM & HadGEM2 × RCP4.5/8.5 × 2040-69/2070-99):
#      a. Mask future climate to CDL walnut cropland pixels.
#      b. GF-predict current (donor) and future (recipient) climate.
#      c. Compute average recipient offset (mean Euclidean distance in GF space).
#      d. Write offset raster to GeoTIFF.
# =============================================================================

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


CDL_ALLcrops <- vect("CDL_Cali_ALLcrops_Polygon.shp")
CDL_ALLcrops <- project(CDL_ALLcrops, crs(ras_8110))

#TO C852070
ALLcropsmask_C852070 <- mask(crop(ras_C85_2070, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_C852070_points <- as.data.frame(ALLcropsmask_C852070,xy=TRUE,cells = TRUE)
ALLcropsmask_C852070_points <- na.omit(ALLcropsmask_C852070_points)
recOff_xy<-data.frame(ALLcropsmask_C852070_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_C85_2070 <- predict(gf, ALLcropsmask_C852070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_C85_2070), cut(1:nrow(pred_C85_2070), 100, labels = FALSE))
Offset_C85_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_C85_2070[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_C852070_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_C852070
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C85_2070$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_C85_2070$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_C852070.tif", filetype="GTiff", overwrite=TRUE)



#TO C852040
ALLcropsmask_C852040 <- mask(crop(ras_C85_2040, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_C852040_points <- as.data.frame(ALLcropsmask_C852040,xy=TRUE,cells = TRUE)
ALLcropsmask_C852040_points <- na.omit(ALLcropsmask_C852040_points)
recOff_xy<-data.frame(ALLcropsmask_C852040_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_C85_2040 <- predict(gf, ALLcropsmask_C852040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_C85_2040), cut(1:nrow(pred_C85_2040), 100, labels = FALSE))
Offset_C85_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_C85_2040[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_C852040_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_C852040
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C85_2040$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_C85_2040$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_C852040.tif", filetype="GTiff", overwrite=TRUE)



#TO C452070
ALLcropsmask_C452070 <- mask(crop(ras_C45_2070, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_C452070_points <- as.data.frame(ALLcropsmask_C452070,xy=TRUE,cells = TRUE)
ALLcropsmask_C452070_points <- na.omit(ALLcropsmask_C452070_points)
recOff_xy<-data.frame(ALLcropsmask_C452070_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_C45_2070 <- predict(gf, ALLcropsmask_C452070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_C45_2070), cut(1:nrow(pred_C45_2070), 100, labels = FALSE))
Offset_C45_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_C45_2070[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_C452070_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_C452070
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C45_2070$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_C45_2070$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_C452070.tif", filetype="GTiff", overwrite=TRUE)



#TO C452040
ALLcropsmask_C452040 <- mask(crop(ras_C45_2040, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_C452040_points <- as.data.frame(ALLcropsmask_C452040,xy=TRUE,cells = TRUE)
ALLcropsmask_C452040_points <- na.omit(ALLcropsmask_C452040_points)
recOff_xy<-data.frame(ALLcropsmask_C452040_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_C45_2040 <- predict(gf, ALLcropsmask_C452040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_C45_2040), cut(1:nrow(pred_C45_2040), 100, labels = FALSE))
Offset_C45_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_C45_2040[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_C452040_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_C452040
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_C45_2040$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_C45_2040$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_C852040.tif", filetype="GTiff", overwrite=TRUE)



#TO H852070
ALLcropsmask_H852070 <- mask(crop(ras_H85_2070, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_H852070_points <- as.data.frame(ALLcropsmask_H852070,xy=TRUE,cells = TRUE)
ALLcropsmask_H852070_points <- na.omit(ALLcropsmask_H852070_points)
recOff_xy<-data.frame(ALLcropsmask_H852070_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_H85_2070 <- predict(gf, ALLcropsmask_H852070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_H85_2070), cut(1:nrow(pred_H85_2070), 100, labels = FALSE))
Offset_H85_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_H85_2070[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_H852070_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_H852070
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H85_2070$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_H85_2070$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_H852070.tif", filetype="GTiff", overwrite=TRUE)



#TO H852040
ALLcropsmask_H852040 <- mask(crop(ras_H85_2040, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_H852040_points <- as.data.frame(ALLcropsmask_H852040,xy=TRUE,cells = TRUE)
ALLcropsmask_H852040_points <- na.omit(ALLcropsmask_H852040_points)
recOff_xy<-data.frame(ALLcropsmask_H852040_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_H85_2040 <- predict(gf, ALLcropsmask_H852040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_H85_2040), cut(1:nrow(pred_H85_2040), 100, labels = FALSE))
Offset_H85_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_H85_2040[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_H852040_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_H852040
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H85_2040$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_H85_2040$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_H852040.tif", filetype="GTiff", overwrite=TRUE)



#TO H452070
ALLcropsmask_H452070 <- mask(crop(ras_H45_2070, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_H452070_points <- as.data.frame(ALLcropsmask_H452070,xy=TRUE,cells = TRUE)
ALLcropsmask_H452070_points <- na.omit(ALLcropsmask_H452070_points)
recOff_xy<-data.frame(ALLcropsmask_H452070_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_H45_2070 <- predict(gf, ALLcropsmask_H452070_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_H45_2070), cut(1:nrow(pred_H45_2070), 100, labels = FALSE))
Offset_H45_2070 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_H45_2070[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_H452070_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_H452070
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H45_2070$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_H45_2070$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_H452070.tif", filetype="GTiff", overwrite=TRUE)



#TO H452040
ALLcropsmask_H452040 <- mask(crop(ras_H45_2040, CDL_ALLcrops), CDL_ALLcrops)
ALLcropsmask_H452040_points <- as.data.frame(ALLcropsmask_H452040,xy=TRUE,cells = TRUE)
ALLcropsmask_H452040_points <- na.omit(ALLcropsmask_H452040_points)
recOff_xy<-data.frame(ALLcropsmask_H452040_points[,c("x","y","cell")]) ##definitely check if x,y, and cellnumbers are in there
pred_H45_2040 <- predict(gf, ALLcropsmask_H452040_points[,-(1:3)])  #note the removal of the cell ID column with [,-1])
###AVERAGE OFFSETS#####
########running it#########
cl <- makeCluster(3)
registerDoParallel(cl)
breakItAll <- split(1:nrow(pred_H45_2040), cut(1:nrow(pred_H45_2040), 100, labels = FALSE))
Offset_H45_2040 <- foreach(
  i = seq_along(breakItAll),
  .packages = "fields",
  .combine = rbind
) %dopar% {
  d <- rdist(
    pred_H45_2040[breakItAll[[i]], ],
    currClim_trans
  )
  # Return only what you need
  data.frame(
    cell = ALLcropsmask_H452040_points$cell[breakItAll[[i]]],
    avgRecOff = rowMeans(d, na.rm = TRUE)
  )
}
original_raster <- ALLcropsmask_H452040
extent_orig <- ext(original_raster)
res_orig <- res(original_raster)
ncols_orig <- ncol(original_raster)
nrows_orig <- nrow(original_raster)
cell_values <- Offset_H45_2040$avgRecOff 
new_DIraster <- rast(nrows = nrows_orig, ncols = ncols_orig, 
                     ext = extent_orig, res = res_orig)
values(new_DIraster)[recOff_xy$cell]<-Offset_H45_2040$avgRecOff
crs(new_DIraster) <- crs(original_raster)
writeRaster(new_DIraster, "newloop_JC_AVGrecipientOffset_toALLcrops_H852040.tif", filetype="GTiff", overwrite=TRUE)
