library(raster)
library(terra)
library(RStoolbox)

# set directory
setwd("D:/Remote Sensing/Carbon Magelang/Processing/Imagery/20 meters")
rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"))

# ==== PREPARE DATA ==== 
# load images
AP2 <- stack('AP2dbln_7extd.tif') # HH, HV (db & linear)
AP <- stack('AP2db_7extd.tif')    # HH, HV (db only)
S2 <- stack('S2_20extd.tif')      # 2, 3, 4, 5, 6, 7, 8a, 11, 12 (bands order)

# set layer names
names(AP2) <- c('HH_db','HV_db','HH_ln', 'HV_ln')
names(AP) <- c('HH_dB','HV_dB')
names(S2) <- c('b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8a', 'b11', 'b12')

# resample ALOS
AP2.rsp <- resample(AP2, S2, "bilinear",
                    filename = 'AP2dbln_20extd.tif', overwrite=TRUE)

AP.rsp <- resample(AP, S2, "bilinear",
                    filename = 'AP2db_20extd.tif', overwrite=TRUE)

# check data
S2
AP2.rsp
AP.rsp

# plot data
plot(AP.rsp, col = gray(seq(0, 1, length = 65536)))
hist(AP.rsp)

plot(AP2.rsp)
hist(AP2.rsp)
plot(S2)
hist(S2)

# normalize sentinel-2
S2.rsc <- stack(S2/10000)
plot(S2, col = gray(seq(0, 1, length = 4096)))
hist(S2, xlab = 'spectral reflectance')
S2
S2.rsc
summary(S2)
summary(S2.rsc)

# check normalized sentinel-2 data
S2.rsc
hist(S2.rsc, xlab = 'spectral reflectance')
plot(S2.rsc)
plot(S2.rsc, col = gray(seq(0, 1, length = 4096)))
?plot
?hist

# ==== CREATE VARIABLES ==== 
## create variable from Sentinel-2
EVI <- 2.5 * ((S2.rsc[[7]] - S2.rsc[[3]]) / (S2.rsc[[7]] + (6 * S2.rsc[[3]]) - (7.5 * S2.rsc[[1]]) + 1))
NDVI <- (S2.rsc[[7]] - S2.rsc[[3]]) / (S2.rsc[[7]] + S2.rsc[[3]])
RVI <- S2.rsc[[7]]/ S2.rsc[[3]]
SAVI <- ((S2.rsc[[7]] - S2.rsc[[3]]) / (S2.rsc[[7]] + S2.rsc[[3]] + 0.5)) * (1.5)

summary(EVI)

# set layer names and check indices
names(EVI) <- c('EVI')
names(NDVI) <- c('NDVI')
names(RVI) <- c('RVI')
names(SAVI) <- c('SAVI')
EVI
NDVI
RVI
SAVI

summary(EVI)
summary(NDVI)
summary(RVI)
summary(SAVI)

# plot indices

# Define color palette transitioning from red to green
color_palette <- colorRampPalette(c("red", "orange", "yellow", "green"))

# Plot the raster data with red to green color palette
par(mfrow=c(1,1))
par(mar=c(2,3,2,2))
?par
plot(EVI, main = "EVI", col = gray(seq(0, 1, length = 4096)))
plot(NDVI, main = "NDVI", col = gray(seq(0, 1, length = 4096)))
plot(RVI, main = "RVI", col = gray(seq(0, 1, length = 4096)))
plot(SAVI, main = "SAVI", col = gray(seq(0, 1, length = 4096)))

hist(EVI, main = "EVI", xlab = "pixel value")
hist(NDVI, main = "NDVI", xlab = "pixel value")
hist(RVI, main = "RVI", xlab = "pixel value")
hist(SAVI, main = "SAVI", xlab = "pixel value")

## create variable from ALOS-2 PALSAR-2
# create back scatter ratio with linear unit [3,4]
HHHV <- AP2.rsp[[3]] / AP2.rsp[[4]]
HVHH <- AP2.rsp[[4]] / AP2.rsp[[3]]

# set layer names and check back scatter ratio
names(HHHV) <- c('HH_HV')
names(HVHH) <- c('HV_HH')
HHHV
HVHH

# create back scatter indices with linear unit/ power [3,4]
RFDI <- (AP2.rsp[[3]] - AP2.rsp[[4]]) / (AP2.rsp[[3]] + AP2.rsp[[4]])
RVIhh <- (4 * AP2.rsp[[4]]) / (AP2.rsp[[3]] + AP2.rsp[[4]])

# set layer names and check back scatter indices
names(RFDI) <- c('RFDI')
names(RVIhh) <- c('RVI_hh')
RFDI
RVIhh

# plot SAR ratio and indices
par(mfrow=c(2,2))
plot(HHHV, main = "HH/HV", col = gray(seq(0, 1, length = 4096)))
plot(HVHH, main = "HV/HH", col = gray(seq(0, 1, length = 4096)))
plot(RFDI, main = "RFDI", col = gray(seq(0, 1, length = 4096)))
plot(RVIhh, main = "Radar Vegetation Index", col = gray(seq(0, 1, length = 4096)))

hist(HHHV, main = "HH/HV")
hist(HVHH, main = "HV/HH")
hist(RFDI, main = "RFDI")
hist(RVIhh, main = "RVIhh")


# ==== EXPORT DATA ==== 
# convert to SpatRaster (with terra)
Sen2 <- rast(S2.rsc)
ALOS <- rast(AP.rsp)
idxSen <- rast(raster::brick(EVI, NDVI, RVI, SAVI))
idxALOS <- rast((raster::brick(HHHV, HVHH, RFDI, RVIhh)))
predictor <- c(Sen2, ALOS, idxSen, idxALOS)

# check data
Sen2
ALOS
idxSen
idxALOS
predictor
names(predictor)

# write raster
terra::writeRaster(Sen2, "Sen2_final.tif", overwrite=TRUE)
terra::writeRaster(ALOS, "ALOS_final.tif",overwrite=TRUE)
terra::writeRaster(idxSen, "Sen_Idx_final.tif", overwrite=TRUE)
terra::writeRaster(idxALOS, "ALOS_Idx_final.tif", overwrite=TRUE)
terra::writeRaster(predictor, "predictor50_final.tif", overwrite=TRUE)

# ==== FOR RASTER RESULT ==== 
setwd("D:/Remote Sensing/Carbon Magelang/Result/Export_RF5")
res_gedi <- rast('agc_gedi_5.tif')
res_field <- rast('agc_field_5.tif')
plot(res_gedi)
plot(res_field)

res_merge <- c(res_gedi, res_field)
plot(res_merge)
hist(res_merge)

terra::writeRaster(res_merge, "agc_merge_5.tif", overwrite=TRUE)

