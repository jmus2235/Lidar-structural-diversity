## ----call-libraries, results="hide"--------------------------------------------

############### Packages ################### 
library(lidR)
library(gstat)
library(ggplot2)
library(terra)
library(raster)
############### Set working directory ######
#set the working of the downloaded tutorial folder
wd <- "~/R_Scripts/NEON-Lidar-structural-diversity/" #This will depend on your local machine
setwd(wd)

## USER INPUTS
x <- 295500 
y <- 4101500
easting <- '300000'
northing <- '4101000'
min_elev <- 980
max_elev <- 1200
severity <- "unburned"
  ## Severity class to be evaluated
AOI <- (file=sprintf("./data/Shapefiles/dNBR_%s_severity_%s_%s.shp", severity, easting, northing))

  ## DTM to be used for elevation correction
dtm2021 <- (file=sprintf('./data/DEM/NEON_D17_SOAP_DP3_%s_%s_DTM.tif', easting, northing))
dtm2021 <- rast(dtm2021)
dtm <- raster(dtm2021)

## ----read-in-lidar-data--------------------------------------------------------

############ Read in LiDAR data ###########
#2017 1 km2 tile .laz file type for HARV and SOAP

#Watch out for outlier Z points - this function also allows for the
#ability to filter outlier points well above or below the landscape
#(-drop_z_blow and -drop_z_above). See how we have done this here 
#for you.

# SOAP2019 <- readLAS(paste0(wd,"NEON_D17_SOAP_DP1_299000_4096000_classified_point_cloud_colorized_2019.laz"),
#                     filter = "-drop_z_below 673 -drop_z_above 900")

#SOAP2021 <- readLAS(paste0(wd,"NEON_D17_SOAP_DP1_299000_4096000_classified_point_cloud_colorized_2021.laz"),
#                    filter = "-drop_z_below 673 -drop_z_above 900")

File2019 <- (file=sprintf("./data/Lazfiles/NEON_D17_SOAP_DP1_%s_%s_classified_point_cloud_colorized_2019.laz", easting, northing))

File2021 <- (file=sprintf("./data/Lazfiles/NEON_D17_SOAP_DP1_%s_%s_classified_point_cloud_colorized_2021.laz", easting, northing))

SOAP2019 <- readLAS(paste0(wd,File2019),
                    filter = sprintf("-drop_z_below %s -drop_z_above %s", min_elev, max_elev))

SOAP2021 <- readLAS(paste0(wd,File2021),
                    filter = sprintf("-drop_z_below %s -drop_z_above %s", min_elev, max_elev))

## ----plot-and-summarize-laz-file, eval=F, comment=NA---------------------------
############## Look at data specs ######
#Let's check out the extent, coordinate system, and a 3D plot of each 
#.laz file.
## ----plot-SOAP-1km2-point-cloud, eval=F, comment=NA----------------------------

summary(SOAP2019)
#plot(SOAP2019, color='RGB')

summary(SOAP2021)
#plot(SOAP2021, color='RGB')

# ## ----correct-for-elevation-----------------------------------------------------
# 
# ############## Correct for elevation #####
# #Correct for ground height using a kriging function to interpolate 
# #elevation from ground points in the .laz file.
# #If the function will not run, then you may need to checkfor outliers
# #by adjusting the 'drop_z_' arguments when reading in the .laz files.

# ## -- test crop with shapefile
# f <- system.file("extdata", AOI, package = "lidR")
severity_aoi <- sf::st_read(AOI, quiet = TRUE)
subset2019 <- clip_roi(SOAP2019, severity_aoi)
subset2021 <- clip_roi(SOAP2021, severity_aoi)

#dtm <- grid_terrain(subset2019, 1, kriging(k = 10L))
subset2019 <- normalize_height(subset2019, dtm)

#dtm <- grid_terrain(subset2021, 1, kriging(k = 10L))
subset2021 <- normalize_height(subset2021, dtm)


# #Will often give a warning if not all points could be corrected, 
# #but visually check to see if it corrected for ground height. 
plot(subset2019, color="RGB")
plot(subset2021, color="RGB")

# #There's only a few uncorrected points and we'll fix these in 
# #the next step. 

subset2019@data$Z[subset2019@data$Z <= 0.1] <- 0
subset2021@data$Z[subset2021@data$Z <= 0.1] <- 0  
# #This line filters out all z_vals below .5 m as we are less 
# #interested in shrubs/trees. 
# #You could change it to zero or another height depending on interests. 

#plot(subset2019, color='RGB')
#plot(subset2021, color='RGB')

# 
# ## ----ground classification based on lidR manual ----
# 
# # Three classification techniques
# 
# #data.200m_class <- classify_ground(data.200m, algorithm = pmf(ws = 5, th = 3)) # Progressive Morphological Filter (PMF)
# 
# subset2019_class <- classify_ground(subset2019, algorithm = csf()) # Cloth simulation filtering (CSF)
# subset2021_class <- classify_ground(subset2021, algorithm = csf()) # Cloth simulation filtering (CSF)
# 
# #data.200m_class <- classify_ground(data.200m, mcc(1.5,0.3)) # Multiscale Curvature Classification (MCC)
# 
# plot(subset2019_class, color = "Classification", size = 3, bg = "white")
# plot(subset2021_class, color = "Classification", size = 3, bg = "white")
# 
#las <- subset3_class
# 
# # function to display cross-section
# plot_crossection <- function(las,
#                              p1 = c(min(las@data$X), mean(las@data$Y)),
#                              p2 = c(max(las@data$X), mean(las@data$Y)),
#                              width = 4, colour_by = NULL)
# {
#   colour_by <- enquo(colour_by)
#   data_clip <- clip_transect(las, p1, p2, width)
#   p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + coord_equal() + theme_minimal()
#   
#   if (!is.null(colour_by))
#     p <- p + aes(color = !!colour_by) + labs(color = "")
#   
#   return(p)
# }
# 
# p1 <- c(727400, 4702500)
# p2 <- c(727600, 4702500)
# 
# plot_crossection(las, p1 , p2, colour_by = factor(Classification))

## ----structural-diversity-function---------------------------------------------

#Zip up all the code we previously used and write function to 
#run all 13 metrics in a single function for 2019. 
structural_diversity_metrics_2019 <- function(subset2019) {
  chm <- grid_canopy(subset2019, res = 1, dsmtin()) 
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  rumple <- rumple_index(chm) 
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) 
  deepgap.fraction <- deepgaps/cells 
  cover.fraction <- 1 - deepgap.fraction 
  vert.sd <- cloud_metrics(subset2019, sd(Z, na.rm = TRUE)) 
  sd.1m2 <- grid_metrics(subset2019, sd(Z), 1) 
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  Zs <- subset2019@data$Z
  Zs <- Zs[!is.na(Zs)]
  entro <- entropy(Zs, by = 1) 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  GFP.AOP <- mean(gap_frac$gf) 
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  out.plot <- data.frame(
    matrix(c(easting, northing, mean.max.canopy.ht,max.canopy.ht, 
             rumple,deepgaps, deepgap.fraction, 
             cover.fraction, top.rugosity, vert.sd, 
             sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP),
           ncol = 15)) 
  colnames(out.plot) <- 
    c("easting", "northing", "mean.max.canopy.ht.aop",
      "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      "deepgap.fraction.aop", "cover.fraction.aop",
      "top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      "entropy.aop", "GFP.AOP.aop",
      "VAI.AOP.aop", "VCI.AOP.aop") 
  print(out.plot)
}

#Zip up all the code we previously used and write function to 
#run all 13 metrics in a single function for 2021. 
structural_diversity_metrics_2021 <- function(subset2021) {
  chm <- grid_canopy(subset2021, res = 1, dsmtin()) 
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  rumple <- rumple_index(chm) 
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) 
  deepgap.fraction <- deepgaps/cells 
  cover.fraction <- 1 - deepgap.fraction 
  vert.sd <- cloud_metrics(subset2021, sd(Z, na.rm = TRUE)) 
  sd.1m2 <- grid_metrics(subset2021, sd(Z), 1) 
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  Zs <- subset2021@data$Z
  Zs <- Zs[!is.na(Zs)]
  entro <- entropy(Zs, by = 1) 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  GFP.AOP <- mean(gap_frac$gf) 
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  out.plot <- data.frame(
    matrix(c(easting, northing, mean.max.canopy.ht,max.canopy.ht, 
             rumple,deepgaps, deepgap.fraction, 
             cover.fraction, top.rugosity, vert.sd, 
             sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP),
           ncol = 15)) 
  colnames(out.plot) <- 
    c("easting", "northing", "mean.max.canopy.ht.aop",
      "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      "deepgap.fraction.aop", "cover.fraction.aop",
      "top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      "entropy.aop", "GFP.AOP.aop",
      "VAI.AOP.aop", "VCI.AOP.aop") 
  print(out.plot)
}


SOAP_structural_diversity_2019 <- structural_diversity_metrics_2019(subset2019)

SOAP_structural_diversity_2021 <- structural_diversity_metrics_2021(subset2021)

combined_results_2019_2021=rbind(SOAP_structural_diversity_2019, 
                       SOAP_structural_diversity_2021)
print(combined_results_2019_2021)

write.csv(combined_results_2019_2021, file = (sprintf("~/R_Scripts/NEON-Lidar-structural-diversity/data_out/combined_results_2019_2021_%s_severity_%s_%s.csv", severity, easting, northing)))
