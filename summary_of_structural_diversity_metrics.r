#GENERATE CANOPY HEIGHT MODEL (CHM) (i.e. a 1 m2 raster grid of
#vegetations heights)
#res argument specifies pixel size in meters and dsmtin is
#for raster interpolation
chm <- grid_canopy(subset3, res = 1, dsmtin())

#visualize CHM
plot(chm)

#MEAN OUTER CANOPY HEIGHT (MOCH)
#calculate MOCH, the mean CHM height value
mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE)

#MAX CANOPY HEIGHT
#calculate HMAX, the maximum CHM height value
max.canopy.ht <- max(chm@data@values, na.rm=TRUE)

#RUMPLE
#calculate rumple, a ratio of outer canopy surface area to
#ground surface area (1600 m^2)
rumple <- rumple_index(chm)

#TOP RUGOSITY
#top rugosity, the standard deviation of pixel values in chm and
#is a measure of outer canopy roughness
top.rugosity <- sd(chm@data@values, na.rm = TRUE)

#DEEP GAPS & DEEP GAP FRACTION
#number of cells in raster (also area in m2)
cells <- length(chm@data@values)
chm.0 <- chm
chm.0[is.na(chm.0)] <- 0 #replace NAs with zeros in CHM
#create variable for the number of deep gaps, 1 m^2 canopy gaps
zeros <- which(chm.0@data@values == 0)
deepgaps <- length(zeros) #number of deep gaps
#deep gap fraction, the number of deep gaps in the chm relative
#to total number of chm pixels
deepgap.fraction <- deepgaps/cells

#COVER FRACTION
#cover fraction, the inverse of deep gap fraction
cover.fraction <- 1 - deepgap.fraction

#HEIGHT SD
#height SD, the standard deviation of height values for all points
#in the plot point cloud
vert.sd <- cloud_metrics(subset3, sd(Z, na.rm = TRUE))

#SD of VERTICAL SD of HEIGHT
#rasterize plot point cloud and calculate the standard deviation
#of height values at a resolution of 1 m^2
sd.1m2 <- grid_metrics(subset3, sd(Z), 1)
#standard deviation of the calculated standard deviations
#from previous line
#This is a measure of internal and external canopy complexity
sd.sd <- sd(sd.1m2[,3], na.rm = TRUE)


#some of the next few functions won't handle NAs, so we need
#to filter these out of a vector of Z points
Zs <- subset3@data$Z
Zs <- Zs[!is.na(Zs)]

#ENTROPY
#entropy, quantifies diversity & evenness of point cloud heights
#by = 1 partitions point cloud in 1 m tall horizontal slices
#ranges from 0-1, with 1 being more evenly distributed points
#across the 1 m tall slices
entro <- entropy(Zs, by = 1)

#GAP FRACTION PROFILE
#gap fraction profile, assesses the distribution of gaps in the
#canopy volume
#dz = 1 partitions point cloud in 1 m horizontal slices
#z0 is set to a reasonable height based on the age and height of
#the study sites
gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
#defines gap fraction profile as the average gap fraction in each
#1 m horizontal slice assessed in the previous line
GFP.AOP <- mean(gap_frac$gf)

#VAI
#leaf area density, assesses leaf area in the canopy volume
#k = 0.5 is a standard extinction coefficient for foliage
#dz = 1 partitions point cloud in 1 m horizontal slices
#z0 is set to the same height as gap fraction profile above
LADen<-LAD(Zs, dz = 1, k=0.5, z0=3)
#vegetation area index, sum of leaf area density values for
#all horizontal slices assessed in previous line
VAI.AOP <- sum(LADen$lad, na.rm=TRUE)

#VCI
#vertical complexity index, fixed normalization of entropy
#metric calculated above
#set zmax comofortably above maximum canopy height
#by = 1 assesses the metric based on 1 m horizontal slices in
#the canopy
VCI.AOP <- VCI(Zs, by = 1, zmax=1200)