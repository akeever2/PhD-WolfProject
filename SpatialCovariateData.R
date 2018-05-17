# ##############################################################-

#       Code for getting spatial covariate data

###############################################################-


#### Import data, load packages ####

# First set the working directory
setwd("C:/Users/allison/Documents/Project/WolfData/")


# Function to streamline installing/opening packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


# Load or install these packages:
packages <- c("ks", "lattice", "adehabitatHR", "maptools", 
              "foreign", "rgdal", "sp", "raster","plot3D",
              "rasterVis", "sf", "colorRamps","rgeos", 
              "tidyverse")


# Run function to install packages - e.g., library command
ipak(packages)


# Pull in the different shapfiles needed
forcov <- raster("./wolfPackPts/forest_cover1.tif")
buff600 <- shapefile("./wolfPackPts/WolfPackLocs_Buffer6002Clipped.shp")
buff1700 <- shapefile("./wolfPackPts/WolfPackLocs_Buffer17002Clipped.shp")
rd2wd <- raster("./Montana Wolf POMS/GISfiles_covariates/Final_2WDrds/rddens_2wd1.tif")
rd4wd <- raster("./Montana Wolf POMS/GISfiles_covariates/Final_4WDrds/4wdRas2.tif")
wolflocs<-shapefile("./wolfPackPts/wolfPackLocsACK.shp")


#### Plot the data for visualization ####

# Make a very basic plot of raster
par(mfrow=c(1,1))
plot(forcov)

# Look at the class of the raster
class(forcov)

# Look at basic raster summary
forcov

# Look at raster data
forcov@data

# Look at structure of raster
str(forcov)

# Look at the projection of the raster (note the use of "@" instead of "$")
forcov@crs@projargs
wolflocs@proj4string@projargs
buff600@proj4string@projargs
forcov <- projectRaster(forcov, crs = wolflocs@proj4string@projargs)

# Look at the spatial extent of the raster
extent(forcov)

# Look at the resolution of the raster
res(forcov)

# Create a raster stack, this didn't work because need to be same extent
all.rasters<-stack(forcov, rd2wd, rd4wd, quick=TRUE)

# Check class
class(all.rasters)

spac.cov<-raster::extract(forcov, wolflocs, buffer=13813)
pi<-raster::intersect(buff600, forcov)
pi<-st_intersection(buff600, forcov)
