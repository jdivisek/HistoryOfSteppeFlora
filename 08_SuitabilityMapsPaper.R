###############################################################################
###                       SUITABILITY MAPS FOR PAPER                        ###
###############################################################################

library(rworldmap)
library(raster)
library(berryFunctions)
library(rangeBuilder)
library(viridis)
library(RColorBrewer)
library(colormap)
library(graticule)
library(rgdal)

# if(!require("V8")) install.packages("V8")
# if(!require("devtools")) install.packages("devtools")
# devtools::install_github("bhaskarvk/colormap")

##set periods (in years before present)
tw <- c(500, 2500, 4500, 6500, 8500, 9500, 12500, 14500, 17500, 20500)

##Load ice sheets, Caspian and Aral sea (grids in WGS84)
ice.paths <- paste0("./Grids/Ice sheets/ice.sheet_", tw, "BP.sdat")
names(ice.paths) <- proj.periods
casp.paths <- paste0("./Grids/Caspian_sea/casp.sea_", tw, "BP.sdat")
names(casp.paths) <- proj.periods
aral.paths <- paste0("./Grids/Aral_sea/aral.sea_", tw, "BP.sdat")
names(aral.paths) <- proj.periods

#load continental ice sheets
ice <- stack(ice.paths)
projection(ice) <- CRS("+proj=longlat +datum=WGS84")
ice <- projectRaster(ice, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                     res=c(4900,4900), method="ngb")

#load caspian sea
casp <- stack(casp.paths)
projection(casp) <- CRS("+proj=longlat +datum=WGS84")

#load aral sea
aral <- stack(aral.paths)
projection(aral) <- CRS("+proj=longlat +datum=WGS84")

#coarse map of the World 
m <- readOGR("./Shapes/worldmap/worldmap_coarse_EqEarth.shp")

#graticules
g <- graticule(lons=seq(-20,60,20), lats=c(40,60), xlim=c(-20,80), ylim=c(20,80), proj="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

#create layout
matr <- matrix(data=c(1:10), ncol=2, nrow = 5, byrow=T)

#labels for periods
proj.period.label <- c("present", "2.1 ka BP", "4.2 ka BP", "6.2 ka BP", "8.2 ka BP", "9.9 ka BP", "12.2 ka BP", "14.3 ka BP", "17.1 ka BP", "21 ka BP")

###FIGURE 1-----------------------------------------------------------------
tiff("Fig1.tif", width = 6.2, height = 9.6, units = "in", res=500, compression = "lzw")
layout(matr)
# layout.show(10)
par(mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), tcl=-0.3, mgp=c(3, 0.5, 0))

##All species
mod.paths <- paste0("./Stacked models/All species except Desert steppe/", proj.periods, ".sdat")
rs <- stack(mod.paths)
# max(cellStats(rs, max))

extent(casp) <- extent(rs)
extent(aral) <- extent(rs)

brks <- seq(0,88,8)
nb <- length(brks)-1 
cols <- colormap(colormap = colormaps$salinity, nshades = nb, reverse = T)

for(i in seq(length(proj.periods),1,-1))
{
  ##load model
  mod <- rs[[i]]
  mod <- mask(mod, casp[[i]], maskvalue=NA, inverse=T)
  mod <- mask(mod, aral[[i]], maskvalue=NA, inverse=T)
  mod <- projectRaster(mod, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                       res=c(4900,4900), method="bilinear")
  
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod, add=T, breaks=brks, col=cols, legend=F)
  if(i > 6) {image(ice[[i]], add=T, col=c(NA, "white"), legend=F)}
 
  plot(m, add=T, lwd=0.4)
  box()
  
  text(-2300000, 7500000, labels=proj.period.label[i], pos=4, cex=1.2)
  
  if(i %in% seq(10,2,-2)) {axis(2, at=c(4921020, 6924035), labels = c("40°N", "60°N"))}
  if(i %in% c(1,2)) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  if(i == 7) {text(3730000,  5000000, expression(sum()), pos=3, xpd=NA, cex=1.8)}
}

mtext("Habitat suitability for all species\nexcept those of desert-steppe", side = 3, line=0, outer=T, xpd=NA, font=2, cex=1)

#legend
par(mfrow=c(1,1), mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), new=TRUE)
r <- raster(matrix(0:max(brks)))
addRasterLegend(r, direction = "vertical", location = "right", 
                side=4, cex.axis = 0.8, breaks=brks, ramp = cols, 
                ncolors = length(cols), locs = brks[-1], labelDist = 0.5)

dev.off()

###FIGURE 2-----------------------------------------------------------------
tiff("Fig2.tif", width = 6.2, height = 9.6, units = "in", res=500, compression = "lzw")
layout(matr)
# layout.show(10)
par(mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), tcl=-0.3, mgp=c(3, 0.5, 0))

##Desert steppe
mod.paths <- paste0("./Stacked models/Desert steppe/", proj.periods, ".sdat")
rs <- stack(mod.paths)
# max(cellStats(rs, max))

brks <- seq(0,2,0.2)
nb <- length(brks)-1 
cols <- colormap(colormap = colormaps$salinity, nshades = nb, reverse = T)

for(i in seq(length(proj.periods),1,-1))
{
  ##load model
  mod <- rs[[i]]
  mod <- mask(mod, casp[[i]], maskvalue=NA, inverse=T)
  mod <- mask(mod, aral[[i]], maskvalue=NA, inverse=T)
  mod <- projectRaster(mod, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                       res=c(4900,4900), method="bilinear")
  
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod, add=T, breaks=brks, col=cols, legend=F)
  if(i > 6) {image(ice[[i]], add=T, col=c(NA, "white"), legend=F)}
  # image(casp[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  # image(aral[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  text(-2300000, 7500000, labels=proj.period.label[i], pos=4, cex=1.2)
  
  if(i %in% seq(10,2,-2)) {axis(2, at=c(4921020, 6924035), labels = c("40°N", "60°N"))}
  if(i %in% c(1,2)) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  if(i == 7) {text(3730000,  5000000, expression(sum()), pos=3, xpd=NA, cex=1.8)}
}

mtext("Habitat suitability for desert-steppe species", side = 3, line=0, outer=T, xpd=NA, font=2, cex=1)

#legend
par(mfrow=c(1,1), mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), new=TRUE)
r <- raster(matrix(0:max(brks)))
addRasterLegend(r, direction = "vertical", location = "right", 
                side=4, cex.axis = 0.8, breaks=brks, ramp = cols, 
                ncolors = length(cols), locs = brks[-1], labelDist = 0.5)

dev.off()
