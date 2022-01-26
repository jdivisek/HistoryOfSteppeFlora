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

#load caspian sea
casp <- stack(casp.paths)
projection(casp) <- CRS("+proj=longlat +datum=WGS84")

#load aral sea
aral <- stack(aral.paths)
projection(aral) <- CRS("+proj=longlat +datum=WGS84")

#coarse map of the World 
m <- getMap("coarse")

#create layout
matr <- matrix(data=c(1:10), ncol=2, nrow = 5, byrow=T)

#labels for periods
proj.period.label <- c("present", "2.1 ka BP", "4.2 ka BP", "6.2 ka BP", "8.2 ka BP", "9.9 ka BP", "12.2 ka BP", "14.3 ka BP", "17.1 ka BP", "21 ka BP")

###FIGURE 1-----------------------------------------------------------------
tiff("Fig1.tif", width = 6.2, height = 8.6, units = "in", res=500, compression = "lzw")
layout(matr)
# layout.show(10)
par(mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), tcl=-0.3, mgp=c(3, 0.5, 0))

##All species
mod.paths <- paste0("./Stacked models/All species except Desert steppe/", proj.periods, ".sdat")
rs <- stack(mod.paths)
# max(cellStats(rs, max))

brks <- seq(0,88,8)
nb <- length(brks)-1 
cols <- colormap(colormap = colormaps$salinity, nshades = nb, reverse = T)

for(i in seq(length(proj.periods),1,-1))
{
  ##load model
  mod <- raster(mod.paths[i])
  
  plot(NA, type="n", xlim=c(-11,70), ylim=c(32.5,72), xaxs="i", yaxs="i", axes=F)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  abline(h=c(40,60), v=c(0,20,40,60), lwd=0.5, col="white")
  
  image(mod, add=T, breaks=brks, col=cols, legend=F)
  if(i > 6) {image(ice[[i]], add=T, col="white", legend=F)}
  image(casp[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  image(aral[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  text(-11,68, labels=proj.period.label[i], pos=4, cex=1.2)
  
  if(i %in% seq(10,2,-2)) {axis(2, at=c(40,60), labels = c("40°N", "60°N"))}
  if(i %in% c(1,2)) {axis(1, at=c(0,20,40,60), labels = c("0°", "20°E", "40°E", "60°E"))}
  
  if(i == 7) {text(77, 42, "∑", pos=3, xpd=NA, cex=1.8)}
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
tiff("Fig2.tif", width = 6.2, height = 8.6, units = "in", res=500, compression = "lzw")
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
  mod <- raster(mod.paths[i])
  
  plot(NA, type="n", xlim=c(-11,70), ylim=c(32.5,72), xaxs="i", yaxs="i", axes=F)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  abline(h=c(40,60), v=c(0,20,40,60), lwd=0.5, col="white")
  
  image(mod, add=T, breaks=brks, col=cols, legend=F)
  if(i > 6) {image(ice[[i]], add=T, col="white", legend=F)}
  image(casp[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  image(aral[[i]], add=T, col=rgb(169,211,236, maxColorValue = 255), legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  text(-11,68, labels=proj.period.label[i], pos=4, cex=1.2)
  
  if(i %in% seq(10,2,-2)) {axis(2, at=c(40,60), labels = c("40°N", "60°N"))}
  if(i %in% c(1,2)) {axis(1, at=c(0,20,40,60), labels = c("0°", "20°E", "40°E", "60°E"))}
  
  if(i == 7) {text(77, 42, "∑", pos=3, xpd=NA, cex=1.8)}
}

mtext("Habitat suitability for desert-steppe species", side = 3, line=0.5, outer=T, xpd=NA, font=2, cex=1)

#legend
par(mfrow=c(1,1), mar=c(0.4,0.4,0.4,0.4), oma=c(2,2,3,4), new=TRUE)
r <- raster(matrix(0:max(brks)))
addRasterLegend(r, direction = "vertical", location = "right", 
                side=4, cex.axis = 0.8, breaks=brks, ramp = cols, 
                ncolors = length(cols), locs = brks[-1], labelDist = 0.5)

dev.off()

