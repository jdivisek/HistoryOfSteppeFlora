##########################################################################
###                    MAPS OF SOURCE AREAS FOR PAPER                  ###
##########################################################################

library(rworldmap)
library(raster)
library(berryFunctions)
library(viridis)
library(rangeBuilder)
library(graticule)

##set periods
tw <- c(500, 6500, 20500)

##Load ice sheets, Caspian and Aral sea (grids in WGS84)
ice.paths <- paste0("./Grids/Ice sheets/ice.sheet_", tw, "BP.sdat")
names(ice.paths) <- proj.periods
casp.paths <- paste0("./Grids/Caspian_sea/casp.sea_", tw, "BP.sdat")
names(casp.paths) <- proj.periods
aral.paths <- paste0("./Grids/Aral_sea/aral.sea_", tw, "BP.sdat")
names(aral.paths) <- proj.periods

#Project ice sheets
ice <- stack(ice.paths)
projection(ice) <- CRS("+proj=longlat +datum=WGS84")
ice <- projectRaster(ice, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                     res=c(4900,4900), method="ngb")

#coarse map of the World 
m <- readOGR("./Shapes/worldmap/worldmap_coarse_EqEarth.shp")

#graticules
g <- graticule(lons=seq(-20,60,20), lats=c(40,60), xlim=c(-20,80), ylim=c(20,80), proj="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

##import grid with asigned species ranges
cgrs.range <- read.delim("spec.UTM.WIDE_cgrs.txt", header = T, row.names = 1)
head(cgrs.range)

colnames(cgrs.range) == clas$abbrev
colnames(cgrs.range) <- clas$species

#shapefile
cgrs.wgs84 <- readOGR("./Shapes/CGRS grid/cgrs_sel.shp")
cgrs.EqEarth <- spTransform(cgrs.wgs84, CRS("+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

###FIGURE 3-----------------------------------------------------------------
##LGM and mid-Holocene source areas for ecological groups

#layout
matr <- matrix(data=NA, ncol=14, nrow = 3, byrow=T)
matr[1,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)
matr[2,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+5
matr[3,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+10

##Source areas in the LGM
mod <- stack(c(paste0("./Stacked source areas LGM/", "Rocky steppe", ".sdat"),
               paste0("./Stacked source areas LGM/", "Grass steppe", ".sdat"),
               paste0("./Stacked source areas LGM/", "Meadow steppe", ".sdat")))
mod <- projectRaster(mod, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                     res=c(4900,4900), method="ngb")
# plot(mod)

##Effective refugia in MH
mod2 <- stack(c(paste0("./Stacked source areas mid-Holocene/", "Rocky steppe", ".sdat"),
                paste0("./Stacked source areas mid-Holocene/", "Grass steppe", ".sdat"),
                paste0("./Stacked source areas mid-Holocene/", "Meadow steppe", ".sdat")))
mod2 <- projectRaster(mod2, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                      res=c(4900,4900), method="ngb")
# plot(mod2)

#maximum numbers of overlapping source areas
maxValue(mod)
maxValue(mod2)
max.brk <- c(14, 27, 33)

#plot
tiff("Fig_3.tif", width = 9.8, height = 5.9, units = "in", res=500, compression = "lzw")
layout(matr)
# layout.show(20)
par(mar=c(0.4,0.4,0.4,0.4), oma=c(2,5,4,0.5), tcl=-0.3, mgp=c(3, 0.5, 0))

for(i in 1:nlayers(mod))
{
  if(i == 1){brks <- seq(-2,max.brk[i],2)}
  if(i == 2){brks <- seq(-3,max.brk[i],3)}
  if(i == 3){brks <- seq(-3,max.brk[i],3)}
  nb <- length(brks)-1 
  cols <- c("#F2F2F2", rev(viridis(nb-1)))
  
  ##LGM
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod[[i]], add=T, breaks=brks, col=cols, legend=F)
  image(ice[[3]], add=T, col="white", legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  #add no of species
  rs <- stack(list.files("./Source areas LGM", full.names = T))
  if(i == 1) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$rocky.steppe)])]]}
  if(i == 2) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$grass.steppe)])]]}
  if(i == 3) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$meadow.steppe)])]]}
  maxValue(rs)
  rs <- rs[[which(maxValue(rs) == 1)]]
  n <- nlayers(rs)
  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Source areas\nin the Last Glacial Maximum", side=3, line=1, font=2)
  }
  if(i == 1){ mtext("Rocky steppe", side=2, line=3, font=2)}
  if(i == 2){ mtext("Grass steppe", side=2, line=3, font=2)}
  if(i == 3){ mtext("Meadow steppe", side=2, line=3, font=2)}

  axis(2, at=c(4921020, 6924035), labels = c("40°N", "60°N"))
  if(i == 3) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  ##Mid-Holocene------------------------------------------------
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod2[[i]], add=T, breaks=brks, col=cols, legend=F)
  image(ice[[2]], add=T, col="white", legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  #add no of species
  rs <- stack(list.files("./Source areas mid-Holocene", full.names = T))
  if(i == 1) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$rocky.steppe)])]]}
  if(i == 2) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$grass.steppe)])]]}
  if(i == 3) {rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[!is.na(clas$meadow.steppe)])]]}
  maxValue(rs)
  rs <- rs[[which(maxValue(rs) == 1)]]
  n <- nlayers(rs)
  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Source areas\nin the mid-Holocene", side=3, line=1, font=2)
  }
  if(i == 3) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  ##legend
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  
  ##add legend
  r <- raster(matrix(0:max(brks)))
  addRasterLegend(r, direction = "vertical", location = c(0.2,0.4,0.05,0.95), 
                  side=4, cex.axis = 1, breaks=brks, ramp = cols, 
                  ncolors = length(cols), locs = brks[-1])
  if(i == 1) {mtext("No. sp.", side=3, line=0, cex=0.8)}
  
  ##present distr
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  cgrs.EqEarth$div <- NA
  
  if(i == 1) {div <- rowSums(cgrs.range[, !is.na(clas$rocky.steppe)])}
  if(i == 2) {div <- rowSums(cgrs.range[, !is.na(clas$grass.steppe)])}
  if(i == 3) {div <- rowSums(cgrs.range[, !is.na(clas$meadow.steppe)])}
  
  cgrs.EqEarth$div <- div[cgrs.EqEarth$CGRSNAME]
  cgrs.EqEarth$div[is.na(cgrs.EqEarth$div)] <- 0
  max(div)
  
  if(i == 1) {brks <- seq(-2,16,2)}
  if(i == 2) {brks <- seq(-3,30,3)}
  if(i == 3) {brks <- seq(-4,40,4)}
  # if(i == 4) {brks <- seq(-1,5,1)}
  
  nb <- length(brks)-1 
  cols <- c("#F2F2F2", rev(magma(nb-1)))
  
  cgrs.EqEarth$div.clas <- cut(cgrs.EqEarth$div, breaks=brks, labels=seq(1, length(brks)-1), include.lowest = FALSE, right=TRUE) 
  
  plot(m, add=T, col="#F2F2F2", border=NA)
  plot(cgrs.EqEarth, add=T, col=cols[cgrs.EqEarth$div.clas], border=NA)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  ##add no of species
  if(i ==1) {n <- sum(!is.na(clas$rocky.steppe))}
  if(i ==2) {n <- sum(!is.na(clas$grass.steppe))}
  if(i ==3) {n <- sum(!is.na(clas$meadow.steppe))}

  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Present distribution", side=3, line=1, font=2)
  }
  if(i == 3) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  #legend
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  r <- raster(matrix(0:max(brks)))
  addRasterLegend(r, direction = "vertical", location = c(0.2,0.4,0.05,0.95), 
                  side=4, cex.axis = 1, breaks=brks, ramp = cols, 
                  ncolors = length(cols), locs = brks[-1])
  if(i == 1) {mtext("No. sp.", side=3, line=0, cex=0.8)}
}

dev.off()

###FIGURE 4-----------------------------------------------------------------

matr <- matrix(data=NA, ncol=14, nrow = 5, byrow=T)
matr[1,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)
matr[2,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+5
matr[3,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+10
matr[4,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+15
matr[5,1:14] <- c(1,1,1,1,2,2,2,2,3,4,4,4,4,5)+20

##Source areas in the LGM
mod <- stack(paste0("./Stacked source ares LGM/Chorotype Ward-Expert/Chorotype ", 1:5, ".sdat"))
mod <- projectRaster(mod, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                     res=c(4900,4900), method="ngb")
# plot(mod)

##Source areas in the mid-Holocene
mod2 <- stack(paste0("./Stacked source areas mid-Holocene/Chorotype Ward-Expert/Chorotype ", 1:5, ".sdat"))
mod2 <- projectRaster(mod2, crs="+proj=eqearth +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                      res=c(4900,4900), method="ngb")
# plot(mod2)

maxValue(mod)
maxValue(mod2)
max.brk <- c(20, 20, 24, 22, 3)

#plot
tiff("Fig_4.tif", width = 10, height = 9.5, units = "in", res=500, compression = "lzw")
layout(matr)
# layout.show(20)
par(mar=c(0.4,0.4,0.4,0.4), oma=c(2,5,4,0.5), tcl=-0.3, mgp=c(3, 0.5, 0))

for(i in 1:nlayers(mod))
{
  if(i == 1) {brks <- seq(-2,max.brk[i],2)}
  if(i == 2) {brks <- seq(-2,max.brk[i],2)}
  if(i == 3) {brks <- seq(-2,max.brk[i],2)}
  if(i == 4) {brks <- seq(-2,max.brk[i],2)}
  if(i == 5) {brks <- seq(-1,max.brk[i],1)}
  nb <- length(brks)-1 
  cols <- c("#F2F2F2", rev(viridis(nb-1)))
  
  ##LGM
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod[[i]], add=T, breaks=brks, col=cols, legend=F)
  image(ice[[3]], add=T, col="white", legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  #add no of species
  rs <- stack(list.files("./Source areas LGM", full.names = T))
  rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[clas$chorotype.ward.expert == i])]]
  maxValue(rs)
  rs <- rs[[which(maxValue(rs) == 1)]]
  n <- nlayers(rs)
  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Source areas\nin the Last Glacial Maximum", side=3, line=1, font=2)
  }
  if(i == 1){ mtext("Temperate\ncontinental species", side=2, line=3, font=2, cex=0.7)}
  if(i == 2){ mtext("Subatlantic-Submediterranean\n-continental species", side=2, line=3, font=2, cex=0.7)}
  if(i == 3){ mtext("Central European\nand Submediterranean species", side=2, line=3, font=2, cex=0.7)}
  if(i == 4){ mtext("Broad-range species", side=2, line=3, font=2, cex=0.7)}
  if(i == 5){ mtext("Narrow-range species", side=2, line=3, font=2, cex=0.7)}
  
  axis(2, at=c(4921020, 6924035), labels = c("40°N", "60°N"))
  if(i == 5) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  ##Mid-Holocene--------------------------
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  image(mod2[[i]], add=T, breaks=brks, col=cols, legend=F)
  image(ice[[2]], add=T, col="white", legend=F)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  #add no of species
  rs <- stack(list.files("./Source areas mid-Holocene", full.names = T))
  rs <- rs[[which(gsub("_", " ", names(rs)) %in% clas$species[clas$chorotype.ward.expert == i])]]
  maxValue(rs)
  rs <- rs[[which(maxValue(rs) == 1)]]
  n <- nlayers(rs)
  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Source areas\nin the mid-Holocene", side=3, line=1, font=2)
  }
  if(i == 5) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  ##legend
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  
  ##add legend
  r <- raster(matrix(0:max(brks)))
  addRasterLegend(r, direction = "vertical", location = c(0.2,0.4,0.05,0.95), 
                  side=4, cex.axis = 1, breaks=brks, ramp = cols, 
                  ncolors = length(cols), locs = brks[-1])
  if(i == 1) {mtext("No. sp.", side=3, line=0, cex=0.8)}
  
  ##present distr
  plot(NA, type="n", xlim=c(-2550000,3210000), ylim=c(4042000,7810000), xaxs="i", yaxs="i", axes=F, asp=1)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(169,211,236, maxColorValue = 255), border = NA)
  plot(g, add=T, lwd=0.5, col="white")
  
  cgrs.EqEarth$div <- NA
  
  if(i == 1) {div <- rowSums(cgrs.range[, clas$chorotype.ward.expert == 1])}
  if(i == 2) {div <- rowSums(cgrs.range[, clas$chorotype.ward.expert == 2])}
  if(i == 3) {div <- rowSums(cgrs.range[, clas$chorotype.ward.expert == 3])}
  if(i == 4) {div <- rowSums(cgrs.range[, clas$chorotype.ward.expert == 4])}
  if(i == 5) {div <- rowSums(cgrs.range[, clas$chorotype.ward.expert == 5])}
  
  cgrs.EqEarth$div <- div[cgrs.EqEarth$CGRSNAME]
  cgrs.EqEarth$div[is.na(cgrs.EqEarth$div)] <- 0
  
  
  if(i == 1) {brks <- seq(-2,22,2)}
  if(i == 2) {brks <- seq(-2,24,2)}
  if(i == 3) {brks <- seq(-3,27,3)}
  if(i == 4) {brks <- seq(-2,24,2)}
  if(i == 5) {brks <- seq(-1,9,1)}
  
  nb <- length(brks)-1 
  cols <- c("#F2F2F2", rev(magma(nb-1)))
  
  cgrs.EqEarth$div.clas <- cut(cgrs.EqEarth$div, breaks=brks, labels=seq(1, length(brks)-1), include.lowest = FALSE, right=TRUE) 
  
  plot(m, add=T, col="#F2F2F2", border=NA)
  plot(cgrs.EqEarth, add=T, col=cols[cgrs.EqEarth$div.clas], border=NA)
  
  plot(m, add=T, lwd=0.4)
  box()
  
  ##add no of species
  n <- sum(clas$chorotype.ward.expert == i)
  text(-2300000, 7500000, paste0("N = ", n), pos=4)
  
  if(i == 1)
  {
    mtext("Present distribution", side=3, line=1, font=2)
  }
  if(i == 5) {axis(1, at=c(-1773396, 0, 1773396), labels = c("0°", "20°E", "40°E"))}
  
  #legend
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  r <- raster(matrix(0:max(brks)))
  addRasterLegend(r, direction = "vertical", location = c(0.2,0.4,0.05,0.95), 
                  side=4, cex.axis = 1, breaks=brks, ramp = cols, 
                  ncolors = length(cols), locs = brks[-1])
  if(i == 1) {mtext("No. sp.", side=3, line=0, cex=0.8)}
  
}

dev.off()

