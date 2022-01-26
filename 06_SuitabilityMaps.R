###############################################################################
###                    EXPORT HABITAT SUITABILITY MAPS                      ###
###############################################################################

library(raster)
library(rgdal)
library(stringr)

##load shapefiles and grids with the extent of ice sheets, Caspian and Aral sea in selected periods
##all data must be in LAEA projection (+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs)
##set periods (in years before present)
tw <- c(500, 2500, 4500, 6500, 8500, 9500, 12500, 14500, 17500, 20500)

ice.paths <- paste0("./Shapes LAEA/Ice sheets/ice.sheet_", tw, "BP.shp")#shapes
names(ice.paths) <- proj.periods
casp.paths <- paste0("./Grids LAEA/Caspian_sea/casp.sea_", tw, "BP.sdat")#grids
names(casp.paths) <- proj.periods
aral.paths <- paste0("./Grids LAEA/Aral_sea/aral.sea_", tw, "BP.sdat")#grids
names(aral.paths) <- proj.periods

#load continental ice sheets
ice <- list()
for(i in proj.periods){ice[[i]] <- readOGR(ice.paths[i])}

#load caspian sea
casp <- stack(casp.paths)
projection(casp) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(casp)[!is.na(values(casp))] <- 1

#load aral sea
aral <- stack(aral.paths)
projection(aral) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(aral)[!is.na(values(aral))] <- 1

##load countries and graticules
country <- readOGR("./Shapes LAEA/Country/worldmap.shp")
grat <- readOGR("./Shapes LAEA/Graticules/grat.shp")
label.coords <- read.delim("./Shapes LAEA/Graticules/coord_labels.txt", header=T, row.names = 1)

#set breaks for habitat suitability classification
brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
cols <- rev(terrain.colors(nb))

for(q in names(me))
{
  print(paste0("********************", q, "********************"))
  
  pdf(file = paste0("./Suitability maps/", q, ".pdf"), width=11, height = 7.5, onefile = T)
  par(fig = c(0,1,0,1), mar=c(0,0.5,1,3))
  
  mod.paths <- paste0("./Projected models LAEA/", proj.periods, "/", q, ".sdat")
  
  for(i in seq(length(proj.periods),1,-1))
  {
    print(proj.periods[i])
    
    ###load model and mask Capian and Aral sea
    mod <- raster(mod.paths[i])
    mod <- mask(mod, mask=casp[[i]], maskvalue=1)
    mod <- mask(mod, mask=aral[[i]], maskvalue=1)
    
    ##plot model
    plot(frame, col=rgb(169,211,236, maxColorValue = 255), border=NA, add=F)
    plot(grat, lwd=0.3, col="white", add=T)
    plot(mod, breaks=brks, col=cols, lab.breaks=brks, asp=1, axes=F, box=F, add=T)
    plot(ice[[i]], lwd=0.5, col="white", border="white", add=T)
    plot(country, add=T, lwd=0.3)
    
    ##add species name 
    mtext(paste0(q, " | ", proj.periods[i]), side=3, line = -2, outer=T, cex=1.2, font=2)
    plot(frame, col=NA, add=T)
    
    #plot evaluation scores
    text(x=600000, y=7400000, paste0("cvAUC = ", round(me.cv[q, "mean.test.auc"], 3), " ± ", round(me.cv[q, "sd.test.auc"], 3)), pos=4, xpd=T)
    text(x=600000, y=7200000, paste0("overfitting = ", round(me.cv[q, "mean.overfitting"], 3), " ± ", round(me.cv[q, "sd.overfitting"], 3)), pos=4, xpd=T)
    text(x=600000, y=7000000, paste0("cvBoyce = ", round(me.cv[q, "mean.test.boyce"], 3), " ± ", round(me.cv[q, "sd.test.boyce"], 3)), pos=4, xpd=T)

    #plot variable contributions
    contr <- me[[q]]@results[str_detect(rownames(me[[q]]@results), "contribution"),]#selected variables
    contr <- sort(contr, decreasing = T)
    names(contr) <- gsub(".contribution", "", names(contr), fixed = T)
    text(x=7400000, y=7400000, paste(names(contr), collapse="\n"), xpd=T, adj=c(0,1))
    text(x=9100000, y=7400000, paste(paste0(round(contr, 1), "%"), collapse="\n"), xpd=T, adj=c(0,1))
    
    ##coordinate labels
    for(w in 1:10){
      if(w < 5) {text(label.coords[w,1:2], label.coords[w,3], pos=1, xpd=NA, offset=0.5, cex=0.7)}
      if(w == 5) {text(label.coords[w,1:2], paste0("   ", label.coords[w,3]), pos=1, xpd=NA, offset=0.5, cex=0.7)}
      if(w > 5) {text(label.coords[w,1:2], label.coords[w,3], pos=2, xpd=NA, offset=0.3, cex=0.7)}
    }
    
  }
  dev.off()
  
}
removeTmpFiles(0)
