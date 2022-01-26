###############################################################################
###             MAPS OF THE LGM SOURCE AREAS (POTENTIAL REFUGIA)            ###
###############################################################################

library(kissmig)
library(raster)
library(rgdal)
library(sp)
library(RSAGA)
library(rlist)

##land for masking
lgm.land <- raster("./Grids LAEA/LGM land/lgm.land.tif")

##load countries and graticules
country <- readOGR("./Shapes LAEA/Country/worldmap.shp")
grat <- readOGR("./Shapes LAEA/Graticules/grat.shp")
label.coords <- read.delim("./Shapes LAEA/Graticules/coord_labels.txt", header=T, row.names = 1)

##load UTM grid
cgrs <- readOGR("./Shapes LAEA/CGRS grid/cgrs_sel2_LAEA.shp", verbose=T) 
cgrs$species <- 0

##create vector of colors
max(unlist(lapply(acces.test, FUN=function(x){nrow(x)})))
my.cols <- c(rgb(94, 171, 182, maxColorValue = 255),
             rgb(165, 72, 136, maxColorValue = 255),
             rgb(88, 49, 87, maxColorValue = 255),
             rgb(209, 122, 23, maxColorValue = 255),
             rgb(247, 194, 38, maxColorValue = 255),
             "red3", "saddlebrown", "pink3", "dimgrey", "deepskyblue3",
             "aquamarine2", "seagreen", "darkolivegreen3", "yellow", "maroon1")


###plot LGM source areas for each species
for(q in names(acces.test))
{
  print(q)
  
  if(!is.na(acces.test[[q]])[1])
  {
    pdf(file = paste0("./Maps of LGM source areas/", q, ".pdf"), width=11, height = 7.5, onefile = T)
    par(fig = c(0,1,0,1), mar=c(0,0.5,1,3))
    
    mod.paths <- paste0("./Projected models LAEA for simulations/", "21k BP", "/", q, ".sdat")
    r <- suppressWarnings(raster(mod.paths))
    
    ###get threshold
    th <- me[[q]]@results["X10.percentile.training.presence.Cloglog.threshold", 1]

    ###get all potential refugia
    refs <- getRef(r, threshold = th, min.size = 42)
    # plot(refs)
    
    ####PLOT ALL POTENTIAL REFUGIA----------------------------------------------
    grps <- suppressWarnings(split(1:nlayers(refs) ,cut(1:61, c(-Inf, 15, 30, 45, 60, Inf), 1:5)))
    
    ###repeat plotting for "groups" of 15 refugia---------------------------------
    for(g in 1:length(grps))
    {
      if(length(grps[[g]]) > 0)
      {
        plot(frame, col=rgb(169,211,236, maxColorValue = 255), border=NA, add=F)
        plot(grat, lwd=0.3, col="white", add=T)
        plot(lgm.land, add=T, legend=F, col="gray90")
        plot(ice[[1]], add=T, legend=F, col="white")
        
        for(i in grps[[g]])
        {
          plot(refs[[i]], add=T, legend=F, col=c(NA, rep(my.cols, 5)[i]), maxpixels=3000000)
        }
        
        plot(country, add=T, lwd=0.3)
        mtext(paste("All suitable areas in the LGM for", q), side=3, line = -2, outer=T, cex=1.2, font=2)
        plot(frame, col=NA, add=T)
        
        ##coordinate labels
        for(w in 1:10){
          if(w < 5) {text(label.coords[w,1:2], label.coords[w,3], pos=1, xpd=NA, offset=0.5, cex=0.7)}
          if(w == 5) {text(label.coords[w,1:2], paste0("   ", label.coords[w,3]), pos=1, xpd=NA, offset=0.5, cex=0.7)}
          if(w > 5) {text(label.coords[w,1:2], label.coords[w,3], pos=2, xpd=NA, offset=0.3, cex=0.7)}
        }
        
        legend(x=0.8e+06, y=7.2e+06, pch=15, col=rep(my.cols, 5)[grps[[g]]],
               legend=paste0(grps[[g]], "."), bty="n", pt.cex=1.3, xpd=NA, cex=0.8)
        
        if(!is.na(best.mig.rate[q])){
          glm.tab <- as.data.frame(round(acces.test[[q]], 3))
          
          legend(x=7.2e+06, y=7.2e+06, ncol = 4L, legend = c(c("", paste0(grps[[g]], ".")),
                                                             "D.diff", glm.tab[grps[[g]],6],
                                                             'Coefficient', glm.tab[grps[[g]],1],
                                                             'p-value', glm.tab[grps[[g]],4]), x.intersp=0, cex=0.8, xpd=NA, bty="o", bg="white", box.col="white", box.lwd = 0.2)
          
          text(x=8e+06, y=1.5e+06, paste("Migration rate:", round(best.mig.rate[q]*4900/2100, 1), "m/year"), pos=4, xpd=NA, cex=0.8)
        }
        else{
          text(x=8e+06, y=1.5e+06, paste("Migration rate:", "NA"), pos=4, xpd=NA, cex=0.8)
        }
        
      }
    }
    
    ###PLOT SIGNIFICANT SOURCE AREAS-------------------------------------------------
    plot(frame, col=rgb(169,211,236, maxColorValue = 255), border=NA, add=F)
    plot(grat, lwd=0.3, col="white", add=T)
    plot(lgm.land, add=T, legend=F, col="gray90")
    plot(ice[[1]], add=T, legend=F, col="white")
    
    if(!is.na(best.mig.rate[q])){
      eref <- raster(paste0("./Source areas LGM/", q, ".tif"))
      plot(eref, add=T, legend=F, col=c(NA, "#00A600"), maxpixels=3000000)}
    
    cgrs$species <- 0
    cgrs$species[which(cgrs$CGRSNAME %in% rownames(cgrs.range)[cgrs.range[,q] == 1])] <- 1
    cgrs.plot <- cgrs[cgrs$species == 1,]
    
    plot(country, add=T, lwd=0.3, col="gray60")
    plot(cgrs.plot, add=T, density = 15, angle=-45, border=NA, col="black", lwd=0.3)
    
    mtext(paste("LGM source areas for", q), side=3, line = -2, outer=T, cex=1.2, font=2)
    plot(frame, col=NA, add=T)
    
    ##coordinate labels
    for(w in 1:10){
      if(w < 5) {text(label.coords[w,1:2], label.coords[w,3], pos=1, xpd=NA, offset=0.5, cex=0.7)}
      if(w == 5) {text(label.coords[w,1:2], paste0("   ", label.coords[w,3]), pos=1, xpd=NA, offset=0.5, cex=0.7)}
      if(w > 5) {text(label.coords[w,1:2], label.coords[w,3], pos=2, xpd=NA, offset=0.3, cex=0.7)}
    }
    
    if(!is.na(best.mig.rate[q])) {text(x=8e+06, y=1.5e+06, paste("Migration rate:", round(best.mig.rate[q]*4900/2100, 1), "m/year"), pos=4, xpd=NA, cex=0.8)}
    if(is.na(best.mig.rate[q])) {text(x=8e+06, y=1.5e+06, paste("Migration rate:", "NA"), pos=4, xpd=NA, cex=0.8)}
    
    polygon(c(0.8e+06, 1e+06, 1e+06, 0.8e+06),
            c(7.0e+06, 7.0e+06, 7.2e+06, 7.2e+06),
            border=NA, col="#00A600", lwd=0.3)
    text(1e+06, 7.1e+06, "Source areas", cex=0.8, pos=4)
    polygon(c(0.8e+06, 1e+06, 1e+06, 0.8e+06),
            c(6.7e+06, 6.7e+06, 6.9e+06, 6.9e+06),
            density = 15, angle=-45, border=NA, col="black", lwd=0.3)
    text(1e+06, 6.8e+06, "Present distribution", cex=0.8, pos=4)
    
    dev.off()
  }
}


