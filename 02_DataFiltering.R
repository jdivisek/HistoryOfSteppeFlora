###############################################################################
###                            DATA FILTERING                               ###
###############################################################################

library(rworldmap)
library(rgeos)
library(raster)
library(sp)

###Filtering with standardized variables----------------------------------------------------------
m <- getMap(resolution = "coarse", projection = NA)

sort(unlist(lapply(spe.coord, nrow)))#get numbers of occurrence points

spe.coord.f <- list()
filtering.info <- list(stats = list(), used.axes = list())

##plots results in maps
pdf("Filtering.pdf", 8.5, 9, onefile = T)
layout(matrix(c(1,1,1,2,2,2), ncol=1))
par(mar=c(5,5,0,1), oma=c(0,0,2,0))

for(i in names(spe.coord))#
{
  print(paste0("********************* ", i ," *********************"))
  
  if(nrow(spe.coord[[i]])> 100)
  {
    ##transform coordinates
    p <- SpatialPoints(spe.coord[[i]][,1:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    p <- spTransform(p, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
    
    PCA <- as.data.frame(scale(cbind(spe.clim[[i]], coordinates(p))))##No PCA analysis here - raw variables used
    rm(p)
    
    Result <- PCAfilter(OccurenceData=spe.coord[[i]][,1:2], 
                        PCAscores=PCA, 
                        Resolution=0.3, 
                        Seed=1234, 
                        Priority = spe.coord[[i]]$priority,
                        OriginalPoints = TRUE,
                        SupplInfo = TRUE)
    
    spe.coord.f[[i]] <- Result$filtered[,c("X", "Y", "priority")]
    filtering.info$stats[[i]] <- Result$suppl$Records
    filtering.info$sources$Original[[i]] <- table(spe.coord[[i]]$priority)
    filtering.info$sources$Filtered[[i]] <- table(Result$filtered$priority)
    
    ###Plot filtering
    #Selected points
    plot(Result$original[,c("X.1", "Y.1")], type="n", xlab="Longitude", ylab="Latitude", asp=0)
    plot(m, add=T, col="gray85")
    points(Result$original[,c("X.1", "Y.1")])
    points(Result$filtered[,c("X", "Y")], pch=16, col="green", cex=0.5)
    box()
    
    mtext(paste(i, "| Original:", Result$suppl$Records$Original, "| Filtered:",
                Result$suppl$Records$Filtered, "| Omitted:", Result$suppl$Records$OmittedData, sep=" "),
          outer=T, line=0.6, side=3, font=2)
    legend("bottomleft", pch=c(1,16), col=c("black", "green"), legend = c("Original", "Filtered"))
    
    #Omited points
    plot(Result$original[,c("X.1", "Y.1")], type="n", xlab="Longitude", ylab="Latitude", asp=0)
    plot(m, add=T, col="gray85")
    points(Result$original[,c("X.1", "Y.1")])
    points(Result$omitted[,c("X", "Y")], pch=16, col="red", cex=0.5)
    box()
    
    legend("bottomleft", pch=c(1,16), col=c("black", "red"), legend = c("Original", "Omitted"))
    
    ###Plot sources of data
    plot(spe.coord[[i]][,1:2], type="n", xlab="Longitude", ylab="Latitude", asp=0)
    plot(m, add=T, col="gray85")
    points(spe.coord[[i]][,1:2], col=c("red3", "blue", "black")[spe.coord[[i]]$priority])
    box()
    legend("bottomleft", pch=rep(1,3), col=c("red3", "blue", "black"), legend = c("True", "GBIF - PRESERVED_SPECIMEN", "random"))
    
    mtext(paste(i, "| Original:", Result$suppl$Records$Original, "| Filtered:",
                Result$suppl$Records$Filtered, "| Omitted:", Result$suppl$Records$OmittedData, sep=" "),
          outer=T, line=0.6, side=3, font=2)
    
    plot(Result$filtered[,c("X", "Y")], type="n", xlab="Longitude", ylab="Latitude", asp=0)
    plot(m, add=T, col="gray85")
    points(Result$filtered[,c("X", "Y")], col=c("red3", "blue", "black")[Result$filtered$priority])
    box()
    legend("bottomleft", pch=rep(1,3), col=c("red3", "blue", "black"), legend = c("True", "GBIF - PRESERVED_SPECIMEN", "random"))
    
    rm(Result, PCA)
  }
  
  if(nrow(spe.coord[[i]]) <= 100)
  {
    spe.coord.f[[i]] <- as.data.frame(spe.coord[[i]])
    filtering.info$stats[[i]] <- list(Original=nrow(spe.coord[[i]]), 
                                      Filtered=nrow(spe.coord[[i]]),
                                      OmittedData=0)
    filtering.info$used.axes[[i]] <- NA
  }
}

dev.off()

##store filtering statistics
no.records <- do.call(rbind, lapply(filtering.info$stats, unlist))
write.table(no.records, "no.records.txt", sep='\t', dec=".")
