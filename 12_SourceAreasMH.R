###############################################################################
###      SIMULATE MIGRATION FROM SUITABLE AREAS IN THE MID-HOLOCENE         ###
###############################################################################

library(kissmig)
library(raster)
library(rgdal)
library(sp)
library(RSAGA)
library(rlist)

saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env$workspace <- paste0(getwd(),"/temp")

# rsaga.get.libraries(path = saga.env$modules)
# rsaga.get.modules("shapes_grid", env = saga.env)
# rsaga.get.usage("shapes_grid", 2, env = saga.env)

##load grids with the extent of ice sheets, Caspian and Aral sea in selected periods
##all data must be in LAEA projection (+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs)
##set periods (in years before present)
tw <- c(500, 2500, 4500, 6500, 8500, 9500, 12500, 14500, 17500, 20500)

ice.paths <- paste0("./Grids LAEA/Ice sheets/ice.sheet_", tw, "BP.shp")#grids
names(ice.paths) <- proj.periods
casp.paths <- paste0("./Grids LAEA/Caspian_sea/casp.sea_", tw, "BP.sdat")#grids
names(casp.paths) <- proj.periods
aral.paths <- paste0("./Grids LAEA/Aral_sea/aral.sea_", tw, "BP.sdat")#grids
names(aral.paths) <- proj.periods

#load continental ice sheets
ice <- stack(ice.paths)
projection(ice) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#load caspian sea
casp <- stack(casp.paths)
projection(casp) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#load aral sea
aral <- stack(aral.paths)
projection(aral) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257222 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##land for masking
mh.land <- raster("./Grids LAEA/Mid-Holocene land/mh.land.tif")

####SIMULATE MID-HOLOCENE SOURCE AREAS (POTENTIAL REFUGIA) FOR ALL SPECIES-----------------------------------------

acces.test.MH  <- list()

for(q in names(me))
{
  print(q)
  
  mod.paths <- paste0("./Projected models LAEA for simulations/", proj.periods, "/", q, ".tif")
  rs <- stack(rev(mod.paths))
  # plot(rs)
  
  ###get threshold
  th <- me[[q]]@results["X10.percentile.training.presence.Cloglog.threshold", 1]
  
  ###get all climatically suitable areas for LGM
  refs <- getRef(rs[[1]], threshold = th, min.size = 42) #min size in No cells in one spetial cluster. 416.493 corresponds to 10 000 km2 with a raster resolution of 4900 m
  # plot(refs)
  
  ##square suitabilities
  rs0 <- rs[[1:7]]
  for(i in 1:nlayers(rs0)) {
    r <- rs0[[i]]/max(values(rs0[[i]]), na.rm = T)
    rs0[[i]] <- r*r} #rescaled to get max = 1
  
  ##IDENTIFY SOURCE AREAS IN THE MID-HOLOCENE USING BEST MIGRATION RATE
  ##only those that are climatically suitable and accessible from LGM
  if(!is.null(refs) & !is.na(best.mig.rate[q]))
  {
    ##calculate accessibility between LGM and MH
    for(i in seq(1, nlayers(refs)))
    {
      k <- kissmig(refs[[i]], rs0, it=best.mig.rate[q], type='FOC', seed=123)
      k <- kissmigAccess(k, rel = TRUE)
      k[k > 0] <- 1
      suppressWarnings(writeRaster(k, filename=paste0("./temp/", names(refs)[i], ".tif"), format="GTiff", overwrite=T))
    }
    k <- suppressWarnings(stack(list.files("./temp", pattern = ".tif", full.names = T)))
    if(nlayers(k) == 1){k <- k[[1]]}
    if(nlayers(k) > 1){k <- calc(k, fun=sum)}
    k[k > 0] <- 1
    file.remove(list.files("./temp", pattern = ".tif", full.names = T))
    
    ###get all potential source areas for MH
    refs <- getRef(rs[[7]], threshold = th, min.size = 42) #min size in No cells in one spetial cluster. 416.493 corresponds to 10 000 km2 with a raster resolution of 4900 m
    # plot(refs, maxpixels=3000000)
    
    ####update source areas - only those accessible from the LGM considered
    if(!is.null(refs))
    {
      for(i in 1:nlayers(refs))
      {
        refs[[i]] <- refs[[i]]*k
      }
      
      if(any(maxValue(refs) == 0)){
        refs <- dropLayer(refs, which(maxValue(refs) == 0))}
      
      if(nlayers(refs) > 0){
        names(refs) <- paste0("Ref", 1:nlayers(refs))
      }
      if(nlayers(refs) == 0){
        refs <- NULL
      }
    }
  }
  
  if(!is.null(refs))
  {
    refs <- calc(refs, sum)
    refs <- getRef(refs, threshold = 0.5, min.size = 42)
  }
  
  ###TEST SOURCE AREAS IN THE MID-HOLOCENE
  if(is.null(refs) | is.na(best.mig.rate[q])){
    acces.test.MH[[q]] <- NA}
  if(!is.null(refs) & !is.na(best.mig.rate[q]))
  {
    #save MH potential source areas
    if(nlayers(refs) == 1){writeRaster(refs[[1]], filename = paste0("./Accessible areas in mid-Holocene/", q, ".tif"), format="GTiff", overwrite=T)}
    if(nlayers(refs) > 1){calc(refs*c(1:nlayers(refs)), fun=sum, filename = paste0("./Accessible areas in mid-Holocene/", q, ".tif"), format="GTiff", overwrite=T)}
    
    ###TEST MIGRATION RATES---------------------------------
    
    rates <- best.mig.rate[q]
    
    ##square suitabilities
    for(i in 1:nlayers(rs)) {
      r <- rs[[i]]/max(values(rs[[i]]), na.rm = T)
      rs[[i]] <- r*r} #rescaled to get max = 1
    
    print(rates)
    
    ###calculate accessibility from each potential source area
    discart <- rep(FALSE, nlayers(refs))
    for(i in seq(1, nlayers(refs)))
    {
      k <- kissmig(refs[[i]], rs[[7:10]], it=rates, type='FOC', seed=123)
      k <- kissmigAccess(k, rel = FALSE)
      
      if(all(is.na(values(k)))){discart[i] <- TRUE}
      else{suppressWarnings(writeRaster(k, filename=paste0("./temp/", names(refs)[i], ".tif"), format="GTiff", overwrite=T))}
    }
    if(any(discart)) {refs <- dropLayer(refs, which(discart == TRUE))}

    ###extract suitability & accessibility values
    rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(c(mod.paths[1], paste0(names(refs), '.tif')), collapse = "; "),
                                              POLYGONS=paste0(getwd(), "/Shapes LAEA/CGRS grid/cgrs_sel2_LAEA.shp"),
                                              METHOD="0",
                                              COUNT="0", MIN="0", MAX="0", RANGE="0", SUM="0", MEAN="1", VAR="0", STDDEV="0", QUANTILE="0", GINI="0",
                                              RESULT="cgrs_env.shp"), 
                       env=saga.env, show.output.on.console=F)
    
    cgrs <- foreign::read.dbf("./temp/cgrs_env.dbf", as.is = TRUE)
    
    ##prepare data for GLM
    glm.data <- cgrs[,-c(1,2)]
    rownames(glm.data) <- cgrs$CGRSNAME
    colnames(glm.data)[1] <- "suitability"
    colnames(glm.data) <- gsub("..MEAN.", "", colnames(glm.data))
    colnames(glm.data) <- gsub("..MEAN", "", colnames(glm.data))
    
    glm.data <- cbind(glm.data, glm.env.data)
    
    ##join species range
    glm.data$species <- 0
    glm.data[rownames(cgrs.range)[cgrs.range[, q] == 1], "species"] <- 1
    glm.data <- glm.data[complete.cases(glm.data), ]
    
    ###RUN GLM MODELS
    acces.test.MH[[q]] <- glm.test2(data=glm.data, y="species", keep.in=colnames(glm.env.data), 
                                     test=names(refs), family="binomial")
    
    
    list.save(glm.data, paste0("./Source areas mid-Holocene/", q, "_TestData.rdata"))
    write.table(acces.test.MH[[q]], paste0("./Source areas mid-Holocene/", q, "_test.txt"), sep="\t", dec=".")
    file.remove(list.files("./temp", full.names = T))
    
    ##find best migration rate
    test <- acces.test.MH[[q]]
    sum.imp <- sum(test$D.diff[test$Estimate > 0 & test$`Pr(>|z|)` < 0.05], na.rm=T)#
    
    if(any(sum.imp > 0))
    {
      
      pdf(paste0("./Source areas mid-Holocene/", q, ".pdf"), 6, 8)
      par(mfrow=c(2,1), mar=c(0.5,0.5,3,0.5), oma=rep(0,4))
      
      eref <- calc(refs, sum)
      plot(eref, axes=F, box=F, legend=T, main=paste0("All suitable ares in the mid-Holocene for ", q))
      plot(country, add=T, lwd=0.3)
      
      if(sum(c(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05), na.rm = T) > 0)
      {
        eref <- calc(refs[[which(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05)]], sum)
        plot(eref, axes=F, box=F, legend=T, main=paste0("Source areas in the mid-Holocene for ", q))
        plot(country, add=T, lwd=0.3)
      }
      else
      {
        eref[!is.na(eref)] <- 0
        plot(eref, axes=F, box=F, legend=T, main=paste0("Source areas in the mid-Holocene for ", q))
        plot(country, add=T, lwd=0.3)
        
        acces.test.MH[[q]] <- NA
      }
      
      dev.off()
      
      eref <- eref*mh.land
      suppressWarnings(writeRaster(eref, filename=paste0("./Source areas mid-Holocene/", q, ".tif"), format="GTiff", overwrite=T))
      
    }
    else
    {
      eref <- rs[[7]]##7 is Mid-Holocene
      eref[!is.na(eref)] <- 0
      suppressWarnings(writeRaster(eref, filename=paste0("./Source areas mid-Holocene/", q, ".tif"), format="GTiff", overwrite=T))
      
      acces.test.MH[[q]] <- NA
    }
    
    list.save(acces.test.MH, "acces.test.MH.rdata")
    
    rm(rs, refs, eref)
  }
  
}


###Get size and the number of source areas in the mid-Holocene-----------------------------------------
mh.sa.stat <- as.data.frame(matrix(data = NA, nrow=length(me), ncol=4, byrow = F))
rownames(mh.sa.stat) <- names(me)
colnames(mh.sa.stat) <- c('Size.acc', 'No.acc', 'Size.sig', 'No.sig')

mod.paths <- list.files("./Source areas mid-Holocene", pattern="tif", full.names = T)
names(mod.paths) <- gsub(".tif", "", list.files("./Source areas mid-Holocene", pattern="tif", full.names = T), fixed = T)

for(q in names(me))
{
  print(q)
  
  if(q %in% names(mod.paths))
  {
    r <- raster(mod.paths[q])
    
    x <- sum(values(r) == 1, na.rm=T)
    
    mh.sa.stat[q, "Size.sig"] <- x
    
    if(x > 0)
    {
      eref <- raster::clump(r, directions=8, gaps=FALSE)
      mh.sa.stat[q, "No.sig"] <- maxValue(eref)
    }
    else
    {
      mh.sa.stat[q, "No.sig"] <- 0
    }
  }
}

###Get size and the number of accessible areas-----------------------------------------

mod.paths <- list.files("./Accessible areas in mid-Holocene", full.names = T)
names(mod.paths) <- gsub(".tif", "", list.files("./Accessible areas in mid-Holocene", full.names = T), fixed=T)

for(q in names(me))
{
  print(q)
  
  if(q %in% names(mod.paths))
  {
    r <- raster(mod.paths[q])
    
    x <- sum(values(r) > 0, na.rm=T)
    
    mh.sa.stat[q, "Size.acc"] <- x
    
    if(x > 0)
    {
      mh.sa.stat[q, "No.acc"] <- maxValue(r)
    }
    else
    {
      mh.sa.stat[q, "No.acc"] <- 0
    }
  }
}
