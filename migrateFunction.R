#Function that identifies significant source areas using migration simulations and accessibility tests

#@data = list including habitat suitability thresholds (for delineation of suitable areas) and paths to suitability grids for each species
#@lgm.land = grid of LGM land extent (land = 1, ocean = NA)
#@rates = vector of migration rates - seq(10,200,10)
#@temp.folder = path to temporary folder
#@out.folder = path to a folder where grids of significant source areas will be stored
#@grid.path = path to UTM 50x50 grid (shapefile)
#@saga.env = SAGA environment
#@glm.env.data = a data.frame with UTM grids in rows and environmental variables in columns. 
#@cgrs.range = a data.frame with species occurrences (presence/absence) in UTM grid cells
#@country = shapefile of countries

migrate <- function(data, lgm.land, rates, temp.folder, out.folder, grid.path, saga.env, glm.env.data, cgrs.range, country)
{
  require(kissmig)
  require(raster)
  require(rgdal)
  require(sp)
  require(RSAGA)
  require(rlist)
  require(foreign)
  
  rs <- suppressWarnings(stack(rev(data[[2]])))
  
  sn <- names(data[[1]])
  dir.create(paste0(temp.folder, "/", sn))
  
  ###get all climatically suitable areas (potential refugia)
  refs <- getRef(rs[[1]], threshold = data[[1]], min.size = 42) #min size = No cells in one spatial cluster. 42 roughly corresponds to 1000 km2 with a raster resolution of 4900 m
  
  
  if(is.null(refs)){return(NA)}
  if(!is.null(refs))
  {
    ###TEST MIGRATION RATES---------------------------------
    
    test.data <- list()
    acces.test <- list()
    
    ##square suitabilities
    for(i in 1:nlayers(rs)) {
      r <- rs[[i]]/max(values(rs[[i]]), na.rm = T)
      rs[[i]] <- r*r} #rescaled to get max = 1
    
    for(w in 1:length(rates))
    {
      
      ###calculate accessibility from each potential source area
      for(i in seq(1, nlayers(refs)))
      {
        k <- kissmig(refs[[i]], rs, it=rates[w], type='FOC', seed=123)
        k <- kissmigAccess(k, rel = TRUE)
        suppressWarnings(writeRaster(k, filename=paste0(paste0(temp.folder, "/", sn), "/", names(refs)[i], ".tif"), format="GTiff", overwrite=T))
      }
      
      ###extract suitability & accessibility values
      rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(c(data[[2]][1], paste0(temp.folder, "/", sn, "/", names(refs), '.tif')), collapse = "; "),
                                                POLYGONS=grid.path,
                                                METHOD="0",
                                                COUNT="0", MIN="0", MAX="0", RANGE="0", SUM="0", MEAN="1", VAR="0", STDDEV="0", QUANTILE=0, GINI="0",
                                                RESULT=paste0(temp.folder, "/", sn, "/cgrs_env.shp")), 
                         env=saga.env, show.output.on.console=F)
      
      cgrs <- foreign::read.dbf(paste0(temp.folder, "/", sn, "/cgrs_env.dbf"), as.is = TRUE)
      
      ##prepare data for GLM
      glm.data <- cgrs[,-c(1,2)]
      rownames(glm.data) <- cgrs$CGRSNAME
      colnames(glm.data)[1] <- "suitability"
      colnames(glm.data) <- gsub("..MEAN.", "", colnames(glm.data))
      colnames(glm.data) <- gsub("..MEAN", "", colnames(glm.data))
      
      glm.data <- cbind(glm.data, glm.env.data)
      
      ##join species range
      glm.data$species <- 0
      glm.data[rownames(cgrs.range)[cgrs.range[, sn] == 1], "species"] <- 1
      glm.data <- glm.data[complete.cases(glm.data), ]
      
      test.data[[w]] <- glm.data
      
      ###RUN GLM MODELS
      acces.test[[w]] <- glm.test(data=glm.data, y="species", keep.in=colnames(glm.env.data), 
                                  test=names(refs), family="binomial")
      
      file.remove(list.files(paste0(temp.folder, "/", sn), full.names = T))
      
    }
    
    list.save(test.data, paste0(temp.folder, "/", sn, "_TestData.rdata"))
    list.save(acces.test, paste0(temp.folder, "/", sn, "_AccesTests.rdata"))
    
    ##find best migration rate
    sum.imp <- unlist(lapply(acces.test, FUN=function(x){sum(x$D.diff[x$Estimate > 0 & x$`Pr(>|z|)` < 0.05], na.rm=T)}))#
    sum.imp[is.nan(sum.imp)] <- 0
    
    if(any(sum.imp > 0))
    {
      bmr <- rates[which.max(sum.imp)]
      
      ###identify significant source areas
      test <- acces.test[[which.max(sum.imp)]]
      
      pdf(paste0(temp.folder, "/", sn, ".pdf"), 6, 8)
      par(mfrow=c(2,1), mar=c(0.5,0.5,3,0.5), oma=rep(0,4))
      
      eref <- calc(refs, sum)
      plot(eref, axes=F, box=F, legend=T, main=paste0("All climatically suitable areas for ", sn))
      plot(country, add=T, lwd=0.3)
      
      if(sum(c(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05), na.rm = T) > 0)
      {
        eref <- calc(refs[[which(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05)]], sum)
        plot(eref, axes=F, box=F, legend=T, main=paste0("Source areas for ", sn))
        plot(country, add=T, lwd=0.3)
      }
      else
      {
        eref[!is.na(eref)] <- 0
        plot(eref, axes=F, box=F, legend=T, main=paste0("Source areas for ", sn))
        plot(country, add=T, lwd=0.3)
      }
      
      dev.off()
      
      eref <- eref*lgm.land
      suppressWarnings(writeRaster(eref, filename=paste0(out.folder, "/", sn, ".tif"), format="GTiff", overwrite=T))
      test$mig.rate <- bmr
      
      write.table(test, file=paste0(temp.folder, "/", sn, "_test.txt"), sep="\t", dec=".")
      return(test)
      
    }
    else
    {
      bmr <- NA
      
      eref <- rs[[1]]##1 is the Last Glacial Maximum
      eref[!is.na(eref)] <- 0
      suppressWarnings(writeRaster(eref, filename=paste0(out.folder, "/", sn, ".tif"), format="GTiff", overwrite=T))
      
      return(bmr)
    }
    
    rm(rs, refs, eref, k)
  }
  
}