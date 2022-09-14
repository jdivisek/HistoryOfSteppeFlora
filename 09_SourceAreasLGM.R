###############################################################################
###             SIMULATE MIGRATION FROM SUITABLE AREAS IN THE LGM           ###
###############################################################################

##simulates migration from each isolated climatically suitable area in the LGM larger than
##1000 km2 and test the effect of its accessibility on current species distribution
##in Europe

library(pbapply)
library(parallel)
library(kissmig)
library(raster)
library(rgdal)
library(sp)
library(RSAGA)
library(rlist)
library(foreign)

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

##mask suitability maps
for (q in names(me))
{
  print(q)
  
  for(i in 1:length(proj.periods))
  {
    r <- suppressWarnings(raster(paste0(getwd(), "/Projected models LAEA/", proj.periods[i], "/", q, ".sdat")))
    
    r <- mask(r, ice[[i]], inverse=T, maskvalue=NA, updatevalue=0, updateNA=T)
    r <- mask(r, casp[[i]], inverse=T, maskvalue=NA)
    r <- mask(r, aral[[i]], inverse=T, maskvalue=NA)
    
    suppressWarnings(writeRaster(r, file=paste0(getwd(), "/Projected models LAEA for simulations/", proj.periods[i], "/", q, ".tif"), format="GTiff", overwrite=T))
  }
}

##land for masking
lgm.land <- raster("./Grids LAEA/LGM land/lgm.land.tif")

##import grid with asigned species ranges
cgrs.range <- read.delim("spec.UTM.WIDE_cgrs.txt", header = T, row.names = 1)
head(cgrs.range)

colnames(cgrs.range) == clas$abbrev
colnames(cgrs.range) <- clas$species

####SIMULATE MIGRATION----------------------------------------------------------------

mod.paths <- list()
for(q in names(me))
{
  mod.paths[[q]] <- list()
  
  #10th percetile threshold
  mod.paths[[q]][[1]] <- me[[q]]@results["X10.percentile.training.presence.Cloglog.threshold", 1]
  names(mod.paths[[q]][[1]]) <- q
  
  #paths to suitability grids
  mod.paths[[q]][[2]] <- paste0(getwd(), "/Projected models LAEA for simulations/", proj.periods, "/", q, ".tif")
  
}

#set SAGA environment
saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env$workspace <- paste0(getwd(),"/temp")

detectCores()

cl <- makeCluster(12)
clusterExport(cl, list("getRef", "glm.test"))
clusterSetRNGStream(cl,123)

acces.test <- pblapply(mod.paths, FUN = migrate, cl = cl, lgm.land = lgm.land,
                       rates = seq(10, 200, 10),
                       temp.folder = paste0(getwd(), "/temp"),
                       out.folder = paste0(getwd(), "/Source areas LGM"),
                       grid.path = paste0(getwd(), "/Shapes LAEA/CGRS grid/cgrs_sel2_LAEA.shp"),
                       saga.env = saga.env,
                       glm.env.data = glm.env.data, #values of selected climatic variables and TRI for UTM grid cells
                       cgrs.range= cgrs.range, #species presence(1) and absence(0) records in UTM grid cells
                       country = country)



stopCluster(cl)

removeTmpFiles(0)

###exctract best migration rate
best.mig.rate <- unlist(lapply(acces.test, FUN=function(x){ifelse(all(is.na(x)),NA, x$mig.rate[1])}))
summary(best.mig.rate)

###Get size and the number of all climatically suitable areas in the LGM-----------------------------------------

lgm.sa.stat <- as.data.frame(matrix(data = NA, nrow=length(me), ncol=4, byrow = F))
rownames(lgm.sa.stat) <- names(me)
colnames(lgm.sa.stat) <- c('Size.all', 'No.all', 'Size.sig', 'No.sig')
for(q in names(me))
{
  print(q)
  mod.paths <- paste0(getwd(), "/Projected models LAEA for simulations/", "21k BP", "/", q, ".sdat")
  r <- raster(mod.paths)
  
  ###get threshold
  th <- me[[q]]@results["X10.percentile.training.presence.Cloglog.threshold", 1]
  
  ###get all potential refugia
  refs <- getRef(r, threshold = th, min.size = 42) ##42 cells ~ 1000 km2
  
  if(!is.null(refs))
  {
    lgm.sa.stat[q, "No.all"] <- nlayers(refs)
    
    r <- calc(refs, sum)
    lgm.sa.stat[q, "Size.all"] <- sum(values(r) == 1, na.rm=T)
  }
  
}

###Get size and the number of source areas in the LGM-----------------------------------------

mod.paths <- list.files(paste0(getwd(), "/Source areas LGM"), full.names = T)
names(mod.paths) <- gsub(".tif", "", list.files(paste0(getwd(), "/Source areas LGM"), full.names = T), fixed=T)

for(q in names(me))
{
  print(q)
  
  if(q %in% names(mod.paths))
  {
    r <- raster(mod.paths[q])
    
    x <- sum(values(r) == 1, na.rm=T)
    
    lgm.sa.stat[q, "Size.sig"] <- x
    
    if(x > 0)
    {
      eref <- raster::clump(r, directions=8, gaps=FALSE)
      lgm.sa.stat[q, "No.sig"] <- maxValue(eref)
    }
    else
    {
      lgm.sa.stat[q, "No.sig"] <- 0
    }
  }
}

