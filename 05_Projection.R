###############################################################################
###                     PROJECTION OF MAXENT MODELS                         ###
###############################################################################

#load required packages
library(rgdal)
library(raster)
library(sp)
library(dismo)
library(rJava)
library(ENMeval)
library(parallel)
library(maptools)
library(vegan)
library(usdm)
library(stringr)
library(ecospat)
library(pbapply)
library(RSAGA)

###PROJECT MODELS TO PAST CLIMATES-------------------------------------------------
##set paths to palaeoclimatic grids at 2.5 arc-min resolution
proj.periods <- c("present", "2.1k BP", "4.2k BP", "6.2k BP", "8.2k BP", "9.9k BP", "12.2k BP", "14.3k BP", "17.1k BP", "21k BP")

proj.paths <- list()

for(i in proj.periods){
  proj.paths[[i]] <- paste0("path_to_paleoclimatic_grids", i, "/", names(env), ".tif")
}


##export rasters---------------------------------------------------------------------
my.export <- function(model, env, path)
{
  filename <- paste0(path, "/", colnames(model@results), ".tif")
  dismo::predict(model, x=env, filename=filename, format="GTiff", 
                 options="COMPRESS=DEFLATE PREDICTOR=2 ZLEVEL=6")
}


for(i in proj.periods)
{
  print(i)
  
  rs <- stack(proj.paths[[i]])#load env variables for past climates
  
  cl <- makeCluster(12)
  
  pblapply(me, FUN = my.export, cl = cl, env=rs, path="./temp")#project models
  
  stopCluster(cl)
  
  file.rename(from=list.files("./temp", full.names = T),
              to=paste0("./Projected models/", i, "/", list.files("./temp")))
}

###PROJECT GRIDS TO LAEA FOR MAPPING------------------------------------------------------
# rsaga.get.libraries(path = saga.env$modules)
# rsaga.get.modules("pj_proj4", env = saga.env)
# rsaga.get.usage("pj_proj4", 3, env = saga.env)

saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env$workspace <- paste0(getwd(), "/temp")
saga.env

for(i in proj.periods)
{
  print(i)
  
  mod.paths <- list.files(paste0("./Projected models/", i), full.names = T)
  
  rsaga.geoprocessor("pj_proj4", 3, list(CRS_PROJ4="+proj=laea +lon_0=10 +lat_0=52 +x_0=4321000 +y_0=3210000 +units=m +a=6378137 +b=6356752.31414 +lat_1=0 +lat_2=0 +no_defs", 
                                         SOURCE=paste0(mod.paths, collapse = "; "),
                                         GRIDS=paste0(gsub(".tif", ".sgrd", gsub("Projected models", "Projected models LAEA", mod.paths, fixed = T), fixed = T), collapse = "; "),
                                         TARGET_DEFINITION="0",
                                         TARGET_USER_SIZE=4900), 
                     env=saga.env, show.output.on.console=F)
}



