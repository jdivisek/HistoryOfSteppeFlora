###############################################################################
###                             MAXENT MODELS                               ###
###############################################################################

#load required packages
library(rJava)
library(raster)
library(dismo)
library(parallel)
library(pbapply)

jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

rp20 <- read.delim("Palearctic_rp20.txt", header = T, row.names = 1)#randomly sampled background points
rp20 <- rp20[complete.cases(raster::extract(env, rp20)),]
nrow(rp20)

###FIT MAXENT MODELS------------------------------------------------------

my.maxent <- function(coord, env, bg, removeDuplicates, args){
  res <- dismo::maxent(x=env, p=coord, a=bg, removeDuplicates=removeDuplicates, args=args)
  return(res)
}

detectCores()

cl <- makeCluster(12)

me <- pblapply(spe.coord.f, FUN = my.maxent, cl = cl, env=env,
                bg=rp20, removeDuplicates=TRUE,
                args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                       "replicates=1", "replicatetype=crossvalidate"))

stopCluster(cl)

##add species names
for(q in names(me))
{
  colnames(me[[i]]@results) <- q
}

##EXPORT SUITABILITY RASTERS---------------------------------------------
my.export <- function(model, env, path)
{
  filename <- paste0(path, "/", colnames(model@results), ".tif")
  dismo::predict(model, x=env, filename=filename, format="GTiff", 
                 options="COMPRESS=DEFLATE PREDICTOR=2 ZLEVEL=6")
}

cl <- makeCluster(12)

pblapply(me, FUN = my.export, cl = cl, env=env, path="./models")

stopCluster(cl)
