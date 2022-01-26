###############################################################################
###                            CROSS-VALIDATION                             ###
###############################################################################

###VALIDATION----------------------------------
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
library(blockCV)

##Width of species distribution ranges----------------------------------------
range <- list()
for(i in names(spe.coord))
{
  range[[i]] <- max(apply(apply(spe.coord[[i]], 2, range),2, diff))
}

sort(unlist(range))
sort(unlist(lapply(spe.coord.f, nrow)))

##Spatial blocks for cross-validation-----------------------------------------
spe.coord.f.cv <- spe.coord.f

for(i in names(spe.coord.f))
{
  set.seed(1234)
  spe.coord.f.cv[[i]]$folds <- spatialBlock(speciesData = SpatialPoints(spe.coord.f[[i]], CRS("+proj=longlat +datum=WGS84")), # sf or SpatialPoints
                                             species = NULL, # the response column (binomial or multi-class)
                                             rasterLayer = env[[1]], # a raster for background (optional)
                                             theRange = round((range[[i]]*71698)*0.1, 0), # size of the blocks in meters
                                             k = 5, # number of folds
                                             selection = "random",
                                             iteration = 100, # find evenly dispersed folds
                                             degMetre = 71698,
                                             biomod2Format = TRUE,
                                             showBlocks = TRUE,
                                             progress = FALSE,
                                             verbose = FALSE)$foldID
}


###CROSS-VALIDATION--------------------------------------------------
rp20.env <- as.data.frame(raster::extract(env, rp20))
head(rp20.env)

detectCores()

cl <- makeCluster(12)

me.cv <- pblapply(spe.coord.f.cv, FUN = cv.maxent, cl = cl, env=env,
                   bg=rp20, bg.env=rp20.env, removeDuplicates=TRUE,
                   maxent.args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1"),
                   boyce.res=100, print.progress=FALSE)

stopCluster(cl)

me.cv <- do.call(rbind.data.frame, me.cv)

summary(me.cv)

