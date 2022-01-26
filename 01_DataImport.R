###############################################################################
###                              DATA IMPORT                                ###
###############################################################################

library(sp)
library(rgdal)
library(raster)
library(RSAGA)
library(parallel)
library(pbapply)

###LOAD DATA------------------------------------------------------------------------------
###species distribution data
spe.paths <- list.files("path_to_folder_with_shapefiles", pattern="shp", full.names = T)
spe.paths <- spe.paths[seq(1,207, 2)]
spe.paths <- sort(spe.paths)

###species classification
clas.spe <- read.delim("Species_FIN_2021-11-19.txt", header = T, as.is = T)
clas.spe <- clas.spe[order(clas.spe$abbrev),]
head(clas.spe)

names(spe.paths) <- clas$species

###make list of coordinates
spe.coord <- list()

for(i in clas$species)
{
  sp <- readOGR(spe.paths[i])
  spe.coord[[i]] <- coordinates(sp)[,1:2]
  colnames(spe.coord[[i]]) <- c("X", "Y")
}

###load climatic variables for Palearctic
env <- stack(list.files("./Palearctic_current_2.5arcmin/", full.names = T))

##load random background points
rp20 <- read.delim("Palearctic_rp20.txt", header = T, row.names = 1) ##generated in ArcGIS in North Pole Lambert Azimuthal Equal Area Projection 
rp20 <- rp20[complete.cases(raster::extract(env, rp20)),]

#make correlation matrix
cors <- cor(rp20, method="pearson")
cors[1:5,1:5]

#make absolute values
cors <- abs(cors)

#make "distances"
cors <- 1-cors

##hierarchical clustering
cl <- hclust(as.dist(cors), method = "complete")
plot(cl, hang = -1, axes = F, ylab = "Pearson correlation")
axis(2, at=seq(0,1, 0.2), labels=seq(1,0,-0.2))
rect.hclust(cl, h=0.3)

#variables with Pearson correlation less than 0.7
env.names <- c("bio1", "bio2", "bio4", "bio8", "bio10", "bio12", "bio15", "bio16", "aridityIndexThornthwaite", "PETseasonality", "tri")
env <- stack(paste0("./Palearctic_current_2.5arcmin/", env.names, ".tif"))
env
plot(env)

###EXTRACT CLIMATIC VARS------------------------------------------------------------------------
my.extract <- function(coord, env)
{
  ext <- as.data.frame(raster::extract(x=env, y=coord))
  return(ext)
}

detectCores()

cl <- makeCluster(12)

spe.clim <- pblapply(spe.coord, FUN = my.extract, cl = cl, env=env)
stopCluster(cl)

##remove points with NA
for(i in 1:length(spe.coord))
{
  spe.coord[[i]] <- spe.coord[[i]][complete.cases(spe.clim[[i]]),]
  spe.clim[[i]] <- spe.clim[[i]][complete.cases(spe.clim[[i]]),]
}



