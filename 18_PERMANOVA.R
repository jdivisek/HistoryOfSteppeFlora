###############################################################################
###                 DIFFERENCES IN LOCATION OF SOURCE AREAS                 ###
###############################################################################

library(raster)
library(vegan)
library(fuzzySim)
library(RSAGA)
library(foreign)
library(rgdal)
library(ggrepel)
library(rlist)

saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env$workspace <- paste0(getwd(), "/temp")

##species classification into chorotypes
chorotype <- clas$chorotype.ward.expert
names(chorotype) <- clas$species

##species classification into ecological groups
ecotype <- rep(NA, nrow(clas))
names(ecotype) <- clas$species
ecotype[clas$species[!is.na(clas$meadow.steppe)]] <- 1 #"Meadow steppe"
ecotype[clas$species[!is.na(clas$grass.steppe)]] <- 2 #"Grass steppe"
ecotype[clas$species[!is.na(clas$rocky.steppe)]] <- 3 #"Rocky steppe"


##paths to models
mod.paths <- paste0(getwd(),"/Source areas LGM/", clas$species, ".tif")
names(mod.paths) <- clas$species

f <- list.files("./Source areas LGM")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
chorotype <- chorotype[names(chorotype) %in% gsub(".tif", "", f, fixed = T)]
ecotype <- ecotype[names(ecotype) %in% gsub(".tif", "", f, fixed = T)]
cbind(names(chorotype), names(ecotype), names(mod.paths))

###assign source areas to UTM 50x50 km grid
rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(mod.paths, collapse = "; "),
                                          POLYGONS=paste0(getwd(), "/Shapes LAEA/CGRS grid/cgrs_sel2_LAEA_LGM_noIce.shp"),
                                          METHOD="0",
                                          COUNT="0", MIN="0", MAX="1", RANGE="0", SUM="0", MEAN="0", VAR="0", STDDEV="0", QUANTILE=0, GINI="0",
                                          RESULT="cgrs_ext.shp"),
                   env=saga.env, show.output.on.console=F)

cgrs <- readOGR("./temp/cgrs_ext.shp")
cgrs.coord <- as.data.frame(coordinates(cgrs))
colnames(cgrs.coord) <- c("Lon", "Lat")
p <- SpatialPointsDataFrame(cgrs.coord, cgrs.coord, proj4string = crs(cgrs))
p <- spTransform(p, CRS("+proj=longlat +datum=WGS84"))
cgrs.coord <- as.data.frame(round(coordinates(p),3))
colnames(cgrs.coord) <- c("Lon", "Lat")
rownames(cgrs.coord) <- rownames(cgrs)

cgrs <- cgrs@data
rownames(cgrs) <- cgrs$CGRSNAME
head(cgrs)
cgrs <- cgrs[, -c(1:3)]
colnames(cgrs) <- names(mod.paths) 

cgrs.coord <- cgrs.coord[rowSums(cgrs) > 0,]
cgrs <- cgrs[rowSums(cgrs) > 0,]

cgrs <- as.data.frame(t(cgrs))
cgrs <- cgrs[rowSums(cgrs)>0,]
nrow(cgrs)

chorotype <- chorotype[rownames(cgrs)]
ecotype <- ecotype[rownames(cgrs)]
cbind(names(chorotype), names(ecotype), rownames(cgrs))

#Permutational Multivariate Analysis of Variance Using Distance Matrices
set.seed(123)
test <- adonis2(cgrs ~ chorotype, method="bray", sqrt.dist=T, by = NULL, binary=TRUE, permutations=9999)
set.seed(111)
test <- adonis2(cgrs[names(ecotype)[!is.na(ecotype)],] ~ ecotype[!is.na(ecotype)], method="bray", sqrt.dist=T, by = NULL, binary=TRUE, permutations=9999)


##Mid-Holocene-------------------------------------------------------------------

##species classification into chorotypes
chorotype <- clas$chorotype
names(chorotype) <- clas$species

chorotype <- clas$chorotype.ward.expert
names(chorotype) <- clas$species

##species classification into ecological groups
ecotype <- rep(NA, nrow(clas))
names(ecotype) <- clas$species
ecotype[clas$species[!is.na(clas$meadow.steppe)]] <- 1 #"Meadow steppe"
ecotype[clas$species[!is.na(clas$grass.steppe)]] <- 2 #"Grass steppe"
ecotype[clas$species[!is.na(clas$rocky.steppe)]] <- 3 #"Rocky steppe"

##paths to models
mod.paths <- paste0("./Source areas mid-Holocene/", clas$species, ".tif")
names(mod.paths) <- clas$species

f <- list.files("./Source areas mid-Holocene/")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
chorotype <- chorotype[names(chorotype) %in% gsub(".tif", "", f, fixed = T)]
ecotype <- ecotype[names(ecotype) %in% gsub(".tif", "", f, fixed = T)]
cbind(names(chorotype), names(ecotype), names(mod.paths))


###assign source areas to UTM 50x50 km grid
rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(mod.paths, collapse = "; "),
                                          POLYGONS= paste0(getwd(), "/Shapes LAEA/CGRS grid/cgrs_sel2_LAEA_MHsel.shp",
                                                           METHOD="0",
                                                           COUNT="0", MIN="0", MAX="1", RANGE="0", SUM="0", MEAN="0", VAR="0", STDDEV="0", QUANTILE=0, GINI="0",
                                                           RESULT="cgrs_ext.shp"),
                                          env=saga.env, show.output.on.console=F)

cgrs <- readOGR("./temp/cgrs_ext.shp")
cgrs.coord <- as.data.frame(coordinates(cgrs))
colnames(cgrs.coord) <- c("Lon", "Lat")
p <- SpatialPointsDataFrame(cgrs.coord, cgrs.coord, proj4string = crs(cgrs))
p <- spTransform(p, CRS("+proj=longlat +datum=WGS84"))
cgrs.coord <- as.data.frame(round(coordinates(p),3))
colnames(cgrs.coord) <- c("Lon", "Lat")
rownames(cgrs.coord) <- rownames(cgrs)

cgrs <- cgrs@data
rownames(cgrs) <- cgrs$CGRSNAME
head(cgrs)
cgrs <- cgrs[, -c(1:2)]
colnames(cgrs) <- names(mod.paths) 

cgrs.coord <- cgrs.coord[rowSums(cgrs) > 0,]
cgrs <- cgrs[rowSums(cgrs) > 0,]

cgrs <- as.data.frame(t(cgrs))
cgrs[1:5,1:5]
cgrs <- cgrs[rowSums(cgrs)>0,]
nrow(cgrs)

chorotype <- chorotype[rownames(cgrs)]
ecotype <- ecotype[rownames(cgrs)]
cbind(names(chorotype), names(ecotype), rownames(cgrs))

#Permutational Multivariate Analysis of Variance Using Distance Matrices
set.seed(123)
test <- adonis2(cgrs ~ chorotype, method="bray", sqrt.dist=T, by = NULL, binary=TRUE, permutations=9999)
set.seed(111)
test <- adonis2(cgrs[names(ecotype)[!is.na(ecotype)],] ~ ecotype[!is.na(ecotype)], method="bray", sqrt.dist=T, by = NULL, binary=TRUE, permutations=9999)

