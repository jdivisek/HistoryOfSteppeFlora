###############################################################################
###                    STACK MID-HOLOCENE SOURCE AREAS                      ###
###############################################################################

#load required packages
library(raster)
library(stringr)
library(RSAGA)

saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env

# rsaga.get.libraries(path = saga.env$modules)
# rsaga.get.modules("grid_calculus", env = saga.env)
# rsaga.get.usage("grid_calculus", 8, env = saga.env)

head(clas.spe)

##Meadow steppe-------------------------------------
table(clas.spe$meadow.steppe)

mod.paths <- paste0(getwd(), "/Source areas mid-Holocene/", clas.spe$species, ".tif")
names(mod.paths) <- clas.spe$species
mod.paths <- mod.paths[!is.na(clas.spe$meadow.steppe)]

f <- list.files("./Source areas mid-Holocene")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas mid-Holocene/", "Meadow steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)

##Grass steppe-------------------------------------
table(clas.spe$grass.steppe)

mod.paths <- paste0(getwd(), "/Source areas mid-Holocene/", clas.spe$species, ".tif")
names(mod.paths) <- clas.spe$species
mod.paths <- mod.paths[!is.na(clas.spe$grass.steppe)]

f <- list.files("./Source areas mid-Holocene")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas mid-Holocene/", "Grass steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)


##Rocky steppe-------------------------------------
table(clas.spe$rocky.steppe)

mod.paths <- paste0(getwd(), "/Source areas mid-Holocene/", clas.spe$species, ".tif")
names(mod.paths) <- clas.spe$species
mod.paths <- mod.paths[!is.na(clas.spe$rocky.steppe)]

f <- list.files("./Source areas mid-Holocene")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas mid-Holocene/", "Rocky steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)

##Chorotypes-----------------------------------

for(q in 1:5)
{
  print(paste0("*****Chorotype ", q))
  
  mod.paths <- paste0(getwd(), "/Source areas mid-Holocene/", clas.spe$species, ".tif")
  names(mod.paths) <- clas.spe$species
  mod.paths <- mod.paths[clas.spe$chorotype == q]
  
  f <- list.files("./Source areas mid-Holocene")
  mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "./Stacked source areas mid-Holocene/Chorotype/Chorotype ", q, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
  
}