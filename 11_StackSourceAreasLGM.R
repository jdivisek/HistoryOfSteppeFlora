###############################################################################
###                         STACK LGM SOURCE AREAS                          ###
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

head(clas)

##Meadow steppe-------------------------------------
table(clas$meadow.steppe)

mod.paths <- paste0(getwd(), "/Source areas LGM/", clas$species, ".tif")
names(mod.paths) <- clas$species
mod.paths <- mod.paths[!is.na(clas$meadow.steppe)]

f <- list.files("./Source areas LGM/")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas LGM/", "Meadow steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)

##Grass steppe-------------------------------------
table(clas$grass.steppe)

mod.paths <- paste0(getwd(), "/Source areas LGM/", clas$species, ".tif")
names(mod.paths) <- clas$species
mod.paths <- mod.paths[!is.na(clas$grass.steppe)]

f <- list.files("./Source areas LGM/")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas LGM/", "Grass steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)

##Rocky steppe-------------------------------------
table(clas$rocky.steppe)

mod.paths <- paste0(getwd(), "/Source areas LGM/", clas$species, ".tif")
names(mod.paths) <- clas$species
mod.paths <- mod.paths[!is.na(clas$rocky.steppe)]

f <- list.files("./Source areas LGM/")
mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
length(mod.paths)

rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                            RESULT=paste0(getwd(), "/Stacked source areas LGM/", "Rocky steppe", ".sgrd")), 
                   env=saga.env, show.output.on.console=F)

##Chorotypes-----------------------------------

for(q in 1:5)
{
  print(paste0("*****Chorotype ", q))
  
  mod.paths <- paste0(getwd(), "/Source areas LGM/", clas$species, ".tif")
  names(mod.paths) <- clas$species
  mod.paths <- mod.paths[clas$chorotype == q]
  
  f <- list.files("./Source areas LGM/")
  mod.paths <- mod.paths[names(mod.paths) %in% gsub(".tif", "", f, fixed = T)]
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked source areas LGM/Chorotype/Chorotype ", q, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
  
}
