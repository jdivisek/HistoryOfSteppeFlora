###############################################################################
###                                STACK MODELS                             ###
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

###SPECIALISTS OF ECOLOGICAL TYPES-----------------------------------------------------------------
head(clas)

##Meadow steppe----------------------------------------------------------------------
table(clas$meadow.steppe)

for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
  mod.paths <- mod.paths[!is.na(clas$meadow.steppe)]
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models LAEA/Meadow steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

##Grass steppe-----------------------------------------------------------------------
table(clas$grass.steppe)

for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
  mod.paths <- mod.paths[!is.na(clas$grass.steppe)]
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models LAEA/Grass steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

##Rocky steppe-----------------------------------------------------------------------
table(clas$rocky.steppe)

for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
  mod.paths <- mod.paths[!is.na(clas$rocky.steppe)]
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models LAEA/Rocky steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

###CHROROTYPES---------------------------------------------------------------------
table(clas$chorotype)

for(q in 1:5)
{
  print(paste0("*****Chorotype ", q))
  
  for(i in proj.periods)
  {
    print(i)
    mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
    mod.paths <- mod.paths[clas$chorotype == q]
    
    rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                                RESULT=paste0(getwd(), "/Stacked models LAEA/Chorotype/Chorotype ", q, "/", i, ".sgrd")), 
                       env=saga.env, show.output.on.console=F)
  }
}

###ALL SPECIES EXCEPT DESERT STEPPE-------------------------------------------------

for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
  mod.paths <- mod.paths[-c(34,54)]#Ephedra distachya a Krascheninnikovia ceratoides
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models LAEA/All species except Desert steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

###DESERT STEPPE------------------------------------------------
for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models LAEA/", i, "/", clas$species, ".sgrd")
  mod.paths <- mod.paths[c(34,54)]#Ephedra distachya a Krascheninnikovia ceratoides
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models LAEA/Desert steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

####STACK WGS84 MODELS FOR MAIN TEXT----------------------------------------------
###All species except desert steppe-----------------------------------------

for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models/", i, "/", clas$species, ".tif")
  mod.paths <- mod.paths[-c(34,54)]#Ephedra distachya a Krascheninnikovia ceratoides
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models/All species except Desert steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

###Desert steppe--------------------------------------
for(i in proj.periods)
{
  print(i)
  mod.paths <- paste0(getwd(), "/Projected models/", i, "/", clas$species, ".tif")
  mod.paths <- mod.paths[c(34,54)]#Ephedra distachya a Krascheninnikovia ceratoides
  
  rsaga.geoprocessor("grid_calculus", 8, list(GRIDS=paste0(mod.paths, collapse = "; "),
                                              RESULT=paste0(getwd(), "/Stacked models/Desert steppe/", i, ".sgrd")), 
                     env=saga.env, show.output.on.console=F)
}

