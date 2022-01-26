###############################################################################
###                        REGIONS OF SPECIES ORIGIN                        ###
###############################################################################
###WHERE SPECIES CAME FROM?

library(vegan)
library(RSAGA)
library(rgdal)

regs <- readOGR("./Shapes LAEA/Regions/refugial_regions.shp")
plot(regs)

saga.env <- rsaga.env()
saga.env$path <- "C:/saga-6.3.0_x64"
saga.env$modules <- "C:/saga-6.3.0_x64/tools"
saga.env$version <- "6.3.0"
saga.env$cores <- 5
saga.env$workspace <- paste0(getwd(), "/temp")

###LGM---------------------------------------------------------------------

mod.paths <- list.files("./Source areas LGM", full.names = T)
names(mod.paths) <- gsub(".tif", "", list.files("./Source areas LGM", full.names = T), fixed = T)

###assign species occurrence in predefined regions (e.g. Pannonian Basin)
rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(mod.paths, collapse = "; "),
                                          POLYGONS=paste0(getwd(), "/Shapes LAEA/Regions/refugial_regions.shp",
                                          METHOD="0",
                                          COUNT="0", MIN="0", MAX="0", RANGE="0", SUM="1", MEAN="0", VAR="0", STDDEV="0", QUANTILE=0, GINI="0",
                                          RESULT="refugial_regions_FIN.shp"), 
                   env=saga.env, show.output.on.console=T)

regs <- readOGR("./temp/refugial_regions_FIN.shp")

regs <- regs@data[,-c(1:3)]
head(regs)

colnames(regs)[-1] <- names(mod.paths)

rownames(regs) <- regs$Name
regs <- t(regs[,-1])
head(regs)

regs.pa <- decostand(regs, "pa")
regs.pa

###MID-HOLOCENE-------------------------------------------------------------
mod.paths <- list.files("./Source areas mid-Holocene", full.names = T)
names(mod.paths) <- gsub(".tif", "", list.files("./Source areas mid-Holocene", full.names = T), fixed = T)

###assign species occurrence in predefined regions (e.g. Pannonian Basin)
rsaga.geoprocessor("shapes_grid", 2, list(GRIDS=paste(mod.paths, collapse = "; "),
                                          POLYGONS=paste0(getwd(), "/Shapes LAEA/Regions/refugial_regions.shp"),
                                          METHOD="0",
                                          COUNT="0", MIN="0", MAX="0", RANGE="0", SUM="1", MEAN="0", VAR="0", STDDEV="0", QUANTILE=0, GINI="0",
                                          RESULT="refugial_regions_FIN.shp"), 
                   env=saga.env, show.output.on.console=T)

regs.MH <- readOGR("./temp/refugial_regions_FIN.shp")

regs.MH <- regs.MH@data[,-c(1:3)]
head(regs.MH)

colnames(regs.MH)[-1] <- names(mod.paths)

rownames(regs.MH) <- regs.MH$Name
regs.MH <- t(regs.MH[,-1])
head(regs.MH)

regs.MH.pa <- decostand(regs.MH, "pa")