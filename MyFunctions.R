##DATA FILTERING IN GRIDDED ENVIRONMENTAL SPACE----------------------------------------------------
##sensu Varela et al. 2014
##inspired by Míša Vojtěchovských

#@OccurenceData - X (Longitude) and Y (Latitude) coordinates of species occurrence points
#@Scores - standardized environmental variables (mean = 0, sd = 1) or PCA scores
#@Resolution - resolution environmental space (in SD units), the higher value, the more points are deleted
#@Seed - seed for random number generation
#@Priority - numeric vector of priority values for selection: lower number = higher priority
#@OriginalPoints - store original points? Default = TRUE
#@SupplInfo - store supplementary data (No. of original, preserved and deleted points and gridding sekvence)? Default = TRUE

PCAfilter <- function(OccurenceData, 
                      Scores, 
                      Resolution, 
                      Seed, 
                      OriginalPoints = TRUE,
                      SupplInfo = TRUE) {
  
  #load packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load("raster")
  
  points.scores  <- as.data.frame(Scores)
  OccurenceData <- as.data.frame(OccurenceData)
  if (OriginalPoints == TRUE) {
    original.points.scores <- points.scores
  }
  
  ##find a maximum width of environmental space
  min.score <- min(apply(Scores, 2, min))
  max.score <- max(apply(Scores, 2, max)) 
  
  Sekvence <- c(rev(seq(0, min.score-Resolution, by=Resolution*-1)), seq(Resolution, max.score+Resolution, by=Resolution))
  
  ##find points belonging to the same grid cell
  for (i in 1:ncol(Scores)) {
    points.scores[paste0("ID", i)] <- as.numeric(cut(points.scores[,i],
                                                     breaks = Sekvence,
                                                     include.lowest = TRUE))
  }
  
  ##make group IDs
  for (i in 1:ncol(Scores)) {
    points.scores$Group <- paste0(points.scores$Group,
                                  points.scores[,ncol(Scores)+i])
  }
  
  #attach geographical coordinates and priority field
  points.scores$X <- OccurenceData[,1]
  points.scores$Y <- OccurenceData[,2]
  if(!is.null(Priority)){points.scores$priority <- Priority}
  
  #randomize data
  set.seed(Seed)
  points.scores <- points.scores[sample(nrow(points.scores)),]
  if(!is.null(Priority)) {points.scores <- points.scores[order(points.scores$priority, decreasing = FALSE),]}
  
  ##filtering
  NumbersOfColumnsToDrop <- c(1:(ncol(Scores)+1)+ncol(Scores))
  Result <- list()
  Result$filtered <- points.scores[!duplicated(points.scores$Group), -NumbersOfColumnsToDrop]
  Result$omitted <- points.scores[duplicated(points.scores$Group), -NumbersOfColumnsToDrop]
  
  #store original points?
  if (OriginalPoints == TRUE) {
    Result$original <- as.data.frame(append(original.points.scores, OccurenceData))
  }
  
  #store supplementary data?
  if (SupplInfo == TRUE) {
    Result$suppl <- list(Seed = Seed, 
                         Records = list(Original = nrow(OccurenceData), 
                                        Filtered = nrow(Result$filtered), 
                                        OmittedData = nrow(OccurenceData) - nrow(Result$filtered)), 
                         Grids = Sekvence)
  }
  
  return(Result)
}

###CROSS-VALIDATION WITH BOYCE INDEX--------------------------------------------------

#@p = presences (X, Y) an fold IDs in the 3rd column 
#@bg = random background points (X, Y)
#@bg.env = a data frame with values of env variables for background points
#@env = rasterStack of environmental variables
#@folds = vector with k folds pro presence points
#@maxent.args = arguments passed to maxent
#@boyce.res = resolution of boyce index (see help for ecospat::ecospat.boyce)

cv.maxent <- function(p, env, bg, bg.env, removeDuplicates=TRUE, 
                      maxent.args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1"),
                      boyce.res=100, print.progress=TRUE)
{
  require(raster)
  require(ecospat)
  require(dismo)
  require(stringi)
  
  #create temp directory
  if(!dir.exists("./tempMaxent"))
  {
    dir.create("./tempMaxent")
  }
  temp.folder.name <- stri_rand_strings(1, 10)
  dir.create(paste0("./tempMaxent/", temp.folder.name))
  
  #create vectors for storing results
  train.auc <- vector("numeric")
  test.auc <- vector("numeric")
  train.boyce <- vector("numeric")
  test.boyce <- vector("numeric")
  
  #get folds
  folds <- p[,3]
  
  #run cross-validation
  for(i in seq(1,max(folds)))
  {
    if(print.progress){print(paste("Cross-validation:", i ,"out of", max(folds), sep= " "))}
    
    ##split data
    train.p <- as.data.frame(p[folds != i, 1:2]) ##select training presences
    test.p <- as.data.frame(p[folds == i, 1:2]) ##select testing presences
    
    ##create SWD file with test points
    test.swd <- as.data.frame(raster::extract(env, test.p))
    test.swd <- data.frame(Species = "species", test.p, test.swd)
    colnames(test.swd)[2:3] <- c("longitude", "latitude")
    test.swd <- test.swd[complete.cases(test.swd),]#check missing value in predictor vars
    
    write.table(test.swd, file=paste0("./tempMaxent/", temp.folder.name, "/TestSWD.csv"), sep = ",", dec=".", row.names = F, quote = FALSE)
    
    ##run MaxEnt
    maxent.args <- c(maxent.args, paste0("testsamplesfile=", getwd(), "/tempMaxent/", temp.folder.name, "/TestSWD.csv"))
    mod <- dismo::maxent(x=env, p=train.p, a=bg, removeDuplicates=removeDuplicates, args=maxent.args)
    
    #extract aucs
    train.auc[i] <- mod@results["Training.AUC",1]
    test.auc[i] <- mod@results["Test.AUC",1]
    
    #calculate Boyce index
    pr.bg <- dismo::predict(mod, bg.env)
    pr.train <- dismo::predict(mod, as.data.frame(raster::extract(env, train.p)))
    pr.test <- dismo::predict(mod, test.swd[, -c(1:3)])
    train.boyce[i] <- ecospat::ecospat.boyce(pr.bg, pr.train, PEplot=F, res=boyce.res)$Spearman.cor
    test.boyce[i] <- ecospat::ecospat.boyce(pr.bg, pr.test, PEplot=F, res=boyce.res)$Spearman.cor
    
    # removeTmpFiles(0)
    file.remove(paste0("./tempMaxent/", temp.folder.name, "/TestSWD.csv"))
  }
  return(list(mean.train.auc=mean(train.auc, na.rm=T),
              sd.train.auc=sd(train.auc, na.rm=T),
              mean.test.auc=mean(test.auc, na.rm=T),
              sd.test.auc=sd(test.auc, na.rm=T),
              mean.overfitting=mean(train.auc-test.auc, na.rm=T),
              sd.overfitting=sd(train.auc-test.auc, na.rm=T),
              mean.train.boyce=mean(train.boyce, na.rm=T),
              sd.train.boyce=sd(train.boyce, na.rm=T),
              mean.test.boyce=mean(test.boyce, na.rm=T),
              sd.test.boyce=sd(test.boyce, na.rm=T),
              mean.boyce.diff=mean(train.boyce-test.boyce, na.rm=T),
              sd.boyce.diff=sd(train.boyce-test.boyce, na.rm=T),
              No.valid.train.boyce=sum(!is.na(train.boyce)),
              No.valid.test.boyce=sum(!is.na(test.boyce))))
  
  unlink(paste0("./tempMaxent/", temp.folder.name), recursive = TRUE, force = TRUE)
}

##GET ALL POTENTIAL SOURCE AREAS (refugia)-----------------------------------------------------------

#@r = climatic suitability raster produced by Maxent
#@threshold = threshold for defining suitable/unsuitable areas
#@min.size = minimum size (in number of pixels) of isolated climatically suitable areas
#@directions = directions for delineation of isolated areas (see help for raster::clump); 8 (default) = Queen scheme

getRef <- function(r, threshold, min.size, directions=8)
{
  require(raster)
  
  #make binary map
  pa <- r
  pa[pa > threshold] <- 1
  pa[pa != 1] <- 0
  
  if(any(values(pa) > 0, na.rm = T))
  {
    #find spatial clusters
    refs <- raster::clump(pa, directions=directions, gaps=FALSE)
    
    #get No. cells per cluster
    a <- table(values(refs))
    
    for(i in 1:length(a)) {if(a[i] < min.size){refs[refs == i] <- NA}}
    
    ##assign new values to single refugia
    a <- unique(values(refs))[-1]
    
    if(length(a) > 0)
    {
      refs <- reclassify(refs, rcl = matrix(c(a, seq(1, length(a))), ncol=2, byrow = F))
      
      ##assign 0 to NA values
      refs[is.na(refs)] <- 0
      
      #make raterstack with single refugia (origins)
      rs <- stack()
      
      for(i in unique(values(refs))[-1])
      {
        rs <- addLayer(rs, refs)
        
        rs[[i]][rs[[i]] != i] <- 0
        rs[[i]][rs[[i]] == i] <- 1
      }
      names(rs) <- paste0("Ref", 1:maxValue(refs))
      
      return(rs)
    }
    else
    {
      return(NULL)
    } 
  }
  else
  {
    return(NULL)
  }
}

###TEST ACCESSIBILITY OF REFUGIA USING GLM---------------------------------

#@glm.data = a data frame containing dependent and explanatory variables
#@y = name of dependent variable
#@keep.in = name(s) of variable(s) that are always included in the model (climatic vars) 
#@test = name(s) of variable(s) that are tested, i.e. accessibility of individual refugia
#@family = GLM family; default = 'binomial'

glm.test2 <- function(data=glm.data, y="species", keep.in="suitability", test=names(refs), family="binomial")
{
  if(length(test) > 1) {test.sel <- names(which(colSums(data[, test]) != 0))}
  if(length(test) == 1) {
    if(sum(data[, test]) > 0) {test.sel <- test}
    if(sum(data[, test]) == 0) {test.sel <- NULL}
  }
  
  tab <- list()
  
  if(length(test.sel) > 0)
  {
    ##standardize explanatory variables
    data[, c(keep.in, test.sel)] <- scale(data[, c(keep.in, test.sel)])
    
    f <- formula(paste(y, "~", paste0(keep.in, collapse = " + "), sep=" "))
    D <- glm(f, data = data, family = family)$deviance
    
    
    for(i in test.sel)
    {
      f <- formula(paste(y, "~", paste0(keep.in, collapse = " +"), "+", i, sep=" "))
      mod.glm <- suppressWarnings(glm(f, data = data, family = family))
      
      if(!is.na(mod.glm$coefficients[i]))
      {
        tab[[i]] <- summary(mod.glm)$coefficients[i,]
        tab[[i]]["D"] <- D
        tab[[i]]["D.diff"] <- D - mod.glm$deviance
      }
      if(is.na(mod.glm$coefficients[i]))
      {
        tab[[i]] <- rep(NA, 6)
        names(tab[[i]]) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "D", "D.diff")
      }
      
    }
    
    tab <- as.data.frame(do.call(rbind, tab))
    
    if(nrow(tab) < length(test))
    {
      tab[setdiff(test, rownames(tab)), ] <- NA
      tab <- tab[test,]
    }
  }
  
  if(length(test.sel) == 0)
  {
    tab <- as.data.frame(matrix(data=NA, ncol = 6, nrow = length(test)))
    rownames(tab) <- test
    colnames(tab) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "D", "D.diff")
  }
  
  
  return(tab)
}


###MODIFIED addRasterLegend FUNCTION-----------

addRasterLegend2 <- function (r, direction, side, location = "right", nTicks = 2, 
                              shortFrac = 0.02, longFrac = 0.3, axisOffset = 0, border = TRUE, 
                              ramp = "terrain", isInteger = "auto", ncolors = 64, 
                              breaks = NULL, minmax = NULL, locs = NULL, cex.axis = 0.8, 
                              labelDist = 0.7, digits = 2, ...) 
{
  if (!inherits(r, "RasterLayer")) {
    stop("r must be a RasterLayer.")
  }
  if (!hasArg("direction")) {
    direction <- "auto"
  }
  if (!direction %in% c("auto", "vertical", "horizontal")) {
    stop("direction must be auto, vertical or horizontal.")
  }
  if (is.character(location)) {
    if (!location %in% c("bottomleft", "bottomright", 
                         "topleft", "topright", "bottom", 
                         "top", "left", "right")) {
      stop("location is not recognized.")
    }
  }
  if (!isInteger %in% c("auto", TRUE, FALSE)) {
    stop("isInteger must be \"auto\", TRUE or FALSE.")
  }
  if (length(ramp) == 1) {
    if (ramp == "terrain") {
      pal <- rev(terrain.colors(ncolors))
    }
  }
  else if (length(ramp) > 1) {
    pal <- colorRampPalette(ramp)(ncolors)
  }
  if (is.null(minmax)) {
    colorbreaks <- seq(minValue(r), maxValue(r), length.out = (ncolors + 
                                                                 1))
  }
  else {
    colorbreaks <- seq(minmax[1], minmax[2], length.out = (ncolors + 
                                                             1))
  }
  if (!is.null(breaks)) {
    colorbreaks <- breaks
    ncolors <- length(breaks) + 1
  }
  n <- length(colorbreaks)
  minX <- grconvertX(par("fig")[1], from = "ndc", 
                     to = "user")
  maxX <- grconvertX(par("fig")[2], from = "ndc", 
                     to = "user")
  minY <- grconvertY(par("fig")[3], from = "ndc", 
                     to = "user")
  maxY <- grconvertY(par("fig")[4], from = "ndc", 
                     to = "user")
  xrange <- maxX - minX
  yrange <- maxY - minY
  minX <- minX + xrange * 0.05
  maxX <- maxX - xrange * 0.05
  minY <- minY + yrange * 0.05
  maxY <- maxY - yrange * 0.05
  if (is.character(location)) {
    if (location == "topleft" & direction %in% c("auto", 
                                                 "vertical")) {
      location <- vector("numeric", length = 4)
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * shortFrac
      location[3] <- maxY - (maxY - minY) * longFrac
      location[4] <- maxY
    }
    else if (location == "topleft" & direction == "horizontal") {
      location <- vector("numeric", length = 4)
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * longFrac
      location[3] <- maxY - (maxY - minY) * shortFrac
      location[4] <- maxY
    }
    else if (location == "topright" & direction %in% 
             c("auto", "vertical")) {
      location <- vector("numeric", length = 4)
      location[1] <- maxX - (maxX - minX) * shortFrac
      location[2] <- maxX
      location[3] <- maxY - (maxY - minY) * longFrac
      location[4] <- maxY
    }
    else if (location == "topright" & direction == 
             "horizontal") {
      location <- vector("numeric", length = 4)
      location[1] <- maxX - (maxX - minX) * longFrac
      location[2] <- maxX
      location[3] <- maxY - (maxY - minY) * shortFrac
      location[4] <- maxY
    }
    else if (location == "bottomleft" & direction %in% 
             c("auto", "vertical")) {
      location <- vector("numeric", length = 4)
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * shortFrac
      location[3] <- minY
      location[4] <- minY + (maxY - minY) * longFrac
    }
    else if (location == "bottomleft" & direction == 
             "horizontal") {
      location <- vector("numeric", length = 4)
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * longFrac
      location[3] <- minY
      location[4] <- minY + (maxY - minY) * shortFrac
    }
    else if (location == "bottomright" & direction %in% 
             c("auto", "vertical")) {
      location <- vector("numeric", length = 4)
      location[1] <- maxX - (maxX - minX) * shortFrac
      location[2] <- maxX
      location[3] <- minY
      location[4] <- minY + (maxY - minY) * longFrac
    }
    else if (location == "bottomright" & direction == 
             "horizontal") {
      location <- vector("numeric", length = 4)
      location[1] <- maxX - (maxX - minX) * longFrac
      location[2] <- maxX
      location[3] <- minY
      location[4] <- minY + (maxY - minY) * shortFrac
    }
    else if (location == "left") {
      location <- vector("numeric", length = 4)
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * shortFrac
      location[3] <- mean(par("usr")[3:4]) - ((maxY - 
                                                 minY) * longFrac)/2
      location[4] <- mean(par("usr")[3:4]) + ((maxY - 
                                                 minY) * longFrac)/2
      direction <- "vertical"
    }
    else if (location == "right") {
      location <- vector("numeric", length = 4)
      location[1] <- maxX - (maxX - minX) * shortFrac
      location[2] <- maxX
      location[3] <- mean(par("usr")[3:4]) - ((maxY - 
                                                 minY) * longFrac)/2
      location[4] <- mean(par("usr")[3:4]) + ((maxY - 
                                                 minY) * longFrac)/2
      direction <- "vertical"
    }
    else if (location == "top") {
      location <- vector("numeric", length = 4)
      location[1] <- mean(par("usr")[1:2]) - ((maxX - 
                                                 minX) * longFrac)/2
      location[2] <- mean(par("usr")[1:2]) + ((maxX - 
                                                 minX) * longFrac)/2
      location[3] <- maxY - (maxY - minY) * shortFrac
      location[4] <- maxY
      direction <- "horizontal"
    }
    else if (location == "bottom") {
      location <- vector("numeric", length = 4)
      location[1] <- mean(par("usr")[1:2]) - ((maxX - 
                                                 minX) * longFrac)/2
      location[2] <- mean(par("usr")[1:2]) + ((maxX - 
                                                 minX) * longFrac)/2
      location[3] <- minY
      location[4] <- minY + (maxY - minY) * shortFrac
      direction <- "horizontal"
    }
  }
  if (direction == "auto") {
    if (((location[2] - location[1])/(par("usr")[2] - 
                                      par("usr")[1])) >= ((location[4] - location[3])/(par("usr")[4] - 
                                                                                       par("usr")[3]))) {
      direction <- "horizontal"
    }
    else {
      direction <- "vertical"
    }
  }
  if (direction == "horizontal") {
    axisOffset <- axisOffset * (par("usr")[4] - par("usr")[3])
  }
  else if (direction == "vertical") {
    axisOffset <- axisOffset * (par("usr")[2] - par("usr")[1])
  }
  if (!hasArg("side")) {
    if (direction == "vertical") {
      if (mean(location[1:2]) <= mean(par("usr")[1:2])) {
        side <- 4
      }
      else {
        side <- 2
      }
    }
    if (direction == "horizontal") {
      if (mean(location[3:4]) > mean(par("usr")[3:4])) {
        side <- 1
      }
      else {
        side <- 3
      }
    }
  }
  if (direction == "horizontal") {
    x <- seq(from = location[1], to = location[2], length.out = n)
    width <- location[3:4]
  }
  else {
    x <- seq(from = location[3], to = location[4], length.out = n)
    width <- location[1:2]
  }
  x <- rep(x, each = 2)
  x <- x[-c(1, length(x))]
  x <- matrix(x, ncol = 2, byrow = TRUE)
  z <- rep(colorbreaks, each = 2)
  z <- z[-c(1, length(z))]
  z <- matrix(z, ncol = 2, byrow = TRUE)
  if (!is.null(locs)) {
    tol <- 1e-10
    maxValueIncluded <- FALSE
    if (any(abs(maxValue(r) - locs) < tol)) {
      locs <- locs[which(abs(maxValue(r) - locs) >= tol)]
      maxValueIncluded <- TRUE
    }
    tickLocs <- x[sapply(locs, function(x) which((x >= z[, 
                                                         1] & x < z[, 2]) == TRUE)), 1]
    if (maxValueIncluded) {
      tickLocs <- c(tickLocs, max(x[, 2]))
      locs <- c(locs, max(colorbreaks))
    }
    tx <- locs
  }
  else {
    tx <- trunc(seq(from = 1, to = nrow(x), length.out = nTicks + 
                      2))
    tickLocs <- x[tx, 1]
    tx <- z[tx, 1]
    tickLocs[length(tickLocs)] <- max(x[, 2])
    tx[length(tx)] <- max(z[, 2])
  }
  if (isInteger == "auto") {
    randomSample <- sample(values(r)[which(!is.na(values(r)))], 
                           size = 1000, replace = TRUE)
    if (identical(randomSample, trunc(randomSample))) {
      tx <- round(tx, 0)
    }
  }
  else if (isInteger) {
    tx <- round(tx, 0)
  }
  if (direction == "horizontal") {
    rect(xleft = x[, 1], ybottom = width[1], xright = x[, 
                                                        2], ytop = width[2], border = pal, col = pal, xpd = NA)
  }
  else {
    rect(xleft = width[1], ybottom = x[, 1], xright = width[2], 
         ytop = x[, 2], border = pal, col = pal, xpd = NA)
  }
  if (border) {
    rect(location[1], location[3], location[2], location[4], 
         border = "black", xpd = NA, lwd=0.4)
  }
  if (side == 1) {
    axis(side, at = tickLocs, pos = location[3] - axisOffset, 
         labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, 
         mgp = c(3, labelDist, 0), ...)
  }
  if (side == 3) {
    axis(side, at = tickLocs, pos = location[4] + axisOffset, 
         labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, 
         mgp = c(3, labelDist, 0), ...)
  }
  if (side == 2) {
    axis(side, at = tickLocs, pos = location[1] - axisOffset, 
         labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, 
         mgp = c(3, labelDist, 0), ...)
  }
  if (side == 4) {
    axis(side, at = tickLocs, pos = location[2] + axisOffset, 
         labels = round(tx, digits), xpd = NA, las = 1, cex.axis = cex.axis, 
         mgp = c(3, labelDist, 0), ...)
  }
}


