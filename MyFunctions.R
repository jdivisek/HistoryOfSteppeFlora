##DATA FILTERING IN GRIDDED ENVIRONMENTAL SPACE----------------------------------------------------
##sensu Varela et al. 2014
##inspired by Míša Vojtěchovských

#@OccurenceData - X (Longitude) and Y (Latitude) coordinates of species occurrence points
#@Scores - standardized environmental variables (mean = 0, sd = 1) or PCA scores
#@Resolution - resolution environmental space (in SD units), the higher value, the more points are deleted
#@Seed - seed for random number generation
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
  
  #attach geographical coordinates
  points.scores$X <- OccurenceData[,1]
  points.scores$Y <- OccurenceData[,2]
  
  #randomize data
  set.seed(Seed)
  points.scores <- points.scores[sample(nrow(points.scores)),]
  
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

glm.test <- function(data=glm.data, y, keep.in, test, family="binomial")
{
  tab <- list()
  
  f <- formula(paste(y, "~", paste0(keep.in, collapse = " + "), sep=" "))
  D <- glm(f, data = data, family = family)$deviance
  
  for(i in test)
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
      tab[[i]] <- rep(NA, 4)
      names(tab[[i]]) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
      tab[[i]]["D"] <- NA
      tab[[i]]["D.diff"] <- NA
    }
    
  }
  
  return(as.data.frame(do.call(rbind, tab)))
}

