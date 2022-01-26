###############################################################################
###         CONTRIBUTIONS OF ACCESSIBILITY FOR LGM SOURCE AREAS             ###
###############################################################################

###LAST GLACIAL MAXIMUM------------------------------------------------------

dev.res1 <- vector("numeric")
dev.res2 <- vector("numeric")
dev.res3 <- vector("numeric")

dev.null <- vector("numeric")

for(q in names(acces.test))
{
  print(q)
  
  if(any(!is.na(acces.test[[q]])))
  {
    r <- raster(paste0("./Source areas LGM/", q, ".tif"))
    
    if(maxValue(r) > 0)
    {
      test <- acces.test[[q]]
      
      mr <- test$mig.rate[1]
      
      glm.data <- list.load(paste0("./Source areas LGM/", q, "_TestData.rdata"))
      glm.data <- glm.data[[mr/10]]
      
      test <- test[which(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05),]
      
      glm.mod1 <- glm(species ~ ., data=glm.data[,c("species", colnames(glm.env.data))], family = "binomial")
      dev.res1[q] <- glm.mod1$deviance
      
      glm.mod2 <- glm(species ~ ., data=glm.data[,c("species", colnames(glm.env.data), rownames(test))], family = "binomial")
      dev.res2[q] <- glm.mod2$deviance
      
      glm.mod3 <- glm(species ~ ., data=glm.data[,c("species", rownames(test))], family = "binomial")
      dev.res3[q] <- glm.mod3$deviance
      
      dev.null[q] <- glm.mod1$null.deviance
      
    }
    else
    {
      dev.res1[q] <- NA
      dev.res2[q] <- NA
      dev.res3[q] <- NA
      dev.null[q] <- NA
    }
  }
  else
  {
    dev.res1[q] <- NA
    dev.res2[q] <- NA
    dev.res3[q] <- NA
    dev.null[q] <- NA
  }
}

acces.contrib <- as.data.frame(cbind(dev.null, dev.res1, dev.res2, dev.res3))
head(acces.contrib)

acces.contrib.R2 <- as.data.frame(acces.contrib)
acces.contrib.R2$dev.res1 <- (1- (acces.contrib.R2$dev.res1 / acces.contrib.R2$dev.null))*100
acces.contrib.R2$dev.res2 <- (1- (acces.contrib.R2$dev.res2 / acces.contrib.R2$dev.null))*100
acces.contrib.R2$dev.res3 <- (1- (acces.contrib.R2$dev.res3 / acces.contrib.R2$dev.null))*100
head(acces.contrib.R2)
acces.contrib.R2 <- round(acces.contrib.R2[,-1], 2)
colnames(acces.contrib.R2) <- c("climate", "climate.acc", "acc")

acces.contrib.R2$chorotype <- clas.spe[rownames(acces.contrib.R2), "chorotype"]

boxplot(climate ~ chorotype, data=acces.contrib.R2)
boxplot(climate.acc ~ chorotype, data=acces.contrib.R2)
boxplot(acc ~ chorotype, data=acces.contrib.R2)

apply(acces.contrib.R2, 2, mean, na.rm=T)
apply(acces.contrib.R2, 2, sd, na.rm=T)

apply(acces.contrib.R2[acces.contrib.R2$chorotype == 1,], 2, mean, na.rm=T)
apply(acces.contrib.R2[acces.contrib.R2$chorotype == 1, ], 2, sd, na.rm=T)

apply(acces.contrib.R2[acces.contrib.R2$chorotype == 2,], 2, mean, na.rm=T)
apply(acces.contrib.R2[acces.contrib.R2$chorotype == 2, ], 2, sd, na.rm=T)

apply(acces.contrib.R2[acces.contrib.R2$chorotype == 3,], 2, mean, na.rm=T)
apply(acces.contrib.R2[acces.contrib.R2$chorotype == 3, ], 2, sd, na.rm=T)

apply(acces.contrib.R2[acces.contrib.R2$chorotype == 4,], 2, mean, na.rm=T)
apply(acces.contrib.R2[acces.contrib.R2$chorotype == 4, ], 2, sd, na.rm=T)

apply(acces.contrib.R2[acces.contrib.R2$chorotype == 5,], 2, mean, na.rm=T)
apply(acces.contrib.R2[acces.contrib.R2$chorotype == 5, ], 2, sd, na.rm=T)

acces.contrib.R2$meadow.steppe <- clas.spe[rownames(acces.contrib.R2), "meadow.steppe"]
acces.contrib.R2$grass.steppe <- clas.spe[rownames(acces.contrib.R2), "grass.steppe"]
acces.contrib.R2$rocky.steppe <- clas.spe[rownames(acces.contrib.R2), "rocky.steppe"]

apply(acces.contrib.R2[!is.na(acces.contrib.R2$meadow.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.R2[!is.na(acces.contrib.R2$meadow.steppe),], 2, sd, na.rm=T)

apply(acces.contrib.R2[!is.na(acces.contrib.R2$grass.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.R2[!is.na(acces.contrib.R2$grass.steppe),], 2, sd, na.rm=T)

apply(acces.contrib.R2[!is.na(acces.contrib.R2$rocky.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.R2[!is.na(acces.contrib.R2$rocky.steppe),], 2, sd, na.rm=T)

###MID-HOLOCENE-----------------------------------------------------------------

dev.res1 <- vector("numeric")
dev.res2 <- vector("numeric")
dev.res3 <- vector("numeric")

dev.null <- vector("numeric")

for(q in names(acces.test.MH))
{
  print(q)
  
  if(any(!is.na(acces.test.MH[[q]])))
  {
    r <- raster(paste0("./Source areas mid-Holocene/", q, ".tif"))
    
    if(maxValue(r) > 0)
    {
      test <- acces.test.MH[[q]]
      
      glm.data <- list.load(paste0("./Source areas mid-Holocene/", q, "_TestData.rdata"))
      
      test <- test[which(test$D.diff > 0 & test$Estimate > 0 & test$`Pr(>|z|)` < 0.05),]
      
      glm.mod1 <- glm(species ~ ., data=glm.data[,c("species", colnames(glm.env.data))], family = "binomial")
      dev.res1[q] <- glm.mod1$deviance
      
      glm.mod2 <- glm(species ~ ., data=glm.data[,c("species", colnames(glm.env.data), rownames(test))], family = "binomial")
      dev.res2[q] <- glm.mod2$deviance
      
      glm.mod3 <- glm(species ~ ., data=glm.data[,c("species", rownames(test))], family = "binomial")
      dev.res3[q] <- glm.mod3$deviance
      
      dev.null[q] <- glm.mod1$null.deviance
      
    }
    else
    {
      dev.res1[q] <- NA
      dev.res2[q] <- NA
      dev.res3[q] <- NA
      dev.null[q] <- NA
    }
  }
  else
  {
    dev.res1[q] <- NA
    dev.res2[q] <- NA
    dev.res3[q] <- NA
    dev.null[q] <- NA
  }
}

acces.contrib.MH <- as.data.frame(cbind(dev.null, dev.res1, dev.res2, dev.res3))
head(acces.contrib.MH)

acces.contrib.MH.R2 <- as.data.frame(acces.contrib.MH)
acces.contrib.MH.R2$dev.res1 <- (1- (acces.contrib.MH.R2$dev.res1 / acces.contrib.MH.R2$dev.null))*100
acces.contrib.MH.R2$dev.res2 <- (1- (acces.contrib.MH.R2$dev.res2 / acces.contrib.MH.R2$dev.null))*100
acces.contrib.MH.R2$dev.res3 <- (1- (acces.contrib.MH.R2$dev.res3 / acces.contrib.MH.R2$dev.null))*100
head(acces.contrib.MH.R2)
acces.contrib.MH.R2 <- round(acces.contrib.MH.R2[,-1], 2)
colnames(acces.contrib.MH.R2) <- c("climate", "climate.acc", "acc")

boxplot(acces.contrib.MH.R2)
acces.contrib.MH.R2$chorotype <- clas.spe[rownames(acces.contrib.MH.R2), "chorotype"]

boxplot(climate ~ chorotype, data=acces.contrib.MH.R2)
boxplot(climate.acc ~ chorotype, data=acces.contrib.MH.R2)
boxplot(acc ~ chorotype, data=acces.contrib.MH.R2)

apply(acces.contrib.MH.R2, 2, mean, na.rm=T)
apply(acces.contrib.MH.R2, 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 1,], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 1, ], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 2,], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 2, ], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 3,], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 3, ], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 4,], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 4, ], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 5,], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[acces.contrib.MH.R2$chorotype == 5, ], 2, sd, na.rm=T)

acces.contrib.MH.R2$meadow.steppe <- clas.spe[rownames(acces.contrib.MH.R2), "meadow.steppe"]
acces.contrib.MH.R2$grass.steppe <- clas.spe[rownames(acces.contrib.MH.R2), "grass.steppe"]
acces.contrib.MH.R2$rocky.steppe <- clas.spe[rownames(acces.contrib.MH.R2), "rocky.steppe"]

apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$meadow.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$meadow.steppe),], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$grass.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$grass.steppe),], 2, sd, na.rm=T)

apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$rocky.steppe),], 2, mean, na.rm=T)
apply(acces.contrib.MH.R2[!is.na(acces.contrib.MH.R2$rocky.steppe),], 2, sd, na.rm=T)
