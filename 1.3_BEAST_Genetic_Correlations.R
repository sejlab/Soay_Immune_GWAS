#  ------------------------------------------------------------------------
# Beast QG MS
# Animal models for iga, ige, igg across all ages, lambs, adults
# AMS - 07/02/18
#  ------------------------------------------------------------------------


load("data/20181107_SoayGRM.Rdata")
load("data/20181107_BEAST_data_formatted.RData")

library(asreml)
library(plyr)
library(ggplot2)
library(magrittr)
library(tidyr)

source("r/ASReml.EstEffects.R")
source("r/ASReml.ExtractPredictors.R")
source("r/pin.R")

ainv <- asreml.Ainverse(pedigree)$ginv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Make model structures                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

BEASTX$Lamb.IgAmp  <- ifelse(BEASTX$LambAdult == "Lambs" , BEASTX$IgAmp, NA)
BEASTX$Lamb.IgEmp  <- ifelse(BEASTX$LambAdult == "Lambs" , BEASTX$IgEmp, NA)
BEASTX$Lamb.IgGmp  <- ifelse(BEASTX$LambAdult == "Lambs" , BEASTX$IgGmp, NA)
BEASTX$Adult.IgAmp <- ifelse(BEASTX$LambAdult == "Adults", BEASTX$IgAmp, NA)
BEASTX$Adult.IgEmp <- ifelse(BEASTX$LambAdult == "Adults", BEASTX$IgEmp, NA)
BEASTX$Adult.IgGmp <- ifelse(BEASTX$LambAdult == "Adults", BEASTX$IgGmp, NA)

head(BEASTX)


x1 <- data.frame(Trait1 = c("Lamb.IgAmp", "Lamb.IgEmp", "Lamb.IgGmp", "Adult.IgAmp", "Adult.IgEmp", "Adult.IgGmp"))
x2 <- x1
names(x2) <- "Trait2"

x1 <- merge(x1, x2) %>%
  apply(1, sort) %>%
  t %>%
  data.frame %>%
  unique()

x1 <- subset(x1, X1 != X2)  

x1 <- separate(x1, X1, c("LambAdult.1", "Trait.1"), sep = "\\.", remove = F)
x1 <- separate(x1, X2, c("LambAdult.2", "Trait.2"), sep = "\\.", remove = F)

x1$Model <- NA

x1$Model <- ifelse(x1$LambAdult.1 == "Lamb" & x1$LambAdult.2 == "Lamb",
                   "trait+trait:Sex+trait:LambAgeAugust,random=~corgh(trait):ped(ID)",
                   x1$Model)

x1$Model <- ifelse(x1$LambAdult.1 == "Adult" & x1$LambAdult.2 == "Adult",
                   "trait+trait:Sex+trait:Age,random=~corgh(trait):ped(ID)+idh(trait):ide(ID)",
                   x1$Model)


x1$Model <- ifelse(x1$LambAdult.1 != x1$LambAdult.2,
                   "trait+trait:Sex,random=~corgh(trait):ped(ID)",
                   x1$Model)

x1$Model <- paste0("cbind(", x1$X1, ",", x1$X2, ") ~ ", x1$Model)
x1$ModelNo <- 1:nrow(x1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Run the Models                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

restab <- NULL
modlist <- list()

for(i in 1:nrow(x1)){
  
  print(paste0("Running model ", i))
  
  eval(parse(text = paste0("fit1 <- asreml(fixed  = ", x1$Model[i], ",
                                  rcov   = ~ units:idh(trait, init = NA),
                                  data = BEASTX,
                                  ginverse = list(ID = ainv),
                                  workspace = 500e+6, pworkspace = 500e+6,
                                  maxiter = 100, na.method.Y = \"include\", na.method.X = \"include\")")))
  
  modlist[[i]] <- fit1
  
  temp <- cbind(ModelNo = i, summary(fit1)$varcomp)
  temp$Effect <- row.names(temp)
  
  restab <- rbind(restab, temp)
  
  rm(fit1, temp)
}


x1$Model
x1$Model <- gsub("ped", "giv", x1$Model)

#~~ GRM Models

restab.grm <- NULL
BEASTX$ID2 <- as.character(BEASTX$ID)
BEASTX <- subset(BEASTX, ID2 %in% dimnames(grminv)[[1]]) %>% droplevels

for(i in 1:nrow(x1)){
  
  print(paste0("Running model ", i))
  
  eval(parse(text = paste0("fit1 <- asreml(fixed  = ", x1$Model[i], ",
                                  rcov   = ~ units:idh(trait, init = NA),
                                  data = BEASTX,
                                  ginverse = list(ID = grminv),
                                  workspace = 500e+6, pworkspace = 500e+6,
                                  maxiter = 100, na.method.Y = \"include\", na.method.X = \"include\")")))
  
  temp <- cbind(ModelNo = i, summary(fit1)$varcomp)
  temp$Effect <- row.names(temp)
  
  restab.grm <- rbind(restab.grm, temp)
  
  save(fit1, file = paste0("bivar", i, ".RData"))
  
  rm(fit1, temp)
  gc()
}

# 
# 
# fit1 <- asreml(fixed  = cbind(Lamb.IgAmp,Lamb.IgEmp) ~ trait+trait:Sex+trait:LambAgeAugust,
#                random=~corgh(trait):giv(ID),
#                rcov   = ~ units:idh(trait, init = NA),
#                data = BEASTX,
#                ginverse = list(ID = grminv),
#                workspace = 500e+6, pworkspace = 500e+6,
#                maxiter = 100, na.method.Y = "include", na.method.X = "include")
# 
