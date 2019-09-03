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

source("r/ASReml.EstEffects.R")
source("r/ASReml.ExtractPredictors.R")
source("r/pin.R")

ainv <- asreml.Ainverse(pedigree)$ginv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Make model structures & run with pedigree    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

models <- data.frame(Model = c("1+Sex+LambAgeAugust,random=~ped(ID)+MumID+BirthYearF+RunDate+Plate",
                               "1+Sex+Age,random=~ped(ID)+ide(ID)+MumID+BirthYearF+CapYearF+RunDate+Plate"),
                     LambAdult = c("Lambs", "Adults"))

models <- merge(models, data.frame(Response = c("IgAmp", "IgGmp", "IgEmp")))

for(i in 1:3) models[,i] <- as.character(models[,i])

ped.results.list <- list()
ped.results.ranef <- NULL
ped.results.fixef <- NULL

for(i in 1:nrow(models)){
  
  x.vec <- as.character(models$Model[i])
  x.vec <- gsub(",random=~", "+", x.vec)
  x.vec <- strsplit(x.vec, split = "\\+")[[1]]
  x.vec[grep("ped", x.vec)] <- "ID"
  x.vec[grep("ide", x.vec)] <- "ID"
  x.vec <- c(x.vec, as.character(models$Response[i]))
  x.vec <- unique(x.vec[-which(x.vec == 1)])
  
  if(models$LambAdult[i] %in% c("Lambs", "Adults")){
    x.data <- subset(BEASTX, LambAdult == models$LambAdult[i])
  } else {
    x.data <- BEASTX
  }
  

  
  x.data <- na.omit(x.data[,x.vec])
  x.data <- droplevels(subset(x.data, ID %in% dimnames(grminv)[[1]]))
  
  if(models$LambAdult[i] == "Lambs"){
    eval(parse(text = paste0("x.data <- subset(x.data, ", models$Response[i], "!= 0)")))
    x.data[,models$Response[i]] <- log10(x.data[,models$Response[i]])
  }

  
  eval(parse(text = paste0("fit1 <- asreml(fixed=", models$Response[i], "~",
                           models$Model[i],"
                          , data=x.data, ginverse=list(ID=ainv),
                          workspace = 500e+6, pworkspace = 500e+6,
                          maxiter = 100)")))
  
  ped.results.list[[i]] <- fit1
  
  ped.results.ranef <- rbind(ped.results.ranef,
                             cbind(Trait = models$Response[i],
                                   LambAdult = models$LambAdult[i],
                                   ASReml.EstEffects(fit1)))
  
  x <- data.frame(summary(fit1, all = T)$coef.fixed)
  x$variable <- row.names(x)
  
  ped.results.fixef <- rbind(ped.results.fixef,
                             cbind(Trait = models$Response[i],
                                   LambAdult = models$LambAdult[i],
                                   x))
  
  rm(fit1, x.data, x.vec, x)
  
}


recodetab <- data.frame(variable = c("ped(ID)!ped",
                                     "ide(ID)!id",
                                     "BirthYearF!BirthYearF.var",
                                     "CapYearF!CapYearF.var",
                                     "MumID!MumID.var", 
                                     "RunDate!RunDate.var", 
                                     "Plate!Plate.var",
                                     "R!variance"),
                        new.variable = c("Additive Genetic",
                                         "Permanent Environment",
                                         "Birth Year",
                                         "Capture Year",
                                         "Mother Identity", 
                                         "Run Date", 
                                         "Plate ID",
                                         "Residual"))

recodetab$new.variable <- factor(recodetab$new.variable,
                                 levels = c("Additive Genetic",
                                            "Permanent Environment",
                                            "Birth Year",
                                            "Capture Year",
                                            "Mother Identity", 
                                            "Run Date", 
                                            "Plate ID",
                                            "Residual"))


ped.results.ranef <- join(ped.results.ranef, recodetab)

rm(recodetab)


ggplot(ped.results.ranef, aes(LambAdult, Effect, fill = new.variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Trait) +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "Variable")


ggplot(ped.results.ranef, aes(new.variable, Effect, fill = new.variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0) +
  facet_grid(LambAdult~Trait) +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x  = element_text (angle = 270, hjust = 0, vjust = 0.3),
        legend.position = "none") +
  labs(x = "Variable")


write.table(ped.results.fixef, "results/1_Pedigree_Fixef_Results.txt", row.names = F, sep = "\t", quote = F)
write.table(ped.results.ranef, "results/1_Pedigree_Ranef_Results.txt", row.names = F, sep = "\t", quote = F)

rm(models, i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Make model structures & run with GRM         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

models <- data.frame(Model = c("1+Sex+LambAgeAugust,random=~giv(ID)+MumID+BirthYearF+RunDate+Plate",
                               "1+Sex+Age,random=~giv(ID)+ide(ID)+MumID+BirthYearF+CapYearF+RunDate+Plate"),
                     LambAdult = c("Lambs", "Adults"))

models <- merge(models, data.frame(Response = c("IgAmp", "IgGmp", "IgEmp")))

models$V_raw <- NA
models$VP <- NA
models$VP_SE <- NA
models$N <- NA
models$N_ids <- NA
models$Mean <- NA

for(i in 1:3) models[,i] <- as.character(models[,i])

grm.results.ranef <- NULL
grm.results.fixef <- NULL

Run.GRM.Models <- TRUE

for(i in 1:nrow(models)){
  
  print(paste("Running model", i))
  
  x.vec <- as.character(models$Model[i])
  x.vec <- gsub(",random=~", "+", x.vec)
  x.vec <- strsplit(x.vec, split = "\\+")[[1]]
  x.vec[grep("giv", x.vec)] <- "ID"
  x.vec[grep("ide", x.vec)] <- "ID"
  x.vec <- c(x.vec, as.character(models$Response[i]))
  x.vec <- unique(x.vec[-which(x.vec == 1)])
  
  if(models$LambAdult[i] %in% c("Lambs", "Adults")){
    x.data <- subset(BEASTX, LambAdult == models$LambAdult[i])
  } else {
    x.data <- BEASTX
  }
  x.data <- na.omit(x.data[,x.vec])
  x.data <- droplevels(subset(x.data, ID %in% dimnames(grminv)[[1]]))
  
  if(Run.GRM.Models == TRUE){
    
    
    eval(parse(text = paste0("fit1 <- asreml(fixed=", models$Response[i], "~",
                             models$Model[i],"
                           , data=x.data, ginverse=list(ID=grminv)
                           , workspace = 500e+6, pworkspace = 500e+6,
                           maxiter = 100)")))
    
    save(fit1, file = paste0("animal_models/", models$Response[i], "_", models$LambAdult[i], ".RData"))
    
    
    grm.results.ranef <- rbind(grm.results.ranef,
                               cbind(Trait = models$Response[i],
                                     LambAdult = models$LambAdult[i],
                                     ASReml.EstEffects(fit1)))
    
    x <- data.frame(summary(fit1, all = T)$coef.fixed)
    x$variable <- row.names(x)
    
    grm.results.fixef <- rbind(grm.results.fixef,
                               cbind(Trait = models$Response[i],
                                     LambAdult = models$LambAdult[i],
                                     x))
    
    rm(fit1)
    
    eval(parse(text = paste0("fit1 <- asreml(fixed=", models$Response[i], "~",
                             models$Model[i],"
                           , data=x.data, ginverse=list(ID=grminv),
                           rcov = ~idv(units), workspace = 500e+6, pworkspace = 500e+6,
                           maxiter = 100)")))

    save(fit1, file = paste0("animal_models/", models$Response[i], "_", models$LambAdult[i], "_vp.RData"))

    
    
  } else {
    
    
    load(paste0("animal_models/", models$Response[i], "_", models$LambAdult[i], ".RData"))
    
    
    grm.results.ranef <- rbind(grm.results.ranef,
                               cbind(Trait = models$Response[i],
                                     LambAdult = models$LambAdult[i],
                                     ASReml.EstEffects(fit1)))
    
    x <- data.frame(summary(fit1, all = T)$coef.fixed)
    x$variable <- row.names(x)
    
    grm.results.fixef <- rbind(grm.results.fixef,
                               cbind(Trait = models$Response[i],
                                     LambAdult = models$LambAdult[i],
                                     x))
    
    rm(fit1)
    
    load(paste0("animal_models/", models$Response[i], "_", models$LambAdult[i], "_vp.RData"))
  }
  
  if(models$LambAdult[i] == "Lambs"){
    temp <- pin(fit1, phen.var ~ V1 + V2 + V3 + V4 + V5 + V7)
  } else {
    temp <- pin(fit1, phen.var ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V9)
  }


  models$VP[i] <-   unlist(temp[1])
  models$VP_SE[i] <- unlist(temp[2])

  models$V_raw[i] <- var(x.data[,which(names(x.data) == models$Response[i])])
  models$Mean[i] <- mean(x.data[,which(names(x.data) == models$Response[i])])

  models$N[i] <- nrow(x.data)
  models$N_ids[i] <- length(unique(x.data$ID))

  
  rm(fit1, x.data, x.vec, x, temp)
  
}


recodetab <- data.frame(variable = c("giv(ID).giv",
                                     "ide(ID)!id",
                                     "BirthYearF!BirthYearF.var",
                                     "CapYearF!CapYearF.var",
                                     "MumID!MumID.var", 
                                     "RunDate!RunDate.var", 
                                     "Plate!Plate.var",
                                     "R!variance"),
                        new.variable = c("Additive Genetic",
                                         "Permanent Environment",
                                         "Birth Year",
                                         "Capture Year",
                                         "Mother Identity", 
                                         "Run Date", 
                                         "Plate ID",
                                         "Residual"))

recodetab$new.variable <- factor(recodetab$new.variable,
                                 levels = rev(c("Additive Genetic",
                                            "Permanent Environment",
                                            "Birth Year",
                                            "Capture Year",
                                            "Mother Identity", 
                                            "Run Date", 
                                            "Plate ID",
                                            "Residual")))


grm.results.ranef <- join(grm.results.ranef, recodetab)

rm(recodetab)

rev(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"))

y <- data.frame(Trait = c("IgAmp", "IgGmp", "IgEmp"),
                Trait2 = c("Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgE"))

grm.results.ranef <- join(grm.results.ranef, y)

ggplot(grm.results.ranef, aes(LambAdult, Effect, fill = new.variable)) +
  geom_bar(stat = "identity", col = "grey50") +
  facet_wrap(~Trait2) +
  scale_fill_manual(values = rev(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"))) +
  labs(fill = "", x = "Age Class", y = "Proportion of Phenotypic Variance") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.text = element_text (size = 12))
ggsave("figs/Proportion_of_variance.png", width = 8, height = 5)

ggplot(grm.results.ranef, aes(new.variable, Effect, fill = new.variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0) +
  facet_grid(LambAdult~Trait) +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x  = element_text (angle = 270, hjust = 0, vjust = 0.3),
        legend.position = "none") +
  labs(x = "Variable")


write.table(grm.results.fixef, "results/1_GRM_Fixef_Results.txt", row.names = F, sep = "\t", quote = F)
write.table(grm.results.ranef, "results/1_GRM_Ranef_Results.txt", row.names = F, sep = "\t", quote = F)
write.table(models,            "results/1_GRM_Sample_Sizes_VP.txt", row.names = F, sep = "\t", quote = F)
