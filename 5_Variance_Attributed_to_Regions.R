#
# How much Va is attributed to associated regions?
# SEJ, AMS
# Jan 2019
#
# 

library(asreml)
library(dplyr)
library(ggplot2)
library(plyr)

load("data/20181107_BEAST_data_formatted.RData", verbose = T)
load("data/20181107_SoayGRM.Rdata")
load("data/20181107_Genabel_Data.RData")

source("r/makeGRM.R")
source("r/ASReml.EstEffects.R")

BEASTX$ID2 <- BEASTX$ID

focalsnps <- read.table("results/4_Top_SNPs_for_Paper.txt", header = T, sep = "\t", stringsAsFactors = F)

maptab <- data.frame(SNP.Name = snpnames(genabeldata),
                     Chromosome = chromosome(genabeldata),
                     Position = map(genabeldata), stringsAsFactors = F)

maptab <- arrange(maptab, Chromosome, Position)

runModels <- FALSE

#~~ How much variance does the region explain? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regvar <- NULL

for(i in 1:nrow(focalsnps)){
  
  snplist.name <- paste0("animal_models\\snplist_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_", focalsnps$Chromosome2[i])
  
  
  if(runModels == TRUE){
    
    models <- data.frame(Model = c("1+Sex+LambAgeAugust,random=~giv(ID)+giv(ID2)+MumID+BirthYearF+RunDate+Plate",
                                   "1+Sex+Age,random=~giv(ID)+ide(ID)+giv(ID2)+MumID+BirthYearF+CapYearF+RunDate+Plate"),
                         LambAdult = c("Lambs", "Adults"),
                         Response = focalsnps$Trait[i], stringsAsFactors = F)
    models <- subset(models, LambAdult == focalsnps$LambAdult[i])
    
    focsnp <- focalsnps$SNP.Name[i]
    
    x2 <- subset(maptab, Chromosome == focalsnps$Chromosome[i])
    
    if(!focsnp %in% maptab$SNP.Name){
      x1 <- strsplit(focsnp, split = "_")[[1]]
      focsnp <- x2$SNP.Name[which(x2$Position > as.numeric(x1[3]))[1]]
    } 
    
    snplist <- x2$SNP.Name[(which(x2$SNP.Name == focsnp) - 10):(which(x2$SNP.Name == focsnp) + 10)]
    
    #~~ Deal with SNPs at ends of chromosomes
    
    if(which(x2$SNP.Name == focsnp) > nrow(x2) - 10) snplist <- x2$SNP.Name[(nrow(x2)-20):nrow(x2)]
    if(which(x2$SNP.Name == focsnp) < 11) snplist <- x2$SNP.Name[1:21]
    
    print(snplist)
    
    writeLines(snplist, con = paste0(snplist.name, ".txt"))
    
    system(paste0("gcta64.exe --bfile data\\20180209_SoayPlates1-83 --autosome --autosome-num 26 --keep data\\UniqueIDs.txt --extract ", snplist.name, ".txt --make-grm-gz --out ", snplist.name, "_GRM" ))
    system(paste0("gcta64.exe --grm-gz ", snplist.name, "_GRM --grm-adj 0 --make-grm-gz --out ", snplist.name, "_GRM.adj"))
    
    grm.region <- read.table(paste0(snplist.name, "_GRM.adj.grm.gz"))  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
    ids.region <- read.table(paste0(snplist.name, "_GRM.adj.grm.id"))  # CONTAINS ID LIST
    
    # crop grm to only IDs I need
    
    grmregion <- makeGRM(grm.region, ids.region, id.vector = BEASTX$ID)
    
    #~~ Make data for model
    
    x.vec <- as.character(models$Model[1])
    x.vec <- gsub(",random=~", "+", x.vec)
    x.vec <- strsplit(x.vec, split = "\\+")[[1]]
    x.vec[grep("giv", x.vec)] <- "ID"
    x.vec[grep("ide", x.vec)] <- "ID"
    x.vec <- c(x.vec, as.character(focalsnps$Trait[i]), "ID2")
    x.vec <- unique(x.vec[-which(x.vec == 1)])
    
    if(models$LambAdult[1] %in% c("Lambs", "Adults")){
      x.data <- subset(BEASTX, LambAdult == focalsnps$LambAdult[i])
    } else {
      x.data <- BEASTX
    }
    
    x.data <- droplevels(na.omit(x.data[,x.vec]))
    
    x.data <- droplevels(subset(x.data, ID %in% dimnames(grminv)[[1]]))
    
    
    try({
      eval(parse(text = paste0("fit1 <- asreml(fixed=", focalsnps$Trait[i], "~",
                               models$Model[1],"
                             , data=x.data, ginverse=list(ID=grminv, ID2 = grmregion),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)")))
      
      save(fit1, file = paste0(snplist.name, ".RData"))
      
      
      
    })
    
    rm(x.data, x.vec, grm.region, grmregion, ids.region, snplist)
    
  } else {
    
    load(file = paste0(snplist.name, ".RData"))
    
  }
  
  x <- cbind(Trait = focalsnps$Trait[i],
             LambAdult = focalsnps$LambAdult[i],
             Chromsome = focalsnps$Chromosome2[i],
             print(ASReml.EstEffects(fit1)))
  
  
  regvar <- rbind(regvar, x)
  
  rm(snplist.name, fit1, x)
}

regvar$Model <- paste0(regvar$Trait, "_", regvar$LambAdult, "_", regvar$Chromsome)


recodetab <- data.frame(variable = c("giv(ID).giv",
                                     "giv(ID2).giv",
                                     "ide(ID)!id",
                                     "BirthYearF!BirthYearF.var",
                                     "CapYearF!CapYearF.var",
                                     "MumID!MumID.var", 
                                     "RunDate!RunDate.var", 
                                     "Plate!Plate.var",
                                     "R!variance"),
                        new.variable = c("Additive Genetic",
                                         "Regional Genetic",
                                         "Permanent Environment",
                                         "Birth Year",
                                         "Capture Year",
                                         "Mother Identity", 
                                         "Run Date", 
                                         "Plate ID",
                                         "Residual"))

regvar <- join(regvar, recodetab)

regvar$new.variable2 <- factor(regvar$new.variable,
                               levels = rev(c("Additive Genetic", "Regional Genetic", "Permanent Environment","Birth Year","Capture Year",
                                              "Mother Identity", "Run Date", "Plate ID","Residual")))

ggplot(regvar, aes(Model, Effect, fill = new.variable2)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(fill = "Variable") +
  theme(axis.text.x  = element_text (size = 12, angle = 270),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14))


x <- reshape2::dcast(subset(regvar, select = c(Model, new.variable, Effect)), Model ~ new.variable)
x$Vp <- rowSums(x[,2:ncol(x)], na.rm = T)
x$GenVar <- x$`Regional Genetic`/(x$`Regional Genetic` + x$`Additive Genetic`)
x$GenVar <- round(x$GenVar, digits = 3)
x$`Regional Genetic` <- round(x$`Regional Genetic`, digits = 4)

subset(x, select = c(Model,  `Regional Genetic`, GenVar))


write.table(regvar, "results/5_Variance_Attributed_by_associated_Regions.txt", row.names = F, sep = "\t", quote = F)

