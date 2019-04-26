#
# Genotype imputation - parse
# SEJ, AMS
# Nov 2018
#
#


library(ggplot2)
library(plyr)
library(reshape)

sig.regions <- read.table("results/2_Regions_for_Imputation_log_lambs.txt", header = T, stringsAsFactors = F)

sig.regions <- data.frame(Min = tapply(sig.regions$Position, sig.regions$Chromosome, min)  - 2e6,
                          Max = tapply(sig.regions$Position, sig.regions$Chromosome, max) + 2e6)
sig.regions$Chromosome <- row.names(sig.regions)


full.imputedgenos <- NULL
full.mapfile <- NULL

for(h in 1:nrow(sig.regions)){
  
  print(paste("Running chromosome", sig.regions$Chromosome[h]))
  
  setwd(paste0("alphaimpute/chr", sig.regions$Chromosome[h]))
  #~~ Load imputed data
  
  imputedgenos <- read.table("Results/ImputeGenotypes.txt")
  imputedgenos[1:5, 1:5]
  
  mapfile <- read.table(paste0("chr", sig.regions$Chromosome[h], "merge.bim"), stringsAsFactors = F)
  head(mapfile)
  
  #~~ Look at the success of imputation
  
  mapfile$ImputeSuccess <- NA
  
  for(i in 1:nrow(mapfile)){
    mapfile$ImputeSuccess[i] <- 1 - (length(which(imputedgenos[,i+1] == 9))/nrow(imputedgenos))
  }
  
  ggplot(mapfile, aes(ImputeSuccess)) + geom_histogram(binwidth = 0.01)
  ggplot(mapfile, aes(1:nrow(mapfile), ImputeSuccess)) + geom_point()
  
  #~~ Convert to genotypes
  
  names(imputedgenos)[2:ncol(imputedgenos)] <- mapfile$V2
  
  
  for(i in 2:ncol(imputedgenos)){
    if(i %in% seq(1, nrow(imputedgenos),  50)) print(paste("Running column", i, "of", ncol(imputedgenos)-1))
    temp <- data.frame(Value = imputedgenos[,i])
    temp2 <- data.frame(Value    = c(0, 1, 2, 9),
                        NewValue = c(paste0(mapfile$V5[i-1], "/", mapfile$V5[i-1]),
                                     paste0(mapfile$V5[i-1], "/", mapfile$V6[i-1]),
                                     paste0(mapfile$V6[i-1], "/", mapfile$V6[i-1]),
                                     NA))
    
    suppressMessages(temp <- join(temp, temp2))
    imputedgenos[,i] <- temp$NewValue
    rm(temp, temp2)
  }
  
  imputedgenos <- melt(imputedgenos, id.vars = "V1")
  
  imputedgenos$variable <- as.character(imputedgenos$variable)
  imputedgenos$value <- as.character(imputedgenos$value)
  
  full.imputedgenos <- rbind(full.imputedgenos, imputedgenos)
  full.mapfile <- rbind(full.mapfile, mapfile)
  
  rm(imputedgenos, mapfile)
  
  setwd("../..")  
}

library(reshape)

full.imputedgenos <- cast(full.imputedgenos, formula = V1 ~ variable)
full.imputedgenos.log.lambs <- full.imputedgenos

save(full.imputedgenos.log.lambs, full.mapfile, file = "alphaimpute/Imputed_GWAS_log_Lambs.RData")
