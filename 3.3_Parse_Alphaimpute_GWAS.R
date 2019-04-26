#
# Parse the Imputed GWAS Results
# Dec 2018
# SEJ
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up workspace                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Load libraries and data

library(GenABEL)
library(ggplot2)
library(plyr)
library(reshape2)
library(biomaRt)
library(tidyr)
library(dplyr)

source("r/multiplot.R")

setwd("alphaimpute/")

load("Imputed_GWAS.RData", verbose = T)

gene.list <- read.table("../data/Ovis_aries.Oar_v3.1.94.genes.txt", header = T, stringsAsFactors = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Read in and tabulate the results       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- dir()
x <- x[grep("RData", x)]
x <- x[grep("GWAS_", x)]


full.fixef <- list()
full.ranef <- list()
full.wald  <- list()
full.n     <- list()

for(i in 1:length(x)){
  
  print(paste("Running line", i, "of", length(x)))
  
  load(x[i])
  
  full.fixef[[i]] <- restab.fixef
  full.ranef[[i]] <- restab.ranef
  full.wald[[i]]  <- restab.wald
  full.n[[i]]     <- restab.n
  
  rm(restab.ranef, restab.fixef, restab.wald, restab.n)
  
}


suppressWarnings({
  full.fixef <- bind_rows(full.fixef)
  full.ranef <- bind_rows(full.ranef)
  full.wald  <- bind_rows(full.wald)
  full.n <- bind_rows(full.n)
})

#~~ Merge with allele frequency information

freqinfo <- full.mapfile
names(freqinfo) <- c("Chromosome", "SNP.Name", "cM", "Position", "A1", "A2", "ImputeSuccess")
freqinfo$Q.2 <- NA
freqinfo$No.Measured <- NA

for(i in 1:nrow(freqinfo)){
  
  if(i %in% seq(1, nrow(freqinfo), 100)) print(paste("Running line", i, "of", nrow(freqinfo)))
  x <- data.frame(table(full.imputedgenos[,which(names(full.imputedgenos) == freqinfo$SNP.Name[i])]))
  x$Var1 <- as.character(x$Var1)
  
  x$Order <- 2
  x$Order[which(x$Var1 == paste0(freqinfo$A1[i], "/", freqinfo$A1[i]))] <- 1
  x$Order[which(x$Var1 == paste0(freqinfo$A2[i], "/", freqinfo$A2[i]))] <- 3
  
  x <- arrange(x, Order)
  
  freqinfo$Q.2[i] <- 1 - (x$Freq[1]+(0.5*x$Freq[2]))/sum(x$Freq)
  freqinfo$NoMeasured[i] <- sum(x$Freq)
  rm(x)
}

freqinfo$Order <- 1:nrow(freqinfo)

snp.wald <- filter(full.wald, variable == "SNP")
head(snp.wald)

snp.wald <- left_join(snp.wald, freqinfo)


#~~ Calculate the P values for P < 2e-16

snp.wald$Pr.Chisq.[which(snp.wald$Pr.Chisq. == 0)] <- pchisq(snp.wald$Wald.statistic[which(snp.wald$Pr.Chisq. == 0)], 
                                                             snp.wald$Df[which(snp.wald$Pr.Chisq. == 0)], 
                                                             lower.tail = F)
table(snp.wald$Df)


#~~ Did you get association for everything?

nrow(snp.wald)/6 == nrow(freqinfo)

#~~ missing SNPs

missing.snps <- freqinfo$SNP.Name[which(!freqinfo$SNP.Name %in% snp.wald$SNP.Name)]
# which(snpnames(genabeldata) %in% missing.snps)
# 
# missing.snps <- unique(snp.wald[which(snp.wald$Df == 0),"SNP.Name"])
# missing.snps <- which(snpnames(genabeldata) %in% missing.snps)
# 
# writeLines(paste0("GWAS_", missing.snps[which(missing.snps %in% seq(1, nsnps(genabeldata), 20))], ".sh"), "test.txt")


#~~ quick look

log10(0.05/(39176+nrow(full.mapfile)))

ggplot(snp.wald, aes(Order, -log10(Pr.Chisq.), col = factor(Chromosome))) +
  geom_hline(yintercept = -log10(0.05/(39176+nrow(full.mapfile))), linetype = "dashed") +
  geom_point() +
  facet_grid(Trait ~ LambAdult, scales = "free_y")

ggplot(snp.wald, aes(Order, -log10(Pr.Chisq.), col = ImputeSuccess)) +
  geom_hline(yintercept = -log10(0.05/(39176+nrow(full.mapfile))), linetype = "dashed") +
  geom_point() +
  facet_grid(Trait ~ LambAdult, scales = "free_y")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create information for plotting data   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get Chromosome information for GWAS X-axis

snp.wald$Model <- paste0(snp.wald$Trait, "_", snp.wald$LambAdult)
snp.wald$Chromosome <- as.character(snp.wald$Chromosome)
snp.wald$Chromosome[which(snp.wald$Chromosome == "X")] <- 27
snp.wald$Chromosome <- as.numeric(as.character(snp.wald$Chromosome ))

# 1. Get the effect sizes into the snp.wald.corrected table

head(full.fixef)

full.fixef2 <- full.fixef[grep("SNP", full.fixef$variable),c("Trait", "LambAdult", "variable", "SNP.Name", "solution")]

full.fixef2 <- dcast(full.fixef2, Trait + LambAdult + SNP.Name ~ variable, value.var = "solution")
head(full.fixef2)

snp.wald <- join(snp.wald, full.fixef2)
head(snp.wald)

#~~ Create three columns of effects

snp.wald$EffectAA <- NA
snp.wald$EffectAB <- NA
snp.wald$EffectBB <- NA

snp.wald$A1 <- as.character(snp.wald$A1)
snp.wald$A2 <- as.character(snp.wald$A2)

snp.wald$EffectAA[which(snp.wald$A1 == "A")] <- snp.wald$`SNP_A/A`[which(snp.wald$A1 == "A")]
snp.wald$EffectAA[which(snp.wald$A1 == "C")] <- snp.wald$`SNP_C/C`[which(snp.wald$A1 == "C")]
snp.wald$EffectAA[which(snp.wald$A1 == "G")] <- snp.wald$`SNP_G/G`[which(snp.wald$A1 == "G")]
snp.wald$EffectAA[which(snp.wald$A1 == "T")] <- snp.wald$`SNP_T/T`[which(snp.wald$A1 == "T")]

snp.wald$EffectBB[which(snp.wald$A2 == "A")] <- snp.wald$`SNP_A/A`[which(snp.wald$A2 == "A")]
snp.wald$EffectBB[which(snp.wald$A2 == "C")] <- snp.wald$`SNP_C/C`[which(snp.wald$A2 == "C")]
snp.wald$EffectBB[which(snp.wald$A2 == "G")] <- snp.wald$`SNP_G/G`[which(snp.wald$A2 == "G")]
snp.wald$EffectBB[which(snp.wald$A2 == "T")] <- snp.wald$`SNP_T/T`[which(snp.wald$A2 == "T")]

temp <- which(snp.wald$A1 == "A" &  snp.wald$A2 == "C")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_A/C`[temp]
temp <- which(snp.wald$A1 == "C" &  snp.wald$A2 == "A")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_C/A`[temp]
temp <- which(snp.wald$A1 == "A" &  snp.wald$A2 == "G")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_A/G`[temp]
temp <- which(snp.wald$A1 == "G" &  snp.wald$A2 == "A")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_G/A`[temp]
temp <- which(snp.wald$A1 == "A" &  snp.wald$A2 == "T")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_A/T`[temp]
temp <- which(snp.wald$A1 == "T" &  snp.wald$A2 == "A")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_T/A`[temp]
temp <- which(snp.wald$A1 == "C" &  snp.wald$A2 == "G")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_C/G`[temp]
temp <- which(snp.wald$A1 == "G" &  snp.wald$A2 == "C")
snp.wald$EffectAB[temp] <- snp.wald$`SNP_G/C`[temp]

snp.wald <- snp.wald[,-grep("SNP", names(snp.wald))[-1]]
head(snp.wald)

#~~ Save the raw output from the SNPs

write.table(snp.wald, "../results/3_Imputed_GWAS_Results.txt", row.names = F, sep = "\t", quote = F)

full.fixef.hd <- full.fixef
full.n.hd     <- full.n
full.ranef.hd <- full.ranef
full.wald.hd  <- full.wald


save(full.fixef.hd, full.n.hd, full.ranef.hd, full.wald.hd, file = "../results/3_HD_GWAS_Full_Results.RData")

setwd("..")