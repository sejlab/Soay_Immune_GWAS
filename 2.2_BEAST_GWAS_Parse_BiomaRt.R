#
# Parse the GWAS Results
# Nov 2018
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
source("r/biomaRt_func.R")

load("data/20181107_Genabel_Data.RData")

setwd("brute_GWAS/")

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
  
  #print(paste("Running line", i, "of", length(x)))
  
  load(x[i])
  
  full.fixef[[i]] <- restab.fixef
  full.ranef[[i]] <- restab.ranef
  full.wald[[i]]  <- restab.wald
  full.n[[i]]     <- restab.n
  
  rm(restab.ranef, restab.fixef, restab.wald, restab.n)
  
}

setwd("..")

suppressWarnings({
  full.fixef <- bind_rows(full.fixef)
  full.ranef <- bind_rows(full.ranef)
  full.wald  <- bind_rows(full.wald)
  full.n <- bind_rows(full.n)
})


#~~ Merge with allele frequency information

freqinfo <- summary.snp.data(gtdata(genabeldata))
head(freqinfo)
freqinfo$SNP.Name <- row.names(freqinfo)
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

nrow(snp.wald)/9 == nsnps(genabeldata)

#~~ missing SNPs

# missing.snps <- snp.names(genabeldata)[which(!snp.names(genabeldata) %in% snp.wald$SNP.Name)]
# which(snpnames(genabeldata) %in% missing.snps)
# 
# missing.snps <- unique(snp.wald[which(snp.wald$Df == 0),"SNP.Name"])
# missing.snps <- which(snpnames(genabeldata) %in% missing.snps)
# 
# writeLines(paste0("GWAS_", missing.snps[which(missing.snps %in% seq(1, nsnps(genabeldata), 20))], ".sh"), "test.txt")


#~~ quick look

ggplot(snp.wald, aes(Order, -log10(Pr.Chisq.), col = Chromosome)) +
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

x <- subset(snp.wald, Model == unique(snp.wald$Model)[1])
x <- arrange(x, Chromosome, Position)

x$Diff <- c(0, diff(x$Position))
x$Diff <- ifelse(x$Diff < 0, 1, x$Diff)
x$Cumu <- cumsum(x$Diff)

chrinfo <- NULL
for(j in unique(x$Chromosome)){
  temp1 <- subset(x, Chromosome == j)
  temp2 <- data.frame(Chromosome = j,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  chrinfo <- rbind(chrinfo, temp2)
}
chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$Chromosome[28] <- "X"

rm(x, temp1, temp2, j)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Set up BiomaRt information for gene names, GO terms etc.   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Calculate lambda, extract top hits, look at genes and make plots

snp.wald.corrected   <- NULL
genes.in.sig.regions <- NULL
genes.GO.terms       <- NULL


for(i in unique(snp.wald$Model)){
  
  x <- subset(snp.wald, Model == i)
  x <- arrange(x, Chromosome, Position)
  x$Diff <- c(0, diff(x$Position))
  x$Diff <- ifelse(x$Diff < 0, 1, x$Diff)
  x$Cumu <- cumsum(x$Diff)
  
  x <- arrange(x, Pr.Chisq.)
  
  x$Exp.P <- seq(1/nrow(x),1,1/nrow(x))
  
  table(x$Df)
  
  mod.lambda <- median(x$Wald.statistic, na.rm = T)/median(qchisq(x$Exp.P, 2, lower.tail = F))
  
  #mod.lambda <- median(subset(x, df == 2)$Wald, na.rm = T)/median(qchisq(subset(x, df == 2)$Exp.P, 2, lower.tail = F))
  
  #mod.lambda <- median(x$Pr.Chisq., na.rm = T)/median(qchisq(x$Exp.P, x$Df, lower.tail = F))
  
  print(paste("Lambda for ", i, "is", mod.lambda))
  
  x$Wald.P.corrected <- pchisq(x$Wald.statistic/mod.lambda, 2, lower.tail = F)
  
  x <- arrange(x, Wald.P.corrected)
  
  x$Exp.P <- seq(1/nrow(x),1,1/nrow(x))
  
  #~~ Association plot
  
  p1 <- ggplot(x, aes(Cumu, -log10(Wald.P.corrected), col = factor(Chromosome %% 2))) + 
    geom_hline(yintercept = -log10(2.245e-6), linetype = "dashed") +
    geom_point() +
    scale_colour_brewer(palette = "Set1") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chromosome) +
    labs(x = "Chromosome", y = "-log10 P") +
    ggtitle(i)
  
  #~~ PP Plot
  
  p2 <- ggplot(x, aes(-log10(Exp.P), -log10(Wald.P.corrected)))+ 
    geom_hline(yintercept = -log10(2.245e-6), linetype = "dashed") +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    scale_colour_brewer(palette = "Set1") +
    theme(legend.position = "none") +
    ggtitle(i)
  
  print(p1)
  print(p2)
  
  png(paste0("figs/2_GWAS_", i, ".png"), width = 9, height = 4, units = "in", res = 300)
  multiplot(p1, p2, layout = matrix(c(1,1,2), nrow=1))
  dev.off()
  
  
  #~~ SNPs above the threshold
  
  y <- filter(x, Wald.P.corrected < 2.245e-6) %>% select(-Diff, -Cumu, -Exp.P) %>% arrange(Wald.P.corrected)
  
  print(y)
  
  #~~ Genes in region
  
  if(nrow(y) > 0){
    
    gene.temp <- NULL
    
    for(j in 1:nrow(y)){
      
      gene.temp <- rbind(gene.temp,
                         subset(gene.list, Chromosome == y$Chromosome[j] & Stop > y$Position[j] - 1e6 & Start < y$Position[j] + 1e6))
      
    }
    
    gene.temp <- unique(gene.temp)
    gene.temp <- arrange(gene.temp, Chromosome, Start)
    gene.temp$Model <- i
    
    gene.info <- genedesc(gene.temp$gene_name)
    
    genes.GO.terms <- rbind(genes.GO.terms, gene.info)
    
    gene.info <- unique(subset(gene.info, select = c(Species, external_gene_name, description)))
    gene.info <- separate(gene.info, description, into = c("description", "delete"), remove = T, sep = "\\[")
    gene.info$delete <- NULL
    
    print(gene.temp)
    print(gene.info)
    genes.in.sig.regions <- rbind(genes.in.sig.regions, gene.temp)
    
  }
  #~~ Add results to table
  
  snp.wald.corrected <- rbind(snp.wald.corrected, x)
  
  rm(x, mod.lambda, gene.temp, gene.info, y, j)
}

# 1. Get the effect sizes into the snp.wald.corrected table

head(full.fixef)

full.fixef2 <- full.fixef[grep("SNP", full.fixef$variable),c("Trait", "LambAdult", "variable", "SNP.Name", "solution")]

full.fixef2 <- dcast(full.fixef2, Trait + LambAdult + SNP.Name ~ variable, value.var = "solution")
head(full.fixef2)

snp.wald.corrected <- join(snp.wald.corrected, full.fixef2)
head(snp.wald.corrected)

#~~ Create three columns of effects

snp.wald.corrected$EffectAA <- NA
snp.wald.corrected$EffectAB <- NA
snp.wald.corrected$EffectBB <- NA

snp.wald.corrected$A1 <- as.character(snp.wald.corrected$A1)
snp.wald.corrected$A2 <- as.character(snp.wald.corrected$A2)

snp.wald.corrected$EffectAA[which(snp.wald.corrected$A1 == "A")] <- snp.wald.corrected$`SNP_A/A`[which(snp.wald.corrected$A1 == "A")]
snp.wald.corrected$EffectAA[which(snp.wald.corrected$A1 == "C")] <- snp.wald.corrected$`SNP_C/C`[which(snp.wald.corrected$A1 == "C")]
snp.wald.corrected$EffectAA[which(snp.wald.corrected$A1 == "G")] <- snp.wald.corrected$`SNP_G/G`[which(snp.wald.corrected$A1 == "G")]
snp.wald.corrected$EffectAA[which(snp.wald.corrected$A1 == "T")] <- snp.wald.corrected$`SNP_T/T`[which(snp.wald.corrected$A1 == "T")]

snp.wald.corrected$EffectBB[which(snp.wald.corrected$A2 == "A")] <- snp.wald.corrected$`SNP_A/A`[which(snp.wald.corrected$A2 == "A")]
snp.wald.corrected$EffectBB[which(snp.wald.corrected$A2 == "C")] <- snp.wald.corrected$`SNP_C/C`[which(snp.wald.corrected$A2 == "C")]
snp.wald.corrected$EffectBB[which(snp.wald.corrected$A2 == "G")] <- snp.wald.corrected$`SNP_G/G`[which(snp.wald.corrected$A2 == "G")]
snp.wald.corrected$EffectBB[which(snp.wald.corrected$A2 == "T")] <- snp.wald.corrected$`SNP_T/T`[which(snp.wald.corrected$A2 == "T")]

temp <- which(snp.wald.corrected$A1 == "A" &  snp.wald.corrected$A2 == "C")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_A/C`[temp]
temp <- which(snp.wald.corrected$A1 == "C" &  snp.wald.corrected$A2 == "A")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_C/A`[temp]
temp <- which(snp.wald.corrected$A1 == "A" &  snp.wald.corrected$A2 == "G")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_A/G`[temp]
temp <- which(snp.wald.corrected$A1 == "G" &  snp.wald.corrected$A2 == "A")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_G/A`[temp]
temp <- which(snp.wald.corrected$A1 == "A" &  snp.wald.corrected$A2 == "T")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_A/T`[temp]
temp <- which(snp.wald.corrected$A1 == "T" &  snp.wald.corrected$A2 == "A")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_T/A`[temp]
temp <- which(snp.wald.corrected$A1 == "C" &  snp.wald.corrected$A2 == "G")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_C/G`[temp]
temp <- which(snp.wald.corrected$A1 == "G" &  snp.wald.corrected$A2 == "C")
snp.wald.corrected$EffectAB[temp] <- snp.wald.corrected$`SNP_G/C`[temp]

snp.wald.corrected <- snp.wald.corrected[,-grep("SNP", names(snp.wald.corrected))[-1]]
head(snp.wald.corrected)

#~~ Save the raw output from the SNPs

write.table(snp.wald.corrected, "results/2_GWAS_Results.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Format the gene tables                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sig.snps <- subset(snp.wald.corrected, Pr.Chisq. < 2.245e-6)


#~~ Get the significant regions for imputation

sig.snps <- unique(subset(sig.snps, select = c(Chromosome, Position)))
sig.snps <- arrange(sig.snps, Chromosome, Position)
sig.snps <- subset(sig.snps, Chromosome != 0)

write.table(sig.snps, "results/2_Regions_for_Imputation.txt", row.names = F, sep = "\t", quote = F)

save(full.fixef, full.n, full.ranef, full.wald, snp.wald.corrected, sig.snps, file  = "results/2_GWAS_Full_Results.RData")
