#
# How much Va is attributed to associated regions?
# SEJ, AMS
# Jan 2019
#
# 

library(asreml)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)

source("r/makeGRM.R")
source("r/ASReml.EstEffects.R")

load("brute_gwas/BEAST_GWAS.RData")
load("alphaimpute/Imputed_GWAS.RData", verbose = T)
load("alphaimpute/Imputed_GWAS_log_Lambs.RData", verbose = T)

BEASTX$logIgEmp <- log10(BEASTX$IgEmp)
BEASTX$logIgEmp[which(BEASTX$logIgEmp == -Inf)] <- NA

names(full.imputedgenos.log.lambs)[1] <- "ID"

ainv <- asreml.Ainverse(pedigree)$ginv

focalsnps <- read.table("results/4_Top_SNPs_for_Paper.txt", header = T, sep = "\t", stringsAsFactors = F)
focalsnps <- focalsnps[,c("Trait", "LambAdult", "SNP.Name")]

models <- data.frame(Model = c("1+Sex+SNP+LambAgeAugust,random=~ped(ID)+MumID+BirthYearF",
                               "1+Sex+SNP+Age,random=~ped(ID)+ide(ID)+MumID+BirthYearF+CapYearF"),
                     LambAdult = c("Lambs", "Adults"))

focalsnps <- join(focalsnps, models)
focalsnps$Trait[which(focalsnps$Trait == "IgEmp" & focalsnps$LambAdult == "Lambs")] <- "logIgEmp"

full.imputedgenos <- full.imputedgenos[,c(1, which(names(full.imputedgenos) %in% focalsnps$SNP.Name))]
full.imputedgenos.log.lambs <- full.imputedgenos.log.lambs[,c(1, which(names(full.imputedgenos.log.lambs) %in% focalsnps$SNP.Name))]

full.imputedgenos <- join(full.imputedgenos, full.imputedgenos.log.lambs)

rm(full.imputedgenos.log.lambs, full.mapfile, models)


ped.results.fixef <- NULL
ped.results.wald <- NULL
ped.results.n <- NULL

for(i in 1:nrow(focalsnps)){
  
  x.vec <- as.character(focalsnps$Model[i])
  x.vec <- gsub(",random=~", "+", x.vec)
  x.vec <- strsplit(x.vec, split = "\\+")[[1]]
  x.vec[grep("ped", x.vec)] <- "ID"
  x.vec[grep("ide", x.vec)] <- "ID"
  x.vec <- c(x.vec, as.character(focalsnps$Trait[i]))
  x.vec <- unique(x.vec[-which(x.vec == 1)])
  
  
  if(focalsnps$LambAdult[i] %in% c("Lambs", "Adults")){
    x.data <- subset(BEASTX, LambAdult == focalsnps$LambAdult[i])
  } else {
    x.data <- BEASTX
  }
  
  genotab <- full.imputedgenos[,c("ID", focalsnps$SNP.Name[i])]
  names(genotab)[2] <- "SNP"
  
  x.data <- join(x.data, genotab)
  
  x.data <- droplevels(na.omit(x.data[,x.vec]))
  
  x.data <- droplevels(subset(x.data, ID %in% pedigree$ID))
  
  if(focalsnps$SNP.Name[i] == "oar3_OAR10_10333145"){
    x.data <- subset(x.data, SNP != "A/A") %>% droplevels
  }
  if(focalsnps$SNP.Name[i] == "oar3_OAR16_12632988"){
    x.data <- subset(x.data, SNP != "G/G") %>% droplevels
  }
  
  newmod <- as.character(focalsnps$Model[i])
  
  newmod <- gsub("Sex+SNP", "Sex:SNP", newmod, fixed = T)
  
  eval(parse(text = paste0("fit1 <- asreml(fixed=", focalsnps$Trait[i], "~",
                           newmod,"
                             , data=x.data, ginverse=list(ID=ainv),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)")))
  
  x <- data.frame(summary(fit1, all = T)$coef.fixed)
  x$variable <- row.names(x)
  
  ped.results.fixef <- rbind(ped.results.fixef,
                             cbind(Trait = focalsnps$Trait[i],
                                   LambAdult = focalsnps$LambAdult[i],
                                   SNP.Name = focalsnps$SNP.Name[i],
                                   x))
  rm(x)
  
  x <- data.frame(wald.asreml(fit1))
  x$variable <- row.names(x)
  
  ped.results.wald <- rbind(ped.results.wald,
                            cbind(Trait = focalsnps$Trait[i],
                                  LambAdult = focalsnps$LambAdult[i],
                                  SNP.Name = focalsnps$SNP.Name[i],
                                  x))
  
  ped.results.n <- rbind(ped.results.n, 
                         data.frame(Trait = focalsnps$Trait[i],
                                    LambAdult = focalsnps$LambAdult[i],
                                    SNP.Name = focalsnps$SNP.Name[i],
                                    table(x.data$SNP, x.data$Sex)))
  
  
}

ped.results.fixef.hold <- ped.results.fixef
ped.results.n.hold <- ped.results.n

ped.results.fixef <- ped.results.fixef[grep("SNP", ped.results.fixef$variable),]
ped.results.fixef <- separate(ped.results.fixef, variable, c("Sex", "SNP.Genotype"), sep = ":")
ped.results.fixef$SNP.Genotype[is.na(ped.results.fixef$SNP.Genotype)] <- ped.results.fixef$Sex[is.na(ped.results.fixef$SNP.Genotype)]
ped.results.fixef$Sex <- gsub("Sex_", "", ped.results.fixef$Sex)
ped.results.fixef$SNP.Genotype <- gsub("SNP_", "", ped.results.fixef$SNP.Genotype)
ped.results.fixef$Model <- paste0(gsub("s", "", ped.results.fixef$LambAdult), " ", 
                                  gsub("mp", "", ped.results.fixef$Trait), "\n", 
                                  ped.results.fixef$SNP.Name)
ped.results.fixef$Sex[grep("SNP", ped.results.fixef$Sex)] <- "Both"

ped.results.fixef$SNP.Genotype[which(ped.results.fixef$SNP.Genotype == "G/A")] <- "A/G"
ped.results.fixef$SNP.Genotype.num <- ped.results.fixef$SNP.Genotype
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$SNP.Genotype.num == "A/A")] <- 1
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$SNP.Genotype.num == "A/G")] <- 2
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$SNP.Genotype.num == "G/A")] <- 2
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$SNP.Genotype.num == "G/G")] <- 3
ped.results.fixef$SNP.Genotype.num <- as.numeric(ped.results.fixef$SNP.Genotype.num)
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$Sex == "F")] <- ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$Sex == "F")]-0.05
ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$Sex == "M")] <- ped.results.fixef$SNP.Genotype.num[which(ped.results.fixef$Sex == "M")]+0.05
ped.results.fixef$std.error[is.na(ped.results.fixef$std.error)] <- 0

ped.results.n$SNP.Genotype.num <- as.character(ped.results.n$Var1)
ped.results.n$SNP.Genotype.num[which(ped.results.n$SNP.Genotype.num == "A/A")] <- 1
ped.results.n$SNP.Genotype.num[which(ped.results.n$SNP.Genotype.num == "A/G")] <- 2
ped.results.n$SNP.Genotype.num[which(ped.results.n$SNP.Genotype.num == "G/A")] <- 2
ped.results.n$SNP.Genotype.num[which(ped.results.n$SNP.Genotype.num == "G/G")] <- 3

ped.results.n$Model <- paste0(gsub("s", "", ped.results.n$LambAdult), " ", 
                                  gsub("mp", "", ped.results.n$Trait), "\n", 
                              ped.results.n$SNP.Name)
ped.results.n$Label <- paste0(ped.results.n$Var2, ": ", ped.results.n$Freq)

library(reshape2)
ped.results.n <- subset(ped.results.n, select = c(Model, SNP.Genotype.num, Var2, Label))
ped.results.n <- dcast(ped.results.n, Model + SNP.Genotype.num ~ Var2)
ped.results.n$Label <- paste0(ped.results.n$F, "\n", ped.results.n$M)

x <- ped.results.fixef %>% group_by(Model) %>% summarise(ymax = max(solution + std.error, na.rm= T),
                                                         ymin = min(solution - std.error, na.rm= T))
ped.results.n <- join(ped.results.n, x)
ped.results.n$Range = ped.results.n$ymax - ped.results.n$ymin

str(ped.results.n)
ped.results.n$SNP.Genotype.num <- as.numeric(ped.results.n$SNP.Genotype.num)
#ped.results.n$y1 <- ped.results.n$y1 - 0.1

ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(data = ped.results.fixef.sex, aes(SNP.Genotype.num, solution, col = Sex, shape = Sex), size = 2) +
  geom_line(data = ped.results.fixef.sex, aes(SNP.Genotype.num, solution, col = Sex)) +
  geom_errorbar(data = ped.results.fixef.sex, aes(x = SNP.Genotype.num, col = Sex, ymin = solution - std.error, ymax = solution + std.error), width = 0) +
  geom_text(data = ped.results.n, aes(SNP.Genotype.num, ymax + 0.2*Range, label = Label)) +
  geom_point(data = ped.results.n, aes(SNP.Genotype.num,ymax + 0.4*Range), alpha = 0) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  scale_x_continuous(breaks = 1:3, labels = c("A/A", "A/G", "G/G")) +
  facet_wrap(~Model, scales = "free", ncol = 2) +
  labs(x = "SNP Genotype", y = "Effect Size (Relative to Model Intercept)") +
  coord_cartesian(xlim = c(0.5, 3.5)) +
  theme(legend.position = "top") 

ggsave("figs/7_Sex_Specific_SNP_Effects.png", width = 6, height = 10)

write.table(ped.results.wald , "results/7_Allele_effects_sex_Specific_wald.txt", row.names = F, sep = "\t", quote = F)
write.table(ped.results.fixef, "results/7_Sex_specific_fixef.txt", row.names = F, sep = "\t", quote = F)
