library(asreml)
library(GenABEL)

load("data/20181107_BEAST_data_formatted.RData")
load("data/20181107_Genabel_Data.RData")

source("r/ASReml.EstEffects.R")

ainv <- asreml.Ainverse(pedigree)$ginv

extract.snp <- function(snpmodel, snpid){
  
  x <- data.frame(Model = as.character(bquote(snpmodel)),
                  SNP.Name = snpid,
                  Wald = wald.asreml(snpmodel)["SNP", "Wald statistic"],
                  Wald.P = wald.asreml(snpmodel)["SNP", "Pr(Chisq)"])
  
  x1 <- data.frame(summary(snpmodel, all = T)$coef.fixed[grep("SNP", row.names(summary(snpmodel, all = T)$coef.fixed)),])
  x1$Genotype <- row.names(x1)
  
  x <- cbind(x, x1)
  x
}


models <- data.frame(Model = c("1+Sex+LambAgeAugust+SNP,random=~ped(ID)+MumID+BirthYearF",
                               "1+Sex+Age+SNP,random=~ped(ID)+ide(ID)+MumID+BirthYearF+CapYearF"),
                     LambAdult = c("Lambs", "Adults"))

models <- merge(models, data.frame(Response = c("IgAmp", "IgGmp", "IgEmp")))

for(i in 1:3) models[,i] <- as.character(models[,i])

save(list = ls(), file = "brute_GWAS/BEAST_GWAS.RData")





