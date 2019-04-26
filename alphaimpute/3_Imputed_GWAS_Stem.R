

Sys.setenv("ASREML_LICENSE_FILE" = "/exports/igmm/software/pkg/el7/apps/asreml/3.0gm/bin/asreml.lic") 

library(asreml)
library(reshape)
library(GenABEL)
library(plyr)

load("BEAST_GWAS.RData")
load("Imputed_GWAS.RData")

rm(genabeldata)

full.imputedgenos <- full.imputedgenos[,c(1, which(names(full.imputedgenos) %in% full.mapfile$V2[startsnp:(startsnp+99)]))]

gc()

#~~ Prepare data frame for results

restab.ranef <- NULL
restab.fixef <- NULL
restab.wald <- NULL
restab.n <- NULL

for(h in startsnp:(startsnp+99)){

  print(paste("Running SNP", h))
  
  genotab <- full.imputedgenos[,c(1, which(names(full.imputedgenos) == full.mapfile$V2[h]))]
  names(genotab)[2] <- "SNP"
  
  ped.results.ranef <- NULL
  ped.results.fixef <- NULL
  ped.results.wald <- NULL
  ped.results.n <- NULL
  
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
    
    x.data <- join(x.data, genotab)
    x.data <- na.omit(x.data[,x.vec])
    #x.data <- droplevels(subset(x.data, ID %in% dimnames(grminv)[[1]]))
    
    eval(parse(text = paste0("fit1 <- asreml(fixed=", models$Response[i], "~",
                             models$Model[i],"
                             , data=x.data, ginverse=list(ID=ainv),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)")))
    
    
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
    rm(x)
    
    x <- data.frame(wald.asreml(fit1))
    x$variable <- row.names(x)
    
    ped.results.wald <- rbind(ped.results.wald,
                               cbind(Trait = models$Response[i],
                                     LambAdult = models$LambAdult[i],
                                     x))
    
    ped.results.n <- rbind(ped.results.n, 
                         data.frame(Trait = models$Response[i],
                                    LambAdult = models$LambAdult[i],
                                    table(x.data$SNP), 
                                    table(unique(subset(x.data, select = c(ID, SNP)))$SNP)))
    
    
    
    
    rm(fit1, x.data, x.vec, x)
    
  }
  

  ped.results.ranef$SNP.Name <- full.mapfile$V2[h]
  ped.results.fixef$SNP.Name <- full.mapfile$V2[h]
  ped.results.wald$SNP.Name <- full.mapfile$V2[h]
  ped.results.n$SNP.Name <- full.mapfile$V2[h]
  
  
  restab.ranef <- rbind(restab.ranef, ped.results.ranef)
  restab.fixef <- rbind(restab.fixef, ped.results.fixef)
  restab.wald <- rbind(restab.wald, ped.results.wald)
  restab.n <- rbind(restab.n, ped.results.n)
  
  rm(ped.results.ranef, ped.results.fixef, ped.results.wald, ped.results.n, genotab, i)
  
  save(restab.ranef, restab.fixef, restab.wald, restab.n, file = paste0("GWAS_", startsnp, ".RData"))
  
}




