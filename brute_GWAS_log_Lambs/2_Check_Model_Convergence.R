setwd("brute_GWAS/")

load("BEAST_GWAS.RData")

converge.tab <- NULL

for(i in seq(1, nsnps(genabeldata), 100)){
  
  x <- readLines(paste0("BEAST_GWAS_outputs/GWAS_", i, ".Rout"))
  x <- x[grep("Converged", x, ignore.case = T)]
  x <- data.frame(x = x, Order = 1:length(x), Script = i)

  converge.tab <- rbind(converge.tab, x)
}
  