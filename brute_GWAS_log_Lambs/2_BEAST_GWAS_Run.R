load("BEAST_GWAS.RData")

library(GenABEL)

for(i in seq(1, nsnps(genabeldata), 100)){
  
  x <- readLines("2_BEAST_GWAS_Stem.R")
  writeLines(c(paste0("startsnp = ", i), x), paste0("GWAS_", i, ".R"))
  
  writeLines(paste0("#!/bin/sh
#$ -cwd
#$ -l h_rt=01:00:00
#$ -V
#$ -l h_vmem=5200M
                    
. /etc/profile.d/modules.sh
module load R
                    
R CMD BATCH GWAS_", i, ".R"),
paste0("GWAS_", i, ".sh"))
  
  system(paste0("qsub GWAS_", i, ".sh"))
}