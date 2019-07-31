for(i in 1:15){
  x <- readLines("1.3_Gen_Corr_Stem.R")
  writeLines(c(paste0("i = ", i), x), paste0("Gen_Corr_", i, ".R"))
}

system("qsub -t 1-15 Gen_Corr_Script.sh")
