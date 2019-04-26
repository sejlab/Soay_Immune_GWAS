#
# Genotype imputation
# SEJ, AMS
# Nov 2018
#
#

library(plyr)

prepareRun <- FALSE

sig.regions <- read.table("results/2_Regions_for_Imputation_log_lambs.txt", header = T, stringsAsFactors = F)
sig.regions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Prepare the data                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(prepareRun == TRUE){
  
  #~~ Format pedigree
  
  pedigree <- read.table("data/Pedigree_2018-02-12.ped", header = T, stringsAsFactors = F)
  head(pedigree)
  
  pedigree[is.na(pedigree)] <- 0
  pedigree$MOTHER[which(pedigree$ID == 46)] <- 0
  
  write.table(pedigree[,c("ID", "FATHER", "MOTHER")], "alphaimpute/pedigree.txt", row.names = F, col.names = F, quote = F)
  
  #~~ Format the HD file
  
  system("plink --file data/20140214_SheepHD_QC1_Polym --sheep --make-just-fam --out data/20140214_SheepHD_QC1_Polym")
  
  HDfam <- read.table("data/20140214_SheepHD_QC1_Polym.fam")
  names(HDfam)[c(1, 3:6)] <- paste0(names(HDfam)[c(1, 3:6)], ".HD")
  
  oldfam <- read.table("data/20180209_SoayPlates1-83.fam")
  oldfam <- join(oldfam, HDfam)
  
  update.sex <- oldfam[which(oldfam$V5 != oldfam$V5.HD),c("V1", "V2", "V5")]
  update.ids <- cbind(HDfam[,1:2], 1, HDfam[,2])
  
  write.table(update.sex, "data/HD_UpdateSex.txt", row.names = F, col.names = F, quote = F)
  write.table(update.ids, "data/HD_UpdateID.txt", row.names = F, col.names = F, quote = F)
  
  system("plink --file data/20140214_SheepHD_QC1_Polym --sheep --update-ids data/HD_UpdateID.txt --recode --out data/20140214_SheepHD_QC1_Polym")
  system("plink --file data/20140214_SheepHD_QC1_Polym --sheep --update-sex data/HD_UpdateSex.txt --recode --out data/20140214_SheepHD_QC1_Polym")
  
  system("rm data\\HD_UpdateSex.txt")
  system("rm data\\HD_UpdateID.txt")
  
  rm(pedigree, oldfam, HDfam, update.ids, update.sex)
  
}

pedigree <- read.table("alphaimpute/pedigree.txt", header = F, stringsAsFactors = F)

#~~ Get region information


sig.regions <- data.frame(Min = tapply(sig.regions$Position, sig.regions$Chromosome, min),
                          Max = tapply(sig.regions$Position, sig.regions$Chromosome, max))
sig.regions$Chromosome <- row.names(sig.regions)

sig.regions$Min <- sig.regions$Min - 2e6
sig.regions$Max <- sig.regions$Max + 2e6


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Run AlphaImpute                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for(h in 1:nrow(sig.regions)){
  system(paste0("mkdir alphaimpute/chr", sig.regions$Chromosome[h]))
  
  #~~ Make new files to merge
  
  dirinfo <- paste0("alphaimpute/chr", sig.regions$Chromosome[h], "/chr", sig.regions$Chromosome[h])
  
  system(paste0("plink --bfile data/20180209_SoayPlates1-83 --sheep --recode --chr ", sig.regions$Chromosome[h], " --from-bp ", sig.regions$Min[h], " --to-bp ", sig.regions$Max[h], " --out ", dirinfo))
  system(paste0("plink --file data/20140214_SheepHD_QC1_Polym --sheep --recode --chr ", sig.regions$Chromosome[h], " --from-bp ", sig.regions$Min[h], " --to-bp ", sig.regions$Max[h], " --out ", dirinfo, "HD"))
  
  system(paste0("plink --file ", dirinfo, "  --sheep --merge ", dirinfo, "HD --out ", dirinfo, "merge"))
  system(paste0("plink --bfile ", dirinfo, "merge --sheep --recode A --output-missing-genotype 9 --out ", dirinfo))
  
  #~~ remove header, family/parent/sex/pheno, recode NA's to 9
  
  setwd(paste0("alphaimpute/chr", sig.regions$Chromosome[h]))
  
  system(paste0("sed -i \"1d\" chr", sig.regions$Chromosome[h], ".raw")) # remove header
  system("cmd", input = paste0("cut -d ' ' -f2,7- chr", sig.regions$Chromosome[h], ".raw > chr", sig.regions$Chromosome[h], ".2.raw"))
  system("cmd", input = paste0("sed \"s/NA/9/g\" chr", sig.regions$Chromosome[h], ".2.raw > chr", sig.regions$Chromosome[h], ".raw"))
  system(paste0("rm chr", sig.regions$Chromosome[h], ".2.raw"))
  
  x <- dir()[grep(paste0("chr", sig.regions$Chromosome[h]), dir())]
  x <- x[-grep("raw", x)]
  x <- x[-grep("merge.bim", x)]
  for(i in x) system(paste0("rm ", i))
  
  #~~ Make spec file
  
  InternalIterations <- 20 # should be about 5, but lower values make it go faster
  
  nsnps <- readLines(paste0("chr", sig.regions$Chromosome[h], ".raw"), n = 1)
  nsnps <- length(strsplit(nsnps, split = " ")[[1]][-1])
  
  CoreAndTailLengths <- "300,350,400,450"
  CoreLengths        <- "250,300,350,400"
  
  if(nsnps < 450){
    x <- seq(50, nsnps, 50)
    if(length(x) >= 5){
      CoreAndTailLengths <- paste(x[(length(x)-3):length(x)], collapse = ",")
      CoreLengths <- paste(x[(length(x)-3):length(x)]-50, collapse = ",")
    } else {
      stop("Window too narrow")
    }
  }
  
  write.table(pedigree, "pedigree.txt", row.names = F, col.names = F, quote = F)
  
  writeLines(paste0("
                    = BOX 1: Input Files ==========================================================
                    PedigreeFile                ,pedigree.txt
                    GenotypeFile                ,chr", sig.regions$Chromosome[h], ".raw
                    SexChrom                    ,No
                    = BOX 2: SNPs ==================================================================
                    NumberSnp                   ,", nsnps, "
                    MultipleHDPanels            ,No
                    = BOX 3: Filtering =======================================================
                    InternalEdit                ,No
                    EditingParameters           ,0.0,0.0,0.0,AllSnpOut
                    = BOX 4: Phasing ================================================================
                    HDAnimalsThreshold          ,90.0
                    NumberPhasingRuns           ,4
                    CoreAndTailLengths          ,", CoreAndTailLengths, "
                    CoreLengths                 ,", CoreLengths, "
                    PedigreeFreePhasing         ,No
                    GenotypeError               ,0.0
                    = BOX 5: Imputation =========================================================
                    InternalIterations          ,", InternalIterations, "
                    ConservativeHaplotypeLibraryUse     ,No
                    WellPhasedThreshold         ,99.0
                    = BOX 6: Hidden Markov Model ================================================
                    HMMOption                   ,No
                    TemplateHaplotypes          ,200
                    BurnInRounds                ,5
                    Rounds                      ,20
                    ParallelProcessors          ,4
                    Seed                        ,-123456789
                    ThresholdForMissingAlleles  ,50.0
                    ThresholdImputed            ,90.0
                    = BOX 7: Running options ====================================================
                    ParallelProcessors          ,2 
                    PreprocessDataOnly          ,No
                    PhasingOnly                 ,No
                    UserDefinedAlphaPhaseAnimalsFile    ,None
                    PrePhasedFile               ,None
                    RestartOption               ,0
                    OutputOnlyGenotypedAnimals  , Yes"
                    
  ), con = "spec.txt"
  )
  
  system("cmd", input = "..\\AlphaImpute.exe spec.txt")

  setwd("../..")
  rm(dirinfo, nsnps, x, i)
}


