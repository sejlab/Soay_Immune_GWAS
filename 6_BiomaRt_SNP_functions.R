
topsnps <- read.table("results/4_Information_on_Top_SNPs.txt", header = T, sep = "\t", stringsAsFactors = F)
head(topsnps)

library(biomaRt)

mart <- useEnsembl("snp")
dataset <- useDataset("oaries_snp", mart=mart)
filters <- listFilters(dataset)

attributes <- listAttributes(dataset)

chr.vec <- sort(unique(topsnps$Chromosome))
chr.vec <- chr.vec[-1]

snp.hold <- NULL

for(i in chr.vec){
  
  print(paste("Running Chromosome", i))
  
  snps.chr <- getBM(attributes=c('refsnp_id', 'chrom_start', 'consequence_type_tv', 'distance_to_transcript'),
                                  filters='chr_name',
                                  values=i,
                                  mart=dataset)
  
  snp.hold <- rbind(snp.hold,
                    #snps.chr[which(snps.chr$chrom_start %in% topsnps$Position),])
                    snps.chr)
  
  save(snp.hold, file = "temp.RData")
  
  rm(snps.chr)
  gc()
  
  
}

# system.time(snps.chr24 <- getBM(attributes=c('refsnp_id', 'chrom_start', 'consequence_type_tv', 'distance_to_transcript'),
#                     filters=c('chr_name', 'start'),
#                     values=list(24, top.24$Position),
#                     mart=dataset))



