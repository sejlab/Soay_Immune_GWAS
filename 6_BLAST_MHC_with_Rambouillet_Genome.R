# BLAST of MHC region SNPs

library(plyr)
library(dplyr)
library(magrittr)
library(seqinr)
library(ggplot2)

#~~ Get MHC Snps and write to file

mhc.snps <- read.table("results/3_Imputed_GWAS_Results.txt", header = T, stringsAsFactors = F)
head(mhc.snps)

mhc.snps <- subset(mhc.snps, Chromosome == 20) %>%
  subset(select = c(SNP.Name, Chromosome, Position)) %>%
  unique() %>%
  arrange(Position)

head(mhc.snps)

setwd("../Genome Sequences/Sheep Genome Rambouillet/")

#~~ Extract flanking sequences of those SNPs

x <- read.fasta("SNPChip_Flanking_Seqs.fasta")
x <- x[which(names(x) %in% mhc.snps$SNP.Name)]

write.fasta(x, names(x), file.out = "MhcSNPs.fasta")

#~~ Make blast DBs and BLAST it

system("cmd", input = "makeblastdb -in MhcSNPs.fasta  -dbtype nucl")

system("cmd", input = "blastn -db GCA_002742125.1_Oar_rambouillet_v1.0_genomic.fna -query MhcSNPs.fasta -out MhcSNPs.out -outfmt 6")

#~~ Load results

mhc.blast <- read.table("MhcSNPs.out", stringsAsFactors = F)
str(mhc.blast)

ggplot(mhc.blast, aes(V3)) + geom_histogram(binwidth = 0.5)
ggplot(mhc.blast, aes(V4)) + geom_histogram(binwidth = 1)

table(table(mhc.blast$V1))

mhc.blast <- subset(mhc.blast, V4 > 150 & V3 > 95 & V2 == "CM008491.1")

table(table(mhc.blast$V1))
table(mhc.blast$V2)

mhc.blast <- join(mhc.blast, plyr::count(mhc.blast, "V1"))

head(mhc.blast)

#~~ Look at positions 

names(mhc.blast)[1] <- c("SNP.Name")
mhc.blast <- join(mhc.blast, mhc.snps)


ggplot(mhc.blast, aes(Position/1e6, V9/1e6)) + 
  geom_point(alpha = 0.2) +
  theme_bw() +
  labs(x = "Oar_v3.1 Chr 20 Position", y = "Rambouillet Chr 20 Position")
  
  


