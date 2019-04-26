library(biomaRt)
library(plyr)

source("r/biomaRt_func.R")

snp.wald.corrected <- read.table("results/2_GWAS_Results.txt", header = T, stringsAsFactors = F)

#~~ Deal with log IgE lamb data

snp.ige <- read.table("results/2_GWAS_Results_log_Lambs.txt", header = T, sep = "\t", stringsAsFactors = F)
snp.ige <- subset(snp.ige, Model == "IgEmp_Lambs")

snp.wald.corrected <- subset(snp.wald.corrected, Model != "IgEmp_Lambs")
snp.wald.corrected <- rbind(snp.wald.corrected, snp.ige)

rm(snp.ige)

#~~ Genes

gene.list <- read.table("data/Ovis_aries.Oar_v3.1.94.genes.txt", header = T, stringsAsFactors = F)

snp.wald.corrected <- subset(snp.wald.corrected, Pr.Chisq. < 2.245e-6)

#~~ extract genes in significant regions

genes.in.sig.regions <- NULL

for(i in 1:nrow(snp.wald.corrected)){
  
  genes.in.sig.regions <- rbind(genes.in.sig.regions,
                                subset(gene.list, Chromosome == snp.wald.corrected$Chromosome[i] & 
                                         Stop > snp.wald.corrected$Position[i] - 1e6 & 
                                         Start < snp.wald.corrected$Position[i] + 1e6)
  )
  
}

genes.in.sig.regions <- unique(genes.in.sig.regions)

genes.in.sig.regions <- arrange(genes.in.sig.regions, Chromosome, Start)


#~~ find orthologues for genes as many unidentified...

orthotab <- data.frame(sheep_gene_id = subset(genes.in.sig.regions, is.na(gene_name))$gene_id)

gene.orthos <- NULL

# mice
x <- getLDS(attributes=c("ensembl_gene_id"),
            filters="ensembl_gene_id", values=orthotab$sheep_gene_id, mart=ensembl_oaries,
            attributesL=c("ensembl_gene_id"), martL=ensembl_mmusculus)

names(x) <- c("sheep_gene_id", "ensembl_gene_id")
x$Species <- "mmmusculus"

x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = x$ensembl_gene_id, mart = ensembl_mmusculus))

gene.orthos <- rbind(gene.orthos, x)
rm(x)

# humans
x <- getLDS(attributes=c("ensembl_gene_id"),
            filters="ensembl_gene_id", values=orthotab$sheep_gene_id, mart=ensembl_oaries,
            attributesL=c("ensembl_gene_id"), martL=ensembl_hsapiens)

names(x) <- c("sheep_gene_id", "ensembl_gene_id")
x$Species <- "hsapiens"

x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = x$ensembl_gene_id, mart = ensembl_hsapiens))

gene.orthos <- rbind(gene.orthos, x)
rm(x)

# cattle
x <- getLDS(attributes=c("ensembl_gene_id"),
            filters="ensembl_gene_id", values=orthotab$sheep_gene_id, mart=ensembl_oaries,
            attributesL=c("ensembl_gene_id"), martL=ensembl_btaurus)

names(x) <- c("sheep_gene_id", "ensembl_gene_id")
x$Species <- "btaurus"

x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = x$ensembl_gene_id, mart = ensembl_btaurus))

gene.orthos <- rbind(gene.orthos, x)
rm(x)

#~~ Format the table:

x <- gene.orthos[,c("sheep_gene_id", "external_gene_name")]
x$external_gene_name <- toupper(x$external_gene_name)
x <- unique(x)
x <- subset(x, external_gene_name != "")
x$external_gene_name

x1 <- data.frame(table(x$sheep_gene_id))
names(x1)[1] <- "sheep_gene_id"


x <- join(x, x1)
rm(x1)

names(x) <- c("gene_id", "orthologue", "orthologue_count")


genes.in.sig.regions <- join(genes.in.sig.regions, x)
genes.in.sig.regions$consensus_locus <- ifelse(is.na(genes.in.sig.regions$gene_name), genes.in.sig.regions$orthologue, genes.in.sig.regions$gene_name)


rm(x, i, orthotab, gene.list)


#~~ Get the GO terms

GO.tab <- genedesc(unique(genes.in.sig.regions$consensus_locus))
head(GO.tab)

temp <- sort(unique(c(grep("immun", GO.tab$description),
                      grep("immun", GO.tab$phenotype_description),
                      grep("immun", GO.tab$name_1006),
                      grep("immun", GO.tab$definition_1006),
                      grep("antibod", GO.tab$description),
                      grep("antibod", GO.tab$phenotype_description),
                      grep("antibod", GO.tab$name_1006),
                      grep("antibod", GO.tab$definition_1006))))

GO.tab <- GO.tab[temp,]

head(GO.tab)


GO.tab$consensus_locus <- GO.tab$external_gene_name
GO.tab$consensus_locus <- toupper(GO.tab$consensus_locus)
tail(GO.tab)

GO.tab <- join(GO.tab, genes.in.sig.regions)

GO.tab <- GO.tab[-grep("nuclear speck", GO.tab$name_1006),]
GO.tab <- GO.tab[-grep("nuclear body", GO.tab$name_1006),]


write.table(genes.in.sig.regions, "results/2_GeneList_in_associated_regions.txt", row.names = F, sep = "\t", quote = F)
write.table(GO.tab, "results/2_Immune_GO_terms_in_associated_regions.txt", row.names = F, sep = "\t", quote = F)
