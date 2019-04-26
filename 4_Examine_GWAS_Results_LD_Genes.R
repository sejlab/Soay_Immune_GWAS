#
# Parse the GWAS Results
# Nov 2018
# SEJ
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up workspace                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Load libraries and data

library(GenABEL)
library(ggplot2)
library(plyr)
library(reshape)
library(tidyr)
library(dplyr)
library(LDheatmap)


source("r/multiplot.R")
colfunc <- colorRampPalette(c("red", "white"))
source("r/biomaRt_func.R")


#~~ Load results and genomic information

load("results/2_GWAS_Full_Results.RData",  verbose = T)
load("results/3_HD_GWAS_Full_Results.RData", verbose = T)

snp.wald.imputed <- read.table("results/3_Imputed_GWAS_Results.txt", header = T, stringsAsFactors = F)

#~~ Load Gene information

gene.list <- read.table("data/Ovis_aries.Oar_v3.1.94.genes.txt", header = T, stringsAsFactors = F)
immune.GO.in.sig.regions <- read.table("results/2_Immune_GO_terms_in_associated_regions.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "", quote = "")
genes.in.sig.regions <- read.table("results/2_GeneList_in_associated_regions.txt", header = T, stringsAsFactors = F, sep = "\t")

#~~ Get sample size info 

samples <- read.table("results/1_GRM_Sample_Sizes_VP.txt", header = T, stringsAsFactors = F, sep = "\t")
samples <- subset(samples, select = c(Response, LambAdult, N, N_ids))
names(samples)[1] <- "Trait"

y <- data.frame(Trait = c("IgAmp", "IgGmp", "IgEmp"),
                Trait2 = c("Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgE"))

samples <- join(samples, y)

#~~ remove stuff you don't need.

rm(full.n, full.ranef, full.wald, full.n.hd, full.ranef.hd, full.wald.hd)

#~~ Make GenABEL object for LD measures

load("data/20181107_Genabel_Data.RData")

HDabel <- load.gwaa.data(phenofile = "data/20140214_SheepHD_QC1_GenABELpheno.txt", 
                         genofile = "data/20140214_SheepHD_QC1_GenABEL.txt")

HDabel <- HDabel[,!snpnames(HDabel) %in% snp.names(genabeldata)]

HDabel <- merge.gwaa.data(HDabel, genabeldata)

rm(genabeldata)

#~~ IgE Logged - add this data in for parsing

snp.ige <- read.table("results/2_GWAS_Results_log_Lambs.txt", header = T, sep = "\t", stringsAsFactors = F)
snp.ige <- subset(snp.ige, Model == "IgEmp_Lambs")

snp.ige.imputed <- read.table("results/3_Imputed_GWAS_Results_log_Lambs.txt", header = T, sep = "\t", stringsAsFactors = F)
snp.ige.imputed <- subset(snp.ige.imputed, Model == "IgEmp_Lambs")


snp.wald.corrected <- subset(snp.wald.corrected, Model != "IgEmp_Lambs")
snp.wald.corrected <- rbind(snp.wald.corrected, snp.ige)

snp.wald.imputed <- subset(snp.wald.imputed, Model != "IgEmp_Lambs")
snp.wald.imputed <- rbind(snp.wald.imputed, snp.ige.imputed)

# snp.ige <- read.table("results/2_GWAS_Results_log_Lambs.txt", header = T, sep = "\t", stringsAsFactors = F)
# snp.ige <- subset(snp.ige, Model %in% c("IgEmp_Lambs", "IgGmp_Lambs"))
# 
# snp.wald.corrected <- subset(snp.wald.corrected, !Model %in% c("IgEmp_Lambs", "IgGmp_Lambs"))
# snp.wald.corrected <- rbind(snp.wald.corrected, snp.ige)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create GWAS Plots for Paper               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

snp.wald.corrected <- join(snp.wald.corrected, y)
head(snp.wald.corrected)

#~~ Make X axis positions

x <- subset(snp.wald.corrected, Model == unique(snp.wald.corrected$Model)[1])
x <- arrange(x, Chromosome, Position)

x$Diff <- c(0, diff(x$Position))
x$Diff <- ifelse(x$Diff < 0, 1, x$Diff)
x$Cumu <- cumsum(x$Diff)

chrinfo <- NULL
for(j in unique(x$Chromosome)){
  temp1 <- subset(x, Chromosome == j)
  temp2 <- data.frame(Chromosome = j,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  chrinfo <- rbind(chrinfo, temp2)
}
chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$Chromosome[28] <- "X"

chrinfo$Chromosome2 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "", 
                         "12", "", "14", "", "16", "", "18", "", "20", "", "", 
                         "", "24", "", "", "X")

rm(x, temp1, temp2, j)

samples$Ntext <- paste0(samples$N, " (", samples$N_ids, ")")
samples

x <- melt(tapply(-log10(snp.wald.corrected$Wald.P.corrected), list(snp.wald.corrected$Trait2), max))
names(x)[1] <- c("Trait2")

samples <- join(samples, x)
samples$Ntext[which(samples$LambAdult == "Lambs")] <- samples$N[which(samples$LambAdult == "Lambs")]

samples$LambAdult <- factor(samples$LambAdult, levels = c("Lambs", "Adults"))

rm(x)

#~~ plot lambs, adults

ggplot() + 
  geom_hline(yintercept = -log10(2.245e-6), linetype = "dashed") +
  geom_point(data = snp.wald.corrected, aes(Cumu, -log10(Wald.P.corrected), col = factor(Chromosome %% 2)), alpha = 0.7) +
  scale_colour_brewer(palette = "Set1") +
  geom_point(data = samples, aes(x=1, y = value*(1.15)), col = "white", alpha = 0) +
  geom_text(data = samples,  aes(x=1, y = value*(1.09), label = Ntext), vjust=0, hjust=0, size = 3.5) +
  facet_grid(Trait2 ~ LambAdult, scales = "free_y") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chromosome2) +
  labs(x = "Chromosome", y = "-log10 P") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 8, hjust = 0.35),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.position = "none")

ggsave("figs/4_GWAS_Plot_Lambs_Adults.png", width = 10, height = 10)

ggplot(data = snp.wald.corrected, aes(-log10(Exp.P), -log10(Wald.P.corrected))) + 
  geom_hline(yintercept = -log10(2.245e-6), linetype = "dashed") +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(Trait2 ~ LambAdult, scales = "free_y") +
  theme(legend.position = "none") +
  labs(x = "Expected -log10 P", y = "Observed -log10 P") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12, hjust = 0.35),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.position = "none")

ggsave("figs/4_PP_Plot_Lambs_Adults.png", width = 8, height = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Extract the Top SNPs and candidate information #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

topsnps <- subset(snp.wald.corrected, Pr.Chisq. < 2.245e-6)
topsnps <- subset(topsnps, select = c(Trait, LambAdult, Df, SNP.Name, Chromosome, Position, A1, A2, NoMeasured, CallRate, Q.2, Model, Pr.Chisq., EffectAA, EffectAB, EffectBB))
topsnps$Chip <- "SNP50"

#~~ Get the HD regions

topsnpsHD <- subset(snp.wald.imputed, Pr.Chisq. < 2.245e-6)
topsnpsHD <- subset(topsnpsHD, select = c(Trait, LambAdult, Df, SNP.Name, Chromosome, Position, A1, A2, NoMeasured, ImputeSuccess, Q.2, Model, Pr.Chisq., EffectAA, EffectAB, EffectBB))
names(topsnpsHD)[which(names(topsnpsHD) == "ImputeSuccess")] <- "CallRate"
topsnpsHD$Chip <- "SNPHD"

topsnps <- rbind(topsnps, topsnpsHD)

table(topsnps$Trait, topsnps$Chromosome)

topsnps <- join(topsnps, y)

x1 <- data.frame(table(topsnps$SNP.Name, topsnps$Model))
names(x1)[1:2] <- c("SNP.Name", "Model")

topsnps <- join(topsnps, x1)
topsnps <- arrange(topsnps, SNP.Name)
head(topsnps)

topsnps <- topsnps[-which(topsnps$Freq == 2 & topsnps$Chip == "SNP50"),]

rm(topsnpsHD)


#~~ Get the most adjacent genes

topsnps$GeneLeftID <- NA
topsnps$GeneLeftName <- NA
topsnps$GeneLeftStart <- NA
topsnps$GeneLeftStop <- NA
topsnps$GeneLeftStrand <- NA

topsnps$GeneRightID <- NA
topsnps$GeneRightName <- NA
topsnps$GeneRightStart <- NA
topsnps$GeneRightStop <- NA
topsnps$GeneRightStrand <- NA

topsnps$InGene <- NA
topsnps$Candidates <- NA

for(i in 1:nrow(topsnps)){
  
  if(topsnps$Chromosome[i] != "0"){
    genetemp <- subset(genes.in.sig.regions, Chromosome == topsnps$Chromosome[i] & Start < topsnps$Position[i])
    genetemp <- genetemp[nrow(genetemp),]
    
    if(nrow(genetemp) > 0){
      
      topsnps$GeneLeftName[i] <- genetemp$consensus_locus
      topsnps$GeneLeftID[i] <- genetemp$gene_id
      topsnps$GeneLeftStart[i] <- genetemp$Start
      topsnps$GeneLeftStop[i] <- genetemp$Stop
      topsnps$GeneLeftStrand[i] <- genetemp$Strand
      
      
      
      if(topsnps$Position[i] < genetemp$Stop){
        topsnps$InGene[i] <- "Left"
      } else {
        genetemp <- subset(genes.in.sig.regions, Chromosome == topsnps$Chromosome[i] & Start > topsnps$Position[i])
        genetemp <- genetemp[1,]
        
        topsnps$GeneRightName[i] <- genetemp$consensus_locus
        topsnps$GeneRightID[i] <- genetemp$gene_id
        topsnps$GeneRightStart[i] <- genetemp$Start
        topsnps$GeneRightStop[i] <- genetemp$Stop
        topsnps$GeneRightStrand[i] <- genetemp$Strand
        
      }
      
    }
    cands <- subset(immune.GO.in.sig.regions, Chromosome == topsnps$Chromosome[i])
    topsnps$Candidates[i] <- paste0(unique(cands$consensus_locus), collapse = ", ")
    
    rm(cands, genetemp)
  }
}

topsnps <- subset(topsnps, CallRate >= 0.95)

write.table(topsnps, "results/4_Information_on_Top_SNPs.txt", row.names = F, sep = "\t", quote = F)

rm(chrinfo, i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Examine the associated regions            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get the SNPs at the centre of each region and how much either side you want to plot

head(topsnps)
topsnps$Chromosome2 <- topsnps$Chromosome

focalsnps <- melt(tapply(topsnps$Pr.Chisq., list(topsnps$Model, topsnps$Chromosome2), min))
names(focalsnps) <- c("Model", "Chromosome2", "Pr.Chisq.")
focalsnps <- na.omit(focalsnps)

focalsnps <- join(focalsnps, subset(topsnps, select = c(Model, Chromosome2, Pr.Chisq., SNP.Name, Position)))

x1 <- na.omit(melt(tapply(topsnps$Position, list(topsnps$Model, topsnps$Chromosome2), min)))
x2 <- na.omit(melt(tapply(topsnps$Position, list(topsnps$Model, topsnps$Chromosome2), max)))

names(x1) <- c("Model", "Chromosome2", "Min.Window")
names(x2) <- c("Model", "Chromosome2", "Max.Window")

focalsnps <- join(focalsnps, x1)
focalsnps <- join(focalsnps, x2)

rm(x1, x2)

table(focalsnps$Model, focalsnps$Chromosome2)

x1 <- na.omit(melt(tapply(focalsnps$Min.Window, list(focalsnps$Model, focalsnps$Chromosome2), min)))
names(x1) <- c("Model", "Chromosome2", "Min.Window.all")
x2 <- na.omit(melt(tapply(focalsnps$Max.Window, list(focalsnps$Model, focalsnps$Chromosome2), max)))
names(x2) <- c("Model", "Chromosome2", "Max.Window.all")

focalsnps <- join(focalsnps, x1)
focalsnps <- join(focalsnps, x2)

rm(x1, x2)

focalsnps$Order <- 1
for(i in 2:nrow(focalsnps)) focalsnps$Order[i] <- ifelse(focalsnps$Min.Window.all[i] == focalsnps$Min.Window.all[i-1],
                                                         focalsnps$Order[i-1] + 1, 1)

focalsnps <- subset(focalsnps, Order == 1)

focalsnps <- subset(focalsnps, select = c(Model, Chromosome2, Pr.Chisq., SNP.Name, Position, Min.Window.all, Max.Window.all))
focalsnps$Chromosome <- focalsnps$Chromosome2
focalsnps$Chromosome[which(focalsnps$Chromosome %in% c("2a", "2b"))] <- 2
focalsnps$Chromosome <- as.numeric(as.character(focalsnps$Chromosome))

focalsnps$Min.Window.all <- focalsnps$Min.Window.all - 1e6
focalsnps$Max.Window.all <- focalsnps$Max.Window.all + 1e6

focalsnps <- separate(focalsnps, Model, c("Trait", "LambAdult"), sep = "_", remove = T)
focalsnps <- subset(focalsnps, Chromosome != 0)

#~~ make the IgA Lambs match the IgA Adults

focalsnps$Min.Window.all[which(focalsnps$Trait == "IgAmp" & 
                                 focalsnps$LambAdult == "Lambs" &
                                 focalsnps$Chromosome == 24)] <- 
  focalsnps$Min.Window.all[which(focalsnps$Trait == "IgAmp" & 
                                   focalsnps$LambAdult == "Adults" &
                                   focalsnps$Chromosome == 24)]

focalsnps$Max.Window.all[which(focalsnps$Trait == "IgAmp" & 
                                 focalsnps$LambAdult == "Lambs" &
                                 focalsnps$Chromosome == 24)] <- 
  focalsnps$Max.Window.all[which(focalsnps$Trait == "IgAmp" & 
                                   focalsnps$LambAdult == "Adults" &
                                   focalsnps$Chromosome == 24)]

focalsnps

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Plot the associated regions               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

snp.wald.imputed$SNP50 <- ifelse(snp.wald.imputed$SNP.Name %in% snp.wald.corrected$SNP.Name, "50K", "Imputed")
snp.wald.imputed$SNP50 <- factor(snp.wald.imputed$SNP50, levels = c("Imputed", "50K"))

LD.info <- NULL

LD.with.focalsnps <- NULL

for(i in 1:nrow(focalsnps)){
  
  #~~ Get gene and association information for the region
  
  regtab <- subset(snp.wald.imputed, Trait == focalsnps$Trait[i] & 
                     LambAdult == focalsnps$LambAdult[i] &
                     Chromosome == focalsnps$Chromosome[i] & 
                     Position > focalsnps$Min.Window.all[i] & 
                     Position < focalsnps$Max.Window.all[i])
  
  regtab <- arrange(regtab, Position)
  
  regtab <- subset(regtab, ImputeSuccess > 0.95)
  
  reggenes <- subset(gene.list, Chromosome == focalsnps$Chromosome[i] & 
                       Stop > focalsnps$Min.Window.all[i] & 
                       Start < focalsnps$Max.Window.all[i])
  
  reggenes <- unique(subset(reggenes, select = c(gene_id, Start, Stop)))
  reggenes$gene_name2 <- 1:nrow(reggenes)
  
  #~~ Do some formatting of genes
  
  temp.p <- -log10(focalsnps$Pr.Chisq.[i])
  temp.p.divider <- ((temp.p *(109/100))-temp.p)/3
  temp.p <- temp.p + (4*temp.p.divider)
  
  reggenes$positions <- rep(seq(temp.p + temp.p.divider, temp.p + (6*temp.p.divider), temp.p.divider), length = nrow(reggenes))
  reggenes <- melt(reggenes, id.vars = c("gene_id", "gene_name2", "positions"))
  
  #~~ Make and save the plots
  
  p1 <- ggplot(regtab, aes(Position/1e6, -log10(Pr.Chisq.))) +
    geom_hline(yintercept = -log10(2.245e-6), linetype = "dashed") +
    geom_point(aes(colour = SNP50, shape = SNP50), alpha = 0.7, size = 1.5) +
    geom_rect(xmin = (focalsnps$Min.Window.all[i]/1e6)-1,
              xmax = (focalsnps$Max.Window.all[i]/1e6)+1,
              ymin = temp.p-temp.p.divider,
              ymax = temp.p + (8*temp.p.divider), 
              fill = "grey95") +
    geom_hline(yintercept = temp.p-temp.p.divider) +
    geom_line(data = reggenes, aes(value/1e6, positions,
                                   group = gene_name2,
                                   colour = gene_id %in% immune.GO.in.sig.regions$gene_id), size = 2) +
    #geom_point(data = subset(reggenes, gene_id %in% immune.GO.in.sig.regions$gene_id), aes(value/1e6, positions), col = "green", alpha = 0.2) +
    scale_colour_manual(values = c("red", "black", "black", "red")) +
    theme_bw() +
    labs(x = paste0("Chromosome ", focalsnps$Chromosome[i]," Position (MB)"), y = "-log10(P)", col = "", shape = "") +
    theme(axis.text.x  = element_text (size = 12),
          axis.text.y  = element_text (size = 12),
          strip.text.x = element_text (size = 12),
          strip.text.y = element_text (size = 12),
          axis.title.y = element_text (size = 14, angle = 90),
          axis.title.x = element_text (size = 14))
  
  p1
  ggsave(paste0("figs/4_Assoc_Region_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], ".png"), width = 8, height = 4)
  
  pdf(paste0("figs/Assoc_Region_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], ".pdf"), width = 8, height = 4)
  p1
  dev.off()
  
  
  #~~ Effect sizes at the top SNP
  
  # if(focalsnps$SNP.Name[i] %in% full.fixef$SNP.Name){
  #   x <- subset(full.fixef, SNP.Name == focalsnps$SNP.Name[i] & 
  #                 Trait == focalsnps$Trait[i] &
  #                 LambAdult == focalsnps$LambAdult[i])
  # } else {
  #   x <- subset(full.fixef.hd, SNP.Name == focalsnps$SNP.Name[i] & 
  #                 Trait == focalsnps$Trait[i] &
  #                 LambAdult == focalsnps$LambAdult[i])
  # }
  # 
  # x <- join(x, y)
  # x <- x[grep("SNP", x$variable),]
  # x$variable <- gsub("SNP_", "", x$variable)
  # 
  # ggplot(x, aes(variable, solution)) +
  #   geom_hline(yintercept = 0) +
  #   geom_point(size = 2) +
  #   geom_errorbar(aes(ymin = solution - std.error, ymax = solution + std.error), width = 0) +
  #   theme_bw() +
  #   theme(axis.text.x  = element_text (size = 12),
  #         axis.text.y  = element_text (size = 12),
  #         axis.title.y = element_text (size = 12, angle = 90),
  #         axis.title.x = element_text (size = 12),
  #         strip.text.x = element_text (size = 12, hjust = 0.02)) +
  #   labs(x ="Genotype", y = "Effect Size") +
  #   ggtitle(paste(focalsnps$Trait[i], focalsnps$LambAdult[i], focalsnps$SNP.Name[i], "MAF = ", formatC(subset(snp.wald.imputed, SNP.Name == x$SNP.Name[1])$Q.2[1], digits = 3)))
  # 
  # ggsave(paste0("figs/4_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], "_Effect_Size.png"), width = 6, height = 6)
  # 
  # rm(x)
  
  #~~ Plot LD of these regions
  
  regtab <- subset(regtab, SNP.Name %in% snpnames(HDabel))
  x <- HDabel[,regtab$SNP.Name]
  x1 <- r2fast(x)

  x2 <- LDheatmap(x1, genetic.distances = map(x), distances = "physical", color = colfunc(20), flip = T, title = paste(focalsnps$Trait[i], focalsnps$LambAdult[i], "Chromosome", focalsnps$Chromosome2))

  x3 <- melt(x1)
  x3 <- na.omit(subset(x3, value < 1.01))
  

  tempmap <- subset(regtab, select = c(SNP.Name, Position))
  names(tempmap) <- c("X1", "X1pos")
  x3 <- join(x3, tempmap)
  tempmap <- subset(regtab, select = c(SNP.Name, Position))
  names(tempmap) <- c("X2", "X2pos")
  x3 <- join(x3, tempmap)
  
  x3$Distance <- x3$X2pos-x3$X1pos
  
  x3 <- subset(x3, Distance > 0)
  
  x4 <- x3[c(which(x3$X1 == focalsnps$SNP.Name[i]), which(x3$X2 == focalsnps$SNP.Name[i])),]

  
  ggplot(x3, aes(Distance, value)) +
    geom_point(alpha = 0.1) +
    stat_smooth() +
    theme_bw() +
    labs(x = paste0("Distance (bp)"), y = "LD") +
    theme(axis.text.x  = element_text (size = 12),
          axis.text.y  = element_text (size = 12),
          strip.text.x = element_text (size = 12),
          strip.text.y = element_text (size = 12),
          axis.title.y = element_text (size = 14, angle = 90),
          axis.title.x = element_text (size = 14)) +
    ggtitle(paste(focalsnps$Trait[i], focalsnps$LambAdult[i], "Chromosome", focalsnps$Chromosome[i]))
  ggsave(paste0("figs/4_LD_Decay_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], ".png"), width = 6, height = 6)
  
  
  LD.info <- rbind(LD.info,
                   data.frame(Nsnps = nrow(regtab),
                              MeanLD = mean(x3$value),
                              RegionSize = max(regtab$Position) - min(regtab$Position),
                              MeanLDoverDistance = mean(x3$value/x3$Distance),
                              SNP.Name = focalsnps$SNP.Name[i],
                              Trait = focalsnps$Trait[i],
                              LambAdult = focalsnps$LambAdult[i],
                              Chromosome = focalsnps$Chromosome))
  
  pdf(paste0("figs/4_LD_Heatmap_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], ".pdf"), width = 10, height = 10)
  LDheatmap(x2)
  dev.off()

  png(paste0("figs/4_LD_Heatmap_", focalsnps$Trait[i], "_", focalsnps$LambAdult[i], "_Chr_", focalsnps$Chromosome2[i], ".png"), width = 10, height = 10, units = "in", res = 300)
  LDheatmap(x2)
  dev.off()
  
  LD.with.focalsnps <- rbind(LD.with.focalsnps, x4)
  
  rm(p1, temp.p, temp.p.divider, reggenes, regtab, x, x1, x2, x3, x4)
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Look at genes that you are interested in         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

genedesc("SNX29")

x <- join(focalsnps, topsnps)

x1 <- data.frame(table(x$SNP.Name))
names(x1)[1] <- "SNP.Name"

# x <- join(x, x1)
# x$Freq <- ifelse(x$Freq == 2 & x$Chip == "SNPHD", 2, 1)
# x <- subset(x, Freq == 1)

rm(x1)

write.table(x, "results/4_Top_SNPs_for_Paper.txt", row.names = F, sep = "\t", quote = F)

write.table(LD.with.focalsnps, "results/4_LD_with_Top_SNPs.txt", row.names = F, sep = "\t", quote = F)


gc()



