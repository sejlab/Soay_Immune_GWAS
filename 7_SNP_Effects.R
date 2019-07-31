
library(asreml)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)

focalsnps <- read.table("results/4_Top_SNPs_for_Paper.txt", header = T, sep = "\t", stringsAsFactors = F)
focalsnps <- focalsnps[,c("Trait", "LambAdult", "SNP.Name")]

load("results/2_GWAS_Full_Results.RData")

focalsnps1 <- join(focalsnps, full.fixef)
focalsnps2 <- join(focalsnps, full.wald)
focalsnps3 <- join(focalsnps, full.n)


load("results/3_HD_GWAS_Full_Results.RData")

full.fixef.hd <- subset(full.fixef.hd, !SNP.Name %in% unique(subset(focalsnps1, !is.na(solution))$SNP.Name))
full.wald.hd <- subset(full.wald.hd, !SNP.Name %in% unique(subset(focalsnps1, !is.na(solution))$SNP.Name))
full.n.hd <- subset(full.n.hd, !SNP.Name %in% unique(subset(focalsnps1, !is.na(solution))$SNP.Name))

focalsnps4 <- join(focalsnps, full.fixef.hd)
focalsnps5 <- join(focalsnps, full.wald.hd)
focalsnps6 <- join(focalsnps, full.n.hd)

load("results/3_HD_GWAS_Full_Results_log_Lambs.RData")

focalsnps7 <- join(focalsnps, full.fixef.hd)
focalsnps8 <- join(focalsnps, full.wald.hd)
focalsnps9 <- join(focalsnps, full.n.hd)

rm(list = grep("full", ls(), value = T))
rm(snp.wald.corrected, sig.snps)

names(focalsnps1) == names(focalsnps4)
names(focalsnps1) == names(focalsnps7)

focalsnps1 <- rbind(focalsnps1, focalsnps4, focalsnps7)
focalsnps1 <- subset(focalsnps1, !is.na(variable))

focalsnps1 <- focalsnps1[grep("SNP", focalsnps1$variable),]
focalsnps1$variable <- gsub("SNP_", "", focalsnps1$variable)


focalsnps2 <- rbind(focalsnps2, focalsnps5, focalsnps8)
focalsnps2 <- subset(focalsnps2, variable == "SNP")
focalsnps2 <- subset(focalsnps2, select = -variable)

focalsnps1 <- join(focalsnps1, focalsnps2)

focalsnps3 <- rbind(focalsnps3, focalsnps6, focalsnps9)
focalsnps3 <- subset(focalsnps3, !is.na(Var1))

focalsnps1 <- join(focalsnps1, focalsnps3)

rm(list = paste0("focalsnps", 2:9))

#~~ Make plots

focalsnps1
str(focalsnps1)

focalsnps1 <- subset(focalsnps1, variable == Var1) %>% unique
focalsnps1$std.error[which(is.na(focalsnps1$std.error))] <- 0

x <- focalsnps1 %>%
  group_by(SNP.Name) %>%
  summarize(ymax = max(solution + std.error, na.rm = T),
            ymin = min(solution - std.error, na.rm = T))

x$Range <- x$ymax - x$ymin

focalsnps <- join(focalsnps1, x)


focalsnps$SNP.Genotype.num <- as.character(focalsnps$Var1)
focalsnps$SNP.Genotype.num[which(focalsnps$SNP.Genotype.num == "A/A")] <- 1
focalsnps$SNP.Genotype.num[which(focalsnps$SNP.Genotype.num == "A/G")] <- 2
focalsnps$SNP.Genotype.num[which(focalsnps$SNP.Genotype.num == "G/A")] <- 2
focalsnps$SNP.Genotype.num[which(focalsnps$SNP.Genotype.num == "G/G")] <- 3

focalsnps$SNP.Genotype.num <- as.numeric(focalsnps$SNP.Genotype.num)

focalsnps$Model <- paste0(gsub("s", "", focalsnps$LambAdult), " ", 
                          gsub("mp", "", focalsnps$Trait), "\n", 
                          focalsnps$SNP.Name)

ggplot(focalsnps) +
  geom_hline(yintercept = 0) +
  geom_point(aes(SNP.Genotype.num, solution), size = 2) +
  geom_line(aes(SNP.Genotype.num, solution)) +
  geom_errorbar(aes(x = SNP.Genotype.num, ymin = solution - std.error, ymax = solution + std.error), width = 0) +
  geom_text(aes(SNP.Genotype.num, ymax + 0.3*Range, label = Freq)) +
  geom_point(aes(SNP.Genotype.num,ymax + 0.4*Range), alpha = 0) +
  theme_bw() +
  scale_x_continuous(breaks = 1:3, labels = c("A/A", "A/G", "G/G")) +
  facet_wrap(~Model, scales = "free", ncol = 2) +
  labs(x = "SNP Genotype", y = "Effect Size (Relative to Model Intercept)") +
  coord_cartesian(xlim = c(0.5, 3.5)) +
  theme(legend.position = "top") 

ggsave("figs/7_SNP_Effects.png", width = 6, height = 10)

write.table(focalsnps[,-20], "results/7_Allele_Effects.txt", row.names = F, sep = "\t", quote = F)
