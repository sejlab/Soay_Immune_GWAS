#
# Parse genetic correlations
# SEJ, AMS
# July 2019
#

library(ggplot2)
library(plyr)

load("gen_corr/1.3_Gen_Corr_Data.RData")

x1

restab <- NULL

for(i in 1:nrow(x1)){
  
  xped <- read.table(paste0("gen_corr/Ped_Res_", i, ".txt"), header = T, sep = "\t")
  xgrm <- read.table(paste0("gen_corr/GRM_Res_", i, ".txt"), header = T, sep = "\t")
  
  xped$Matrix <- "PED"
  xgrm$Matrix <- "GRM"
  
  xped$ModelNo <- i
  xgrm$ModelNo <- i
  
  restab <- rbind(restab, xped, xgrm)
  
}

head(restab)

cortab <- restab[grep("cor", restab$Effect),]

cortab

cortab <- join(cortab, x1)
cortab$Model <- NULL
cortab$Effect <- NULL
cortab$constraint <- NULL
cortab$gamma     <- NULL

ggplot(cortab, aes(X2, X1, fill = component)) +
  geom_tile() +
  geom_text(aes(label = signif(component, digits = 3))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  facet_wrap(~Matrix)

write.table(cortab, "results/1_GRM_Correlations.txt", row.names = F, sep = "\t", quote = F)

cortab <- subset(cortab, Matrix == "GRM")

cortab$label <- paste0(signif(cortab$component, digits = 3), "\n(", signif(cortab$std.error, digits = 2), ")")
cortab$label <- paste0(cortab$label, ifelse(cortab$z.ratio > 3.291, "***",
                                            ifelse(cortab$z.ratio > 2.576, "**",
                                                   ifelse(cortab$z.ratio > 1.96, "*", ""))))

cortab$Model1 <- paste0(cortab$LambAdult.1, "\n", cortab$Trait.1)
cortab$Model2 <- paste0(cortab$LambAdult.2, "\n", cortab$Trait.2, " ")

cortab$Model1 <- gsub("mp", "", cortab$Model1)
cortab$Model2 <- gsub("mp", "", cortab$Model2)

cortab$Model2 <- factor(cortab$Model2, levels = rev(unique(sort(cortab$Model2))))



ggplot(cortab, aes(Model1, Model2, fill = component)) +
  geom_tile(colour = "black", alpha = 0.5) +
  geom_text(aes(label = label)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = "", y = "", fill = "Correlation") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, hjust = 1),
        axis.ticks = element_blank())
ggsave("figs/1_Genetic_Correlations.pdf", width = 7.54, height = 5.90, units = "in", dpi = 300)

subset(cortab, LambAdult.1 != LambAdult.2 & Trait.1 != Trait.2)
