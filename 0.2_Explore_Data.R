#
# Basic summaries and figures
# AMS & SEJ
# Nov 2018
#

library(ggplot2)
library(reshape)
library(plyr)
library(GenABEL)

source("r/multiplot.R")

load("data/20181107_BEAST_data_formatted.RData")
load("data/20181107_Genabel_Data.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNP Data set                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

nsnps(genabeldata)
nids(genabeldata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Distributions of antibody levels              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ format data

BEASTX$LambAdult <- ifelse(BEASTX$Age == 0, "Lambs", "Adults")

Ig.Dist <- subset(BEASTX, select = c(LambAdult, IgAmp, IgEmp, IgGmp))
Ig.Dist <- melt(Ig.Dist, id.vars = "LambAdult")
Ig.Dist <- na.omit(Ig.Dist)

Ig.Dist$variable <- as.character(Ig.Dist$variable)

x <- data.frame(variable = c("IgAmp", "IgGmp", "IgEmp"),
                new.variable = c("Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgE"))

Ig.Dist <- join(Ig.Dist, x)

Ig.Dist$LambAdult <- factor(Ig.Dist$LambAdult, levels = c("Lambs", "Adults"))

ggplot(Ig.Dist, aes(value)) + 
  geom_histogram(col = "grey", binwidth = 0.1) +
  facet_grid(new.variable ~ LambAdult, scales = "free") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14)) +
  labs(x = "Antibody Level", y = "count")

ggsave("figs/0_Distn_of_Antibody_Measures.png", width = 7, height = 9)


#~~ Log Lambs

Ig.Dist <- subset(Ig.Dist, LambAdult == "Lambs")
Ig.Dist$value2 <- log(Ig.Dist$value)

ggplot(Ig.Dist, aes(value2)) + 
  geom_histogram(col = "grey", binwidth = 0.1) +
  facet_grid(new.variable ~ LambAdult, scales = "free") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14)) +
  labs(x = "Antibody Level", y = "count")

ggsave("figs/0_Distn_of_Antibody_Measures_log_Lambs.png", width = 7, height = 9)


rm(x, Ig.Dist)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Correlations between different levels         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

comparisons <- data.frame(x = c("IgAmp", "IgEmp", "IgEmp", "IgAmp", "IgEmp", "IgEmp"),
                          y = c("IgGmp", "IgGmp", "IgAmp", "IgGmp", "IgGmp", "IgAmp"),
                          x.full = c("Anti-Tc IgA", "Anti-Tc IgE", "Anti-Tc IgE", "Anti-Tc IgA", "Anti-Tc IgE", "Anti-Tc IgE"),
                          y.full = c("Anti-Tc IgG", "Anti-Tc IgG", "Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgG", "Anti-Tc IgA"), 
                          z = c("A", "B", "C", "D", "E", "F"),
                          age = c("Lambs", "Lambs", "Lambs", "Adults", "Adults", "Adults"),
                          slope = NA,
                          intercept = NA,
                          R2 = NA,
                          Adj.R2 = NA,
                          P = NA,
                          stringsAsFactors = F)

plot.list <- list()
model.res <- list()

for(i in 1:nrow(comparisons)){
  
  x <- BEASTX[,c("LambAdult", comparisons$x[i], comparisons$y[i])]
  x <- x[which(x$LambAdult == comparisons$age[i]),]
  names(x)[2:3] <- c("x", "y")
  x$Facet <- comparisons$z[i]
  
  model.res[[i]] <- summary(lm(x$y ~ x$x))
  
  plot.list[[i]] <- ggplot(x, aes(x, y)) + 
    geom_point(alpha = 0.2) +
    stat_smooth(method = "lm") +
    theme_bw() +
    facet_wrap(~Facet) +
    theme(axis.text.x  = element_text (size = 12),
          axis.text.y  = element_text (size = 12),
          axis.title.y = element_text (size = 14, angle = 90),
          axis.title.x = element_text (size = 14),
          strip.text.x = element_text (size = 14, hjust = 0.02),
          strip.background = element_blank()) +
    labs(x = comparisons$x.full[i], y = comparisons$y.full[i])
  
  comparisons$slope[i]     <- model.res[[i]]$coefficients[2,1]
  comparisons$intercept[i] <- model.res[[i]]$coefficients[1,1]
  comparisons$Adj.R2[i]    <- model.res[[i]]$adj.r.squared
  comparisons$R2[i]        <- model.res[[i]]$r.squared
  comparisons$P[i]         <- model.res[[i]]$coefficients[2,4]
  
  
  
}

plot.list <- plot.list[c(1, 4, 2, 5, 3, 6)]

multiplot(plotlist = plot.list, cols = 3)


png(filename = "figs/0_Corrs_of_Antibody_Measures.png", width = 9, height = 7, units = "in", res = 300)
multiplot(plotlist = plot.list, cols = 3)
dev.off()

comparisons <- subset(comparisons, select = -c(x, y, z, R2))
names(comparisons)[1:2] <- c("Measure X", "Measure Y")

comparisons$slope     <- formatC(comparisons$slope, digits = 3, format = "f")
comparisons$intercept <- formatC(comparisons$intercept, digits = 3, format = "f")
comparisons$Adj.R2    <- formatC(comparisons$Adj.R2, digits = 3, format = "f")
comparisons$P         <- formatC(comparisons$P, digits = 2, format = "e")

write.table(comparisons, "figs/Table_S1.txt", row.names = F, sep = "\t", quote = F)

#plot(BEASTX$IgEmp, BEASTX$IgGmp)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Association with Age                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- subset(BEASTX, select = c(AgeF, IgAmp, IgEmp, IgGmp, Sex))
x <- subset(x, !is.na(AgeF) & !is.na(Sex))

x <- melt(x, id.vars = c("AgeF", "Sex"))
head(x)

y <- data.frame(variable = c("IgAmp", "IgGmp", "IgEmp"),
                new.variable = c("Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgE"))

x <- join(x,y)
x$Age <- as.numeric(as.character(x$Age))

ggplot(x, aes(AgeF, value, fill = Sex)) + 
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~new.variable, ncol = 1, scales = "free_y") +
  scale_fill_brewer(palette =  "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12)) +
  labs(x = "Age", y = "Antibody Level")
  
ggsave("figs/0_Distn_of_Antibody_Measures_with_Age.png", width = 9, height = 9)

x$LambAdult = "Adults"
p1 <- ggplot(subset(x, Age != 0), aes(Age, value)) + 
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm") +
  facet_grid(new.variable~LambAdult, scales = "free_y") +
  scale_colour_brewer(palette =  "Set1") +
  scale_x_continuous(breaks = seq(0, 16, 1)) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12, hjust = 0.02),
        strip.text.y = element_text (size = 12)) +
  labs(x = "Age (Years)", y = "Antibody Level")



x1 <- subset(BEASTX, select = c(LambAgeAugust, IgAmp, IgEmp, IgGmp, Sex, AgeF))
x1 <- subset(x1, !is.na(LambAgeAugust) & !is.na(Sex) & AgeF == 0)
x1$AgeF <- NULL

x1 <- melt(x1, id.vars = c("LambAgeAugust", "Sex"))
head(x1)

x1 <- join(x1,y)
x1$LambAgeAugust <- as.numeric(as.character(x1$LambAgeAugust))

x1$LambAdult = "Lambs"
p2 <- ggplot(x1, aes(LambAgeAugust, value)) + 
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm") +
  facet_grid(new.variable~LambAdult, scales = "free_y") +
  scale_colour_brewer(palette =  "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12, hjust = 0.02),
        strip.text.y = element_text (size = 12)) +
  labs(x = "Age (Days)", y = "Antibody Level")


multiplot(p2, p1, cols = 2)


png(filename = "figs/0_Antibody_corrs_within_Age.png", width = 9, height = 9, units = "in", res = 300)
multiplot(p2, p1, cols = 2)
dev.off()

rm(x, x1, y, p1, p2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Association with Sex                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- subset(BEASTX, select = c(LambAdult, IgAmp, IgEmp, IgGmp, Sex))
x <- subset(x, !is.na(LambAdult) & !is.na(Sex))

x <- melt(x, id.vars = c("LambAdult", "Sex"))
head(x)

y <- data.frame(variable = c("IgAmp", "IgGmp", "IgEmp"),
                new.variable = c("Anti-Tc IgA", "Anti-Tc IgG", "Anti-Tc IgE"))

x <- join(x,y)

x$LambAdult <- factor(x$LambAdult, levels = c("Lambs", "Adults"))

ggplot(x, aes(Sex, value, fill = Sex)) + 
  geom_boxplot(outlier.alpha = 0.2, notch = T, width = 0.6) +
  facet_grid(new.variable~LambAdult, scales = "free_y") +
  scale_fill_brewer(palette =  "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        legend.position = "none") +
  labs(x = "Sex", y = "Antibody Level")

ggsave("figs/0_Distn_of_Antibody_Measures_with_Age_and_Sex.png", width = 5, height = 7)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Correlation from one year to the next         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(BEASTX)

Adults <- subset(BEASTX, Age != 0)
Adults <- subset(Adults, select = c(ID, CapYear, IgAmp, IgEmp, IgGmp))
head(Adults)

Adults <- melt(Adults, id.vars = c("ID", "CapYear"))

AdultsPlus1 <- Adults
AdultsPlus1$CapYear <- AdultsPlus1$CapYear - 1
names(AdultsPlus1)[4] <- "value.t1"

Adults <- join(Adults, AdultsPlus1)
rm(AdultsPlus1)
head(Adults)

Adults <- join(Adults, y)

Adults$Difference <- Adults$value.t1 - Adults$value

p1 <- ggplot(Adults, aes(value, value.t1)) + 
  geom_point(alpha = 0.4) + 
  stat_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~new.variable, scales = "free") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12)) +
  labs(x = "Antibody Level at Time t", y = "Antibody Level at Time t+1")

p2 <- ggplot(Adults, aes(Difference)) + 
  geom_histogram(binwidth = 0.1, col = "darkgrey") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~new.variable, scales = "free") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12)) +
  labs(x = "Change in Antibody Level from t to t+1")


png(filename = "figs/0_Corrs_of_Antibody_Measures_t_t1.png", width = 9, height = 7, units = "in", res = 300)
multiplot(p1, p2)
dev.off()

multiplot(p1, p2)


#~~ Correlations


comparisons <- data.frame(x = c("IgAmp", "IgEmp", "IgGmp"),
                          x.full = c("Anti-Tc IgA", "Anti-Tc IgE", "Anti-Tc IgG"),
                          slope = NA,
                          intercept = NA,
                          Adj.R2 = NA,
                          P = NA,
                          stringsAsFactors = F)

for(i in 1:nrow(comparisons)){
  
  x <- subset(Adults, variable == comparisons$x[i])
  
  y <- summary(lm(x$value.t1 ~ x$value))
  
  comparisons$slope[i]     <- y$coefficients[2,1]
  comparisons$intercept[i] <- y$coefficients[1,1]
  comparisons$Adj.R2[i]    <- y$adj.r.squared
  comparisons$P[i]         <- y$coefficients[2,4]
  
  
  
}


comparisons$slope     <- formatC(comparisons$slope, digits = 3, format = "f")
comparisons$intercept <- formatC(comparisons$intercept, digits = 3, format = "f")
comparisons$Adj.R2    <- formatC(comparisons$Adj.R2, digits = 3, format = "f")
comparisons$P         <- formatC(comparisons$P, digits = 2, format = "e")


