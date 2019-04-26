#
# Format animal model results
# AMS & SEJ
# Nov 2018
#
#

load("data/20181107_BEAST_data_formatted.RData")

library(asreml)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape)

source("r/ASReml.EstEffects.R")
source("r/pin.R")


#~~ Load fixed effects results tables and add meaningful definitions

grm.results.fixef <- read.table("results/1_GRM_Fixef_Results.txt", header = T, stringsAsFactors = F)

grm.results.fixef$new.variable <- NULL

grm.results.ranef <- read.table("results/1_GRM_Ranef_Results.txt", header = T, stringsAsFactors = F, sep = "\t")


recodetab <- data.frame(Trait = c("IgAmp", "IgEmp", "IgGmp"),
                        Trait2 = c("Anti-Tc IgA", "Anti-Tc IgE", "Anti-Tc IgG"), stringsAsFactors = F)


grm.results.fixef <- left_join(grm.results.fixef, recodetab)
grm.results.ranef <- left_join(grm.results.ranef, recodetab)


models <- read.table("results/1_GRM_Sample_Sizes_VP.txt", header = T, stringsAsFactors = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Load the GRM models to get VP, sample sizes, wald tests from  #
#    the GRM model outputs.                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Make table for the model information

models$Trait <- models$Response
models$ID <- paste0(models$Trait, "_", models$LambAdult)


#~~ Make something to save the random effects and wald tests

grm.results.wald  <- NULL


#~~ Load models and extract information

for(i in 1:nrow(models)){
  
  load(paste0("animal_models/", models$ID[i], "_vp.RData"), verbose = T) # loads fit1
  
  grm.results.wald <- rbind(grm.results.wald,
                             cbind(Trait = models$Trait[i], LambAdult = models$LambAdult[i], wald.asreml(fit1), variable = row.names(wald.asreml(fit1))))
  
  rm(fit1)
}

#~~ Tidy up the model information

models <- join(models, recodetab)

models$V_raw   <- formatC(models$V_raw, digits = 4, format = "f")
models$VP      <- formatC(models$VP, digits = 4, format = "f")
models$VP_SE   <- formatC(models$VP_SE, digits = 4, format = "f")
models$mean   <- formatC(models$mean, digits = 4, format = "f")


#~~ Make the definitions clearer for the random effects table

grm.results.wald$Trait <- as.character(grm.results.wald$Trait)
grm.results.wald <- join(grm.results.wald, recodetab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Deal with the fixed effects                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grm.results.fixef <- grm.results.fixef[,c("Trait2", "LambAdult", "variable", "solution", "std.error", "z.ratio")]
grm.results.wald  <- grm.results.wald [,c("Trait2", "LambAdult", "variable", "Df", "Wald statistic", "Pr(Chisq)")]


#~~ Format the fixed effect results

grm.results.fixef <- separate(grm.results.fixef, variable, c("variable", "level"), remove = T)
grm.results.fixef$variable[grm.results.fixef$variable == ""] <- "(Intercept)"
grm.results.fixef$level[is.na(grm.results.fixef$level)] <- ""
grm.results.fixef$level[grm.results.fixef$level == "Intercept"] <- ""
grm.results.fixef$LambAdult <- factor(grm.results.fixef$LambAdult, levels = c("Lambs", "Adults"))
#grm.results.fixef$level[which(grm.results.fixef$level %in% 0:9)] <- paste0("0", grm.results.fixef$level[which(grm.results.fixef$level %in% 0:9)])


grm.results.fixef <- arrange(grm.results.fixef, Trait2, LambAdult, variable, level)

# grm.results.fixef$solution  <- formatC(grm.results.fixef$solution, digits = 4, format = "f")
# grm.results.fixef$std.error <- formatC(grm.results.fixef$std.error, digits = 4, format = "f")
# grm.results.fixef$z.ratio   <- formatC(grm.results.fixef$z.ratio, digits = 4, format = "f")

grm.results.fixef <- join(grm.results.fixef, subset(models, select = c(Trait2, LambAdult, N, N_ids)))

names(grm.results.fixef) <- c("Trait", "Age", "Effect", "Effect.Level", "Solution", "Standard.Error", "Z.Ratio", "N_measures", "N_uniqueIDs")
write.table(grm.results.fixef, "results/1_GRM_Fixef_Results_FORMATTED.txt", row.names = F, sep = "\t", quote = F)

grm.results.wald <- subset(grm.results.wald, !variable %in% c("residual (MS)"))
grm.results.wald <- arrange(grm.results.wald, Trait2)

write.table(grm.results.wald, "results/1_GRM_Fixef_Wald.txt", row.names = F, sep = "\t", quote = F)

grm.results.wald$`Wald statistic`  <- formatC(grm.results.wald$`Wald statistic`, digits = 2, format = "f")
grm.results.wald$`Pr(Chisq)`  <- formatC(grm.results.wald$`Pr(Chisq)` , digits = 2, format = "e")

write.table(grm.results.wald, "results/1_GRM_Fixef_Wald_FORMATTED.txt", row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Deal with the random effects                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ reorder everything random effects

grm.results.ranef <- grm.results.ranef[,c("Trait2", "LambAdult", "new.variable", "gamma", "component", "std.error", 
                                          "z.ratio", "constraint", "Effect", "SE")]



models.hold <- models

models <- models.hold
models$mean <- NULL

x <- subset(grm.results.ranef, select = c(Trait2, LambAdult, new.variable, Effect))
x$Effect  <- formatC(x$Effect, digits = 4, format = "f")
x <- cast(x, formula = Trait2 + LambAdult ~ new.variable)
head(x)

models <- join(models, x)

x <- subset(grm.results.ranef, select = c(Trait2, LambAdult, new.variable, SE))
x$SE  <- formatC(x$SE, digits = 4, format = "f")
x$SE <- paste0("(", x$SE, ")")
x <- cast(x, formula = Trait2 + LambAdult ~ new.variable)
head(x)
x <- cbind(x, data.frame(Trait = NA, ID = NA, V_raw = NA, VP = NA, VP_SE = NA,
                          N = NA, N_ids = NA, Model = NA, Response = NA, Mean = NA))


models <- rbind(models, x)
head(models)



models <- arrange(models, Trait2, LambAdult)

models <- models[,c("Trait2", "LambAdult", "V_raw",  "N", "N_ids", "Mean", "VP", "VP_SE",
                    "Additive Genetic", "Birth Year", "Capture Year", "Mother Identity",
                    "Permanent Environment", "Plate ID", "Run Date", "Residual")]

models$VP[seq(2, nrow(models), 2)] <- paste0("(", models$VP_SE[seq(1, nrow(models), 2)], ")")

models$VP_SE <- NULL

models

write.table(models, "results/1_GRM_Ranef_Results_FORMATTED_For_Paper.txt", row.names = F, sep = "\t", quote = F)

#~~ Make nice figures

head(grm.results.ranef)

grm.results.ranef$new.variable2 <- factor(grm.results.ranef$new.variable,
                                  levels = rev(c("Additive Genetic","Permanent Environment","Birth Year","Capture Year",
                                             "Mother Identity", "Run Date", "Plate ID","Residual")))


grm.results.ranef$LambAdult <- relevel(as.factor(grm.results.ranef$LambAdult) , ref = "Lambs")

library(RColorBrewer)

alex.scheme <- rev(c("#E69F02","#F0E442","#CC79A7","#5BB2E8","#0470B0","#90EE8E","#009F73","#999A96"))

ggplot(grm.results.ranef, aes(LambAdult, Effect, fill = new.variable2)) +
  geom_bar(stat = "identity", colour = "grey50") +
  facet_wrap(~Trait2) +
  scale_fill_manual(values = alex.scheme) +
  labs(fill = "", x = "Age Class", y = "Proportion of Phenotypic Variance") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14))

ggsave("figs/1_Stacked_Variance_Components.png", height = 7, width = 10)

ggplot(grm.results.ranef, aes(new.variable2, Effect, fill = new.variable2)) +
  geom_bar(stat = "identity", colour = "grey50") +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0) +
  facet_grid(LambAdult~Trait2) +
  scale_fill_manual(values = alex.scheme) +
  labs(x = "Variable", y = "Proportion of VP", fill = "Variable") +
  theme_bw() +
  theme(axis.text.x  = element_text (angle = 270, hjust = 0, vjust = 0.3, size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.position = "none")

ggsave("figs/1_Variance_Components.png", height = 10, width = 10)
