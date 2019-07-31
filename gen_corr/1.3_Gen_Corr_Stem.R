Sys.setenv("ASREML_LICENSE_FILE" = "/exports/igmm/software/pkg/el7/apps/asreml/3.0gm/bin/asreml.lic") 

library(asreml)
library(magrittr)

load("1.3_Gen_Corr_Data.RData")

print(paste0("Running model ", i))

eval(parse(text = paste0("fit1 <- asreml(fixed  = ", x1$Model[i], ",
                                  rcov   = ~ units:idh(trait, init = NA),
                                  data = BEASTX,
                                  ginverse = list(ID = ainv),
                                  workspace = 500e+6, pworkspace = 500e+6,
                                  maxiter = 100, na.method.Y = \"include\", na.method.X = \"include\")")))

save(fit1, file = paste0("GenCorr_Ped_", i, ".RData"))

restab <- cbind(ModelNo = i, summary(fit1)$varcomp)
restab$Effect <- row.names(restab)

restab

write.table(restab, paste0("Ped_Res_", i, ".txt"), row.names = F, sep = "\t", quote = F)

x1$Model
x1$Model <- gsub("ped", "giv", x1$Model)

#~~ GRM Models

BEASTX$ID2 <- as.character(BEASTX$ID)
BEASTX <- subset(BEASTX, ID2 %in% dimnames(grminv)[[1]]) %>% droplevels


eval(parse(text = paste0("fit1 <- asreml(fixed  = ", x1$Model[i], ",
                                  rcov   = ~ units:idh(trait, init = NA),
                                  data = BEASTX,
                                  ginverse = list(ID = grminv),
                                  workspace = 500e+6, pworkspace = 500e+6,
                                  maxiter = 100, na.method.Y = \"include\", na.method.X = \"include\")")))

restab.grm <- cbind(ModelNo = i, summary(fit1)$varcomp)
restab.grm$Effect <- row.names(restab.grm)

restab.grm

write.table(restab.grm, paste0("GRM_Res_", i, ".txt"), row.names = F, sep = "\t", quote = F)


save(fit1, file = paste0("GenCorr_GRM_", i, ".RData"))



