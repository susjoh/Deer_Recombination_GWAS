
load("../recsumm_a.Rdata")

library(coxme)

source("../makeGRM.R")

# ALL IDS

grm.rest <- read.table(paste0(i,"_wo_GRM_adj.grm.gz"))
ids.rest <- read.table(paste0(i,"_wo_GRM_adj.grm.id"))

grm.reg <- read.table(paste0(i,"_GRM_adj.grm.gz"))
ids.reg <- read.table(paste0(i,"_GRM_adj.grm.id"))

reginv  <- makeGRM(grm.reg, ids.reg, id.vector = unique(recsumm$RRID))
restinv <- makeGRM(grm.rest, ids.rest, id.vector = unique(recsumm$RRID))

model1 <- lmekin(TotalRecombCount ~ RRID.Sex + RRID.Fhat3 + (1|RRID) + (1|RRID), data=recsumm, varlist=restinv)

model2 <- lmekin(TotalRecombCount ~ RRID.Sex + RRID.Fhat3 + (1|RRID), data=recsumm, varlist=list(restinv, reginv))

model1.m <- lmekin(TotalRecombCount ~ RRID.Fhat3 + (1|RRID) + (1|RRID), data=recsumm.m, varlist=restinv)

model2.m <- lmekin(TotalRecombCount ~ RRID.Fhat3 + (1|RRID), data=recsumm.m, varlist=list(restinv, reginv))

model1.f <- lmekin(TotalRecombCount ~ RRID.Fhat3 + (1|RRID) + (1|RRID), data=recsumm.f, varlist=restinv)

model2.f <- lmekin(TotalRecombCount ~ RRID.Fhat3 + (1|RRID), data=recsumm.f, varlist=list(restinv, reginv))


rm(grm.rest, ids.rest, grm.reg, ids.reg)


save.image(file = paste0("model", i, ".RData"))
