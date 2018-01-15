
auto.vec <- 12

AnalysisSuffix <- "e"

mapdata <- read.table("gcta/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)

binwidths <- c(6, 10, 20)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create SNP lists                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.tab <- read.table(paste0("gcta/regh2_",  AnalysisSuffix, "/regh2_master_",  AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

regh2.tab$Analysis.ID <- as.character(regh2.tab$Analysis.ID)


load("gcta/recsumm_a.Rdata")

library(asreml)

source("r/makeGRM.R")

recsumm$RRID2 <- recsumm$RRID
recsumm.f$RRID2 <- recsumm.f$RRID
recsumm.m$RRID2 <- recsumm.m$RRID

source("r/ASReml.EstEffects.R")

regh2.tab <- subset(regh2.tab, CEL.LG != 34)

for(i in regh2.tab$Analysis.ID[1:nrow(regh2.tab)]){
  
  print(i)
  
  try({
    
    grm.rest <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/", i,"_wo_GRM.grm.gz"))
    ids.rest <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_wo_GRM.grm.id"))
    
    grm.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM.grm.gz"))
    ids.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM.grm.id"))
    
    reginv  <- makeGRM(grm.reg, ids.reg, id.vector = unique(recsumm$RRID))
    restinv <- makeGRM(grm.rest, ids.rest, id.vector = unique(recsumm$RRID))
    
    
    model1 <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                     random = ~ giv(RRID) + ide(RRID),
                     data = recsumm,
                     ginverse =  list(RRID = restinv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    model2 <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                     random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                     data = recsumm,
                     ginverse =  list(RRID = restinv, RRID2 = reginv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    model1.m <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.m,
                       ginverse =  list(RRID = restinv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    model2.m <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                       data = recsumm.m,
                       ginverse =  list(RRID = restinv, RRID2 = reginv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    model1.f <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.f,
                       ginverse =  list(RRID = restinv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    model2.f <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                       data = recsumm.f,
                       ginverse =  list(RRID = restinv, RRID2 = reginv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6, trace = F)
    
    x <- list(list(summary(model1, all = T), model1$converge, ASReml.EstEffects(model1)),
              list(summary(model2, all = T), model2$converge, ASReml.EstEffects(model2)),
              list(summary(model1.m, all = T), model1.m$converge, ASReml.EstEffects(model1.m)),
              list(summary(model2.m, all = T), model2.m$converge, ASReml.EstEffects(model2.m)),
              list(summary(model1.f, all = T), model1.f$converge, ASReml.EstEffects(model1.f)),
              list(summary(model2.f, all = T), model2.f$converge, ASReml.EstEffects(model2.f)))
    
    
    
    save(x, file = paste0("gcta/regh2_", AnalysisSuffix, "/", i,".RData"))
    rm(model1, model2, model1.m, model2.m, model1.f, model2.f, x)
  }) 
  
}

