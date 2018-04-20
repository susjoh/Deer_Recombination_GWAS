
AnalysisSuffix <- "b"

mapdata <- read.table("gcta/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)

binwidths <- c(20)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create SNP lists                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.tab <- read.table(paste0("gcta/regh2_",  AnalysisSuffix, "/regh2_master_",  AnalysisSuffix, "_X.txt"), header = T, stringsAsFactors = F)

regh2.tab$Analysis.ID <- as.character(regh2.tab$Analysis.ID)

regh2.tab <- subset(regh2.tab, Window.Size == 20)

load("gcta/recsumm_a.Rdata", verbose = T)

library(asreml)

source("r/makeGRM.R")

recsumm$RRID2 <- recsumm$RRID
recsumm.f$RRID2 <- recsumm.f$RRID
recsumm.m$RRID2 <- recsumm.m$RRID

source("r/ASReml.EstEffects.R")


grm.rest <- read.table(paste0("gcta/gcta/Deer_autoGRM_adj.grm.gz"))
ids.rest <- read.table(paste0("gcta/gcta/Deer_autoGRM_adj.grm.id"))

restinv <- makeGRM(grm.rest, ids.rest, id.vector = unique(recsumm$RRID))


for(i in subset(regh2.tab, CEL.LG == 34)$Analysis.ID){
  
  print(i)
  
  try({
    
    grm.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM.grm.gz"))
    ids.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM.grm.id"))
    
    reginv  <- makeGRM(grm.reg, ids.reg, id.vector = unique(recsumm$RRID), inv = F)

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
    
    x <- list(
      list(summary(model1, all = T), model1$converge, ASReml.EstEffects(model1)),
              list(summary(model2, all = T), model2$converge, ASReml.EstEffects(model2)),
              list(summary(model1.m, all = T), model1.m$converge, ASReml.EstEffects(model1.m)),
              list(summary(model2.m, all = T), model2.m$converge, ASReml.EstEffects(model2.m)),
              list(summary(model1.f, all = T), model1.f$converge, ASReml.EstEffects(model1.f)),
              list(summary(model2.f, all = T), model2.f$converge, ASReml.EstEffects(model2.f)))
    
    
    
    save(x, file = paste0("gcta/regh2_", AnalysisSuffix, "/", i,".RData"))
    rm(
      model1, model2, model1.m, model2.m,
      model1.f, model2.f, x)
  }) 
  
}



for(i in subset(regh2.tab, CEL.LG == 35)$Analysis.ID){
  
  print(i)
  
  try({
    
    grm.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM_adj.grm.gz"))
    ids.reg <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/",i,"_GRM_adj.grm.id"))
    
    reginv  <- makeGRM(grm.reg, ids.reg, id.vector = unique(recsumm$RRID))

    
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
    
    
    x <- list(
      list(summary(model1, all = T), model1$converge, ASReml.EstEffects(model1)),
      list(summary(model2, all = T), model2$converge, ASReml.EstEffects(model2)),
      list(summary(model1.m, all = T), model1.m$converge, ASReml.EstEffects(model1.m)),
      list(summary(model2.m, all = T), model2.m$converge, ASReml.EstEffects(model2.m)),
      list(summary(model1.f, all = T), model1.f$converge, ASReml.EstEffects(model1.f)),
      list(summary(model2.f, all = T), model2.f$converge, ASReml.EstEffects(model2.f)))
    
    
    
    save(x, file = paste0("gcta/regh2_", AnalysisSuffix, "/", i,".RData"))
    rm(
      model1, model2, model1.m, model2.m,
      model1.f, model2.f, x)
  }) 
  
}




#~~~~~~~~~~~~ PARSE RESULTS


mapdata$Cumu.Order <- 1:nrow(mapdata)
mapdata$Diff <- c(mapdata$Estimated.Mb.Position[1],diff(mapdata$Estimated.Mb.Position))
mapdata$Diff <- ifelse(mapdata$Diff < 0, 10000, mapdata$Diff)
mapdata$Cumu.Position <- cumsum(mapdata$Diff)


regh2.tab <- subset(regh2.tab, CEL.LG != 34)
names(regh2.tab)[3] <- "Start.Order"
regh2.tab$Stop.Order <- regh2.tab$Start.Order + regh2.tab$Window.Size - 1

temp <- mapdata[,c("CEL.order", "Cumu.Order", "Cumu.Position", "CEL.LG", "Estimated.Mb.Position", "BTA.Position")]

library(plyr)

names(temp) <- c("Start.Order", "Start.Cumu.Order", "Start.Cumu.Position", "CEL.LG", "Start.Position", "Start.BTA.Position")
regh2.tab <- join(regh2.tab, temp)
names(temp) <- c("Stop.Order", "Stop.Cumu.Order", "Stop.Cumu.Position", "CEL.LG", "Stop.Position", "Stop.BTA.Position")
regh2.tab <- join(regh2.tab, temp)


regh2.tab$Window.Bp <- regh2.tab$Stop.Position - regh2.tab$Start.Position

tapply(regh2.tab$Window.Bp, regh2.tab$Window.Size, mean)
tapply(regh2.tab$Window.Bp, regh2.tab$Window.Size, sd)


head(regh2.tab)
head(mapdata)

rm(temp)

tail(regh2.tab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse the likelihood results                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

restab <- list()
vartab <- list()

for(j in 1:nrow(regh2.tab)){
  
  if(j %in% seq(1, nrow(regh2.tab), 50)) print(paste("Running row", j, "of", nrow(regh2.tab)))
  
  if(file.exists(paste0("gcta/regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[j], ".RData"))){
    
    load(paste0("gcta/regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[j], ".RData"))
    
    model1 <- x[[1]][1][[1]]
    model2 <- x[[2]][1][[1]]
    model1.m <- x[[3]][1][[1]]
    model2.m <- x[[4]][1][[1]]
    model1.f <- x[[5]][1][[1]]
    model2.f <- x[[6]][1][[1]]
    
    
    
    restab[[j]] <- data.frame(Analysis.ID = regh2.tab$Analysis.ID[j],
                              RRID.Sex     = c("Both", "Male", "Female"),
                              model1.Li    = c(model1$loglik, model1.m$loglik, model1.f$loglik),
                              model2.Li    = c(model2$loglik, model2.m$loglik, model2.f$loglik))
    
    vartab[[j]] <- rbind(cbind(x[[1]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Both"  , Model = "wo.region", Effect = row.names(x[[1]][3][[1]])),
                         cbind(x[[2]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Both"  , Model = "wi.region", Effect = row.names(x[[2]][3][[1]])),
                         cbind(x[[3]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Male"  , Model = "wo.region", Effect = row.names(x[[3]][3][[1]])),
                         cbind(x[[4]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Male"  , Model = "wi.region", Effect = row.names(x[[4]][3][[1]])),
                         cbind(x[[5]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Female", Model = "wo.region", Effect = row.names(x[[5]][3][[1]])),
                         cbind(x[[6]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Female", Model = "wi.region", Effect = row.names(x[[6]][3][[1]])))
    
    
    rm(model1, model1.f, model1.m, model2, model2.f, model2.m, x)
    
  }
  
}

rm(j)

#~~ Convert lists to frames

restab <- data.table::rbindlist(restab)

gc()

vartab <- data.table::rbindlist(vartab)

#~~ Join with region information

restab <- join(restab, regh2.tab)
vartab <- join(vartab, regh2.tab)

#~~ Get test statistics from likelihood ratio test, chidist with 1df

restab$Chi2 <- 2*(restab$model2.Li - restab$model1.Li)
restab$P <- 1- pchisq(restab$Chi2, df = 1)

head(restab)

