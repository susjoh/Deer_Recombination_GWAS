#
# Run Animal Models and Genome-Wide Association Studies
# Susan Johnston
# April 2017
#


library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)
library(RepeatABEL)
library(asreml)


#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- subset(mapdata, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position))

#~~ Load ASReml and GenABEL Functions

source("r/ASReml.EstEffects.R")
source("r/ASReml.ExtractPredictors.R")
source("r/makeGRM.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Explore recombination rate variation    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(paste0("gcta/recsumm_", AnalysisSuffix, ".Rdata"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Is recombination rate heritable? GRM             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

grm.auto <- read.table("gcta/Deer31.recoded.v2.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("gcta/Deer31.recoded.v2.grm.id")  # CONTAINS ID LIST

grminv <- makeGRM(grm.auto, ids.auto, recsumm$RRID)

set.seed(568)

recsumm <- subset(recsumm, select = c(TotalRecombCount, RRID.Fhat3, RRID, RRID.Sex))
recsumm$RRID2 <- recsumm$RRID
nmales <- length(which(recsumm$RRID.Sex == "M"))

runSampling <- FALSE

if(runSampling == TRUE){
  
  for(i in 1:100){
    
    try({
      print(paste("Running iteration", i))
      
      # load(paste0("sampling/MF_iteration_", i, ".RData"))
      # recsumm.m$RRID2 <- recsumm.m$RRID
      # recsumm.f$RRID2 <- recsumm.f$RRID
      
      
      recsumm.m <- recsumm[sample(which(recsumm$RRID.Sex == "M"), size = nmales, replace = T),]
      recsumm.f <- recsumm[sample(which(recsumm$RRID.Sex == "F"), size = nmales, replace = T),]


      #~~ Run GRM models for Total Recombination Count

      grm.summ.RR.m <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                              random = ~ giv(RRID) + ide(RRID),
                              data = recsumm.m,
                              ginverse =  list(RRID = grminv),
                              na.method.X = "omit", na.omit.Y = "na.omit",
                              workspace = 500e+6, pworkspace = 500e+6, trace = F)

      grm.summ.RR.f <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                              random = ~ giv(RRID) + ide(RRID),
                              data = recsumm.f,
                              ginverse =  list(RRID = grminv),
                              na.method.X = "omit", na.omit.Y = "na.omit",
                              workspace = 500e+6, pworkspace = 500e+6, trace = F)

      asreml.results.m <- ASReml.EstEffects(grm.summ.RR.m)
      asreml.results.f <- ASReml.EstEffects(grm.summ.RR.f)


      gwas.acc.m   <- rGLS(TotalRecombCount ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.m, id.name = "RRID", GRM = grminv)
      gwas.acc.f   <- rGLS(TotalRecombCount ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.f, id.name = "RRID", GRM = grminv)

      lambda.m <- lambda(gwas.acc.m)$estimate
      lambda.f <- lambda(gwas.acc.f)$estimate

      lambda.m <- ifelse(lambda.m < 1, 1, lambda.m)
      lambda.f <- ifelse(lambda.m < 1, 1, lambda.f)

      gwas.acc.m <- results(gwas.acc.m)
      gwas.acc.f <- results(gwas.acc.f)

      gwas.acc.m$SNP.Name <- row.names(gwas.acc.m)
      gwas.acc.f$SNP.Name <- row.names(gwas.acc.f)

      gwas.acc.m$chi2.1df.adj <- qchisq(p = 1-gwas.acc.m$P1df, df =  1)/lambda.m
      gwas.acc.f$chi2.1df.adj <- qchisq(p = 1-gwas.acc.f$P1df, df =  1)/lambda.f

      gwas.acc.m$Pc1df <- 1-pchisq(gwas.acc.m$chi2.1df.adj, 1)
      gwas.acc.f$Pc1df <- 1-pchisq(gwas.acc.f$chi2.1df.adj, 1)

      
      #########################
      
      
      grm.rest <- read.table(paste0("gcta/regh2_c/LG12.c.290.10_wo_GRM_adj.grm.gz"))
      ids.rest <- read.table(paste0("gcta/regh2_c/LG12.c.290.10_wo_GRM_adj.grm.id"))
      
      grm.reg <- read.table(paste0("gcta/regh2_c/LG12.c.290.10_GRM_adj.grm.gz"))
      ids.reg <- read.table(paste0("gcta/regh2_c/LG12.c.290.10_GRM_adj.grm.id"))
      
      reginv  <- makeGRM(grm.reg, ids.reg, id.vector = unique(recsumm$RRID))
      restinv <- makeGRM(grm.rest, ids.rest, id.vector = unique(recsumm$RRID))
      
    
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
      
      x <- list(list(summary(model1.m, all = T), model1.m$converge, ASReml.EstEffects(model1.m)),
                list(summary(model2.m, all = T), model2.m$converge, ASReml.EstEffects(model2.m)),
                list(summary(model1.f, all = T), model1.f$converge, ASReml.EstEffects(model1.f)),
                list(summary(model2.f, all = T), model2.f$converge, ASReml.EstEffects(model2.f)))
      
      
      
      save(x, recsumm.m, recsumm.f, asreml.results.m, asreml.results.f, gwas.acc.m, gwas.acc.f, file = paste0("sampling/MF_iteration_", i, ".RData"))
      
      rm(x, recsumm.m, recsumm.f, asreml.results.m, asreml.results.f, gwas.acc.m, gwas.acc.f, lambda.m, lambda.f)
    })
  }
}


#~~ Parse the data

h2res <- NULL
gwasres <- NULL
topres <- NULL
restab <- list()
vartab <- list()

for(i in 1:100){
  
  print(i)
  load(paste0("sampling/MF_iteration_", i, ".RData"))
  
  
  h2res <- rbind(h2res,
                 cbind(asreml.results.f, Iteration = i, Sex = "F", EffectName = row.names(asreml.results.f)),
                 cbind(asreml.results.m, Iteration = i, Sex = "M", EffectName = row.names(asreml.results.m)))
  
  
  
  topres <- rbind(topres,
                  cbind(gwas.acc.m[which(gwas.acc.m$SNP.Name %in% c("cela1_red_10_25661750", "cela1_red_10_26005249")),], Iteration = i, Sex = "M"),
                  cbind(gwas.acc.f[which(gwas.acc.f$SNP.Name %in% c("cela1_red_10_25661750", "cela1_red_10_26005249")),], Iteration = i, Sex = "F"))
  
  
  gwas.acc.m <- arrange(gwas.acc.m, P1df)[1:1000,]
  gwas.acc.f <- arrange(gwas.acc.f, P1df)[1:1000,]
  
  
  gwasres <- rbind(gwasres,
                   cbind(gwas.acc.m, Iteration = i, Sex = "M"),
                   cbind(gwas.acc.f, Iteration = i, Sex = "F"))

  

  model1.m <- x[[1]][1][[1]]
  model2.m <- x[[2]][1][[1]]
  model1.f <- x[[3]][1][[1]]
  model2.f <- x[[4]][1][[1]]
  
  
  
  restab[[i]] <- data.frame(Analysis.ID = i,
                            RRID.Sex     = c("Male", "Female"),
                            model1.Li    = c(model1.m$loglik, model1.f$loglik),
                            model2.Li    = c(model2.m$loglik, model2.f$loglik))
  
  vartab[[i]] <- rbind(cbind(x[[1]][3][[1]], Analysis.ID = i, RRID.Sex = "Male"  , Model = "wo.region", Effect = row.names(x[[1]][3][[1]])),
                       cbind(x[[2]][3][[1]], Analysis.ID = i, RRID.Sex = "Male"  , Model = "wi.region", Effect = row.names(x[[2]][3][[1]])),
                       cbind(x[[3]][3][[1]], Analysis.ID = i, RRID.Sex = "Female", Model = "wo.region", Effect = row.names(x[[3]][3][[1]])),
                       cbind(x[[4]][3][[1]], Analysis.ID = i, RRID.Sex = "Female", Model = "wi.region", Effect = row.names(x[[4]][3][[1]])))
  
  
  rm(model1.f, model1.m, model2.f, model2.m, x)
  
  rm(recsumm.m, recsumm.f, asreml.results.m, asreml.results.f, gwas.acc.m, gwas.acc.f)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Examine heritability estimates         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)

obsdata <- data.frame(Sex = c("F", "M"), Effect = c(0.11, 0.07), SE = c( 0.06, 0.11))
obsdata$Jitter <- c(1, 1.4)

h2res$Jitter <- NA

for(i in 1:100){
  h2res$Jitter[which(h2res$Iteration == i & h2res$Sex == "F")] <- runif(1, min = 0.9, max = 1.1)
  h2res$Jitter[which(h2res$Iteration == i & h2res$Sex == "M")] <- runif(1, min = 1.3, max = 1.5)
}

  
head(h2res)


ggplot(subset(h2res, EffectName == "giv(RRID).giv"), aes(Jitter, Effect)) + 
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE), width = 0, alpha = 0.7) +
  geom_point(data = obsdata, aes(Jitter, Effect), colour = "red") +
  geom_errorbar(data = obsdata, aes(Jitter, ymin = Effect - SE, ymax = Effect + SE), width = 0, colour = "red") +
  coord_cartesian(xlim = c(0.7, 1.7))


ggplot(subset(h2res, EffectName == "giv(RRID).giv"), aes(Sex, Effect)) + 
  geom_jitter(width = 0.2) +
  geom_boxplot(notch = T, alpha = 0) +
  geom_point(data = obsdata, aes(as.numeric(Sex)+0.05, Effect), colour = "red") +
  geom_errorbar(data = obsdata, aes(as.numeric(Sex)+0.05, ymin = Effect - SE, ymax = Effect + SE), width = 0, colour = "red") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
labs(y = "Heritability estimate")

ggsave("figs/5_Male_Female_Sampling_heritability.png", width = 6, height = 6)


t.test(subset(h2res, EffectName == "giv(RRID).giv")$Effect ~ subset(h2res, EffectName == "giv(RRID).giv")$Sex, var.equal = F, paired = T)  # cela1_red_10_25661750


?t.test

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Examine GWAS estimates                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Top GWAS hits
topres$RRID.Sex <- NA
topres$RRID.Sex[which(topres$Sex == "F")] <- "Female"
topres$RRID.Sex[which(topres$Sex == "M")] <- "Male"


ggplot(topres, aes(RRID.Sex, -log10(Pc1df))) +
  geom_jitter(width = 0.2) +
  geom_boxplot(notch = T, alpha = 0) +
  facet_wrap(~SNP.Name) +
  theme_bw() +
  labs(y = "-log10(P)") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())

ggsave("figs/5_Male_Female_Sampling_GWAS.png", width = 8, height = 6)

topres.1 <- subset(topres, SNP.Name == "cela1_red_10_25661750")
topres.2 <- subset(topres, SNP.Name == "cela1_red_10_26005249")


t.test(-log10(topres.1$Pc1df) ~ topres.1$Sex, var.equal = F)  # cela1_red_10_25661750
t.test(-log10(topres.2$Pc1df) ~ topres.1$Sex, var.equal = F)  # cela1_red_10_26005249

#~~ top regh2 hits

restab <- do.call("rbind", restab)
head(restab)

restab$Chi2 <- 2*(restab$model2.Li - restab$model1.Li)
restab$P <- 1- pchisq(restab$Chi2, df = 1)

ggplot(restab, aes(RRID.Sex, -log10(P))) +
  geom_jitter(width = 0.2) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw() +
  labs(x = "Sex", y = "-log10(P)") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())

ggsave("figs/5_Male_Female_Sampling_reg_h2.png", width = 6, height = 6)


t.test(-log10(restab$P) ~ restab$RRID.Sex, var.equal = F)


#~~ Prepare to merge with map information

mapdata$Diff <- c(1, diff(mapdata$Estimated.Mb.Position))
mapdata$Diff <- ifelse(mapdata$Diff < 0, 1e5, mapdata$Diff)
mapdata$Cumu <- cumsum(mapdata$Diff)

str(mapdata)
str(gwasres)

gwasres <- join(gwasres, mapdata)

ggplot(gwasres, aes(Cumu, -log10(Pc1df))) + 
  geom_point() +
  facet_wrap(~Sex)

head(gwasres)

gwasres$Order <- NA

for(i in 1:100) gwasres$Order[which(gwasres$Iteration == i)] <- c(1:1000, 1:1000)

table(subset(gwasres, Order %in% 1)$CEL.LG, subset(gwasres, Order %in% 1)$Sex)
