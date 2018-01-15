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

#~~ Read in pedigree files

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)
pedigree <- pedigree[,c(1, 3, 2)]
for(i in 1:3) pedigree[which(pedigree[,i] == 0),i] <- NA
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])


famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ Read in phenotype file

recsumm <- read.table(paste0("results/3_Recombination_Phenotype_Data_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

rectab  <- read.table(paste0("results/2_Per_Chromosome_Recomb_Final_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)


#~~ Get IDs that are in recsumm

famtab <- read.table("gcta/Deer31.recoded.v2.fam")
head(famtab)
famtab <- subset(famtab, V2 %in% recsumm$RRID)

write.table(famtab[,1:2], "gcta/idlist.txt", row.names = F, col.names = F, quote = F)

#~~ read in SNP positions and PAR SNPs

chrom.map <- read.table(paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- subset(mapdata, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position))

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Load ASReml and GenABEL Functions

source("r/ASReml.EstEffects.R")
source("r/ASReml.ExtractPredictors.R")
source("r/makeGRM.R")
source("r/multiplot.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Explore recombination rate variation    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm <- subset(recsumm, TotalChrCount == 33)
rectab <- subset(rectab, RRID %in% recsumm$RRID)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1, col = "white") + labs(x = "Total Crossover Count")
ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1) + facet_wrap(~RRID.Sex, ncol = 1)

ggplot(recsumm, aes(RRID.Sex, TotalRecombCount, fill = RRID.Sex)) +
  geom_boxplot(notch = T, width = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank())

tapply(recsumm$TotalRecombCount, recsumm$RRID.Sex, mean)

ggplot(recsumm, aes(TotalRecombCount, TotalInfLoci)) + geom_point() + stat_smooth(method = "lm")

ggplot(recsumm, aes(TotalRecombCount, RRID.Fhat3)) + geom_point() + stat_smooth(method = "lm") +
  labs(x = "Autosomal Crossover Count", y = "Inbreeding Coefficient (Fhat3)")
ggsave("figs/4_Inbreeding_vs_CrossoverCount.png", width = 6, height = 6)

summary(lm(RRID.Fhat3 ~ TotalRecombCount, data = recsumm))

ggplot(recsumm, aes(TotalInfLoci, RRID.Fhat3)) + geom_point() + stat_smooth(method = "lm")


ggplot(recsumm, aes(RRID.Sex, TotalRecombCount)) + geom_boxplot(width = 0.5, notch = T)
ggplot(recsumm, aes(RRID.Age, TotalRecombCount, colour = RRID.Sex)) + geom_point() + stat_smooth()

ggplot(recsumm, aes(RRID.Sex, AcroCentroPC.20/TotalRecombCount)) + geom_boxplot(width = 0.5, notch = T)
ggplot(recsumm, aes(RRID.Age, AcroCentroPC.20/TotalRecombCount, colour = RRID.Sex)) + geom_point() + stat_smooth()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Is recombination rate heritable? PEDIGREE        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Pedigree Models

names(pedigree) <- c("RRID", "RRID.FATHER", "RRID.MOTHER")


ainv <- asreml.Ainverse(pedigree)$ginv

recsumm$RRID <- as.factor(recsumm$RRID)
recsumm$RRID.Sex <- as.factor(recsumm$RRID.Sex)
recsumm$RRID.MOTHER <- as.factor(recsumm$RRID.MOTHER)
recsumm$RRID.FATHER <- as.factor(recsumm$RRID.FATHER)

recsumm$AcroCentroPC.20.Adjusted <- recsumm$AcroCentroPC.20/recsumm$TotalRecombCount
recsumm$AcroTeloPC.80.Adjusted <- 1-recsumm$AcroCentroPC.20.Adjusted
recsumm$AcroTeloPC.80 <- recsumm$TotalRecombCount-recsumm$AcroCentroPC.20


recsumm.m <- droplevels(recsumm[which(recsumm$RRID.Sex == "M"),])
recsumm.f <- droplevels(recsumm[which(recsumm$RRID.Sex == "F"),])

save(recsumm, recsumm.m, recsumm.f, file = paste0("gcta/recsumm_", AnalysisSuffix, ".Rdata"))

#~~ Run Pedigree models for Total Recombination Count

ped.summ.RR <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3 + RRID.Age,
                      random = ~ ped(RRID) + ide(RRID) + RRID.MOTHER + RRID.BirthYear,
                      data = recsumm,
                      ginverse =  list(RRID = ainv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(ped.summ.RR)
wald.asreml(ped.summ.RR)


ped.summ.RR.m <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3 + RRID.Age,
                        random = ~ ped(RRID) + ide(RRID) + RRID.MOTHER + RRID.BirthYear,
                        data = recsumm.m,
                        ginverse =  list(RRID = ainv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(ped.summ.RR.m)
wald.asreml(ped.summ.RR.m)


ped.summ.RR.f <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3 + RRID.Age,
                        random = ~ ped(RRID) + ide(RRID) + RRID.MOTHER + RRID.BirthYear,
                        data = recsumm.f,
                        ginverse =  list(RRID = ainv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(ped.summ.RR.f)
wald.asreml(ped.summ.RR.f)

#~~ Run Pedigree models for peri-centromeric recombination proportion

ped.summ.PCR <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Sex + RRID.Fhat3,
                       random = ~ ped(RRID) + ide(RRID),
                       data = recsumm,
                       family = asreml.binomial(),
                       ginverse =  list(RRID = ainv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.PCR)
wald.asreml(ped.summ.PCR)


ped.summ.PCR.m <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.m,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.PCR.m)
wald.asreml(ped.summ.PCR.m)


ped.summ.PCR.f <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.f,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.PCR.f)
wald.asreml(ped.summ.PCR.f)



#~~ Run Pedigree models for peri-centromeric recombination proportion

ped.summ.TER <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Sex + RRID.Fhat3,
                       random = ~ ped(RRID) + ide(RRID),
                       data = recsumm,
                       family = asreml.binomial(),
                       ginverse =  list(RRID = ainv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.TER)
wald.asreml(ped.summ.TER)


ped.summ.TER.m <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.m,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.TER.m)
wald.asreml(ped.summ.TER.m)


ped.summ.TER.f <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.f,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.TER.f)
wald.asreml(ped.summ.TER.f)


#~~ Run Pedigree models for mean distance to first crossover

recsumm$Mean.CentroToCO <- recsumm$Mean.CentroToCO/1e6

ped.summ.CCO <- asreml(fixed = Mean.CentroToCO ~ RRID.Sex + RRID.Fhat3,
                       random = ~ ped(RRID) + ide(RRID),
                       data = recsumm,
                       ginverse =  list(RRID = ainv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.CCO)
wald.asreml(ped.summ.CCO)


ped.summ.CCO.m <- asreml(fixed = Mean.CentroToCO/1e6 ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.m,
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.CCO.m)
wald.asreml(ped.summ.CCO.m)


ped.summ.CCO.f <- asreml(fixed = Mean.CentroToCO/1e6 ~ RRID.Fhat3,
                         random = ~ ped(RRID) + ide(RRID),
                         data = recsumm.f,
                         ginverse =  list(RRID = ainv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(ped.summ.CCO.f)
wald.asreml(ped.summ.CCO.f)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Is recombination rate heritable? GRM             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


grm.auto <- read.table("gcta/Deer31.recoded.v2.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("gcta/Deer31.recoded.v2.grm.id")  # CONTAINS ID LIST

grminv <- makeGRM(grm.auto, ids.auto, recsumm$RRID) # vector of IDs from the datasset that you use for the asreml model
dim(grminv)


#~~ Run GRM models for Total Recombination Count

grm.summ.RR <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                      random = ~ giv(RRID) + ide(RRID),
                      data = recsumm,
                      ginverse =  list(RRID = grminv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(grm.summ.RR)
wald.asreml(grm.summ.RR)

grm.summ.RR.m <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.m,
                        ginverse =  list(RRID = grminv), 
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(grm.summ.RR.m)
wald.asreml(grm.summ.RR.m)


grm.summ.RR.f <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.f,
                        ginverse =  list(RRID = grminv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

ASReml.EstEffects(grm.summ.RR.f)
wald.asreml(grm.summ.RR.f)


#~~ Run without heritabiltiy for significance testing


grm.summ.RR.woh2 <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                      random = ~ ide(RRID),
                      data = recsumm,
                      ginverse =  list(RRID = grminv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

grm.summ.RR.m.woh2 <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ ide(RRID),
                        data = recsumm.m,
                        ginverse =  list(RRID = grminv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

grm.summ.RR.f.woh2 <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ ide(RRID),
                        data = recsumm.f,
                        ginverse =  list(RRID = grminv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)


grm.summ.RR.forvp  <- asreml(fixed = TotalRecombCount ~ RRID.Sex + RRID.Fhat3,
                      random = ~ giv(RRID) + ide(RRID),
                      data = recsumm,
                      ginverse =  list(RRID = grminv), rcov = ~idv(units),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

grm.summ.RR.m.forvp <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.m,
                        ginverse =  list(RRID = grminv), rcov = ~idv(units),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)


grm.summ.RR.f.forvp <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.f,
                        ginverse =  list(RRID = grminv), rcov = ~idv(units),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)







#~~ Run GRM models for peri-centromeric recombination proportion

grm.summ.PCR <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm,
                       family = asreml.binomial(),
                       ginverse =  list(RRID = grminv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.PCR)
wald.asreml(grm.summ.PCR)


grm.summ.PCR.m <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Fhat3,
                         random = ~ giv(RRID) + ide(RRID),
                         data = recsumm.m,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = grminv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.PCR.m)
wald.asreml(grm.summ.PCR.m)


grm.summ.PCR.f <- asreml(fixed = AcroCentroPC.20.Adjusted ~ RRID.Fhat3,
                         random = ~ giv(RRID) + ide(RRID),
                         data = recsumm.f,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = grminv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.PCR.f)
wald.asreml(grm.summ.PCR.f)



#~~ Run GRM models for peri-centromeric recombination proportion

grm.summ.TER <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Sex + RRID.Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm,
                       family = asreml.binomial(),
                       ginverse =  list(RRID = grminv),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.TER)
wald.asreml(grm.summ.TER)


grm.summ.TER.m <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Fhat3,
                         random = ~ giv(RRID) + ide(RRID),
                         data = recsumm.m,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = grminv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.TER.m)
wald.asreml(grm.summ.TER.m)


grm.summ.TER.f <- asreml(fixed = AcroTeloPC.80.Adjusted ~ RRID.Fhat3,
                         random = ~ giv(RRID) + ide(RRID),
                         data = recsumm.f,
                         family = asreml.binomial(),
                         ginverse =  list(RRID = grminv),
                         na.method.X = "omit", na.omit.Y = "na.omit",
                         workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.TER.f)
wald.asreml(grm.summ.TER.f)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Bivariate Animal Models                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create female and male specific columns to allow the model to run:

recsumm$Male.RR   <- ifelse(recsumm$RRID.Sex == "M"  , recsumm$TotalRecombCount, NA)
recsumm$Female.RR <- ifelse(recsumm$RRID.Sex == "F", recsumm$TotalRecombCount, NA)

#~~ Run the raw model.

grm.summ.bivar.raw <- asreml(fixed  = cbind(Female.RR, Male.RR) ~ trait + trait:RRID.Fhat3, 
                             random = ~ corgh(trait, init = c(0.5, 1, 1)):giv(RRID)
                             + idh(trait):ide(RRID),
                             rcov   = ~ units:idh(trait, init = NA),
                             data = recsumm,
                             ginverse = list(RRID = grminv),
                             workspace = 500e+6, pworkspace = 500e+6,
                             maxiter = 100)

#~~ Constrain rA to be zero

grm.summ.bivar.r0 <- asreml(fixed     = cbind(Female.RR, Male.RR) ~ trait+ trait:RRID.Fhat3, 
                            random    = ~ idh(trait):giv(RRID) + idh(trait):ide(RRID),
                            rcov      = ~ units:idh(trait, init = NA),
                            data      = recsumm,
                            ginverse  = list(RRID = grminv),
                            workspace = 500e+6, pworkspace = 500e+6)

#~~ Constrain rA to be one (or very close to one)


myinit.1        <- c(0.99999, 1, 1)
names(myinit.1) <- c("F", "P", "P")

grm.summ.bivar.r1 <- asreml(fixed     = cbind(Female.RR, Male.RR) ~ trait+ trait:RRID.Fhat3, 
                            random    = ~ corgh(trait, init = myinit.1):giv(RRID) + idh(trait):ide(RRID),
                            rcov      = ~ units:idh(trait, init = NA),
                            data      = recsumm,
                            ginverse  = list(RRID = grminv),
                            workspace = 500e+6, pworkspace = 500e+6)



#~~ Constrain variances to be the same, estaimate rA


grm.summ.bivar.cons <- asreml(fixed     = cbind(Female.RR, Male.RR) ~ trait+ trait:RRID.Fhat3, 
                              random    = ~ corgv(trait):giv(RRID) + idh(trait):ide(RRID),
                              rcov      = ~ units:idh(trait, init = NA),
                              data      = recsumm,
                              ginverse  = list(RRID = grminv),
                              G.param   = grm.summ.bivar.raw$G.param,
                              R.param   = grm.summ.bivar.raw$R.param,
                              workspace = 500e+6, pworkspace = 500e+6)


#~~ Constrain variances to be the same and rA =1

grm.summ.bivar.cons.r1 <- asreml(fixed     = cbind(Female.RR, Male.RR) ~ trait+ trait:RRID.Fhat3, 
                                 random    = ~ giv(RRID, var=T) + idh(trait):ide(RRID),
                                 rcov      = ~ units:idh(trait, init = NA),
                                 data      = recsumm,
                                 ginverse  = list(RRID = grminv),
                                 workspace = 500e+6, pworkspace = 500e+6)

#~~ Constrain variances to be the same and rA =0

grm.summ.bivar.cons.r0 <- asreml(fixed     = cbind(Female.RR, Male.RR) ~ trait+ trait:RRID.Fhat3, 
                                 random    = ~ idv(trait):giv(RRID) + idh(trait):ide(RRID),
                                 rcov      = ~ units:idh(trait, init = NA),
                                 data      = recsumm,
                                 ginverse  = list(RRID = grminv),
                                 workspace = 500e+6, pworkspace = 500e+6)

#~~ Constrain residual variances to be the same

grm.summ.bivar.cons.resid <- asreml(fixed  = cbind(Female.RR, Male.RR) ~ trait + trait:RRID.Fhat3, 
                                    random = ~ corgh(trait, init = c(0.5, 1, 1)):giv(RRID)
                                    + idh(trait):ide(RRID),
                                    rcov   = ~ units:idv(trait, init = NA),
                                    data = recsumm,
                                    ginverse = list(RRID = grminv),
                                    workspace = 500e+6, pworkspace = 500e+6,
                                    maxiter = 100)


#~~ Summary shows the additive genetic correlation in the first row, individual variances in rest of rows:


heritability.results <- list()

heritability.results[["bivar.fixef"]] <- summary(grm.summ.bivar.raw, all = T)$coef.fixed
heritability.results[["bivar.ranef"]] <- summary(grm.summ.bivar.raw)$varcomp

# Significantly different from zero?
heritability.results[["bivar.sig.diff.from.zero"]] <- 1 - pchisq(2*(summary(grm.summ.bivar.raw)$loglik - summary(grm.summ.bivar.r0)$loglik), df = 1)

# Significantly different from one?
heritability.results[["bivar.sig.diff.from.one"]] <- 1 - pchisq(2*(summary(grm.summ.bivar.raw)$loglik - summary(grm.summ.bivar.r1)$loglik), df = 1)

# Significantly different additive genetic variance?
heritability.results[["bivar.sig.diff.Va"]] <- 1 - pchisq(2*(summary(grm.summ.bivar.raw)$loglik - summary(grm.summ.bivar.cons)$loglik), df = 1)

# Significantly different residual variance?
heritability.results[["bivar.sig.diff.Vr"]] <- 1 - pchisq(2*(summary(grm.summ.bivar.raw)$loglik - summary(grm.summ.bivar.cons.resid)$loglik), df = 1)


heritability.results

save(heritability.results, ped.summ.RR,  ped.summ.RR.m , ped.summ.RR.f , ped.summ.PCR ,
     ped.summ.PCR.m, ped.summ.PCR.f, ped.summ.TER, ped.summ.TER.m, ped.summ.TER.f,
     grm.summ.RR, grm.summ.RR.m, grm.summ.RR.f, grm.summ.RR.f.woh2, grm.summ.RR.m.woh2, grm.summ.RR.woh2,
     grm.summ.RR.forvp, grm.summ.RR.m.forvp, grm.summ.RR.f.forvp,
     grm.summ.PCR, grm.summ.PCR.m, grm.summ.PCR.f, grm.summ.TER, grm.summ.TER.m, grm.summ.TER.f, file = "results/4_Heritability_Results.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. GWAS of Recombination Rate                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


grm.auto <- read.table("gcta/Deer31.recoded.v2.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("gcta/Deer31.recoded.v2.grm.id")  # CONTAINS ID LIST

grminv <- makeGRM(grm.auto, ids.auto, recsumm$RRID) # vector of IDs from the datasset that you use for the asreml model
grminv.m <- makeGRM(grm.auto, ids.auto, recsumm.m$RRID) # vector of IDs from the datasset that you use for the asreml model
grminv.f <- makeGRM(grm.auto, ids.auto, recsumm.f$RRID) # vector of IDs from the datasset that you use for the asreml model


gwas.acc.all <- rGLS(TotalRecombCount ~ RRID.Sex + RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm, id.name = "RRID", GRM = grminv)
gwas.acc.m   <- rGLS(TotalRecombCount ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.m, id.name = "RRID", GRM = grminv)
gwas.acc.f   <- rGLS(TotalRecombCount ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.f, id.name = "RRID", GRM = grminv)


# gwas.pcr.all <- rGLS(AcroCentroPC.20 ~ 1, genabel.data = abeldata, phenotype.data = recsumm, id.name = "RRID")
# gwas.pcr.m   <- rGLS(AcroCentroPC.20 ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.m, id.name = "RRID")
# gwas.pcr.f   <- rGLS(AcroCentroPC.20 ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.f, id.name = "RRID")
# 
# gwas.ter.all <- rGLS(AcroTeloPC.80 ~ 1, genabel.data = abeldata, phenotype.data = recsumm, id.name = "RRID")
# gwas.ter.m   <- rGLS(AcroTeloPC.80 ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.m, id.name = "RRID")
# gwas.ter.f   <- rGLS(AcroTeloPC.80 ~ RRID.Fhat3, genabel.data = abeldata, phenotype.data = recsumm.f, id.name = "RRID")


#~~ Process output and add deer positions

res.vec <- ls()[grep("gwas", ls())]

restab <- NULL

for(i in res.vec){
  
  eval(parse(text = paste0("x <- results(", i, ")")))
  x$SNP.Name <- row.names(x)
  x$Analysis <- i
  x <- arrange(x, P1df)
  x$ExpP <- seq(1/nrow(x), 1, 1/nrow(x))
  
  x$chi2.1df.adj <- qchisq(p = 1-x$P1df, df =  1)/eval(parse(text = paste0("lambda(", i, ")$estimate")))
  x$Pc1df <- 1-pchisq(x$chi2.1df.adj, 1)
  restab <- rbind(restab, x)
  
  rm(x)
}

#~~ restab analysis type

restab$Phenotype <- sapply(restab$Analysis, function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
restab$Sex <- sapply(restab$Analysis, function(x) strsplit(x, split = ".", fixed = T)[[1]][3])

head(restab)

#~~ Get Cumulative positions

unique(restab$SNP.Name[which(!restab$SNP.Name %in% mapdata$SNP.Name)])

mapdata <- rbind(mapdata,
                 data.frame(SNP.Name = unique(restab$SNP.Name[which(!restab$SNP.Name %in% mapdata$SNP.Name)]),
                            CEL.LG = 0,
                            Estimated.Mb.Position = 0))


head(mapdata)
tail(mapdata)


mapdata <- arrange(mapdata, CEL.LG, Estimated.Mb.Position)
mapdata$Diff <- c(diff(mapdata$Estimated.Mb.Position), 1e6)
mapdata$Diff[which(mapdata$Diff < 0)] <- 1e6
mapdata$Cumu <- c(0, cumsum(mapdata$Diff)[-nrow(mapdata)])

restab <- join(restab, mapdata)

write.table(restab, paste0("results/4_GWAS_results_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Create allele and Manhattan Plots           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Look at genotype effects

head(restab)

x <- data.frame(as.character.gwaa.data(abeldata[,unique(restab$SNP.Name)[1:10]]))
x$RRID <- row.names(x)
x <- melt(x, id.vars = "RRID")

y <- recsumm[,c("RRID", "RRID.Sex", "TotalRecombCount")]
y$RRID <- as.character(y$RRID)

x <- join(x, y)
x <- subset(x, !is.na(TotalRecombCount))

rm(y)

ggplot(x, aes(value, TotalRecombCount, fill = RRID.Sex)) +
  geom_boxplot(notch = T) + 
  facet_wrap(~variable, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Genotype", y = "Crossover Count")


#~~ Get Chromosome Positions

chrinfo <- NULL

for(i in na.omit(unique(mapdata$CEL.LG))){
  
  temp1 <- arrange(subset(mapdata, CEL.LG == i), Estimated.Mb.Position)
  
  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$CEL.LG)))
colourscale <- colourscale[1:(length(colourscale)/2)]

bonf = 0.05/35263.64

chrinfo$CEL.LG[35] <- "X"
chrinfo$CEL.LG2 <- c("", 1, "", "", 4, "", "", 7, "", "", 10, "", "", 13, "", "", 16, "", "", 19, "", "", 22, "", "", "", "", 27, "", "", "", 31, "", "", "X")


#~~ PLOT

png(paste0("figs/4_GWAS_Manhattan_Plot_", AnalysisSuffix, ".png"), width = 12, height = 12, units = "in", res = 300)
multiplot(
  
  ggplot(subset(restab, Phenotype == "acc"), aes(Cumu,-log10(Pc1df), col = factor(CEL.LG))) +
    geom_point(size = 2, alpha = 0.4) +
    geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
    scale_colour_manual(values = colourscale) +
    theme(legend.position="none") +
    theme(axis.text.x  = element_text (size = 16),
          axis.text.y  = element_text (size = 14),
          strip.text.x = element_text (size = 16),
          axis.title.y = element_text (size = 16, angle = 90),
          axis.title.x = element_text (size = 16),
          strip.background = element_blank()) +
    scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG2) +
    labs(x ="Chromosome", y = "-log10 P") +
    facet_wrap(~Sex, ncol = 1)
  ,
  ggplot(subset(restab, Phenotype == "acc"), aes(x= -log10(ExpP), y = -log10(Pc1df))) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          strip.background = element_blank()) +
    labs(x="Expected -log10 P",y="Observed -log10 P")+
    facet_wrap(~Sex, ncol = 1)
  , cols = 2, layout = matrix(c(1,1,2), 1))

dev.off()

#~~ Examine associated region in more detail

png(paste0("figs/4_GWAS_CEL12_Plot_", AnalysisSuffix, ".png"), width = 8, height = 12, units = "in", res = 300)

ggplot(subset(restab, Phenotype == "acc" & CEL.LG == 12), aes(Estimated.Mb.Position,-log10(P1df))) +
  geom_point(size = 2, alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.05) +
  geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P") +
  facet_wrap(~Sex, ncol = 1)

dev.off()

png(paste0("figs/4_GWAS_CEL12_Zoom_Plot_", AnalysisSuffix, ".png"), width = 8, height = 12, units = "in", res = 300)

ggplot(subset(restab, Phenotype == "acc" & 
                CEL.LG == 12 & 
                Estimated.Mb.Position > 15e6 & 
                Estimated.Mb.Position < 35e6), aes(Estimated.Mb.Position,-log10(P1df))) +
  geom_point(size = 2, alpha = 0.5) +
  geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P") +
  facet_wrap(~Sex, ncol = 1)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Examine variation in recombination rate in trans   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

runTransAnalysis <- FALSE

if(runTransAnalysis == TRUE){
  transres <- NULL
  
  for(lg in 1:34){
    
    print(paste("Running chromosome", lg, "of", 34))
    
    #~~ Create phenotypic dataset
    
    x <- subset(rectab, !CEL.LG %in% c(lg, 34))
    y <- data.frame(TransRecombCount = tapply(x$RecombCount, x$Family, sum))
    y$Family <- row.names(y)
    
    x <- join(recsumm, y)
    x.m <- subset(x, RRID.Sex == "M")
    x.f <- subset(x, RRID.Sex == "F")
    
    #~~ Create genabel subset
    tempmap <- subset(mapdata, CEL.LG == lg)
    
    tempabel <- abeldata[,as.character(tempmap$SNP.Name)]
    
    gwas.acc.all.trans <- rGLS(TransRecombCount ~ RRID.Sex + RRID.Fhat3, genabel.data = tempabel, phenotype.data = x, id.name = "RRID", GRM = grminv)
    gwas.acc.m.trans   <- rGLS(TransRecombCount ~ RRID.Fhat3, genabel.data = tempabel, phenotype.data = x.m, id.name = "RRID", GRM = grminv)
    gwas.acc.f.trans   <- rGLS(TransRecombCount ~ RRID.Fhat3, genabel.data = tempabel, phenotype.data = x.f, id.name = "RRID", GRM = grminv)
    
    head(results(gwas.acc.all.trans))
    
    res.vec <- ls()[grep("gwas", ls())]
    res.vec <- res.vec[grep("trans", res.vec)]
    
    for(i in res.vec){
      
      eval(parse(text = paste0("x <- results(", i, ")")))
      x$SNP.Name <- row.names(x)
      x$Analysis <- i
      transres <- rbind(transres, x)
      
      rm(x)
    }
    
    rm(gwas.acc.all.trans, gwas.acc.m.trans, gwas.acc.f.trans, tempabel, x, x.m, x.f, y)
    
  }
  
  write.table(transres, "temp.txt", row.names = F, sep = "\t", quote = F)
  
  
  #~~ Format table
  
  
  transres$Phenotype <- sapply(transres$Analysis, function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
  transres$Sex <- sapply(transres$Analysis, function(x) strsplit(x, split = ".", fixed = T)[[1]][3])
  
  head(transres)
  
  transres <- join(transres, mapdata)
  
  
  #~~ Deal with multiple testing
  
  transres <- arrange(transres, Analysis, P1df)
  transres$ExpP <- rep(seq(1/(nrow(transres)/3), 1, 1/(nrow(transres)/3)), times = 3)
  transres$chi2.1df <- qchisq(p = 1-transres$P1df, df =  1)
  transres$null.chi2.1df <- qchisq(p = 1-transres$ExpP, df =  1)
  
  #~~ Lambda is calculated as the median observed chi2 value divided by the median expected chi2 value
  
  lambdas <- tapply(transres$chi2.1df, transres$Analysis, median, na.rm = T)/
    tapply(transres$null.chi2.1df, transres$Analysis, median, na.rm = T)
  
  lambdas
  
  transres$chi2.1df.adj <- ifelse(transres$Analysis == names(lambdas)[1], transres$chi2.1df/lambdas[1],
                                  ifelse(transres$Analysis == names(lambdas)[2], transres$chi2.1df/lambdas[2],
                                         ifelse(transres$Analysis == names(lambdas)[3], transres$chi2.1df/lambdas[3],NA)))
  
  transres$Pc1df <- 1-pchisq(transres$chi2.1df.adj, 1)
  
  transres$null.chi2.1df <- NULL
  
  write.table(transres, paste0("results/4_GWAS_results_trans_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
} else {
  transres <- read.table( paste0("results/4_GWAS_results_trans_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F, sep = "\t")
}
#~~ Look at the top hits

transres <- arrange(transres, Pc1df)
transres[1:10,]


#~~ Make the plot!

png(paste0("figs/4_GWAS_Manhattan_TRANS_", AnalysisSuffix, ".png"), width = 12, height = 12, units = "in", res = 300)

multiplot(
  
  ggplot(subset(transres, Phenotype == "acc"), aes(Cumu,-log10(Pc1df), col = factor(CEL.LG))) +
    geom_point(size = 2, alpha = 0.4) +
    geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
    scale_colour_manual(values = colourscale) +
    theme(legend.position="none") +
    theme(axis.text.x  = element_text (size = 16),
          axis.text.y  = element_text (size = 14),
          strip.text.x = element_text (size = 16),
          axis.title.y = element_text (size = 16, angle = 90),
          axis.title.x = element_text (size = 16),
          strip.background = element_blank()) +
    scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG2) +
    labs(x ="Chromosome", y = "-log10 P") +
    facet_wrap(~Sex, ncol = 1)
  ,
  ggplot(subset(transres, Phenotype == "acc"), aes(x= -log10(ExpP), y = -log10(Pc1df))) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          strip.background = element_blank()) +
    labs(x="Expected -log10 P",y="Observed -log10 P")+
    facet_wrap(~Sex, ncol = 1)
  , cols = 2, layout = matrix(c(1,1,2), 1))

dev.off()



#~~ Examine associated region in more detail

png(paste0("figs/4_GWAS_CEL12_Plot_TRANS_", AnalysisSuffix, ".png"), width = 8, height = 12, units = "in", res = 300)

ggplot(subset(transres, Phenotype == "acc" & CEL.LG == 12), aes(Estimated.Mb.Position,-log10(P1df))) +
  geom_point(size = 2, alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.05) +
  geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P") +
  facet_wrap(~Sex, ncol = 1)

dev.off()

png(paste0("figs/4_GWAS_CEL12_Zoom_Plot_TRANS_", AnalysisSuffix, ".png"), width = 8, height = 12, units = "in", res = 300)

ggplot(subset(transres, Phenotype == "acc" & 
                CEL.LG == 12 & 
                Estimated.Mb.Position > 15e6 & 
                Estimated.Mb.Position < 35e6), aes(Estimated.Mb.Position,-log10(P1df))) +
  geom_point(size = 2, alpha = 0.5) +
  geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P") +
  facet_wrap(~Sex, ncol = 1)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Chromosome partitioning           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm$RRID2 <- recsumm$RRID
recsumm.f$RRID2 <- recsumm.f$RRID
recsumm.m$RRID2 <- recsumm.m$RRID


runChrPartitionPrep <- FALSE

if(runChrPartitionPrep == TRUE){
  setwd(paste0("gcta/chr_part_", AnalysisSuffix))
  
  
    for(i in 1:33){
    
    print(paste("Running chromosome", i)) 
    snplist <- subset(mapdata, CEL.LG == i)$SNP.Name
    write.table(snplist, paste0("snplist.", i, ".txt"), row.names = F, col.names = F, quote = F)
    
    snplist <- subset(mapdata, !CEL.LG %in% c(i, 34))$SNP.Name
    write.table(snplist, paste0("snplist.wo.", i, ".txt"), row.names = F, col.names = F, quote = F)
    
    
    system(paste0("..\\gcta.exe --bfile ..\\Deer31.recoded.v2 --autosome-num 33 --keep ..\\idlist.txt --extract snplist.",i, ".txt  --make-grm --out chr", i, "_GRM"))
    system(paste0("..\\gcta.exe --grm chr", i, "_GRM --grm-adj 0 --make-grm --out chr", i, "_GRM_adj"))
    system(paste0("..\\gcta.exe --bfile ..\\Deer31.recoded.v2 --autosome-num 33 --keep ..\\idlist.txt  --extract snplist.wo.", i,".txt --make-grm --out chr", i, "_wo_GRM"))
    system(paste0("..\\gcta.exe --grm chr", i, "_wo_GRM --grm-adj 0 --make-grm --out chr", i, "_wo_GRM_adj"))
    
    
    grm.rest <- read.table(paste0("chr", i,"_wo_GRM_adj.grm.gz"))
    ids.rest <- read.table(paste0("chr",i,"_wo_GRM_adj.grm.id"))
    
    grm.reg <- read.table(paste0("chr",i,"_GRM_adj.grm.gz"))
    ids.reg <- read.table(paste0("chr",i,"_GRM_adj.grm.id"))
    
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
    
    
    
    save(x, file = paste0("chr", i,".RData"))
    rm(model1, model2, model1.m, model2.m, model1.f, model2.f, x)
  }
setwd("../..")
}

# 
# chr.res <- NULL
# 
# for(i in 1:33){
#   
#   load(paste0("gcta/chr_part_", AnalysisSuffix, "/chr", i, ".RData"), verbose = T)
#   
#   chr.res <- rbind(chr.res,
#                    data.frame(Analysis.ID = i,
#                               RRID.Sex     = c("Both", "Male", "Female"),
#                               model1.Li    = c(x[[1]][2], x[[3]][2], x[[5]][2]),
#                               model2.Li    = c(x[[2]][2]$loglik, x[[4]][2]$loglik, x[[6]][2]$loglik),
#                               model1.Va2   = c(x[[1]]$vcoef[[1]], x[[3]]$vcoef[[1]], x[[5]]$vcoef[[1]]),
#                               model2.Va2   = c(x[[2]]$vcoef[[1]][1], x[[4]]$vcoef[[1]][1], x[[6]]$vcoef[[1]][1]),
#                               model2.Vreg2 = c(x[[2]]$vcoef[[1]][2], x[[4]]$vcoef[[1]][2], x[[6]]$vcoef[[1]][2]),
#                               model2.Vr = c(x[[2]]$sigma^2, x[[4]]$sigma^2, x[[6]]$sigma^2)))
#   
#   rm(x)
#   
# }
# 
# chr.res$LogLi.P <- 1- pchisq(2*(chr.res$model1.Li - chr.res$model2.Li), df = 1)







