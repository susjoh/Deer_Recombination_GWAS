#
# Create phenotype files & GRMs for use in analysis
# Susan Johnston
# April 2017
#
#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

auto.vec <- 1:33
sex.vec <- 34
chr.vec <- 1:34
acro.vec <- c(1:4, 6:33)

runGCTA = FALSE

#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

#~~ Read in pedigree files

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)
famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ read in SNP positions and PAR SNPs

chrom.map <- read.table(paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- subset(mapdata, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position))

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Load window variation from previous paper

bin.tab <- read.table("data/TableS4_Recombination_Landscape_Info.txt", header = T, stringsAsFactors = F)

#~~ Load per chromosome rectab information

rectab <- read.table(paste0("results/2_Per_Chromosome_Recomb_Final_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ Read in individual deer data

phenotab <- read.table("data/20160725_Deer_BasicData.txt", header = T, sep = "\t")
head(phenotab)

annfit <- read.table("data/AnnFitData_Recoded.txt", header = T)
head(annfit)

recoded.ids <- read.table("data/DeerRecodedIDs.txt", header= T)
names(recoded.ids) <- c("Code", "ID")

phenotab <- join(phenotab, recoded.ids)

phenotab$Sex <- ifelse(phenotab$Sex == 1, "F",
                       ifelse(phenotab$Sex == 2, "M", "Unk"))

rm(recoded.ids)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Bring together basic phenotype information    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Remove the sex chromosomes

rectab <- subset(rectab, CEL.LG != sex.vec)

#~~ Create the basic recombination phenotypes

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalChrCount    = tapply(rectab$RecombCount, rectab$Family, length),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

temp.tab <- unique(subset(rectab, select = c(Family, RRID, ANIMAL)))

recsumm$Family <- row.names(recsumm)                    
recsumm$RRID.Sex    <- NA
recsumm$RRID.Sex[grep("FATHER", recsumm$Family)] <- "M"
recsumm$RRID.Sex[grep("MOTHER", recsumm$Family)] <- "F"
recsumm <- join(recsumm, temp.tab)

rm(temp.tab)

head(recsumm) #1345

#~~ Plot

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)
ggplot(recsumm, aes(TotalInfLoci, TotalRecombCount)) + geom_point(alpha = 0.2) + stat_smooth(method = "lm")

#~~ Add RRID basic phenotypic information to recsumm

head(phenotab)

x <- subset(phenotab, select = c(ID, BirthMonth, BirthYear, CaptureWt, HindLeg))
names(x)[2:5] <- paste0("RRID.", names(x)[2:5])
names(x)[1] <- "RRID"

recsumm <- join(recsumm, x)
head(recsumm)

#~~ Add ANIMAL basic phenotypic information to recsumm

head(phenotab)

x <- subset(phenotab, select = c(ID, Sex, BirthMonth, BirthYear, CaptureWt, HindLeg))
names(x)[2:6] <- paste0("ANIMAL.", names(x)[2:6])
names(x)[1] <- "ANIMAL"

recsumm <- join(recsumm, x)
head(recsumm)

#~~ Add further phenotypes

recsumm$RRID.Age <- recsumm$ANIMAL.BirthYear - recsumm$RRID.BirthYear
hist(recsumm$RRID.Age)


#~~ Calculate Inbreeding coefficients and make GRM while you are at it

if(runGCTA == TRUE){
  
  setwd("gcta")
  
  system("cmd", input = "plink --file Deer31.recoded.v2 --make-bed --chr-set 33 no-y no-xy no-mt --out Deer31.recoded.v2")
  
  system("cmd", input = "gcta.exe --bfile Deer31.recoded.v2  --autosome  --ibc --autosome-num  33 --out Deer31.recoded.v2")
  
  system("cmd", input = "gcta.exe --bfile Deer31.recoded.v2  --autosome  --make-grm --autosome-num  33 --out Deer31.recoded.v2")
  
  setwd("..")
  
}

deer.ibc <- read.table("gcta/Deer31.recoded.v2.ibc", header = T)

deer.ibc <- subset(deer.ibc, select = c(IID, Fhat3))
names(deer.ibc) <- c("RRID", "RRID.Fhat3")

recsumm <- join(recsumm, deer.ibc)

names(deer.ibc) <- c("ANIMAL", "ANIMAL.Fhat3")

recsumm <- join(recsumm, deer.ibc)

head(recsumm)

#~~ Add parent information

head(pedigree)
names(pedigree) <- c("ANIMAL", "ANIMAL.MOTHER", "ANIMAL.FATHER")

recsumm <- join(recsumm, pedigree)

names(pedigree) <- c("RRID", "RRID.MOTHER", "RRID.FATHER")

recsumm <- join(recsumm, pedigree)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Quantifying variation in recombination rate   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create temporary data frame to concatenate everything

rectemp <- data.frame(Family = recsumm$Family)


#~~ Proportion of chromosome from the centromere on acrocentric chromosomes

prop.vec <- c(0.2, 0.4, 0.6, 0.8)

for (i in prop.vec){
  print(paste("Running proportion", i))
  
  y <- data.frame(RR = sapply(rectab$data, function (x) {
    x <- strsplit(x, split = "")[[1]]
    x <- x[1:(length(x)*i)]
    x <- x[which(x %in% c(0, 1))]
    
    x <- length(which(c(-999, x) != c(x, -999))) -2
    if(x < 0) x <- NA
    x
  }
  ),
  Family = rectab$Family,
  CEL.LG = rectab$CEL.LG)
  
  head(y)
  
  y <- subset(y,  CEL.LG %in% acro.vec)
  
  x <- data.frame(RR = tapply(y$RR, y$Family, sum, na.rm = T))
  x$Family <- row.names(x)
  head(x)
  
  names(x)[1] <- paste0("AcroCentroPC.", i*100)
  head(x)
  
  rectemp <- join(rectemp, x)
  
  rm(x, y)
}


#~~ Examine the variation in recombination landscape

ggplot(bin.tab, aes(Window, cM.Female)) + stat_smooth(method = "loess", span = 0.15)
ggplot(bin.tab, aes(Window, cM.Female, group = CEL.LG)) + stat_smooth(method = "loess", span = 0.15, se = F)


#~~ Extract recombination counts for the first 20Mb of the chromosome in acros

head(mapdata)

chrom.map <- join(chrom.map, mapdata)
chrom.map.temp <- subset(chrom.map, Estimated.Mb.Position < 20e6)

max.temp <- data.frame(Order.20 = tapply(chrom.map.temp$Order, chrom.map.temp$CEL.LG, max))
max.temp$CEL.LG <- row.names(max.temp)

rectab <- join(rectab, max.temp)

foo <- function (x, y) {
  x <- strsplit(x, split = "")[[1]]
  x <- x[1:y]
  x <- x[which(x %in% c(0, 1))]
  x <- length(which(c(-999, x) != c(x, -999))) -2
  if(x < 0) x <- NA
  x
  }

y <- data.frame(RR = mapply(foo, rectab$data, rectab$Order.20),
                Family = rectab$Family,
                CEL.LG = rectab$CEL.LG)

y <- subset(y,  CEL.LG %in% acro.vec)

x <- data.frame(RR = tapply(y$RR, y$Family, sum, na.rm = T))
x$Family <- row.names(x)
head(x)

names(x)[1] <- paste0("AcroCentroMb.20")
head(x)

rectemp <- join(rectemp, x)

rm(x, y, i, foo)


head(recsumm)
recsumm <- join(recsumm, rectemp)




#~~ Distance from centromere


doub.xovers <- NULL


for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  x <- parse_crossovers(chrompicfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp"),
                        familyPedigree = famped)

  doub.xovers <- rbind(doub.xovers, check_double_crossovers(x))
  
  rm(x)
  
}

head(doub.xovers)

doub.xovers.hold <- doub.xovers

doub.xovers <- doub.xovers.hold

#~~ Keep acrocentric autosomes

doub.xovers <- subset(doub.xovers, !analysisID %in% c("5a", "34a"))

#~~ Table the First segements, and get the position mid-way between the start and stop span.


doub.xovers <- subset(doub.xovers, Family %in% recsumm$Family)
doub.xovers <- subset(doub.xovers, Type == "First")

doub.xovers$Order <- round(apply(doub.xovers[,c("StopPos", "StopSpan")], 1, mean), digits = 0)
doub.xovers$CEL.LG <- gsub("a", "", doub.xovers$analysisID)

doub.xovers <- join(doub.xovers, subset(chrom.map, select = c(CEL.LG, Order, Estimated.Mb.Position)))

centrotab <- data.frame(Mean.CentroToCO = tapply(doub.xovers$Estimated.Mb.Position, doub.xovers$Family, mean),
                        Median.CentroToCO = tapply(doub.xovers$Estimated.Mb.Position, doub.xovers$Family, median))
centrotab$Family <- row.names(centrotab)
head(centrotab)

recsumm <- join(recsumm, centrotab)

write.table(recsumm, paste0("results/3_Recombination_Phenotype_Data_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)



