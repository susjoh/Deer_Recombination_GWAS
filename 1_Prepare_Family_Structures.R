#
# Prepare family structures and run Cri-MAP to extract crossover positions
# Susan Johnston
# April 2017
#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

auto.vec <- 1:33    # 1:29
sex.vec <- 34
chr.vec <- 1:34


#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Read in pedigree file

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)

#~~ read in SNP positions and PAR SNPs

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)
mapdata <- arrange(mapdata, CEL.LG, CEL.order)

head(mapdata)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Define thresholds for pedigree construction

mend.err.locus.threshold <- 0.01
mend.err.pair.threshold  <- 0.001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Determine Working Pedigrees for Crimap    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove non-genotyped parents from the pedigree

pedigree    <- subset(pedigree, ANIMAL %in% idnames(abeldata))
pedigree$MOTHER[which(!pedigree$MOTHER %in% idnames(abeldata))] <- 0
pedigree$FATHER[which(!pedigree$FATHER %in% idnames(abeldata))] <- 0

head(pedigree)

#~~ Determine all P-O pairs where the parent has both parents and mate are genotyped.

po.pairs <- melt(pedigree, id.vars = "ANIMAL")
po.pairs <- subset(po.pairs, value != 0)

head(po.pairs)
po.pairs$Selected <- NA

famped <- NULL

for(i in 1:nrow(po.pairs)){

  if(i %in% seq(1, nrow(po.pairs), 500)) print(paste("Running line", i, "of", nrow(po.pairs)))
  
  #~~ Extract the focal parent
  
  x1 <- pedigree[which(pedigree$ANIMAL == po.pairs$value[i]),]
  
  #~~ Extract the offspring
  
  x2 <- pedigree[which(pedigree$ANIMAL == po.pairs$ANIMAL[i]),]
  
  #~~ Extract focal parent's parents and mate
  
  x3 <- data.frame(ANIMAL = c(unlist(x1[1,2:3]), ifelse(po.pairs$variable[i] == "MOTHER", x2$FATHER[1], x2$MOTHER[1])), 
                   MOTHER = 0, FATHER = 0)
  
  x4 <- rbind(x3, x1, x2)
  x4$Family <- paste("Offspring", po.pairs$ANIMAL[i], po.pairs$variable[i], po.pairs$value[i], sep = "_")
  
  x4 <- subset(x4, ANIMAL != 0)
  
  if(nrow(x4) == 5){
    po.pairs$Selected[i] <- "yes"
    famped <- rbind(famped, x4)
  } else {
    po.pairs$Selected[i] <- "no"
  }
  
  rm(x1, x2, x3, x4)
  
}

head(famped)
rm(i)

nrow(famped)/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create crimap files & Run 1st instance    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

AnalysisSuffix <- "a"

system.time({
  
  for(i in chr.vec){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    create_crimap_input(gwaa.data = abeldata, 
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        snplist = subset(mapdata, CEL.LG == i)$SNP.Name, 
                        is.X = ifelse(i == sex.vec, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   is.X = ifelse(i == sex.vec, TRUE, FALSE),
                   pseudoautoSNPs = pseudoautoSNPs,
                   genabel.phdata = phdata(abeldata),
                   save.mendfile = TRUE)
    
  }
  
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Deal with Mendelian Errors                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


menderr <- NULL

for(i in chr.vec){
  
  temperr <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mndverbose"),
                        header = T, sep = "\t", stringsAsFactors = F) 
  
  if(i == sex.vec) temperr <- subset(temperr, select = -sex)
  
  menderr <- rbind(menderr,
                   cbind(temperr, analysisID = paste0(i, AnalysisSuffix)))
  
  rm(temperr)
  
}

#~~ Do particular pairs perform badly?

menderr <- subset(menderr, SNP.Name %in% snpnames(abeldata))

poor.ids <- subset(menderr, select = c(ANIMAL, MOTHER, Mat.Mismatch, Family))
poor.ids2 <- subset(menderr, select = c(ANIMAL, FATHER, Pat.Mismatch, Family))

names(poor.ids) <- c("ANIMAL", "PARENT", "Mismatch", "Family")
names(poor.ids2) <- c("ANIMAL", "PARENT", "Mismatch", "Family")

poor.ids$PARENT.Type <- "MOTHER"
poor.ids2$PARENT.Type <- "FATHER"

poor.ids <- rbind(poor.ids, poor.ids2)

rm(poor.ids2)

poor.ids <- droplevels(subset(poor.ids, Mismatch == "yes"))
poor.ids$ID.Parent.Family <- paste(poor.ids$ANIMAL, poor.ids$PARENT, poor.ids$Family, sep = "_")

poor.ids <- data.frame(table(poor.ids$ID.Parent))
names(poor.ids)[1] <- "ID.Parent"

ggplot(poor.ids, aes(Freq)) + geom_histogram(binwidth = 1, col = "grey")

table(poor.ids$Freq)

poor.ids <- poor.ids[which(poor.ids$Freq > 50),]

#~~ which families have these bad ids?

poor.ids$Family <- sapply(as.character(poor.ids$ID.Parent), function (x) paste(strsplit(x, split = "_")[[1]][3:5], collapse = "_"))

famped <- subset(famped, !Family %in% poor.ids$Family) #18 families removed

write.table(poor.ids, paste0("results/1_ID_pairs_with_high_Mend_errors_", AnalysisSuffix, ".txt"),
            row.names = F, quote = F)


write.table(famped, 
            paste("results/1_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)









