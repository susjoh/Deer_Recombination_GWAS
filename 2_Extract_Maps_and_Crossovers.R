#
# Extract linkage maps and crossover information
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

auto.vec <- 1:33
sex.vec <- 34
chr.vec <- 1:34

#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

firstRun <- FALSE

#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Read in pedigree files

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)
famped <- read.table(paste0("results/1_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ read in SNP positions and PAR SNPs

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)
mapdata <- arrange(mapdata, CEL.LG, CEL.order)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create CRI-Map Files and run maps/chrompic  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun == TRUE){
  
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
                          clear.existing.analysisID = TRUE,
                          use.mnd = TRUE)
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
      run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
    }
    
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse and check sex-averaged maps           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chrom.map <- NULL

for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  chrom.map <- rbind(chrom.map,
                     parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp")))
  
  
}

chrom.map$CEL.LG <- as.numeric(gsub(AnalysisSuffix, "", chrom.map$analysisID))

maxvals <- data.frame(CEL.LG = chr.vec,
                      Estimated.Length.Mb = tapply(mapdata$Estimated.Mb.Position, mapdata$CEL.LG, max),
                      max.cM =              tapply(chrom.map$cMPosition, chrom.map$CEL.LG, max))

ggplot(chrom.map, aes(Order, cMPosition)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free")

ggplot(maxvals, aes(Estimated.Length.Mb, max.cM)) +
  geom_text(aes(label = CEL.LG)) +
  stat_smooth(method = "lm")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Parse crossovers and double-recombinants       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- NULL
doub.xovers <- NULL


for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  x <- parse_crossovers(chrompicfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp"),
                        familyPedigree = famped)
  rectab <- rbind(rectab, x)
  
  doub.xovers <- rbind(doub.xovers, check_double_crossovers(x))
  
  rm(x)
  
}

#~~ Merge with the span positions

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)

write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_raw_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_raw_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_raw_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. QC and deal with singleton double xovers        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(doub.xovers)
head(chrom.map)

#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)

#~~ SAVE FIGURE

ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + 
  geom_histogram(binwidth = 1) + 
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Span Distance (cM)",
       y = "Count")
ggsave("figs/2_Double_Crossover_Histogram.png", width = 6, height = 4.5)


#~~ Remove all singletons phased runs shorter than 10cM

bad.doubs <- rbind(subset(doub.xovers, Singleton == "yes"))
head(bad.doubs)

x <- subset(chrom.map, select = c(analysisID, Order, SNP.Name))
names(x) <- c("analysisID", "StartPos", "SNP.Name")

bad.doubs <- join(bad.doubs, x)

head(bad.doubs)
bad.doubs <- subset(bad.doubs, select = c(SNP.Name, Family))
bad.doubs$RRID <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][4])
bad.doubs$ANIMAL <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][2])
bad.doubs$Family <- NULL

bad.doubs <- melt(bad.doubs, id.vars = "SNP.Name")
bad.doubs$variable <- NULL
names(bad.doubs) <- c("SNP.Name", "ANIMAL")

#~~ Make a master .mnd file and add bad doubles (singletons)

mnd.master <- NULL

for(i in chr.vec){
  mnd.file <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mnd"), header = T)
  mnd.master <- rbind(mnd.master, mnd.file)
  rm(mnd.file)
}

mnd.master <- rbind(mnd.master, bad.doubs)
mnd.master <- unique(mnd.master)

write.table(mnd.master, paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Rerun without Singletons                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun == TRUE){
  
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
                          clear.existing.analysisID = TRUE,
                          use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"))
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
      run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
    }
    
  })
}



chrom.map <- NULL

for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  chrom.map <- rbind(chrom.map,
                     parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp")))
  
  
}

chrom.map$CEL.LG <- as.numeric(gsub(AnalysisSuffix, "", chrom.map$analysisID))

maxvals <- data.frame(CEL.LG = chr.vec,
                      Estimated.Length.Mb = tapply(mapdata$Estimated.Mb.Position, mapdata$CEL.LG, max),
                      max.cM =              tapply(chrom.map$cMPosition, chrom.map$CEL.LG, max))

ggplot(chrom.map, aes(Order, cMPosition)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free")

ggplot(maxvals, aes(Estimated.Length.Mb, max.cM)) +
  geom_text(aes(label = CEL.LG)) +
  stat_smooth(method = "lm")


rectab <- NULL
doub.xovers <- NULL


for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  x <- parse_crossovers(chrompicfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp"),
                        familyPedigree = famped)
  rectab <- rbind(rectab, x)
  
  doub.xovers <- rbind(doub.xovers, check_double_crossovers(x))
  
  rm(x)
  
}

#~~ Merge with the span positions

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)

write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_QC1_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_QC1_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_QC1_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. QC and deal with singleton double xovers        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(doub.xovers)
head(chrom.map)

#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)

#~~ SAVE FIGURE

ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + 
  geom_histogram(binwidth = 1) + 
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Span Distance (cM)",
       y = "Count")
ggsave("figs/2_Double_Crossover_Histogram_QC1.png", width = 6, height = 4.5)


#~~ Remove all singletons phased runs shorter than 10cM

bad.doubs <- rbind(subset(doub.xovers, Singleton == "yes"),
                   subset(doub.xovers, SpancMDist < 10))
head(bad.doubs)

x <- subset(chrom.map, select = c(analysisID, Order, SNP.Name))
names(x) <- c("analysisID", "StartPos", "SNP.Name")

bad.doubs <- join(bad.doubs, x)

head(bad.doubs)
bad.doubs <- subset(bad.doubs, select = c(SNP.Name, Family))
bad.doubs$RRID <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][4])
bad.doubs$ANIMAL <- sapply(bad.doubs$Family, function(x) strsplit(x, split = "_")[[1]][2])
bad.doubs$Family <- NULL

bad.doubs <- melt(bad.doubs, id.vars = "SNP.Name")
bad.doubs$variable <- NULL
names(bad.doubs) <- c("SNP.Name", "ANIMAL")


mnd.master <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), header = T, stringsAsFactors = F)
mnd.master <- unique(rbind(mnd.master, doub.xovers))


write.table(mnd.master, paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"), row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Rerun without Singletons and short segments, FINAL  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun == TRUE){
  
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
                          clear.existing.analysisID = TRUE,
                          use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr_", AnalysisSuffix, "_MASTER.mnd"))
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
      run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
      
    }
    
  })
}



chrom.map <- NULL

for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  chrom.map <- rbind(chrom.map,
                     parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp")))
  
  
}

chrom.map$CEL.LG <- as.numeric(gsub(AnalysisSuffix, "", chrom.map$analysisID))

maxvals <- data.frame(CEL.LG = chr.vec,
                      Estimated.Length.Mb = tapply(mapdata$Estimated.Mb.Position, mapdata$CEL.LG, max),
                      max.cM =              tapply(chrom.map$cMPosition, chrom.map$CEL.LG, max))

ggplot(chrom.map, aes(Order, cMPosition)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free")

ggplot(maxvals, aes(Estimated.Length.Mb, max.cM)) +
  geom_text(aes(label = CEL.LG)) +
  stat_smooth(method = "lm")


rectab <- NULL
doub.xovers <- NULL


for(i in chr.vec){
  
  print(paste("Running chromosome", i, "of", length(chr.vec)))
  
  x <- parse_crossovers(chrompicfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp"),
                        familyPedigree = famped)
  rectab <- rbind(rectab, x)
  
  doub.xovers <- rbind(doub.xovers, check_double_crossovers(x))
  
  rm(x)
  
}

#~~ Merge with the span positions

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StartSpan", "StartSpan.cM")

doub.xovers <- join(doub.xovers, x)

x <- subset(chrom.map, select = c(analysisID, Order, cMPosition))
names(x) <- c("analysisID", "StopSpan", "StopSpan.cM")

doub.xovers <- join(doub.xovers, x)

rm(x)

write.table(rectab     , paste0("results/2_Per_Chromosome_Recomb_QC2_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, paste0("results/2_Double_Xovers_QC2_", AnalysisSuffix, ".txt")        , row.names = F, sep = "\t", quote = F)
write.table(chrom.map  , paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt")               , row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. Extract and format crossover data                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab      <- read.table(paste0("results/2_Per_Chromosome_Recomb_QC2_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)
doub.xovers <- read.table(paste0("results/2_Double_Xovers_QC2_", AnalysisSuffix, ".txt")         , header = T, stringsAsFactors = F)
chrom.map <- read.table(paste0("results/2_cM_Map_QC2_", AnalysisSuffix, ".txt")               , header = T, stringsAsFactors = F)



#~~ Are there particular problematic individuals?

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, Family %in% row.names(recsumm))
doub.xovers <- subset(doub.xovers, Family %in% row.names(recsumm))

#~~ Determine span distance

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers$SpancMDist    <- doub.xovers$StopSpan.cM - doub.xovers$StartSpan.cM

doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpancMDist, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpancMDist)) + geom_point(alpha = 0.1)



#~~ Create new famped!

famped <- subset(famped, Family %in% rectab$Family)


#~~ remove the remaining short crossovers directly from the rectab object

doub.xovers <- subset(doub.xovers, SpancMDist < 10)
head(doub.xovers)

for(i in 1:nrow(doub.xovers)){
  
  row.n <- which(rectab$UniqueID == doub.xovers$UniqueID[i])
  if(length(row.n) > 0){
    x <- rectab[row.n,"data"]
    x <- strsplit(x, split = "")[[1]]
    x[doub.xovers$StartPos[i]:doub.xovers$StopPos[i]] <- "-"
    
    #~~ revise crossover number
    
    x1 <- x[which(x %in% c(0, 1))]
    rectab$RecombCount[row.n] <- length(which(c(-999, x1) != c(x1, -999)))-2
    
    rectab$data[row.n] <- paste(x, collapse = "")
    
    rm(x, x1, row.n)
  }
}

head(rectab)

rectab$No.Inf.Loci <- sapply(rectab$data, function(x) nchar(gsub("-", "", x)))

rectab$CEL.LG <- gsub(AnalysisSuffix, "", rectab$analysisID)


write.table(rectab, paste0("results/2_Per_Chromosome_Recomb_Final_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)
write.table(famped, paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)


