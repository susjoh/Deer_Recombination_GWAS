#
# SHAPEIT Phasing
# Susan Johnston
# September 2017
#
#

library(ggplot2)
library(plyr)
library(dplyr)
library(GenABEL)
library(reshape)

source("r/Shapeit_Functions.R")
source("r/Beagle_Functions.R")


#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

#~~ Define the window sizes and focal SNPs


snp.window.vec <- c(2, 3, 4, 5, 7, 10, 20, 50)
focal.snp.vec <- c("cela1_red_10_26005249", "cela1_red_10_25661750", "cela1_red_10_21092495")

#~~ Read in map information

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- subset(mapdata, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position, CEL.order))

FirstRun <- F

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Prepare input files              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(FirstRun == T){
  
  #~~ Get rid of IDs with low genotyping success and create new object for chr 12
  
  system("cmd", input = paste0("plink --bfile gcta/Deer31.recoded.v2 --recode --chr 12 --chr-set 33 no-y no-xy no-mt --mind 0.05 --make-bed --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  
  #~~ Recode the families
  
  tempfam <- read.table("shapeit/Deer31.recoded.v2.chr12.QCed.fam")
  head(tempfam)
  tempfam <- tempfam[,1:2]
  tempfam$v3 <- 1
  tempfam$V4 <- tempfam$V2
  
  write.table(tempfam, "shapeit/PlinkFamRecode.txt", row.names = F, col.names = F, quote = F)
  
  rm(tempfam)
  
  system("cmd", input = paste0("plink --bfile shapeit/Deer31.recoded.v2.chr12.QCed --recode --chr-set 33 no-y no-xy no-mt --update-ids shapeit/PlinkFamRecode.txt --make-bed --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  
  #~~ Read Pedigree file, format and write to beagle directory
  
  pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header= T)
  head(pedigree)
  pedigree <- cbind(Family = 1, pedigree)
  pedigree <- pedigree[,c(1, 2, 4, 3)]
  
  pedigree <- subset(pedigree, ANIMAL %in% tempfam$V2)
  
  pedigree$FATHER[which(!pedigree$FATHER %in% tempfam$V2)] <- 0
  pedigree$MOTHER[which(!pedigree$MOTHER %in% tempfam$V2)] <- 0
  
  # pedigree$FATHER[which(pedigree$FATHER != 0 & pedigree$MOTHER == 0)] <- 0
  # pedigree$MOTHER[which(pedigree$MOTHER != 0 & pedigree$FATHER == 0)] <- 0
  
  write.table(pedigree, "shapeit/DeerPedigree.txt", row.names = F, col.names = F, quote = F, sep = "\t")
  
  #~~ Deal with sex issues
  
  tempfam <- read.table("shapeit/Deer31.recoded.v2.chr12.QCed.fam")
  head(tempfam)
  
  tempfam$Parent.Status <- ifelse(tempfam$V2 %in% pedigree$FATHER, "FATHER", ifelse(tempfam$V2 %in% pedigree$MOTHER, "MOTHER", NA))
  table(tempfam$Parent.Status, tempfam$V5)
  
  tempfam <- subset(tempfam, !is.na(Parent.Status) & V5 == 0)
  tempfam$V5 <- ifelse(tempfam$Parent.Status == "FATHER", 1, 2)
  
  tempfam <- tempfam[,c(1, 2, 5)]
  
  write.table(tempfam, "shapeit/SexRecode.txt", row.names = F, col.names = F, quote = F, sep = "\t")
  
  rm(tempfam)
  
  
  system("cmd", input = paste0("plink --bfile shapeit/Deer31.recoded.v2.chr12.QCed --recode --update-sex shapeit/SexRecode.txt --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  system("cmd", input = paste0("plink --file shapeit/Deer31.recoded.v2.chr12.QCed --recode --make-bed --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  
  
  #~~ Add pedigree information and make bed
  
  system("cmd", input = paste0("plink --bfile shapeit/Deer31.recoded.v2.chr12.QCed --recode --update-parents shapeit/DeerPedigree.txt --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  system("cmd", input = paste0("plink --file shapeit/Deer31.recoded.v2.chr12.QCed --recode --make-bed --out shapeit/Deer31.recoded.v2.chr12.QCed"))
  
  

  
  #~~ tidy up files
  
  dir("shapeit")
  system("rm shapeit/*~")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Prepare files for SHAPEIT                   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  writeLines("#!/bin/bash", con = "shapeit/Run_Deer_Haplotyping.sh")
  
  for(i in 1:length(snp.window.vec)){
    for(j in 1:length(focal.snp.vec)){
      
      #~~ take the top SNP and SNPs either side.
      
      snp.window  <- snp.window.vec[i]
      focal.snp   <- focal.snp.vec[j]
      window.snps <- as.character(mapdata[(which(mapdata$SNP.Name == focal.snp) - snp.window):(which(mapdata$SNP.Name == focal.snp) + snp.window),"SNP.Name"])
      
      out.prefix <- paste0("Run_", focal.snp, "_span_", snp.window, "_ped")
      
      writeLines(window.snps, paste0("shapeit/", out.prefix, ".snplist"))
      
      #~~ Subset the PLINK data
      
      system("cmd", input = paste0("plink --bfile shapeit/Deer31.recoded.v2.chr12.QCed --recode --chr-set 33 no-y no-xy no-mt --extract shapeit/", out.prefix, ".snplist --make-bed --out shapeit/", out.prefix))
      
      #~~ Make a genetic map (currently flawed as don't have population estimates)
      
      x <- read.table(paste0("shapeit/",out.prefix, ".map"))
      x$V5 <- c(diff(x$V3), 0)
      x$V5[which(x$V5 == 0)] <- 0.001
      x <- x[,c(4, 5, 3)]
      names(x) <- c("pposition", "rrate", "gposition")
      
      write.table(x, paste0("shapeit/", out.prefix, ".genmap"), row.names = F, quote = F)
      
      write.table(paste0("./shapeit --input-bed ", out.prefix, ".bed ", out.prefix, ".bim ", out.prefix, ".fam --input-map ", 
                         out.prefix, ".genmap --output-max ", out.prefix, ".haps ", out.prefix, ".sample --duohmm"),
                  "shapeit/Run_Deer_Haplotyping.sh", row.names = F, col.names = F, quote = F, append = T)
      
      
      write.table(paste0("./shapeit -convert --input-haps ", out.prefix, " --output-vcf ", out.prefix, ".vcf"),
                  "shapeit/Run_Deer_Haplotyping.sh", row.names = F, col.names = F, quote = F, append = T)
      
      
      
    }
  }
  
  
  #~~ Write the whole chromosome to file

  window.snps <- as.character(mapdata[which(mapdata$CEL.LG == 12),"SNP.Name"])
  window.snps <- window.snps[-which(window.snps == "cela1_red_10_57228626")]  
  
  out.prefix <- paste0("Run_Chr_12_ped")
  
  writeLines(window.snps, paste0("shapeit/", out.prefix, ".snplist"))
  

  #~~ Subset the PLINK data
  
  system("cmd", input = paste0("plink --bfile shapeit/Deer31.recoded.v2.chr12.QCed --recode --chr-set 33 no-y no-xy no-mt --extract shapeit/", out.prefix, ".snplist --make-bed --out shapeit/", out.prefix))
  
  #~~ Make a genetic map (currently flawed as don't have population estimates)
  
  x <- read.table(paste0("shapeit/",out.prefix, ".map"))
  x$V5 <- c(diff(x$V3), 0)
  x$V5[which(x$V5 == 0)] <- 0.001
  x <- x[,c(4, 5, 3)]
  names(x) <- c("pposition", "rrate", "gposition")
  
  write.table(x, paste0("shapeit/", out.prefix, ".genmap"), row.names = F, quote = F)
  
  write.table(paste0("./shapeit --input-bed ", out.prefix, ".bed ", out.prefix, ".bim ", out.prefix, ".fam --input-map ", 
                     out.prefix, ".genmap --output-max ", out.prefix, ".haps ", out.prefix, ".sample --duohmm"),
              "shapeit/Run_Chr_12.sh", row.names = F, col.names = F, quote = F)
  
  
  write.table(paste0("./shapeit -convert --input-haps ", out.prefix, " --output-vcf ", out.prefix, ".vcf"),
              "shapeit/Run_Chr_12.sh", row.names = F, col.names = F, quote = F, append = T)
  
  
  
  
    
  
  system("rm shapeit/*.ped")
  system("rm shapeit/*.map")
  system("rm shapeit/*.snplist")
  system("rm shapeit/*.nosex")
  system("rm shapeit/*.log")
  
  #~~ NEED TO RUN THE SCRIPT ON EDDIE
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse SHAPEIT File                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

haplo.res <- NULL

for(i in 1:length(snp.window.vec)){
  for(j in 1:length(focal.snp.vec)){
    
    #~~ take the top SNP and SNPs either side.
    
    snp.window  <- snp.window.vec[i]
    focal.snp   <- focal.snp.vec[j]
    window.snps <- as.character(mapdata[(which(mapdata$SNP.Name == focal.snp) - snp.window):(which(mapdata$SNP.Name == focal.snp) + snp.window),"SNP.Name"])
    
    out.prefix <- paste0("Run_", focal.snp, "_span_", snp.window, "_ped")
    
    #~~ Read in files
    
    x <- ExtractHaploSHAPEIT(paste0("shapeit/", out.prefix, ".vcf"))$haplotypes
    
    x$snp.window <- snp.window
    x$focal.snp  <- focal.snp
    x$out.prefix <- out.prefix
    
    haplo.res <- rbind(haplo.res, x)
    
  }
}

rm(x, i, j, out.prefix, snp.window, focal.snp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Check consistency between estimates       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- subset(haplo.res, focal.snp == "cela1_red_10_26005249")
x <- subset(x, snp.window %in% c(20, 50))
x1 <- subset(x, snp.window == 50)
x <- subset(x, snp.window == 20)

head(x)
x <- subset(x, select = -c(focal.snp, out.prefix, snp.window))
x1 <- subset(x1, select = -c(focal.snp, out.prefix, snp.window))
names(x1)[3] <- "Haplo.Long"

x <- join(x, x1)
head(x)
x$Haplo.Sub <- substr(x$Haplo.Long, 31, 71)
table(x$Haplo == x$Haplo.Sub)


