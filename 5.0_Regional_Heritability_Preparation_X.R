
auto.vec <- 1:33
sex.vec <- 34
chr.vec <- 1:34

AnalysisSuffix <- "b"

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)

binwidths <- c(20)

pseudoautosomalSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

mapdata$CEL.LG[which(mapdata$SNP.Name %in% pseudoautosomalSNPs)] <- sex.vec + 1

recsumm <- read.table("results/3_Recombination_Phenotype_Data_a.txt", header = T, sep = "\t")
head(recsumm)
recsumm <- subset(recsumm, select = c(RRID, RRID.Sex))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create SNP lists                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.tab <- NULL

for(i in binwidths){
  
  for(j in sex.vec:(sex.vec+1)){
    
    x <- subset(mapdata, CEL.LG == j)
    
    snplist.start.vec <- sort(c(seq(1, nrow(x), i/2), nrow(x)-i+1))
    snplist.start.vec <- snplist.start.vec[-which(snplist.start.vec > nrow(x)-i+1)]
    
    for(k in snplist.start.vec){
      
      snplist <- x$SNP.Name[k:((k+i)-1)]
      
      analysis.id <- paste0("LG", j, ".", AnalysisSuffix, ".", k, ".", i)
      
      regh2.tab <- rbind(regh2.tab,
                         data.frame(Analysis.ID = analysis.id,
                                    CEL.LG = j,
                                    Start.Pos = k,
                                    Window.Size = i))
      
      write.table(snplist,
                  paste0("gcta/regh2_", AnalysisSuffix, "/snplist.", analysis.id, ".txt"),
                  row.names = F,
                  col.names = F,
                  quote = F)
    }
    
    
  }
  
}


regh2.tab <- unique(regh2.tab)

write.table(regh2.tab, paste0("gcta/regh2_", AnalysisSuffix, "/regh2_master_", AnalysisSuffix, "_X.txt"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Write scripts                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd("gcta")

#~~ update sex in the deer data:

idlist <- read.table("idlist.txt")
names(idlist) <- c("Fam", "RRID")

library(plyr)
idlist <- join(idlist, recsumm)
head(idlist)

idlist$RRID.Sex <- as.character(idlist$RRID.Sex)

idlist$RRID.Sex[which(idlist$RRID.Sex == "M")] <- 1
idlist$RRID.Sex[which(idlist$RRID.Sex == "F")] <- 2

idlist = unique(idlist)

write.table(idlist, "sexupdate.txt", row.names =  F, col.names = F, quote = F)

head(idlist)
idlist.f <- droplevels(subset(idlist, RRID.Sex == 2))

write.table(idlist.f[,1:2], "idlist.f.txt", row.names =  F, col.names = F, quote = F)

system("plink --bfile Deer31.recoded.v2 --autosome-num 33 --update-sex sexupdate.txt --make-bed --out Deer31.recoded.v3")

setwd("regh2_b/")

for(i in 1:nrow(regh2.tab)){
  if(regh2.tab$CEL.LG[i] == 34){
    system("cmd", input = paste0("gcta64.exe --bfile ..\\Deer31.recoded.v3 --autosome-num 33 --keep ..\\idlist.txt --extract snplist.",regh2.tab$Analysis.ID[i], ".txt --update-sex ..\\sexupdate.txt  --make-grm-xchr-gz --dc 1 --out ", regh2.tab$Analysis.ID[i], "_GRM"))
  } else {
    system("cmd", input = paste0("gcta64.exe --bfile ..\\Deer31.recoded.v3 --autosome-num 35 --keep idlist.txt --extract snplist.",regh2.tab$Analysis.ID[i], ".txt  --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_GRM"))
    system("cmd", input = paste0("gcta64.exe --grm-gz ", regh2.tab$Analysis.ID[i], "_GRM --grm-adj 0 --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_GRM_adj"))
  }
}

setwd("../..")



