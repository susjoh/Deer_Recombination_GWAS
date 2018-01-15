
auto.vec <- 1:33
sex.vec <- 34
chr.vec <- 1:34

AnalysisSuffix <- "b"

mapdata <- read.table("TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)

binwidths <- c(20)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create SNP lists                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.tab <- NULL

for(i in binwidths){
  
  for(j in auto.vec){
    
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
                  paste0("regh2_", AnalysisSuffix, "/snplist.", analysis.id, ".txt"),
                  row.names = F,
                  col.names = F,
                  quote = F)
    }        
    
    
  }
  
}


regh2.tab <- unique(regh2.tab)

write.table(regh2.tab, paste0("regh2_", AnalysisSuffix, "/regh2_master_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Write scripts                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for(i in 1:nrow(regh2.tab)){
#for(i in 5001:10000){
  
  write.table(paste0("i = \"", regh2.tab$Analysis.ID[i], "\""),
              paste0("regh2_", AnalysisSuffix, "/",regh2.tab$Analysis.ID[i], ".R"),
              row.names = F, col.names = F, quote = F)
  
  system(paste0("cat 5.1_Regional_Heritability_Template.R >> regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[i], ".R"))
  
  writeLines(paste0("#!/bin/sh
#$ -cwd
#$ -l h_rt=00:30:00
#$ -V
#$ -l h_vmem=5200M

. /etc/profile.d/modules.sh
module load R
../gcta64 --bfile ../Deer31.recoded.v2 --autosome-num 33 --keep ../idlist.txt --extract snplist.",regh2.tab$Analysis.ID[i], ".txt  --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_GRM
../gcta64 --grm-gz ", regh2.tab$Analysis.ID[i], "_GRM --grm-adj 0 --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_GRM_adj
../gcta64 --bfile ../Deer31.recoded.v2 --autosome-num 33 --autosome --keep ../idlist.txt  --exclude snplist.", regh2.tab$Analysis.ID[i],".txt --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_wo_GRM
../gcta64 --grm-gz ", regh2.tab$Analysis.ID[i], "_wo_GRM --grm-adj 0 --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_wo_GRM_adj
R CMD BATCH ", regh2.tab$Analysis.ID[i], ".R"),
             paste0("regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[i], ".sh"))

setwd(paste0("regh2_", AnalysisSuffix))
  system(paste0("qsub ", regh2.tab$Analysis.ID[i], ".sh"))
  
 setwd("..")
}






