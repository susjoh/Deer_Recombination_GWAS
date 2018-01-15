

AnalysisSuffix <- "c"

mapdata <- read.table("TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)
mapdata <- subset(mapdata, Q.2 >= 0.05)


binwidths <- c(6, 10, 20)

auto.vec <- 12
startSNP <- 241
stopSNP <- 401

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create SNP lists                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.tab <- NULL

for(i in binwidths){
  

    x <- subset(mapdata, CEL.LG == auto.vec)
    
    snplist.start.vec <- startSNP:stopSNP
    
    for(k in snplist.start.vec){
      
      snplist <- x$SNP.Name[k:((k+i)-1)]
      
      analysis.id <- paste0("LG", auto.vec, ".", AnalysisSuffix, ".", k, ".", i)
      
      regh2.tab <- rbind(regh2.tab,
                         data.frame(Analysis.ID = analysis.id,
                                    CEL.LG = auto.vec,
                                    Start.Pos = k,
                                    Window.Size = i))
      
      write.table(snplist,
                  paste0("regh2_", AnalysisSuffix, "/snplist.", analysis.id, ".txt"),
                  row.names = F,
                  col.names = F,
                  quote = F)
    }        

}


regh2.tab <- unique(regh2.tab)

write.table(regh2.tab, paste0("regh2_", AnalysisSuffix, "/regh2_master_", AnalysisSuffix, ".txt"), row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Write scripts                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for(i in 1:nrow(regh2.tab)){

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
../gcta64 --grm-gz ", regh2.tab$Analysis.ID[i], "_wo_GRM --grm-adj 0 --make-grm-gz --out ", regh2.tab$Analysis.ID[i], "_wo_GRM_adj"),
paste0("regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[i], ".sh"))
  
  setwd(paste0("regh2_", AnalysisSuffix))
  system(paste0("qsub ", regh2.tab$Analysis.ID[i], ".sh"))
  
  setwd("..")
}






