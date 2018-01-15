ldPlot <- function(gwaa.genabel, snp.subset, mapdata){
  
  if(!all(names(mapdata) == c("SNP.Name", "CEL.LG", "Estimated.Mb.Position", "CEL.order", "BTA.Position"))) stop ("mapdata object must have names SNP.Name CEL.LG Estimated.Mb.Position CEL.order BTA.Position")
  
  system.time(ldtab <- r2fast(gwaa.genabel, snpsubset = snp.subset))
  ldtab[1:5, 1:5]
  ldtab <- melt(ldtab)
  ldtab$X1 <- as.character(ldtab$X1)
  ldtab$X2 <- as.character(ldtab$X2)
  head(ldtab)
  ldtab <- subset(ldtab, value < 1.01)
  
  
  
  names(mapdata) <- c("X1", "X1.LG", "X1.Pos", "X1.order", "X1.BTA.Position")
  ldtab <- join(ldtab, mapdata)
  names(mapdata) <- c("X2", "X2.LG", "X2.Pos", "X2.order", "X2.BTA.Position")
  ldtab <- join(ldtab, mapdata)
  names(mapdata) <- c("SNP.Name", "CEL.LG", "Estimated.Mb.Position", "CEL.order", "BTA.Position")
  head(ldtab)
  
  library(ggplot2)
  
  print(ggplot(ldtab, aes(X2.order, X1.order, fill = value)) + 
    geom_tile() + 
    scale_fill_continuous(low = "white", high = "red"))
  
  ldtab
}