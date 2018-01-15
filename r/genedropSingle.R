#~~ Genedrop function for a single locus

genedropSingle <- function(ped, allele.ids, allele.founder.freqs, verbose = T){
  
  require(plyr)
  require(reshape2)
  require(kinship2)
  
  
  #~~ melt pedfile to get a unique row for each gamete transfer
  
  transped <- melt(ped[,c("ID", "MOTHER", "FATHER")], id.vars = "ID")
  transped$variable <- as.character(transped$variable)
  
  #~~ assign pedfile IDs to cohort and merge with transped
  
  cohorts <- data.frame(ID = ped[,1],
                        Cohort = kindepth(ped[,1], ped[,3], ped[,2]))
  
  suppressMessages(transped <- join(transped, cohorts))
  
  #~~ Redefine columns
  
  names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort")
  transped$Key    <- paste(transped$Parent.ID, transped$Offspring.ID, sep = "_")
  
  #~~ Recode founder gametes to 0 e.g. if one parent is unknown
  
  transped$Cohort[which(is.na(transped$Parent.ID))] <- 0
  transped$Cohort[which(transped$Parent.ID == 0)] <- 0
  
  
  #~~ Create founder haplotypes by sampling allele frequencies
  
  head(transped)
  
  #~~ Create a list to sample haplotypes
  
  haplo.list <- list()
  haplo.list[1:length(unique(transped$Offspring.ID))] <- list(list(MOTHER = NA, FATHER = NA))
  names(haplo.list) <- unique(as.character(transped$Offspring.ID))
  
  #~~ Sample the founder haplotypes
  
  system.time(for(i in which(transped$Cohort == 0)){
    
    if(transped$Parent.ID.SEX[i] == "MOTHER") {
      
      haplo.list[as.character(transped$Offspring.ID[i])][[1]]$MOTHER <- sample(allele.ids,
                                                                               size = 1,
                                                                               replace = T,
                                                                               prob = allele.founder.freqs)
      
    }
    
    if(transped$Parent.ID.SEX[i] == "FATHER") {
      
      haplo.list[as.character(transped$Offspring.ID[i])][[1]]$FATHER <- sample(allele.ids,
                                                                               size = 1,
                                                                               replace = T,
                                                                               prob = allele.founder.freqs)
      
    }
    
  })
  
  
  #~~ Sample the non-founder haplotypes by cohort. This loops through cohorts sequentially as
  #   parental haplotypes must exist before sampling.  
  
  for(cohort in 1:max(transped$Cohort)){
    
    if(verbose == TRUE) print(paste("Simulating haplotypes for cohort", cohort, "of", max(transped$Cohort)))
    
    for(j in which(transped$Cohort == cohort)){
      
      #~~ Sample one of the parental alleles at random
      
      haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[as.character(transped$Parent.ID[j])][[1]][sample.int(2, 1)]
      
    }
  }
  
  #~~ Condense haplotypes into genotypes
  
  genotype.list <- lapply(haplo.list, function(x) data.frame(Allele1 = x[1][[1]],
                                                             Allele2 = x[2][[1]]))
  
  geno.table <- do.call(rbind, genotype.list)
  
  geno.table$ID <- row.names(geno.table)
  
  geno.table <- join(geno.table, cohorts)
  
  geno.table
  
}
