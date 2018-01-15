#
# Investigate the genes and LD in the associated regions
# Susan Johnston
# July 2017
#
#


library(ggplot2)
library(plyr)
library(dplyr)
library(GenABEL)
library(reshape)

#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

#~~ Read in phenotype file

recsumm <- read.table(paste0("results/3_Recombination_Phenotype_Data_", AnalysisSuffix, ".txt"), header = T, stringsAsFactors = F)

#~~ Read in map file

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- unique(subset(mapdata, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position, CEL.order, BTA.Position)))

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Load GWAS Results

gwastab <- read.table(paste0("results/4_GWAS_results_", AnalysisSuffix, ".txt"), header = T, sep = "\t", stringsAsFactors = F)
gwastab <- subset(gwastab, Phenotype == "acc")

gwastab$CisTrans <- "cistrans"

transgwastab <- read.table(paste0("results/4_GWAS_results_trans_", AnalysisSuffix, ".txt"), header = T, sep = "\t", stringsAsFactors = F)
transgwastab <- subset(transgwastab, Phenotype == "acc")

transgwastab$CisTrans <- "trans"

gwastab <- rbind(gwastab, transgwastab)
rm(transgwastab)

#~~ Load Regional heritability results
AnalysisSuffix <- "b"


regtab <- read.table(paste0("results/5_Regional_Heritability_Results_", AnalysisSuffix, ".txt"), header = T, sep = "\t", stringsAsFactors = F)

#~~ Load annotation

bos.genes <- read.table("data/Bos_taurus.UMD3.1.89.chr.gff3.geneinfo", header = T, sep = "\t", quote = "\"", stringsAsFactors = F)

#~~ Add SNP summary information to gwastab

snp.summ <- summary.snp.data(gtdata(abeldata))
snp.summ$SNP.Name <- row.names(snp.summ)
snp.summ <- subset(snp.summ, select = c(SNP.Name, Q.2, P.11, P.12, P.22))
head(snp.summ)

gwastab <- join(gwastab, snp.summ)

gc()

source("r/ldPlot.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Revisit Significant GWAS Regions            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gwastab <- arrange(gwastab, Pc1df)

#~~ Look at genotype effects

x <- data.frame(as.character.gwaa.data(abeldata[,unique(gwastab$SNP.Name)[1:10]]))
x$RRID <- row.names(x)
x <- melt(x, id.vars = "RRID")

y <- recsumm[,c("RRID", "RRID.Sex", "TotalRecombCount")]
y$RRID <- as.character(y$RRID)

x <- join(x, y)
x <- subset(x, !is.na(TotalRecombCount))

x$value <- factor(x$value, levels = sort(unique(as.character(x$value))))
x <- na.omit(x)

ggplot(x, aes(value, TotalRecombCount, fill = RRID.Sex)) +
  geom_boxplot(notch = T) + 
  facet_wrap(~variable, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Genotype", y = "Crossover Count")
ggsave(paste0("figs/4_Top_genotype_hits_pheno_", AnalysisSuffix, ".png"), width = 10, height = 10)


#~~ Examine associated region in more detail

gwastab <- join(gwastab, mapdata)

ggplot(subset(gwastab, Phenotype == "acc" & CEL.LG == 12), aes(BTA.Position,-log10(Pc1df))) +
  geom_point(size = 2, alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.05) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="BTA_3.1 Position (bp)", y = "-log10 P") +
  facet_wrap(~Sex, ncol = 1)

rm(x, y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Annotation Profiles                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Candidate genes

cand.genes <- c("RNF212", "RNF212B", "REC8", "PIGG", "CPLX1", "HEI10", "PRDM9", "GCLM", "MSH4",
                "REC114", "PABPN1", "FMN1", "NEK9", "CEP55", "SMC3", "CCNB1IP1", "C14orf39", "SMEK1",
                "RAD21", "CCDC43")

cand.tab <- unique(bos.genes[grep(paste(cand.genes, collapse = "|"), bos.genes$Name, ignore.case = T),])
cand.tab <- subset(cand.tab, select = c(Chromosome, Type, RefStart, RefStop, Name))


#~~ Which the candidate genes not annotated?

cand.genes[which( cand.genes %in% bos.genes$Name)]

paste(cand.genes[which(!cand.genes %in% bos.genes$Name)], collapse = " or ") 

#(RNF212 or CPLX1 or HEI10 or PRDM9 or REC114 or FMN1 or SMEK1) AND "Bos taurus"[porgn:__txid9913]
# HEI10 is the only one missing, but actually it is a orthologue of CCNB1IP1


cand.tab <- rbind(cand.tab, data.frame(Chromosome = c(10, 10, 10, 10, 21, 1, "X"),
                                       Type = "BLAST",
                                       RefStart = c(109178914, 109070252, 20022794, 29962320, 56834182, 45033572, 143415608),
                                       RefStop  = c(109170517, 109071838, 20023190, 29970835, 56834937, 450345522, 143416557),
                                       Name = c("RNF212", "CPLX1", "REC114", "FMN1", "SMEK1", "PRDM9", "PRDM9")))


#~~~~ Look at your data

gwastab[which(gwastab$Analysis == "gwas.acc.all"),][1:20,]

# There are several regions of interest: CEL.LG 12 (RNF212B and REC8 are here),
# CEL16 (although may be spurious rare allele), CEL31, CEL11, CEL5 & CEL32

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. CEL12                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.gwas <- subset(gwastab, Analysis == "gwas.acc.all" & CEL.LG == 12 & CisTrans == "cistrans")
arrange(x.gwas, Pc1df)[1:20,]

x.regh2 <- subset(regtab, CEL.LG == 12 & RRID.Sex == "Both")
arrange(x.regh2, P)[1:20,]

ggplot() +
  geom_point(data = x.gwas , aes(BTA.Position, -log10(Pc1df))) +
  geom_line (data = x.regh2, aes(Mid.BTA.Position         , -log10(P    )), colour = "red")



x <- subset(x.gwas, BTA.Position > 10e6 & BTA.Position < 35e6)
x <- arrange(x, BTA.Position)
head(x)

subset(bos.genes, Chromosome == x$Chromosome[1] &  RefStart > x$Position[1] & RefStop < x$Position[nrow(x)])[,c(1, 2, 3, 4, 7, 9)]


#~~ Examine LD

ldtab <- ldPlot(abeldata, x$SNP.Name, mapdata = mapdata) 

y <- subset(gwastab, SNP.Name %in% x$SNP.Name & Analysis == "gwas.acc.all")
ldtab <- subset(ldtab, X1 == arrange(x, Pc1df)$SNP.Name[1] | X2 == arrange(x, Pc1df)$SNP.Name[1])
head(ldtab)
ldtab$SNP.Name <- ifelse(ldtab$X1 == arrange(x, Pc1df)$SNP.Name[1], ldtab$X2, ldtab$X1)

y <- join(y, subset(ldtab, select = c(SNP.Name, value)))
y$value[1] <- 1

ggplot(y, aes(BTA.Position,-log10(Pc1df), colour = value)) +
  geom_point(size = 2, alpha = 1) +
  scale_colour_gradient(low = "red", high = "black") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P", col = "LD")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CEL32                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.gwas <- subset(gwastab, Analysis == "gwas.acc.all" & CEL.LG == 32 & CisTrans == "cistrans")
arrange(x.gwas, Pc1df)[1:20,]

x.regh2 <- subset(regtab, CEL.LG == 32 & RRID.Sex == "Both")
arrange(x.regh2, P)[1:20,]

ggplot() +
  geom_point(data = x.gwas , aes(Estimated.Mb.Position, -log10(Pc1df))) +
  geom_line (data = x.regh2, aes(Mid.Position         , -log10(P    )), colour = "red")



x <- subset(x.gwas, Estimated.Mb.Position > 35e6 & Estimated.Mb.Position < 47e6)
x <- arrange(x, Estimated.Mb.Position)
head(x)

arrange(x, Pc1df)$SNP.Name[1]

subset(bos.genes, Chromosome == x$Chromosome[1] &  RefStart > x$Position[1] & RefStop < x$Position[nrow(x)])[,c(1, 2, 3, 4, 7, 9)]

#~~ Examine LD

ldtab <- ldPlot(abeldata, x$SNP.Name, mapdata = mapdata) 

y <- subset(gwastab, SNP.Name %in% x$SNP.Name & Analysis == "gwas.acc.all")
ldtab <- subset(ldtab, X1 == arrange(x, Pc1df)$SNP.Name[1] | X2 == arrange(x, Pc1df)$SNP.Name[1])
head(ldtab)
ldtab$SNP.Name <- ifelse(ldtab$X1 == arrange(x, Pc1df)$SNP.Name[1], ldtab$X2, ldtab$X1)

y <- join(y, subset(ldtab, select = c(SNP.Name, value)))
y$value[1] <- 1

ggplot(y, aes(Estimated.Mb.Position,-log10(Pc1df), colour = value)) +
  geom_point(size = 2, alpha = 1) +
  scale_colour_gradient(low = "red", high = "black") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P", col = "LD")

rm(x, x.gwas, x.regh2, y)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CEL16                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.gwas <- subset(gwastab, Analysis == "gwas.acc.all" & CEL.LG == 16 & CisTrans == "cistrans")
arrange(x.gwas, Pc1df)[1:20,]

x.regh2 <- subset(regtab, CEL.LG == 16 & RRID.Sex == "Both")
arrange(x.regh2, P)[1:20,]

ggplot() +
  geom_point(data = x.gwas , aes(Estimated.Mb.Position, -log10(Pc1df))) +
  geom_line (data = x.regh2, aes(Mid.Position         , -log10(P    )), colour = "red")



x <- subset(x.gwas, Estimated.Mb.Position > 30e6 & Estimated.Mb.Position < 38e6)
x <- arrange(x, Estimated.Mb.Position)
head(x)

arrange(x, Pc1df)$SNP.Name[1]

subset(bos.genes, Chromosome == x$Chromosome[1] &  RefStart > x$Position[1] & RefStop < x$Position[nrow(x)])[,c(1, 2, 3, 4, 7, 9)]


#~~ Examine LD

ldtab <- ldPlot(abeldata, x$SNP.Name, mapdata = mapdata)

y <- subset(gwastab, SNP.Name %in% x$SNP.Name & Analysis == "gwas.acc.all")
ldtab <- subset(ldtab, X1 == arrange(x, Pc1df)$SNP.Name[1] | X2 == arrange(x, Pc1df)$SNP.Name[1])
head(ldtab)
ldtab$SNP.Name <- ifelse(ldtab$X1 == arrange(x, Pc1df)$SNP.Name[1], ldtab$X2, ldtab$X1)

y <- join(y, subset(ldtab, select = c(SNP.Name, value)))
y$value[1] <- 1

ggplot(y, aes(Estimated.Mb.Position,-log10(Pc1df), colour = value)) +
  geom_point(size = 2, alpha = 1) +
  scale_colour_gradient(low = "red", high = "black") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P", col = "LD")

rm(x, x.gwas, x.regh2, y)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CEL30                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x.gwas <- subset(gwastab, Analysis == "gwas.acc.all" & CEL.LG == 30 & CisTrans == "cistrans")
arrange(x.gwas, Pc1df)[1:20,]

x.regh2 <- subset(regtab, CEL.LG == 16 & RRID.Sex == "Both")
arrange(x.regh2, P)[1:20,]

ggplot() +
  geom_point(data = x.gwas , aes(Estimated.Mb.Position, -log10(Pc1df))) +
  geom_line (data = x.regh2, aes(Mid.Position         , -log10(P    )), colour = "red")



x <- subset(x.gwas, Estimated.Mb.Position > 35e6 & Estimated.Mb.Position < 45e6)
x <- arrange(x, Estimated.Mb.Position)
head(x)

arrange(x, Pc1df)$SNP.Name[1]

subset(bos.genes, Chromosome == x$Chromosome[1] &  RefStart > x$Position[1] & RefStop < x$Position[nrow(x)])[,c(1, 2, 3, 4, 7, 9)]


#~~ Examine LD

ldtab <- ldPlot(abeldata, x$SNP.Name, mapdata = mapdata)

y <- subset(gwastab, SNP.Name %in% x$SNP.Name & Analysis == "gwas.acc.all")
ldtab <- subset(ldtab, X1 == arrange(x, Pc1df)$SNP.Name[1] | X2 == arrange(x, Pc1df)$SNP.Name[1])
head(ldtab)
ldtab$SNP.Name <- ifelse(ldtab$X1 == arrange(x, Pc1df)$SNP.Name[1], ldtab$X2, ldtab$X1)

y <- join(y, subset(ldtab, select = c(SNP.Name, value)))
y$value[1] <- 1

ggplot(y, aes(Estimated.Mb.Position,-log10(Pc1df), colour = value)) +
  geom_point(size = 2, alpha = 1) +
  scale_colour_gradient(low = "red", high = "black") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x ="Estimated Position (bp)", y = "-log10 P", col = "LD")

rm(x, x.gwas, x.regh2, y)


