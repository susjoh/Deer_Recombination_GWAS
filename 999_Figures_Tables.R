#
#
# Figures and Tables for Publication
# Susan Johnston
# July 2017
#
#
#

library(ggplot2)
library(asreml)
library(reshape)
library(plyr)

source("r/ASReml.EstEffects.R")
source("r/pin.R")
source("r/multiplot.R")

AnalysisSuffix <- "a"

load("gcta/recsumm_a.Rdata")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SEX DIFFERENCES IN RR                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


recsumm$RRID.Sex2 <- ifelse(recsumm$RRID.Sex == "F", "Females", "Males")

pdf("figs/SexDiffsInACC.pdf", width = 5, height = 5, useDingbats = F)

ggplot(recsumm, aes(RRID.Sex2, TotalRecombCount)) +
  geom_boxplot(notch = T, width = 0.4, fill = "lightgrey") +
  labs(x = "Sex", y = "Autosomal Crossover Count") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 16, colour = "black"),
        axis.text.y  = element_text (size = 14, colour = "black"),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16,
                                     angle = 90,
                                     margin=margin(0,20,0,0)),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  scale_y_continuous(breaks = seq(10, 50, 5))

dev.off()  




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HERITABILITY                          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

model.vec <- load("results/4_Heritability_Results.RData", verbose = T)


#~~ Fixed effect results  

fixwald <- rbind(cbind(data.frame(wald.asreml(grm.summ.RR))  , Analysis = "Both"  , Effect = row.names(data.frame(wald.asreml(grm.summ.RR)))),
                 cbind(data.frame(wald.asreml(grm.summ.RR.f)), Analysis = "Female", Effect = row.names(data.frame(wald.asreml(grm.summ.RR.f)))),
                 cbind(data.frame(wald.asreml(grm.summ.RR.m)), Analysis = "Male"  , Effect = row.names(data.frame(wald.asreml(grm.summ.RR.m)))))

fixwald <- fixwald[-grep("residual", fixwald$Effect),]
fixwald <- arrange(fixwald, Analysis, Effect)


fixres <- rbind(cbind(data.frame(summary(grm.summ.RR, all = T)$coef.fixed)  , Analysis = "Both"  , Effect = row.names(summary(grm.summ.RR, all = T)$coef.fixed)),
                cbind(data.frame(summary(grm.summ.RR.f, all = T)$coef.fixed), Analysis = "Female", Effect = row.names(summary(grm.summ.RR.f, all = T)$coef.fixed)),
                cbind(data.frame(summary(grm.summ.RR.m, all = T)$coef.fixed), Analysis = "Male"  , Effect = row.names(summary(grm.summ.RR.m, all = T)$coef.fixed)))

fixres <- fixres[-grep("_F", fixres$Effect),]
fixres <- arrange(fixres, Analysis, Effect)

fixres
fixwald$Effect <- gsub("RRID.Sex", "RRID.Sex_M", fixwald$Effect)

fixres <- join(fixres, fixwald)
fixres <- subset(fixres, select = c(Analysis, Effect, solution, std.error, z.ratio, Pr.Chisq.))

fixres$solution <- round(fixres$solution, digits = 3)
fixres$std.error <- round(fixres$std.error, digits = 3)
fixres$Pr.Chisq. <- round(fixres$Pr.Chisq., digits = 3)
fixres$Pr.Chisq. <- ifelse(fixres$Pr.Chisq. < 0.00001, "<2e-16", fixres$Pr.Chisq.)

fixres$Effect <- c("Intercept", "Fhat3", "Male", "Intercept", "Fhat3", "Intercept", "Fhat3")
names(fixres) <- c("Analysis", "Effect", "Estimate", "S.E.", "Z", "P (Wald Test)")

write.table(fixres, "doc/tables/4_Fixed_Effects_Table.txt", row.names = F, sep = " & ", quote = F, eol = " \\\\\n")

rm(fixwald)

#~~ tabulate random effect results

h2res <- rbind(cbind(Analysis = "Both"   , ASReml.EstEffects(grm.summ.RR), Comp = c("VA", "VPE", "VR")),
               cbind(Analysis = "Females", ASReml.EstEffects(grm.summ.RR.f), Comp = c("VA", "VPE", "VR")),
               cbind(Analysis = "Males"  , ASReml.EstEffects(grm.summ.RR.m), Comp = c("VA", "VPE", "VR")))


#~~ deal with effect sizes

h2res.effects  <- subset(h2res, select = c(Analysis, Effect, SE, Comp))

h2res.effects$Effect.Text <- round(h2res.effects$Effect, digits=2)
h2res.effects$SE.Text <- round(h2res.effects$SE, digits=2)
h2res.effects$Res.Text <- paste0(h2res.effects$Effect.Text, " (", h2res.effects$SE.Text, ")")
h2res.effects <- subset(h2res.effects, select = c(Analysis, Comp, Res.Text))
h2res.effects <- cast(h2res.effects, Analysis ~ Comp)
names(h2res.effects)[2:4] <- c("h2", "pe2", "e2")

#~~ deal with variance components sizes

h2res.variance <- subset(h2res, select = c(Analysis, component, std.error, Comp))

h2res.variance$Effect.Text <- round(h2res.variance$component, digits=2)
h2res.variance$SE.Text <- round(h2res.variance$std.error, digits=2)
h2res.variance$Res.Text <- paste0(h2res.variance$Effect.Text, " (", h2res.variance$SE.Text, ")")
h2res.variance <- subset(h2res.variance, select = c(Analysis, Comp, Res.Text))
h2res.variance <- cast(h2res.variance, Analysis ~ Comp)

h2res <- join(h2res.variance, h2res.effects)

#~~ Get VP


x <- rbind(round(pin(grm.summ.RR.forvp, phen.var ~ V1 + V2 + V4), digits = 2),
           round(pin(grm.summ.RR.f.forvp, phen.var ~ V1 + V2 + V4), digits = 2),
           round(pin(grm.summ.RR.m.forvp, phen.var ~ V1 + V2 + V4), digits = 2))

h2res$VP <- paste0(x$Estimate, " (", x$SE, ")")


#~~ Get significance of heritabilities

h2res$P_h2 <- c(1 - pchisq(2*(grm.summ.RR$loglik - grm.summ.RR.woh2$loglik), 1),
                1 - pchisq(2*(grm.summ.RR.f$loglik - grm.summ.RR.f.woh2$loglik), 1),
                1 - pchisq(2*(grm.summ.RR.m$loglik - grm.summ.RR.m.woh2$loglik), 1))


#~~ Get means and sample sizes

h2res$Mean <- paste0(round(c(mean(recsumm$TotalRecombCount), tapply(recsumm$TotalRecombCount, recsumm$RRID.Sex, mean)), digits = 2), " (",
                     round(c(sd(recsumm$TotalRecombCount), tapply(recsumm$TotalRecombCount, recsumm$RRID.Sex, sd)), digits = 2), ")")
h2res$Nobs <- round(c(nrow(recsumm), nrow(recsumm.f), nrow(recsumm.m)), digits = 2)
h2res$Nfid <- round(c(length(unique(recsumm$RRID)),
                      length(unique(recsumm.f$RRID)),
                      length(unique(recsumm.m$RRID))), digits = 2)

#~~ Get the number of crossovers

rectab <- read.table("results/2_Per_Chromosome_Recomb_Final_a.txt", header = T, stringsAsFactors = F, sep = "\t")
rectab <- subset(rectab, RRID %in% recsumm$RRID)
h2res$Nxovers <- c(sum(rectab$RecombCount),
                   sum(subset(rectab, RRID %in% recsumm.f$RRID)$RecombCount),
                   sum(subset(rectab, RRID %in% recsumm.m$RRID)$RecombCount))

h2res <- h2res[,c("Analysis", "Nobs", "Nfid", "Mean", "Nxovers", "VP", "VA", "VPE", "VR", "h2", "pe2", "e2", "P_h2")]
write.table(h2res, "doc/tables/4_Heritability_Table.txt", row.names = F, sep = " & ", quote = F, eol = " \\\\\n")



heritability.results

temp <- heritability.results$bivar.ranef

temp$gamma[1] /((temp$gamma[2]* temp$gamma[3] )^0.5) 


#~~ Remove objects you don't need
eval(parse(text = paste("rm(", paste(model.vec, collapse = ", "), ")")))
rm(rectab, x, h2res.effects, h2res.variance, model.vec)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GENOME-WIDE ASSOCIATION               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

restab <- read.table("results/4_GWAS_results_a.txt", header = T, stringsAsFactors = F)

mapdata <- unique(subset(restab, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position, Cumu)))
mapdata <- arrange(mapdata, Cumu)

#~~ Get Chromosome Positions

chrinfo <- NULL

for(i in na.omit(unique(mapdata$CEL.LG))){
  
  temp1 <- arrange(subset(mapdata, CEL.LG == i), Estimated.Mb.Position)
  
  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$CEL.LG)))
colourscale <- colourscale[1:(length(colourscale)/2)]

bonf = 0.05/35263.64

chrinfo$CEL.LG[35] <- "X"
chrinfo$CEL.LG2 <- c("", 1, "", "", 4, "", "6", "", "8", "", 10, "", "", 13, "", "", 16, "", "", 19, "", "", 22, "", "", "", "", 27, "", "", "", 31, "", "", "X")


#~~ PLOT

restab$PlotSex <- ifelse(restab$Sex == "all", " A. Both",
                         ifelse(restab$Sex == "f", "B. Females", "C. Males"))
restab <- droplevels(subset(restab, Phenotype == "acc"))

png(paste0("figs/4_GWAS_Manhattan_Plot_", AnalysisSuffix, "_ForPaper.png"), width = 9, height = 9, units = "in", res = 300)

multiplot(
  
  ggplot(restab, aes(Cumu,-log10(Pc1df), col = factor(CEL.LG))) +
    geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
    geom_point(size = 2, alpha = 0.5) +
    scale_colour_manual(values = colourscale) +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 12),
          axis.text.y  = element_text (size = 14),
          strip.text.x = element_text (size = 16, hjust = 0),
          axis.title.y = element_text (size = 16, angle = 90),
          axis.title.x = element_text (size = 16),
          strip.background = element_blank(),
          legend.position="none") +
    scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG2) +
    scale_y_continuous(breaks = 0:7) +
    labs(x ="CEL Linkage Group", y = "-log10 P") +
    facet_wrap(~PlotSex, ncol = 1) +
    coord_cartesian(ylim = c(0, 7))
  
  ,
  
  ggplot(restab, aes(x= -log10(ExpP), y = -log10(Pc1df))) +
    geom_point(size = 2, alpha = 0.5) +
    geom_abline(intercept=0,slope=1) +
    theme_bw() +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          strip.background = element_blank()) +
    labs(x="Expected -log10 P",y="Observed -log10 P")+
    scale_y_continuous(breaks = 0:7) +
    scale_x_continuous(breaks = 0:7) +
    
    facet_wrap(~PlotSex, ncol = 1)+
    coord_cartesian(ylim = c(0, 7))
  
  
  , cols = 2, layout = matrix(c(1,1,2), 1))

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# Tables

mapdata <- read.table("data/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
mapdata <- subset(mapdata, select = c(SNP.Name, cMPosition.SexAveraged, CEL.LG, Q.2))


restab <- read.table("results/4_GWAS_results_a.txt", header = T, stringsAsFactors = F)
temp <- read.table("results/4_GWAS_results_trans_a.txt", header = T, stringsAsFactors = F)

restab$Type <- "Cis and Trans"
temp$Type <- "Trans only"

restab <- rbind(restab, temp)
rm(temp)

restab <- subset(restab, Phenotype == "acc")

restab <- join(restab, mapdata)

restab <- subset(restab, select = c(SNP.Name, CEL.LG, cMPosition.SexAveraged, Chromosome,
                                    Position, A1, A2, effB, se_effB, P1df, chi2.1df.adj, Pc1df,
                                    Sex, Estimated.Mb.Position, Type, Q.2))
names(restab) <- c("SNP.Name", "CEL.LG", "CEL.cM", "BTA.Chr", 
                   "BTA.Position", "A1", "A2", "effB", "se_effB", "P1df", "chi2.1df.adj", 
                   "Pc1df", "Sex", "CEL.Est.bp.Position", "Cis.Trans", "MAF")

restab <- arrange(restab, Cis.Trans, Sex, CEL.LG, CEL.Est.bp.Position)

write.table(restab, "results/TABLE_ACC_GWAS_Results.txt", row.names = F, sep = ",", quote = F)


tophits <- rbind(arrange(subset(restab, Sex == "all"& Cis.Trans == "Cis and Trans"), Pc1df) [1:10,],
                 arrange(subset(restab, Sex == "f"  & Cis.Trans == "Cis and Trans"), Pc1df) [1:10,],
                 arrange(subset(restab, Sex == "m"  & Cis.Trans == "Cis and Trans"), Pc1df) [1:10,])


tophits$Cis.Trans <- NULL
tophits$CEL.Est.bp.Position <- NULL
tophits$P1df <- NULL

head(tophits)
tophits$CEL.cM <- round(tophits$CEL.cM, digits = 1)
tophits$effB <- round(tophits$effB, digits = 2)
tophits$se_effB <- round(tophits$se_effB, digits = 2)
tophits$chi2.1df.adj <- round(tophits$chi2.1df.adj, digits = 2)
tophits$Pc1df <- format(tophits$Pc1df, digits = 3)
tophits$MAF <- round(tophits$MAF, digits = 2)



tophits$Sex2 <- ifelse(tophits$Sex == "all", " A. Both",
                         ifelse(tophits$Sex == "f", "B. Females", "C. Males"))
tophits$Sex <- NULL
tophits <- cbind(Sex = tophits$Sex2, tophits[,-ncol(tophits)])
tophits$SNP.Name <- gsub("_", "\\_", tophits$SNP.Name, fixed = T)

write.table(tophits, "doc/tables/4_Top_GWAS_Hits.txt", row.names = F, sep = " & ", quote = F, eol = "\\\\\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# Candidate Genes


bos.genes <- read.table("data/Bos_taurus.UMD3.1.89.chr.gff3.geneinfo", header = T, sep = "\t", quote = "\"", stringsAsFactors = F)
head(bos.genes)

snp.vec <- unique(gsub("-201", "", na.omit(subset(bos.genes, Chromosome == 10 & RefStart > 25505249 & RefStart < 26505249)$Name)))
snp.vec <- unique(gsub("-202", "", snp.vec))
snp.vec <- unique(gsub("-203", "", snp.vec))








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# REGIONAL HERITABILITY                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm$RRID.Sex2 <- ifelse(recsumm$RRID.Sex == "F", "Females", "Males")

restab <- read.table("results/4_GWAS_results_a.txt", header = T, stringsAsFactors = F)

mapdata <- unique(subset(restab, select = c(SNP.Name, CEL.LG, Estimated.Mb.Position, Cumu)))
mapdata <- arrange(mapdata, Cumu)

regh2 <- read.table("results/5_Regional_Heritability_Results_b.txt", header = T, stringsAsFactors = F)
regvar <- read.table("results/5_Regional_Heritability_VarComp_b.txt", header = T, stringsAsFactors = F)

head(regh2)
head(regvar)

regvar$z.ratio <- ifelse(regvar$constraint == "Boundary", NA, regvar$z.ratio)

regvar1 <- subset(regvar, Model == "wi.region" & Effect.1 == "giv(RRID).giv")
regvar2 <- subset(regvar, Model == "wi.region" & Effect.1 == "giv(RRID2).giv")

regvar1 <- subset(regvar1, select = c(Analysis.ID, Effect, SE, z.ratio, RRID.Sex, CEL.LG, Start.Order, Window.Size))
regvar2 <- subset(regvar2, select = c(Analysis.ID, Effect, SE, z.ratio, RRID.Sex, CEL.LG, Start.Order, Window.Size))

names(regvar1)[2:4] <- c("Genome.h2", "Genome.se", "Genome.Z")
names(regvar2)[2:4] <- c("Region.h2", "Region.se", "Region.Z")

regvar <- join(regvar1, regvar2)
head(regvar)

regh2 <- join(regh2, regvar)

#~~ Get chromosome info

chrinfo <- NULL

for(i in na.omit(unique(regh2$CEL.LG))){
  
  temp1 <- arrange(subset(regh2, CEL.LG == i), Mid.Cumu.Position)
  
  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Mid.Cumu.Position"],
                      Stop = temp1[nrow(temp1),"Mid.Cumu.Position"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

names(chrinfo) <- c("CEL.LG", "Start", "Stop")

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$CEL.LG2 <- c(1, "", "", 4, "", "", 7, "", "", 10, "", "", 13, "", "", 16, "", "", 19, "", "", 22, "", "", "", "", 27, "", "", "", 31, "", "")



bonf20 = 0.05/(length(unique(subset(regh2, Window.Size == 20)$Analysis.ID))/2)


#~~ Plots

regh2$PlotSex <- ifelse(regh2$RRID.Sex == "Both", "A. Both",
                         ifelse(regh2$RRID.Sex == "Female", "B. Females", "C. Males"))
regh2 <- subset(regh2, Window.Size == 20)

png("figs/5_Regional_Heritability_b_ForPaper.png", width = 7, height = 9, units = "in", res = 300)

ggplot(regh2,aes(Mid.Cumu.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  geom_line() +
  geom_point() +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 16, hjust = 0),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank(),
        legend.position="none") +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(x ="CEL Linkage Group", y = "-log10 P") +
  facet_wrap(~PlotSex, ncol = 1) +
  coord_cartesian(ylim = c(0, 8))

dev.off()

#~~ Full results table

head(regh2)

regh2 <- arrange(regh2, RRID.Sex, CEL.LG, Start.Order)
head(mapdata)
mapdata$Cumu <- NULL

names(mapdata) <- c("Start.SNP.Name", "CEL.LG", "Start.Position")
regh2 <- join(regh2, mapdata)

names(mapdata) <- c("Stop.SNP.Name", "CEL.LG", "Stop.Position")
regh2 <- join(regh2, mapdata)

head(regh2)

reg.final <- subset(regh2, select = c(PlotSex, CEL.LG, Start.Order, Window.Size, Chi2, P, Start.SNP.Name, Stop.SNP.Name, Genome.h2, Genome.se, Genome.Z, Region.h2, Region.se, Region.Z))
head(reg.final)
names(reg.final)[1] <- "Sex"

write.table(reg.final, "results/TABLE_ACC_RegHeritability_Results.txt", sep = "\t", row.names = F, quote = F)

#~~ Top hits for the paper


tophits <- rbind(arrange(subset(reg.final, Sex == "A. Both")  , P) [1:10,],
                 arrange(subset(reg.final, Sex == "B. Females"), P) [1:10,],
                 arrange(subset(reg.final, Sex == "C. Males")  , P) [1:10,])

tophits$Start.SNP.Name <- gsub("_", "\\_", tophits$Start.SNP.Name, fixed = T)
tophits$Stop.SNP.Name <- gsub("_", "\\_", tophits$Stop.SNP.Name, fixed = T)


tophits$Chi2 <- round(tophits$Chi2, digits = 2)
tophits$P <- format(tophits$P, digits = 3)

tophits$Genome.h2 <- round(tophits$Genome.h2, digits = 3)
tophits$Genome.se <- round(tophits$Genome.se, digits = 3)
tophits$Genome.Z <- round(tophits$Genome.Z, digits = 2)
tophits$Region.h2 <- round(tophits$Region.h2, digits = 3)
tophits$Region.se <- round(tophits$Region.se, digits = 3)
tophits$Region.Z <- round(tophits$Region.Z, digits = 2)


write.table(tophits, "doc/tables/4_Top_Regh2_Hits.txt", row.names = F, sep = " & ", quote = F, eol = "\\\\\n")

#~~ Zoom into the region

bonf = 0.05/35263.64

restab$PlotSex <- ifelse(restab$Sex == "all", "A. Both",
                         ifelse(restab$Sex == "f", "B. Females", "C. Males"))


png("figs/5_Regional_Heritability_CEL12_b_ForPaper.png", width = 7, height = 9, units = "in", res = 300)

ggplot() +
  geom_hline(yintercept=-log10(bonf20),linetype=1, alpha = 0.4, size = 1, colour = "red") +
  geom_hline(yintercept=-log10(bonf),linetype=5, alpha = 0.4, size = 1) +
  geom_point(data = subset(restab, CEL.LG == 12), aes(Estimated.Mb.Position/1e6, -log10(Pc1df)), alpha = 0.5) +
  geom_line (data = subset(regh2 , CEL.LG == 12), aes(Mid.Position/1e6, -log10(P), group = CEL.LG), colour = "red") + 
  geom_point (data = subset(regh2 , CEL.LG == 12), aes(Mid.Position/1e6, -log10(P), group = CEL.LG), colour = "red", alpha = 0.5) + 
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 16, hjust = 0),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank(),
        legend.position="none") +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  scale_x_continuous(breaks = seq(0, 110, 10)) +
  labs(x ="CEL.LG Estimated Position (MB)", y = "-log10 P") +
  facet_wrap(~PlotSex, ncol = 1) +
  coord_cartesian(ylim = c(0, 8))

dev.off()

bos.genes <- read.table("data/Bos_taurus.UMD3.1.89.chr.gff3.geneinfo", header = T, sep = "\t", quote = "\"", stringsAsFactors = F)
head(bos.genes)

snp.vec <- subset(bos.genes, Chromosome == 27 & RefStart > 38731584 & RefStart < 41274975)
snp.vec <- subset(snp.vec, select = c(Name, RefStart, RefStop))
snp.vec <- na.omit(snp.vec)
snp.vec$Name <- gsub("-201", "", snp.vec$Name)
snp.vec$Name <- gsub("-202", "", snp.vec$Name)
snp.vec$Name <- gsub("-203", "", snp.vec$Name)
snp.vec$Name <- gsub("-204", "", snp.vec$Name)

snp.vec <- unique(snp.vec)

