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
library(asreml)

source("r/Shapeit_Functions.R")
source("r/Beagle_Functions.R")


#~~ Add CriMAP requirements

AnalysisSuffix <- "a"

#~~ Read in haplotypes

haplo.res <- ExtractHaploSHAPEIT("shapeit/Run_Chr_12_ped.vcf",  make.ref.alleles = T)$haplotypes
str(haplo.res)

haplo.map <- read.table("shapeit/Run_Chr_12_ped.vcf", skip=6, stringsAsFactors = F)[,1:3]
names(haplo.map) <- c("Chr", "Position", "SNP.Name")

head(haplo.map)
haplo.map$BTA.Position <- sapply(haplo.map$SNP.Name, function(x) strsplit(x, split = "_")[[1]][4])
haplo.map$BTA.Position <- as.numeric(haplo.map$BTA.Position)

recsumm <- read.table("results/3_Recombination_Phenotype_Data_a.txt", header = T, sep = "\t")
head(recsumm)

recsumm <- subset(recsumm, select = c(RRID, RRID.Sex, TotalRecombCount))


pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)
pedigree <- pedigree[,c(1, 3, 2)]
for(i in 1:3) pedigree[which(pedigree[,i] == 0),i] <- NA
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])
ainv <- asreml.Ainverse(pedigree)$ginv

recsumm$RRID <- as.character(recsumm$RRID)

source("r/ASReml.EstEffects.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Test each allele against the other in the top haplotype   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

regh2.snp.start <- which(haplo.map$SNP.Name == "cela1_red_10_20476277")

topdata <- haplo.res
topdata$HaploTest <- sapply(topdata$Haplo, function(x) substring(x, regh2.snp.start, regh2.snp.start + 9))
topdata <- topdata[,-3]
head(topdata)

table(topdata$HaploTest)

temp <- vector2Dendrogram(topdata$HaploTest, count.cutoff = 10)
haplo.vec <- temp$haplotype.counts$Haplo


top.fixef <- NULL
top.tab <- data.frame(Haplotype = haplo.vec,
                      LogLi = NA,
                      P = NA)

for(i in 1:length(haplo.vec)){
  
  tempdata <- topdata
  #tempdata <- subset(topdata, ID != 871)
  tempdata$HaploTest <- ifelse(tempdata$HaploTest == haplo.vec[i], "B", "A")
  tempdata <- cast(tempdata, ID ~ HaploID)
  tempdata$Geno1 <- apply(tempdata[,2:3], 1, function(x) sort(x)[1])
  tempdata$Geno2 <- apply(tempdata[,2:3], 1, function(x) sort(x)[2])
  tempdata$Genotype <- paste0(tempdata$Geno1, "/", tempdata$Geno2)
  tempdata <- subset(tempdata, select = c(ID, Genotype))
  names(tempdata)[1] <- "RRID"
  
  head(tempdata)
  
  tempdata <- join(recsumm, tempdata)
  
  tempdata$RRID <- as.factor(tempdata$RRID)
  tempdata$RRID.Sex <- as.factor(tempdata$RRID.Sex)
  tempdata <- subset(tempdata, RRID.Sex == "F")
  
  ped.summ.RR.f <- asreml(fixed = TotalRecombCount ~ Genotype,
                          random = ~ ped(RRID),
                          data = tempdata,
                          ginverse =  list(RRID = ainv),
                          na.method.X = "omit", na.omit.Y = "na.omit",
                          workspace = 500e+6, pworkspace = 500e+6, trace = F)
  
  x <- data.frame(summary(ped.summ.RR.f, all = T)$coef.fixed)
  x$Genotype <- row.names(x)
  x$Ref <- haplo.vec[i]
  
  y <-  data.frame(table(tempdata$Genotype))
  names(y)[1] <- "Genotype"
  y$Genotype <- paste0("Genotype_", y$Genotype)
  
  x <- join(x, y)
  x
  
  top.fixef <- rbind(top.fixef, x)
  top.tab$LogLi[i] <- summary(ped.summ.RR.f, all = T)$loglik
  top.tab$P[i] <- data.frame(wald.asreml(ped.summ.RR.f))["Genotype", "Pr.Chisq."]
  
  rm(x, tempdata, ped.summ.RR.f)
  
}


toptophits <- subset(top.fixef, !Genotype %in% c("(Intercept)", "Genotype", "Genotype_A/A", "RRID.Sex_F", "RRID.Sex_M"))
arrange(toptophits, z.ratio)
x <- subset(top.fixef, Genotype == "Genotype_A/A")
x <- subset(x, select = c(Ref, Freq))
names(x)[2] <- "FreqAA"
head(x)

toptophits <- join(toptophits, x)
toptophits$FreqAA[which(toptophits$Genotype == "Genotype_B/B")] <- NA

toptophits$yval <- ifelse(toptophits$Genotype == "Genotype_A/B", 13, 12)
toptophits <- subset(toptophits, Freq > 5)
toptophits$Sig <- ifelse(toptophits$z.ratio < -1.96, "*",
                         ifelse(toptophits$z.ratio > 1.96, "*", ""))

toptophits$Genotype <- gsub("Genotype_", "", toptophits$Genotype)

ggplot(toptophits, aes(Ref, solution, fill = Genotype, ymin = solution - std.error, ymax = solution + std.error)) +
  geom_text(aes(x = Ref, y = 14, label = FreqAA), size = 4) +
  geom_text(aes(x = Ref, y = yval, label = Freq, colour = factor(yval, levels = c(13, 12))), size = 4, show.legend = F) +
  geom_hline(yintercept = 0) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_text(aes(x = Ref, y = 11, label = Sig), position=position_dodge(width = 0.9)) +
  geom_errorbar(position=position_dodge(width = 0.9), width = 0) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(breaks = seq(0, 16, 2)) +
  labs(x = "Haplotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12, angle = 90),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())

ggsave("figs/6_Haplotype_Effect_Sizes_Females.png", width = 10, height = 6)

#~~ Make Latex table
tophits.save <- toptophits
tophits.save[1:3] <- round(tophits.save[1:3], digits = 3)
tophits.save$FreqAA[which(is.na(tophits.save$FreqAA))] <- tophits.save$FreqAA[which(is.na(tophits.save$FreqAA))-1]
tophits.save$yval <- NULL
tophits.save <- tophits.save[,c(5, 4, 1:3, 6:8)]
names(tophits.save)[1] <- "Haplotype"
tophits.save$`""` <- "\\\\"
write.table(tophits.save, "doc/tables/Fixed_Effect_Sizes_Haplos.txt", sep = " & ", row.names = F, quote = F)



i = which(haplo.vec == "AGGAGAGAAG")

tempdata <- topdata
tempdata$HaploTest <- ifelse(tempdata$HaploTest == haplo.vec[i], "B", "A")
tempdata <- cast(tempdata, ID ~ HaploID)
tempdata$Geno1 <- apply(tempdata[,2:3], 1, function(x) sort(x)[1])
tempdata$Geno2 <- apply(tempdata[,2:3], 1, function(x) sort(x)[2])
tempdata$Genotype <- paste0(tempdata$Geno1, "/", tempdata$Geno2)
tempdata <- subset(tempdata, select = c(ID, Genotype))
names(tempdata)[1] <- "RRID"

head(tempdata)

tempdata <- join(recsumm, tempdata)

tempdata$RRID <- as.factor(tempdata$RRID)

tempdata$RRID.Sex <- as.character(tempdata$RRID.Sex)
tempdata$RRID.Sex[which(tempdata$RRID.Sex == "F")] <- "Female"
tempdata$RRID.Sex[which(tempdata$RRID.Sex == "M")] <- "Male"

x <- subset(tempdata, select = -TotalRecombCount)

x1 <- data.frame(table(x$RRID.Sex, x$Genotype))
x2 <- data.frame(table(unique(x)$RRID.Sex, unique(x)$Genotype))
names(x2)[3] <- "IDFreq"

x1 <- join(x1, x2)
x1$IDFreq <- paste0("(", x1$IDFreq, ")")
names(x1)[1:2] <- c("RRID.Sex", "Genotype")

ggplot(tempdata, aes(Genotype, TotalRecombCount)) +
  geom_text(data = x1, aes(x = Genotype, y = 51, label = Freq)) +
  geom_text(data = x1, aes(x = Genotype, y = 49, label = IDFreq)) +
  geom_boxplot(notch= T) +
  facet_wrap(~RRID.Sex) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "Haplotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())

tapply(tempdata$TotalRecombCount, list(tempdata$Genotype, tempdata$RRID.Sex), mean)
tapply(tempdata$TotalRecombCount, list(tempdata$Genotype, tempdata$RRID.Sex), median)


ggsave("figs/6_Haplotype_Effect_Sizes_AGGAGAGAAG.png", width = 10, height = 6)

x2 <- x1
x2$IDFreq[5] <- "(3)"
x2$Freq[5] <- 7

subset(tempdata, RRID =  871)

ggplot(subset(tempdata, RRID !=  871), aes(Genotype, TotalRecombCount)) +
  geom_text(data = x2, aes(x = Genotype, y = 51, label = Freq)) +
  geom_text(data = x2, aes(x = Genotype, y = 49, label = IDFreq)) +
  geom_boxplot(notch= T) +
  facet_wrap(~RRID.Sex) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "Haplotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())

tapply(tempdata$TotalRecombCount, list(tempdata$Genotype, tempdata$RRID.Sex), mean)
tapply(tempdata$TotalRecombCount, list(tempdata$Genotype, tempdata$RRID.Sex), median)

ggsave("figs/6_Haplotype_Effect_Sizes_AGGAGAGAAG_wo871.png", width = 10, height = 6)


#~~ Get for other haplo

i = which(haplo.vec == "AGAGAAGAGA")

tempdata2 <- topdata
tempdata2$HaploTest <- ifelse(tempdata2$HaploTest == haplo.vec[i], "B", "A")
tempdata2 <- cast(tempdata2, ID ~ HaploID)
tempdata2$Geno1 <- apply(tempdata2[,2:3], 1, function(x) sort(x)[1])
tempdata2$Geno2 <- apply(tempdata2[,2:3], 1, function(x) sort(x)[2])
tempdata2$Genotype <- paste0(tempdata2$Geno1, "/", tempdata2$Geno2)
tempdata2 <- subset(tempdata2, select = c(ID, Genotype))
names(tempdata2)[1] <- "RRID"

head(tempdata2)

tempdata2 <- join(recsumm, tempdata2)

tempdata2$RRID <- as.factor(tempdata2$RRID)

tempdata2$RRID.Sex <- as.character(tempdata2$RRID.Sex)
tempdata2$RRID.Sex[which(tempdata2$RRID.Sex == "F")] <- "Female"
tempdata2$RRID.Sex[which(tempdata2$RRID.Sex == "M")] <- "Male"


#~~ Get for top SNP

top.snp.all <- which(haplo.map$SNP.Name == "cela1_red_10_26005249")

tempdata.snp <- haplo.res
tempdata.snp$HaploTest <- sapply(tempdata.snp$Haplo, function(x) substring(x, top.snp.all, top.snp.all))  #324, 329
tempdata.snp <- tempdata.snp[,-3]
head(tempdata.snp)
tempdata.snp <- cast(tempdata.snp, ID ~ HaploID)
tempdata.snp$Geno1 <- apply(tempdata.snp[,2:3], 1, function(x) sort(x)[1])
tempdata.snp$Geno2 <- apply(tempdata.snp[,2:3], 1, function(x) sort(x)[2])
tempdata.snp$Genotype <- paste0(tempdata.snp$Geno1, "/", tempdata.snp$Geno2)
tempdata.snp <- subset(tempdata.snp, select = c(ID, Genotype))
names(tempdata.snp)[1] <- "RRID"

head(tempdata.snp)

tempdata.snp <- join(recsumm, tempdata.snp)

tempdata.snp$RRID <- as.factor(tempdata.snp$RRID)
tempdata.snp$RRID.Sex <- as.factor(tempdata.snp$RRID.Sex)

x <- subset(tempdata.snp, select = -TotalRecombCount)

x1 <- data.frame(table(x$RRID.Sex, x$Genotype))
x2 <- data.frame(table(unique(x)$RRID.Sex, unique(x)$Genotype))
names(x2)[3] <- "IDFreq"

x1 <- join(x1, x2)
x1$IDFreq <- paste0("(", x1$IDFreq, ")")
names(x1)[1:2] <- c("RRID.Sex", "Genotype")
x1

ggplot(tempdata.snp, aes(Genotype, TotalRecombCount)) +
  geom_text(data = x1, aes(x = Genotype, y = 51, label = Freq)) +
  geom_text(data = x1, aes(x = Genotype, y = 49, label = IDFreq)) +
  geom_boxplot(notch= T) +
  facet_wrap(~RRID.Sex) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "Haplotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())


tapply(tempdata.snp$TotalRecombCount, list(tempdata.snp$Genotype, tempdata.snp$RRID.Sex), mean)
tapply(tempdata.snp$TotalRecombCount, list(tempdata.snp$Genotype, tempdata.snp$RRID.Sex), median)

names(tempdata.snp)[4] <- "SNP.Genotype"
tempdata.snp <- subset(tempdata.snp, select = c(RRID, SNP.Genotype))
tempdata.snp <- unique(tempdata.snp)                       


#~~ Get for top SNP

top.snp.female <- which(haplo.map$SNP.Name == "cela1_red_10_25661750")

tempdata.snp2 <- haplo.res
tempdata.snp2$HaploTest <- sapply(tempdata.snp2$Haplo, function(x) substring(x, top.snp.female, top.snp.female))  #324, 329
tempdata.snp2 <- tempdata.snp2[,-3]
head(tempdata.snp2)
tempdata.snp2 <- cast(tempdata.snp2, ID ~ HaploID)
tempdata.snp2$Geno1 <- apply(tempdata.snp2[,2:3], 1, function(x) sort(x)[1])
tempdata.snp2$Geno2 <- apply(tempdata.snp2[,2:3], 1, function(x) sort(x)[2])
tempdata.snp2$Genotype <- paste0(tempdata.snp2$Geno1, "/", tempdata.snp2$Geno2)
tempdata.snp2 <- subset(tempdata.snp2, select = c(ID, Genotype))
names(tempdata.snp2)[1] <- "RRID"

head(tempdata.snp2)

tempdata.snp2 <- join(recsumm, tempdata.snp2)

tempdata.snp2$RRID <- as.factor(tempdata.snp2$RRID)
tempdata.snp2$RRID.Sex <- as.factor(tempdata.snp2$RRID.Sex)

x <- subset(tempdata.snp2, select = -TotalRecombCount)

x1 <- data.frame(table(x$RRID.Sex, x$Genotype))
x2 <- data.frame(table(unique(x)$RRID.Sex, unique(x)$Genotype))
names(x2)[3] <- "IDFreq"

x1 <- join(x1, x2)
x1$IDFreq <- paste0("(", x1$IDFreq, ")")
names(x1)[1:2] <- c("RRID.Sex", "Genotype")
x1

ggplot(tempdata.snp2, aes(Genotype, TotalRecombCount)) +
  geom_text(data = x1, aes(x = Genotype, y = 51, label = Freq)) +
  geom_text(data = x1, aes(x = Genotype, y = 49, label = IDFreq)) +
  geom_boxplot(notch= T) +
  facet_wrap(~RRID.Sex) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "Haplotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())


tapply(tempdata.snp2$TotalRecombCount, list(tempdata.snp2$Genotype, tempdata.snp2$RRID.Sex), mean)
tapply(tempdata.snp2$TotalRecombCount, list(tempdata.snp2$Genotype, tempdata.snp2$RRID.Sex), median)

names(tempdata.snp2)[4] <- "SNP2.Genotype"
tempdata.snp2 <- subset(tempdata.snp2, select = c(RRID, SNP2.Genotype))
tempdata.snp2 <- unique(tempdata.snp2)                       

names(tempdata2)[4] <- "Haplo2.Genotype"
tempdata2 <- subset(tempdata2, select = c(RRID, Haplo2.Genotype))
tempdata2 <- unique(tempdata2)                       


rm(x, y, x1, x2, temp, i, regh2.snp.start, top.snp.all, top.snp.female)

#~~ Make a combined dataset

combined.data <- join(tempdata, tempdata.snp)

combined.data <- join(combined.data, tempdata.snp2)

combined.data <- join(combined.data, tempdata2)


combined.data$Combined.Genotype  <- paste0(combined.data$Genotype, "_", combined.data$SNP.Genotype)
combined.data$Combined.Genotype2 <- paste0(combined.data$Genotype, "_", combined.data$SNP2.Genotype)

head(combined.data)

table(combined.data$Combined.Genotype)

cor.test(as.numeric(as.factor(combined.data$SNP.Genotype)), as.numeric(as.factor(combined.data$Genotype)))
cor.test(as.numeric(as.factor(combined.data$SNP.Genotype)), as.numeric(as.factor(combined.data$Haplo2.Genotype)))
cor.test(as.numeric(as.factor(combined.data$SNP2.Genotype)), as.numeric(as.factor(combined.data$Genotype)))
cor.test(as.numeric(as.factor(combined.data$SNP2.Genotype)), as.numeric(as.factor(combined.data$Haplo2.Genotype)))



model1 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * Genotype,
                      random = ~ ped(RRID),
                      data = combined.data,
                      ginverse =  list(RRID = ainv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6, trace = F)


model2 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * SNP.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)



model3 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * SNP2.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)



model4 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * Genotype + RRID.Sex * SNP.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)



model5 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * Genotype + RRID.Sex * SNP2.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)



model6 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * Combined.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)



model7 <- asreml(fixed = TotalRecombCount ~ RRID.Sex * Combined.Genotype2,
                 random = ~ ped(RRID),
                 data = combined.data,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit",
                 workspace = 500e+6, pworkspace = 500e+6, trace = F)


final.fixef <- NULL

for(i in 1:7){
  
  eval(parse(text = paste0("x <- data.frame(summary(model", i, ", all = T)$coef.fixed)")))
  x$Genotype <- row.names(x)
  x$Ref <- haplo.vec[i]
  x$Model <- paste0("model", i)
  
  final.fixef <- rbind(final.fixef, x)
  
  rm(x)

}

row.names(final.fixef) <- 1:nrow(final.fixef)

final.sampsizes <- NULL

for(i in names(combined.data)[4:8]){
  
  eval(parse(text = paste0("x <- subset(combined.data, select = c(RRID, RRID.Sex, ", i, "))")))
  eval(parse(text = paste0("y <- data.frame(table(x$RRID.Sex, x$", i, "))")))
  eval(parse(text = paste0("x <- unique(subset(combined.data, select = c(RRID, RRID.Sex, ", i, ")))")))
  eval(parse(text = paste0("y1 <- data.frame(table(x$RRID.Sex, x$", i, "))")))
  names(y1)[3] <- "IDFreq"
  y <- join(y, y1)
  y$Marker <- i
  
  final.sampsizes <- rbind(final.sampsizes, y)

}

final.sampsizes$Genotype <- paste0(final.sampsizes$Marker, "_", final.sampsizes$Var2)
names(final.sampsizes)[1] <- "RRID.Sex"


#~~ Sort out the two tables

head(final.fixef)
write.table(final.fixef, "test.txt", row.names = F, sep = "\t", quote = F)

final.fixef$RRID.Sex <- NA
final.fixef$RRID.Sex[grep("^Genotype", final.fixef$Genotype)] <- "Female"
final.fixef$RRID.Sex[grep("^SNP"     , final.fixef$Genotype)] <- "Female"
final.fixef$RRID.Sex[grep("^Combined", final.fixef$Genotype)] <- "Female"

final.fixef$RRID.Sex[grep("^RRID.Sex_Male:", final.fixef$Genotype)] <- "Male"
final.fixef <- final.fixef[-grep("^RRID.Sex_Female:", final.fixef$Genotype),]

final.fixef$Genotype <- gsub("^RRID.Sex_Male:", "", final.fixef$Genotype)

final.fixef <- join(final.fixef, final.sampsizes)


ggplot(subset(final.fixef, !is.na(RRID.Sex)), aes(Var2, solution)) +
  geom_point() +
  facet_grid(Model~RRID.Sex, scales = "free") +
  theme_bw() +
  labs(x = "SNP Genotype", y = "Effect Size") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank())


head(combined.data)

combined.melt <- melt(combined.data, id.vars = c("RRID", "RRID.Sex", "TotalRecombCount"))

head(combined.melt)

head(final.sampsizes)
names(combined.melt)[c(4, 5)] <- c("Marker", "Var2")

combined.melt <- join(combined.melt, final.sampsizes)

combined.melt$RealMarker <- ifelse(combined.melt$Marker == "Genotype", "Haplotype",
                                   ifelse(combined.melt$Marker == "SNP.Genotype", "cela1_red_10_26005249",
                                          ifelse(combined.melt$Marker == "SNP2.Genotype", "cela1_red_10_25661750",
                                                 ifelse(combined.melt$Marker == "Combined.Genotype", "Haplotype + cela1_red_10_26005249",
                                                        "Haplotype + cela1_red_10_25661750"))))

combined.melt$IDFreq2 <- paste0("(",combined.melt$IDFreq, ")")


x <- subset(combined.melt, Marker %in% c("Genotype"))
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))
ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 51, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 49, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))

x <- subset(combined.melt, Marker %in% c("SNP.Genotype"))
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 51, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 49, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))

x <- subset(combined.melt, Marker %in% c("SNP2.Genotype"))
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 51, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 49, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))


x <- subset(combined.melt, Marker %in% c("Combined.Genotype"))
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 51, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 49, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12, angle = 270),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))

x <- subset(combined.melt, Marker %in% c("Combined.Genotype2"))
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 51, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 49, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12, angle = 270),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))







x <- droplevels(subset(combined.melt, RealMarker %in% c("cela1_red_10_26005249", "Haplotype")))
x$Var2[which(x$Var2 == "A/G")] <- "A/B"
x$Var2[which(x$Var2 == "G/G")] <- "B/B"
x <- droplevels(x)
head(x)
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

head(y)

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 55, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 50, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))

top.fixef


x <- droplevels(subset(combined.melt, Freq >= 20 & RealMarker %in% c("cela1_red_10_25661750", "cela1_red_10_26005249", "Haplotype")))
x$Var2 <- as.character(x$Var2)
x$Var2[which(x$Var2 == "A/G")] <- "A/B"
x$Var2[which(x$Var2 == "G/G")] <- "B/B"
x <- droplevels(x)
head(x)
y <- unique(subset(x, select = c(RRID.Sex, Freq, IDFreq2, Var2, RealMarker)))

head(y)

ggplot(x, aes(Var2, TotalRecombCount)) +
  geom_text(data = y, aes(x = Var2, 55, label = Freq)) +
  geom_text(data = y, aes(x = Var2, 50, label = IDFreq2)) +
  geom_boxplot(notch = T) +
  facet_grid(RealMarker~RRID.Sex, scales = "free_x") +
  theme_bw() +
  labs(x = "Genotype", y = "Autosomal Crossover Count") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12))

x

#~~ fixef...

head(final.fixef)
final.fixef <- subset(final.fixef, Model %in% c("model1", "model2", "model3"))
final.fixef$RealMarker <- NA
final.fixef$RealMarker[which(final.fixef$Model == "model1")] <- "Haplotype"
final.fixef$RealMarker[which(final.fixef$Model == "model2")] <- "cela1_red_10_26005249"
final.fixef$RealMarker[which(final.fixef$Model == "model3")] <- "cela1_red_10_25661750"

final.fixef <- subset(final.fixef, select = -c(Ref, Marker))

x <- NULL

for(i in unique(final.fixef$RealMarker)){
  y <- subset(final.fixef, RealMarker == i)
  y$solution[which(y$RRID.Sex == "Male")] <- y$solution[which(y$RRID.Sex == "Male")] + y$solution[which(y$Genotype == "RRID.Sex_Male")]+ y$solution[which(y$Genotype == "(Intercept)")]
  y$solution[which(y$RRID.Sex == "Female")] <- y$solution[which(y$RRID.Sex == "Female")] + y$solution[which(y$Genotype == "(Intercept)")]
  y <- subset(y, !is.na(RRID.Sex))
  x <- rbind(x, y)
  rm(y)
}

x

ggplot(x, aes(Var2, solution)) +
  geom_point() +
  geom_errorbar(aes(ymax = solution + std.error, ymin = solution - std.error), width = 0) +
  facet_grid(RealMarker ~ RRID.Sex, scales = "free")


#~~ Run sex specific models

combined.data.m <- droplevels(subset(combined.data, RRID.Sex == "Male"))
combined.data.f <- droplevels(subset(combined.data, RRID.Sex == "Female"))

model1.m <- asreml(fixed = TotalRecombCount ~ Genotype,
                 random = ~ ped(RRID),
                 data = combined.data.m,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit")

model2.m <- asreml(fixed = TotalRecombCount ~ SNP.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data.m,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit")

model3.m <- asreml(fixed = TotalRecombCount ~ SNP2.Genotype,
                 random = ~ ped(RRID),
                 data = combined.data.m,
                 ginverse =  list(RRID = ainv),
                 na.method.X = "omit", na.omit.Y = "na.omit")


model4.m <- asreml(fixed = TotalRecombCount ~ Haplo2.Genotype,
                   random = ~ ped(RRID),
                   data = combined.data.m,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit")



model1.f <- asreml(fixed = TotalRecombCount ~ Genotype,
                   random = ~ ped(RRID),
                   data = combined.data.f,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit")

model2.f <- asreml(fixed = TotalRecombCount ~ SNP.Genotype,
                   random = ~ ped(RRID),
                   data = combined.data.f,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit")

model3.f <- asreml(fixed = TotalRecombCount ~ SNP2.Genotype,
                   random = ~ ped(RRID),
                   data = combined.data.f,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit")


model4.f <- asreml(fixed = TotalRecombCount ~ Haplo2.Genotype,
                   random = ~ ped(RRID),
                   data = combined.data.f,
                   ginverse =  list(RRID = ainv),
                   na.method.X = "omit", na.omit.Y = "na.omit")



new.fixef <- NULL

for(i in c("1.m", "1.f", "2.m", "2.f", "3.m", "3.f", "4.m", "4.f")){
  
  eval(parse(text = paste0("x <- data.frame(summary(model", i, ", all = T)$coef.fixed)")))
  x$Genotype <- row.names(x)
  x$Ref <- haplo.vec[i]
  x$Model <- paste0("model", i)
  
  new.fixef <- rbind(new.fixef, x)
  
  rm(x)
  
}

row.names(new.fixef) <- 1:nrow(new.fixef)

new.fixef$RRID.Sex <- NA
new.fixef$RRID.Sex[grep(".m", new.fixef$Model, fixed = F)] <- "Male"
new.fixef$RRID.Sex[grep(".f", new.fixef$Model, fixed = F)] <- "Female"

new.fixef$Var2 <- sapply(new.fixef$Genotype, function(x) strsplit(x, split = "_")[[1]][2])
new.fixef$Var2[which(is.na(new.fixef$Var2))] <- "A/A (Intercept)"
new.fixef <- new.fixef[-which(new.fixef$Var2 == "A/A"),]
new.fixef <- arrange(new.fixef, Model, Var2)

new.fixef$RealMarker <- NA
new.fixef$RealMarker[grep("model1", new.fixef$Model)] <- "Haplotype AGGAGAGAAG"
new.fixef$RealMarker[grep("model2", new.fixef$Model)] <- "cela1_red_10_26005249"
new.fixef$RealMarker[grep("model3", new.fixef$Model)] <- "cela1_red_10_25661750"
new.fixef$RealMarker[grep("model4", new.fixef$Model)] <- "Haplotype AGAGAAGAGA"

new.fixef <- subset(new.fixef, select = -c(Genotype, Ref))

y <- unique(subset(final.fixef, select = c(RealMarker, Var2, RRID.Sex, Freq, IDFreq)))
y <- droplevels(subset(y, !is.na(Var2)))
y$Var2 <- as.character(y$Var2)
y$Var2[which(y$Var2 == "A/A")] <- "A/A (Intercept)"

y$RealMarker[which(y$RealMarker == "Haplotype")] <- "Haplotype AGGAGAGAAG"

temptab <- subset(combined.data, select = c(RRID, RRID.Sex, Haplo2.Genotype))
x <- data.frame(table(temptab$RRID.Sex, temptab$Haplo2.Genotype))
x1 <- data.frame(table(unique(temptab)$RRID.Sex, unique(temptab)$Haplo2.Genotype))
names(x1)[3] <- "IDFreq" 
x <- join(x, x1)
names(x)[1] <- c("RRID.Sex")
x$Var2 <- as.character(x$Var2)
x$Var2[which(x$Var2 == "A/A")] <- "A/A (Intercept)"
x$RealMarker <- "Haplotype AGAGAAGAGA"
y <- rbind(y, x)

new.fixef <- join(new.fixef, y)

new.fixef$Wald.P <- NA

new.fixef$Wald.P[which(new.fixef$Model == "model1.m")] <- wald.asreml(model1.m)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model1.f")] <- wald.asreml(model1.f)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model2.m")] <- wald.asreml(model2.m)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model2.f")] <- wald.asreml(model2.f)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model3.m")] <- wald.asreml(model3.m)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model3.f")] <- wald.asreml(model3.f)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model4.m")] <- wald.asreml(model4.m)$`Pr(Chisq)`[2]
new.fixef$Wald.P[which(new.fixef$Model == "model4.f")] <- wald.asreml(model4.f)$`Pr(Chisq)`[2]


new.fixef <- subset(new.fixef, select = -Model)
new.fixef <- new.fixef[,c("RealMarker", "RRID.Sex", "Var2",  
                          "Freq", "IDFreq", "solution", "std.error", "z.ratio", "Wald.P")]

new.fixef$RealMarker[which(new.fixef$RRID.Sex == "Male")] <- ""
new.fixef$RealMarker[which(new.fixef$RRID.Sex == "Female" & new.fixef$Var2 != "A/A (Intercept)")] <- ""
new.fixef$RRID.Sex[which(new.fixef$Var2 != "A/A (Intercept)")] <- ""
new.fixef$Wald.P <- formatC(new.fixef$Wald.P, format = "e", digits = 3)

new.fixef$Wald.P[which(new.fixef$Var2 != "A/A (Intercept)")] <- ""

#new.fixef$RealMarker[2] <- "(AGGAGAGAAG)"
names(new.fixef)[1:5] <- c("Locus", "Sex", "Genotype", "Count", "ID.Count")

new.fixef[6:8] <- round(new.fixef[6:8], digits = 3)

new.fixef$`""` <- "\\\\"

write.table(new.fixef, file = "doc/tables/Fixed_Effect_Sizes.txt", sep = " & ", quote = F, row.names = F)
