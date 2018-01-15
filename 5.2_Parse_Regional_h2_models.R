#
# Parse the regional heritability models
# Susan Johnston
# July 2017
#
# 

library(plyr)
library(asreml)
library(ggplot2)
library(asreml)

source("r/ASReml.EstEffects.R")


AnalysisSuffix <- "a"

mapdata   <- read.table("gcta/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)

regh2.tab <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/regh2_master_", AnalysisSuffix, ".txt"), header = T)
regh2.tab1 <- read.table(paste0("gcta/regh2_", AnalysisSuffix, "/regh2_master_", AnalysisSuffix, "_1.txt"), header = T)

regh2.tab <- rbind(regh2.tab, regh2.tab1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Format the region file and add map information  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


mapdata$Cumu.Order <- 1:nrow(mapdata)
mapdata$Diff <- c(mapdata$Estimated.Mb.Position[1],diff(mapdata$Estimated.Mb.Position))
mapdata$Diff <- ifelse(mapdata$Diff < 0, 10000, mapdata$Diff)
mapdata$Cumu.Position <- cumsum(mapdata$Diff)


regh2.tab <- subset(regh2.tab, CEL.LG != 34)
names(regh2.tab)[3] <- "Start.Order"
regh2.tab$Stop.Order <- regh2.tab$Start.Order + regh2.tab$Window.Size - 1

temp <- mapdata[,c("CEL.order", "Cumu.Order", "Cumu.Position", "CEL.LG", "Estimated.Mb.Position", "BTA.Position")]

names(temp) <- c("Start.Order", "Start.Cumu.Order", "Start.Cumu.Position", "CEL.LG", "Start.Position", "Start.BTA.Position")
regh2.tab <- join(regh2.tab, temp)
names(temp) <- c("Stop.Order", "Stop.Cumu.Order", "Stop.Cumu.Position", "CEL.LG", "Stop.Position", "Stop.BTA.Position")
regh2.tab <- join(regh2.tab, temp)


regh2.tab$Window.Bp <- regh2.tab$Stop.Position - regh2.tab$Start.Position

tapply(regh2.tab$Window.Bp, regh2.tab$Window.Size, mean)
tapply(regh2.tab$Window.Bp, regh2.tab$Window.Size, sd)


head(regh2.tab)
head(mapdata)

rm(temp)

tail(regh2.tab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse the likelihood results                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

restab <- list()
vartab <- list()

for(j in 1:nrow(regh2.tab)){
  
  if(j %in% seq(1, nrow(regh2.tab), 50)) print(paste("Running row", j, "of", nrow(regh2.tab)))
  
  if(file.exists(paste0("gcta/regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[j], ".RData"))){
    
    load(paste0("gcta/regh2_", AnalysisSuffix, "/", regh2.tab$Analysis.ID[j], ".RData"))
    
    model1 <- x[[1]][1][[1]]
    model2 <- x[[2]][1][[1]]
    model1.m <- x[[3]][1][[1]]
    model2.m <- x[[4]][1][[1]]
    model1.f <- x[[5]][1][[1]]
    model2.f <- x[[6]][1][[1]]
    
    
    
    restab[[j]] <- data.frame(Analysis.ID = regh2.tab$Analysis.ID[j],
                              RRID.Sex     = c("Both", "Male", "Female"),
                              model1.Li    = c(model1$loglik, model1.m$loglik, model1.f$loglik),
                              model2.Li    = c(model2$loglik, model2.m$loglik, model2.f$loglik))
    
    vartab[[j]] <- rbind(cbind(x[[1]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Both"  , Model = "wo.region", Effect = row.names(x[[1]][3][[1]])),
                         cbind(x[[2]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Both"  , Model = "wi.region", Effect = row.names(x[[2]][3][[1]])),
                         cbind(x[[3]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Male"  , Model = "wo.region", Effect = row.names(x[[3]][3][[1]])),
                         cbind(x[[4]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Male"  , Model = "wi.region", Effect = row.names(x[[4]][3][[1]])),
                         cbind(x[[5]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Female", Model = "wo.region", Effect = row.names(x[[5]][3][[1]])),
                         cbind(x[[6]][3][[1]], Analysis.ID = regh2.tab$Analysis.ID[j], RRID.Sex = "Female", Model = "wi.region", Effect = row.names(x[[6]][3][[1]])))


    rm(model1, model1.f, model1.m, model2, model2.f, model2.m, x)

  }
  
}

rm(j)

#~~ Convert lists to frames

restab1 <- data.table::rbindlist(restab[1:5000])
restab2 <- data.table::rbindlist(restab[5001:10000])
restab3 <- data.table::rbindlist(restab[10001:nrow(regh2.tab)])

restab <- rbind(restab1, restab2, restab3)

gc()


vartab1 <- data.table::rbindlist(vartab[1:5000])
vartab2 <- data.table::rbindlist(vartab[5001:10000])
vartab3 <- data.table::rbindlist(vartab[10001:nrow(regh2.tab)])


vartab <- rbind(vartab1, vartab2, vartab3)


rm(restab1, restab2, restab3, vartab1, vartab2, vartab3)

#~~ Join with region information

restab <- join(restab, regh2.tab)
vartab <- join(vartab, regh2.tab)

#~~ Get test statistics from likelihood ratio test, chidist with 1df

restab$Chi2 <- 2*(restab$model2.Li - restab$model1.Li)
restab$P <- 1- pchisq(restab$Chi2, df = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Examine associated regions                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(restab)

restab$Mid.Position <- restab$Start.Position + 0.5*(restab$Stop.Position - restab$Start.Position)
restab$Mid.Cumu.Position <- restab$Start.Cumu.Position + 0.5*(restab$Stop.Cumu.Position - restab$Start.Cumu.Position)
restab$Mid.BTA.Position <- restab$Start.BTA.Position + 0.5*(restab$Stop.BTA.Position - restab$Start.BTA.Position)


#~~ Get chromosome info

chrinfo <- NULL

for(i in na.omit(unique(restab$CEL.LG))){
  
  temp1 <- arrange(subset(restab, CEL.LG == i), Mid.Cumu.Position)
  
  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Mid.Cumu.Position"],
                      Stop = temp1[nrow(temp1),"Mid.Cumu.Position"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

names(chrinfo) <- c("CEL.LG", "Start", "Stop")

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
chrinfo$CEL.LG2 <- c(1, "", "", 4, "", "", 7, "", "", 10, "", "", 13, "", "", 16, "", "", 19, "", "", 22, "", "", "", "", 27, "", "", "", 31, "", "")



bonf50 = 0.05/(length(unique(subset(regh2.tab, Window.Size == 50)$Analysis.ID))/2)
bonf20 = 0.05/(length(unique(subset(regh2.tab, Window.Size == 20)$Analysis.ID))/2)
bonf10 = 0.05/(length(unique(subset(regh2.tab, Window.Size == 10)$Analysis.ID))/2)


#~~ Plots


ggplot(subset(restab, Window.Size == 50),
       aes(Mid.Cumu.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) +
  geom_line() +
  geom_hline(yintercept=-log10(bonf50),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
  labs(x ="Chromosome", y = "-log10 P")

ggsave("figs/5_Regional_Heritability_window_50.png", width = 8, height = 10)

x <- subset(restab, Window.Size == 20 & CEL.LG == 6)
tail(x)

ggplot(subset(restab, Window.Size == 20),
       aes(Mid.Cumu.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) +
  geom_line() +
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
  labs(x ="Chromosome", y = "-log10 P")

ggsave("figs/5_Regional_Heritability_window_20.png", width = 8, height = 10)


ggplot(subset(restab, Window.Size == 10),
       aes(Mid.Cumu.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_line() +
  geom_hline(yintercept=-log10(bonf10),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
  labs(x ="Chromosome", y = "-log10 P")

ggsave("figs/5_Regional_Heritability_window_10.png", width = 8, height = 10)

#~~ CEL12

ggplot(subset(restab, Window.Size == 10 & CEL.LG == 12),
       aes(Mid.BTA.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_line() +
  geom_point() +
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  labs(x ="Chromosome", y = "-log10 P")

ggsave("figs/5_Regional_Heritability_window_10_CEL.12.png", width = 8, height = 10)

subset(vartab, Analysis.ID == arrange(subset(restab, CEL.LG == 12), P)$Analysis.ID[1])[,1:11]

#~~ CEL32

ggplot(subset(restab, Window.Size == 20 & CEL.LG == 32),
       aes(Mid.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_line() +
  geom_point() +
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  labs(x ="Chromosome", y = "-log10 P")
ggsave("figs/5_Regional_Heritability_window_20_CEL.32.png", width = 8, height = 10)

subset(vartab, Analysis.ID == arrange(subset(restab, CEL.LG == 32), P)$Analysis.ID[1])[,1:11]

#~~ CEL16

ggplot(subset(restab, Window.Size == 10 & CEL.LG == 16),
       aes(Mid.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_line() +
  geom_point() +
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  labs(x ="Chromosome", y = "-log10 P")
ggsave("figs/5_Regional_Heritability_window_10_CEL.16.png", width = 8, height = 10)

subset(vartab, Analysis.ID == arrange(subset(restab, CEL.LG == 16), P)$Analysis.ID[1])[,1:11]

#~~ CEL30

ggplot(subset(restab, Window.Size == 20 & CEL.LG == 30),
       aes(Mid.Position, -log10(P), group = CEL.LG, col = factor(CEL.LG %% 2))) + 
  geom_line() +
  geom_point() +
  geom_hline(yintercept=-log10(bonf20),linetype=2, alpha = 0.6, size = 1) +
  facet_wrap(~RRID.Sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 14),
        axis.title.y = element_text (size = 12, angle = 90),
        axis.title.x = element_text (size = 12),
        strip.background = element_blank()) +
  labs(x ="Chromosome", y = "-log10 P")
ggsave("figs/5_Regional_Heritability_window_20_CEL.30.png", width = 8, height = 10)

subset(vartab, Analysis.ID == arrange(subset(restab, CEL.LG == 30), P)$Analysis.ID[1])[,1:11]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Write to file                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

write.table(restab, paste0("results/5_Regional_Heritability_Results_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(vartab, paste0("results/5_Regional_Heritability_VarComp_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")




