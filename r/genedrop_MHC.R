#
# Basic Genedrop Analysis
#
# February 2017
#

#~~ Read in the genedrop function. To run, it requires the libraries plyr, reshape2 and kinship2. This
#   is for a single locus only and its quite inefficient at the moment...

source("genedrop/genedropSingle.R")

#~~ Load in the pedigree. It must be ID, MOTHER, FATHER, missing values are NA, and *not* factors

pedigree <- read.table("genedrop/4_Updated_Pedigree_Feb2017.txt", header = T, stringsAsFactors = F)

str(pedigree)

#~~ If you want to know which ID's are founders (i.e. for calculating the founder allele frequency) then:

library(kinship2)

founders <- data.frame(ID = ped[,1],
                       Cohort = kindepth(ped[,1], ped[,2], ped[,3])) # Those with Cohort 0 are founders!

founders <- subset(founders, Cohort == 0)


#~~ The genedrop is run by running the genedropSingle function with the allele names and the
#   frequencies of the alleles in the founder individuals.

sim01 <- genedropSingle(pedigree,
                        allele.ids = c("A", "B", "C", "D"),
                        allele.founder.freqs = c(0.1, 0.2, 0.3, 0.4))

#~~ A quick look at the output - it gives a genotype for each individual, and also a cohort column which
#   indicates what "generation" the individual is in

head(sim01)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Looking at the true values         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Let's make this a fake dataset to compare other populations to, by sampling 6,000 of the IDs.

example.data <- sim01[sample(1:nrow(sim01), replace = F, size = 6000), c(3, 1, 2)]
head(example.data)

#~~ Merge it with birthyear data (you may want to use this with observation data)

basedata <- read.table("genedrop/BirthYear.txt", header = T, sep = "\t")
head(basedata)

library(plyr)
example.data <- join(example.data, basedata)
head(example.data)

#~~ melt the data to work out the frequencies in each generation.

example.data$Allele1 <- as.character(example.data$Allele1)
example.data$Allele2 <- as.character(example.data$Allele2)

example.melt <- melt(example.data, id.vars = c("ID", "BirthYear"))

#~~ Get the frequency per allele per year and frequency per year to calculate proportions

example.freq <- data.frame(table(example.melt$BirthYear, example.melt$value))
names(example.freq) <- c("BirthYear", "Allele", "Count")

example.year <- data.frame(table(example.melt$BirthYear))
names(example.year) <- c("BirthYear", "TotalCount")

example.freq <- join(example.freq, example.year)
example.freq$Allele.Freq <- example.freq$Count/example.freq$TotalCount

head(example.freq)

#~~ Look at the number of individuals per birth year:

library(ggplot2)

example.freq$BirthYear <- as.numeric(as.character(example.freq$BirthYear)) # Make these numeric
example.year$BirthYear <- as.numeric(as.character(example.year$BirthYear))

ggplot(example.year, aes(BirthYear, TotalCount)) +
  geom_line()

#~~ Lets cut off the data above 1990

example.data <- subset(example.data, BirthYear >= 1990)
head(example.data)

example.freq <- subset(example.freq, BirthYear >= 1990)

ggplot(example.freq, aes(BirthYear, Allele.Freq, colour = Allele)) +
  geom_line() +
  facet_wrap(~Allele)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Doing Gene Drop                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ For each gene drop simulation, we want to calculate the allele frequencies
#   for each locus per year as above.

# Let's say we run for for 100 repetitions (Assuming you have the example.data and use BirthYear as above!!
#  As uses this name below!). It should take about 15 minutes.
#  If you just want to have a look, run it for 10 iterations!
# In reality you have so many alleles I'd run it for at least 1000!

iterations <- 100

simulation.results <- NULL  # create an object to save the results

for(i in 1:iterations){
  
  #~~ Message
  
  print(paste("Running iteration", i, "of", iterations))
  
  #~~ Run genedrop
  
  sim.temp <- genedropSingle(pedigree,
                             allele.ids = c("A", "B", "C", "D"),
                             allele.founder.freqs = c(0.1, 0.2, 0.3, 0.4),
                             verbose = F) # verbose just means whether to print out "running cohort i of 11"
  
  #~~ Make alleles characters and not factors, then melt
  
  sim.temp$Allele1 <- as.character(sim.temp$Allele1)
  sim.temp$Allele2 <- as.character(sim.temp$Allele2)
  
  sim.melt <- melt(sim.temp, id.vars = c("ID", "Cohort"))
  
  #~~ Join it to the real data:
  
  head(example.data)
  
  sim.data <- unique(example.data[,c("ID", "BirthYear")])
  sim.data <- join(sim.data, sim.melt)
  
  #~~ Get the frequency per allele per year and frequency per year to calculate proportions
  
  sim.freq <- data.frame(table(sim.data$BirthYear, sim.data$value))
  names(sim.freq) <- c("BirthYear", "Allele", "Count")
  
  sim.year <- data.frame(table(sim.data$BirthYear))
  names(sim.year) <- c("BirthYear", "TotalCount")
  
  sim.freq <- join(sim.freq, sim.year)
  sim.freq$Allele.Freq <- sim.freq$Count/sim.freq$TotalCount
  
  head(sim.freq)
  
  sim.freq$Iteration <- i
  
  simulation.results <- rbind(simulation.results, sim.freq)
  
  rm(sim.data, sim.freq, sim.melt, sim.year)
  
}

head(simulation.results)


#~~ Plot the results of the simulations

simulation.results$BirthYear <- as.numeric(as.character(simulation.results$BirthYear))

ggplot(simulation.results, aes(BirthYear, Allele.Freq, group = Iteration)) +
  geom_line() +
  facet_wrap(~Allele)

#~~ We can plot the real data on top

example.freq$Iteration <- 0

ggplot(simulation.results, aes(BirthYear, Allele.Freq, group = Iteration)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~Allele) +
  geom_line(data = example.freq, aes(BirthYear, Allele.Freq), colour = "red", size = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Does the frequency change???
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Are allele frequencies changing at a greater rate than expected due to chance?

#~~ IN the true data, there is a slope:

ggplot(example.freq, aes(BirthYear, Allele.Freq)) +
  geom_line() +
  facet_wrap(~Allele) +
  stat_smooth(method = "lm")

true.slope <- NULL

for(i in unique(example.freq$Allele)){
  fit <- lm(Allele.Freq ~ BirthYear, data = subset(example.freq, Allele == i))
  x <- data.frame(Allele = i,
                  Slope = fit$coefficients[2])
  true.slope <- rbind(true.slope, x)
  rm(x)
}

#~~ In the simulated data there are slopes:

sim.slope <- NULL

for(h in 1:iterations){
  
  sim.data <- droplevels(subset(simulation.results, Iteration == h))
  
  for(i in unique(sim.data$Allele)){
    
    fit <- lm(Allele.Freq ~ BirthYear, data = subset(sim.data, Allele == i))
    
    x <- data.frame(Allele = i,
                    Slope = fit$coefficients[2],
                    Iteration = h)
    
    sim.slope <- rbind(sim.slope, x)
    
    rm(x, sim.data)
  }
  
}

#~~ How does each slope relate? Is it in either of the 2.5% tails?

head(sim.slope)

head(true.slope)


ggplot(sim.slope, aes(Slope)) +
  geom_histogram() +
  facet_wrap(~Allele) +
  geom_vline(data = true.slope, aes(xintercept = Slope), colour = "red")

for(i in unique(true.slope$Allele)){
  
  print(paste("Running Allele", i))
  
  x <- droplevels(subset(sim.slope, Allele == i))
  y <- true.slope[which(true.slope$Allele == i), "Slope"] #slope
  
  print("How many slopes were greater than the observed?")
  print(table(x$Slope > y))
  
  
}





