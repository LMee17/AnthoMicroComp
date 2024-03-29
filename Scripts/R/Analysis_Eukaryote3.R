#21st December 2022
#Edited: 4th Jan 2023
#Edited: 9th February 2023: changed sociality category names
#Analysis : Eukaryote Count Data
#See Prokaryote Script for Resource details
set.seed(0908)

dir.create("output/Eukaryote/")
dir.create("output/Eukaryote/Composition_Analysis/")

## Libraries ####
library(tidyverse)
library(vegan)
library(ape)


## Functions ####

#distmat = distance/dissimilarity matrix, 
#title = analysis level for saving (Eukaryote/eukaryote/virus)
#function to assess stress vs dimensionality on this dataset
scree.NMDS <- function(distmat, kingdom, filename){
  strs <- vector(length = 10)
  for(i in 1:10){
    mds <- metaMDS(distmat, distance = "bray", autotransform = F, k = i)
    strs[i] <- mds$stress
    out <- data.frame(NoDimensions = 1:10, Stress = strs)
  }
  out$Dec[out$Stress > 0.2] <- "Not Good"
  out$Dec[out$Stress < 0.2 & out$Stress > 0.1] <- "Passable"
  out$Dec[out$Stress < 0.1 & out$Stress > 0.05] <- "Good"
  out$Dec[out$Stress < 0.05] <- "Excellent"
  out$Dec <- factor(out$Dec, levels = c("Not Good",
                                        "Passable",
                                        "Good", 
                                        "Excellent"))
  ggplot(data = out, aes(y = Stress, x = NoDimensions)) +
    geom_point(aes(colour = Dec)) +
    geom_hline(yintercept = 0.2, colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0.1, colour = "blue", linetype = "dashed") +
    geom_hline(yintercept = 0.05, colour = "green", linetype = "dashed") +
    scale_x_continuous(limits = c(1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_colour_manual(values = c("Excellent" = "green",
                                   "Good" = "blue",
                                   "Passable" = "purple",
                                   "Not Good" = "red")) + 
    labs(y = "Stress", 
         x = "Number of Dimensions",
         colour = "",
         title = "Stress versus Dimensionality")
  ggsave(paste("output/", kingdom, "/Composition_Analysis/",
               filename, "_StressVDimensionScreePlot.pdf", sep = ""))
  return(out)
}

##Load Metadata #####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Oct23.tsv",
                  sep = "\t", header = T)
#microbial metadata
taxkey <- read.table("input/Metadata/Taxa_HitKey_Dec22.tsv",
                     sep ="\t", header = T)
#taxomy
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t", quote = "")

##Step One: Readying Counts and MetaData ####
#count table
cnt <- read.table("input/Counts/Eukaryote_RawCounts.tsv")

#remove any empty taxa rows that appear in less than 5% samples
lt5 <- apply(cnt, 1, function(x) sum(x > 0)) <= ncol(cnt)*0.05
cnt <- cnt[!rownames(cnt) %in% names(lt5)[lt5 == T],]
#remove any samples that have < 100 reads
cnt <- cnt[colSums(cnt)>100]

#write up for ease of later iterations
write.table(cnt, "input/Counts/Eukaryote_Filtered_Raw.tsv",
            col.names = T, row.names = T, quote = F)

#for PCA / negative binomial regression, I will need a cnt matrix of relative abundances
cnt.rel <- apply(cnt,2, FUN=function(x){ x / sum(x)})

#write up
write.table(cnt.rel, "input/Counts/Eukaryote_RelativeAbundance.tsv",
            col.names = T, row.names = T, quote = F)

#when I tried to reduce by tribes, I lost unnacceptable numbers of samples
#so to try and ensure there's some balance in phylogeny (I mean ... as much as I 
#can) i'll reduce by family instead
#reduce the counts by tribes = ie keep only tribes that have at least 4 members
fam.tab <- table(met$Family)
metdf <- met[!met$Family %in% names(fam.tab)[fam.tab < 4],]
#I have more solitary samples than I did with prokaryote analysis!

cnt <- cnt[,names(cnt) %in% metdf$Sample.ID]
cnt.rel <- as.data.frame(cnt.rel)
cnt.rel <- cnt.rel[,names(cnt.rel) %in% metdf$Sample.ID]

write.table(cnt, "input/Counts/Eukaryote_Filtered_FamilyReduced_Raw.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

#sample metadata
#keep only the samples that are in the table
metdf <- met[met$Sample.ID %in% colnames(cnt),]
(table(metdf$Sociality))
#there's only 12 solitary samples! But it is what it is.



##Step Two: Produce Distance Matrix ####
# First step is to calculate a distance matrix. 
# using avgdist from vegan to rarefy counts to the lowest number of the smallest sample
# set by sociality (to try and conserve sample number) with 10,000 iterations
table(metdf$Sociality)
#solitary animals = smallest sociality sample group
#what is the lowest # reads in solitary animals ? 
readCnts <- cnt %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  group_by(Sample) %>%
  pivot_longer(-Sample) %>%
  mutate(Tot = sum(value)) %>%
  select(Sample, Tot) %>%
  unique() 

low_read <- inner_join(metdf, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Solitary") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

#the lowest # reads is 323
str(low_read)
# Using the Bray-Curtis distance metric
#vegan requires microbial hit as columns and samples as rows
cnt.t <- t(cnt)
dist <- avgdist(cnt.t, dmethod="bray", sample=low_read, iterations = 10000)
# The following sampling units were removed because they were below sampling depth: CA_N_007, CA_N_030, SRR11426433, SRR12053407, SRR12053408, SRR12053426, SRR13442822, SRR15496107, SRR3948521, SRR3948533, SRR3948536, SRR3948544, SRR3948562, SRR3948566, SRR3948570, SRR3948575, SRR8754044, SRR8867394

#write up
d <- as.matrix(dist)
write.table(d, "output/Eukaryote/Composition_Analysis/BC_DistanceMatrix.tsv",
            sep = "\t", row.names = T, col.names = T, quote = F)


##Step Three: NonMetric Multidimensional Scaling ####
pro.scree <- scree.NMDS(dist, "Eukaryote", "Eukaryote")
pro.scree

if(any(pro.scree$Dec == "Excellent") == T){
  kval <- pro.scree %>%
    subset(Dec == "Excellent") %>%
    arrange(NoDimensions) %>%
    head(n = 1) %>%
    select(NoDimensions) %>%
    as.numeric()
} else {
  print("Please assess by eye")
  stopifnot(any(pro.scree$Dec == "Excellent"))
}

#below 0.05 is excellent so going for 7 dimensions though I could reduce it to 4 as well 
#as this is still good
nmds <- metaMDS(dist, k = kval, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(nmds)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#save
pdf("output/Eukaryote/Composition_Analysis/NMDS_StressPlot.pdf")
stressplot(nmds)
dev.off()

#extract dimension axis
#"scores" pulls out the NMDS coordinates
ext <- as.data.frame(scores(nmds))
#can also use as.data.frame(pro.nmds$points)
#make plottable
ext$Sample <- rownames(ext)
rownames(ext) <- 1:nrow(ext)
nmds.plot <- merge(ext, metdf, by.x = "Sample", by.y = "Sample.ID")

#plot
ggplot(nmds.plot, aes(x = NMDS1, y = NMDS2))+
  geom_point() 
#check its the same as the vegan version
plot(nmds, display = "sites")

#save
write.table(nmds.plot, "output/Eukaryote/Composition_Analysis/NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#getting centroid plots for these for further exploration
#to calculate centroids I need more than / equal to 4 samples per group
#the rarefaction step in producing the matrix removed samples and so I need to 
#check what I have left for each factor (post distance)
metdf.pd <- metdf[metdf$Sample.ID %in% labels(dist),]

##Step Four: Checking Dispersion ####
#if the dispersion of different groups are very variable I'm likely to have false positive
#adonis / anosim results, whatever I decide to pursue
#first need to make sure the # rows and samples in the metadata match those in the 
#distance matrix 
if (nrow(metdf) != length(labels(dist))){
  metdf <- metdf[metdf$Sample.ID %in% labels(dist),]
}

#run through all variables, run betadisper to compute dispersions from median per
#levels within that variable, run anova to assess whether differences are significant
#control for multiple testing
#set up explanatory variables to loop through
vars <- c("Sociality", "Family", "Tribe", "Continent")
#and better labels to use
varslabs <- c("Sociality", "Host Family", "Host Tribe", "Continent")
disp.df <- data.frame()
for (i in 1:length(vars)){
  print(varslabs[i])
  disp.df[i,1] <- paste(varslabs[i])
  disp <- betadisper(dist, metdf[,vars[i]])
  disp.df[i,2] <- paste(names(disp$group.distances), collapse = ",")
  disp.df[i,3] <- paste(as.numeric(disp$group.distances), collapse = ",")
  a.disp <- anova(disp)
  disp.df[i,4] <- a.disp$Df[1]
  disp.df[i,5] <- a.disp$`F value`[1]
  disp.df[i,6] <- a.disp$`Pr(>F)`[1]
}
names(disp.df) <- c("Variable", "Factor_Levels", "AvgDistanceToMedian", "ANOVA_DF",
                    "ANOVA_Fvalue", "ANOVA_pvalue")
disp.df$adj_pvalue <- p.adjust(disp.df$ANOVA_pvalue, method = "fdr")
disp.df$Sig <- ifelse(disp.df$adj_pvalue < 0.05, "Sig", "Non-Sig")

#what variable have significant differences in dispersions ? 
disp.df %>%
  filter(adj_pvalue < 0.05) %>%
  select(Variable)

#are there any that don't ? 
disp.df %>%
  filter(adj_pvalue > 0.05) %>%
  select(Variable)
#only library selection isn't affected by dispersion

#write up
write.table(disp.df, 
            "output/Eukaryote/Composition_Analysis/VariableDispersion_ANOVA.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

##Step Five: PERMANOVA / Adonis ####
#it won't work with NA values, so I need to remove
torem <- vector()
j <- 1
for(i in 1:ncol(metdf)){
  if(any(is.na(metdf[,i])) == T){
    torem[j] <- print(names(metdf[i]))
    j <- j + 1
  }
}
var.key <- data.frame(Var = vars, VarLab = varslabs)
vars.sim <- var.key[!var.key$Var %in% torem,]
multiAOV <- data.frame()
for (i in 1:length(vars.sim$Var)){
  print(vars.sim$Var[i])
  test <- adonis2(dist ~ metdf[,paste(vars.sim$Var[i])], 
                  data = metdf, permutations = 9999)
  multiAOV[i,1] <- vars.sim$VarLab[i]
  multiAOV[i,2] <- test$Df[1]
  multiAOV[i,3] <- test$R2[1] 
  multiAOV[i,4] <- test$F[1]
  multiAOV[i,5] <- test$`Pr(>F)`[1]
}
names(multiAOV) <- c("Variable", "DF", "R2", "F", "pvalue")
multiAOV$adjP <-p.adjust(multiAOV$pvalue, method = "fdr")

multiAOV
#Quite a few significants but ... nothing is not affected by dispersal other than 
#year collected
#considering just the most significant (alpha = 0.001) and the dispersals...
multiAOV <- inner_join(multiAOV, disp.df, by = "Variable") %>%
  select(Variable, DF, R2, 'F', pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

multiAOV %>%
  filter(adjP < 0.001) %>%
  arrange(-R2)

#Year collected is significant, and not affected by dispersion but I doubt
#its anything to actually do with the year, rather the project

####Pairwise ####
#sociality
socs <- as.character(unique(metdf$Sociality))
cbn <- combn(x=socs, m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Sociality == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Sociality, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.soc <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.soc %>%
  filter(padj < 0.05)
#nothing significantly different in sociality (depsite it being significant above)

write.table(pairwiseAOV.soc, 
            "output/Eukaryote/Composition_Analysis/PairwiseAOV_Sociality.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

#sociality
cbn <- combn(x=unique(metdf$Family), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Family == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Family, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.fam <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.fam %>%
  filter(padj < 0.05)
#no differences between families

write.table(pairwiseAOV.fam, 
            "output/Eukaryote/Composition_Analysis/PairwiseAOV_Family.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

#continent
cbn <- combn(x=unique(metdf$Continent), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Continent == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Continent, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.con <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.con %>%
  filter(padj < 0.05)
#only diff between oceania and europe

write.table(pairwiseAOV.con, 
            "output/Eukaryote/Composition_Analysis/PairwiseAOV_Continent.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

###Modellling #####
p <- 9999

all <- adonis2(dist ~ Sociality + Continent + Family,
               data = metdf,
               permutations = p)
all


##SessionLogs####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Analysis_Eukaryote_Oct23.txt")
