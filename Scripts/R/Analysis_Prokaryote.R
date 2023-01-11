#21st December 2022
#Edited: 4th Jan 2023
#Analysis : Prokaryotic Count Data
#Prokaryotes
#Resources:
#NMDS/PCoA/PCA:
#https://ourcodingclub.github.io/tutorials/ordination/
#https://ourcodingclub.github.io/tutorials/data-clustering/index.html#Distance
#https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/principal-coordinates-analysis
#https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/nmds
#Adonis / PERMAnova / ANOSIM:
#https://rdrr.io/rforge/vegan/man/adonis.html
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html
#https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
#Rarefaction: 
#https://www.youtube.com/watch?v=ht3AX5uZTTQ&list=PLmNrK_nkqBpJuhS93PYC-Xr5oqur7IIWf

set.seed(1549)

dir.create("output/Prokaryote/Composition_Analysis/")
dir.create("output/Prokaryote/Composition_Analysis/AllTribes/")
dir.create("output/Prokaryote/Composition_Analysis/TribeReduced/")

## Libraries ####
library(tidyverse)
library(vegan)
library(ape)
library(gg3D)

## Functions ####

#distmat = distance/dissimilarity matrix, 
#title = analysis level for saving (prokaryote/eukaryote/virus)
#function to assess stress vs dimensionality on this dataset
scree.NMDS <- function(distmat, kingdom, filename, version){
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
  if (version == 1){
    route <- "AllTribes"
  } else {
    route <- "TribeReduced"
  }
  ggsave(paste("output/", kingdom, "/", route, 
               "/", filename, "_StressVDimensionScreePlot.pdf", sep = ""))
  return(out)
}

##Load Metadata #####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_RNAOnly_Dec22.tsv",
                  sep = "\t", header = T)
#microbial metadata
taxkey <- read.table("input/Metadata/Taxa_HitKey_Dec22.tsv",
                     sep ="\t", header = T)
#taxomy
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t", quote = "")

##Step One: Readying Counts and MetaData ####
#count table
cnt <- read.table("input/Counts/Prokaryote_RawCounts.tsv")

#remove any empty taxa rows that appear in less than 5% samples
lt5 <- apply(cnt, 1, function(x) sum(x > 0)) <= ncol(cnt)*0.05
cnt <- cnt[!rownames(cnt) %in% names(lt5)[lt5 == T],]
#remove any samples that have < 100 reads
cnt <- cnt[colSums(cnt)>100]

#write up for ease of later iterations
write.table(cnt, "input/Counts/Prokaryote_Filtered_Raw.tsv",
            col.names = T, row.names = T, quote = F)

#for PCA / negative binomial regression, I will need a cnt matrix of relative abundances
cnt.rel <- apply(cnt,2, FUN=function(x){ x / sum(x)})

#write up
write.table(cnt.rel, "input/Counts/Prokaryote_RelativeAbundance.tsv",
            col.names = T, row.names = T, quote = F)

#sample metadata
#keep only the samples that are in the table
metdf <- met[met$Sample.ID %in% colnames(cnt),]

#add a second socialty field (complex versus primitive eusocial)
#this is just to perhaps balance the sample sets more fairly
metdf$Sociality2 <- metdf$Sociality
metdf$Sociality2[metdf$Tribe == "Apini"] <- "Complex Eusocial"
metdf$Sociality2[metdf$Tribe == "Meliponini"] <- "Complex Eusocial"
metdf$Sociality2[metdf$Tribe == "Bombini"] <- "Primitive Eusocial"

#there will be two version of each count tables going forward:
#1 with everything included and 1 that reduced samples so there were at least 4
#members in each tribe (at least four are needed to compute centroids later)
tri.tab <- table(metdf$Tribe)
metdf2 <- metdf[!metdf$Tribe %in% names(tri.tab)[tri.tab < 4],]

cnt2 <- cnt[,names(cnt) %in% metdf2$Sample.ID]
cnt.rel <- as.data.frame(cnt.rel)
cnt.rel2 <- cnt.rel[,names(cnt.rel) %in% metdf2$Sample.ID]

write.table(cnt2, "input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv",
            sep = "\t", col.names = T, row.names = T)

##Step Two: Principal Components ####
#most of these steps assume that rows are samples / sites, not columns
#using relative abundance data
rel.t <- t(cnt.rel)
rel.t2 <- t(cnt.rel2)
#not using scales as the variables are not on different scales
pca <- rda(rel.t, scale = FALSE)
pca2 <- rda(rel.t2, scale = FALSE)

###Version 1 #####
# Now plot a bar plot of relative eigenvalues. 
#This is the percentage variance explained by each axis
barplot(as.vector(pca$CA$eig)/sum(pca$CA$eig)) 
# over 25% of variance explained by first PC (would probably expect a good pca to be up
# to 60 but there is so much noise in this dataset)

# Calculate the percent of variance explained by first two to four axes
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:2]) # 41.1%
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:4]) # 61.6%

#extract the principal coordinates from the pca object
#prepare plottable dataframe
tmp <- as.data.frame(pca$CA$u)
tmp$Sample <- rownames(tmp)
#add sample metadata
propc.plot <- merge(tmp, metdf, by.x = "Sample", by.y = "Sample.ID")

#manually set palette
#colour-blind friendly palette sourced: 
#https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
myPal <- c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
           "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
           "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
           "#9acd32")

#set up explanatory variables to loop through
vars <- c("Sociality", "Sociality2", "Sex", "YearCollected", "Month", "LibraryLayout",
          "LibrarySelection", "Platform_Spec", "Family", "Tribe", "Continent")
#and better labels to use
varslabs <- c("Sociality", "Sociality ", "Sex", "Year Collected", "Month Collected",
              "Library Layout", "Library Selection", "Platform",
              "Host Family", "Host Tribe", "Continent")

#looking at PC1 and PC2 ... what % of the variance do they represent?
one <- round(sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1])*100, digits = 2)
two <- round(sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[2])*100, digits = 2)

#run through variables and plot
for (i in 1:length(vars)){
  print(vars[i])
  ggplot(data = propc.plot, 
         aes(x = PC1, y = PC2, colour = as.factor(propc.plot[,vars[i]]))) +
    geom_point(size = 2, alpha = 0.75) + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("PC1 (", one, "%)", sep = ""),
         y = paste("PC2 (", two, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    guides(alpha = "none") + 
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/PCA_PC1PC2_", title, ".pdf", sep = ""))
}

#looking at pc3 and pc4
three <- round(sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[3])*100, digits = 2)
four <- round(sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[4])*100, digits = 2)

for (i in 1:length(vars)){
  ggplot(data = propc.plot, 
         aes(x = PC3, y = PC4, colour = as.factor(propc.plot[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("PC3 (", three, "%)", sep = ""),
         y = paste("PC4 (", four, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/PCA_PC3PC4_", title, ".pdf", sep = ""))
}

#Looking at three axes at once
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:3]) # 52.6%
#may be worth doing some 3D plots
for (i in 1:length(vars)){
  ggplot(data = propc.plot, 
         aes(x = PC1, y = PC2, z = PC3,
             colour = as.factor(propc.plot[,vars[i]]))) +
    scale_colour_manual(values = myPal) +
    labs(colour = paste(varslabs[i])) +
    theme_void() +
    axes_3D() +
    stat_3D() +
    labs_3D(labs = c(paste("PC1 (", one, "%)", sep = ""),
                     paste("PC2 (", two, "%)", sep = ""), 
                     paste("PC3 (", three, "%)", sep = "")),
            hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90))
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/PCA_3D_", title, ".pdf", sep = ""))
}

###Version 2#####
barplot(as.vector(pca2$CA$eig)/sum(pca2$CA$eig)) 
# over 25%

sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[1:2]) # 42.7%
sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[1:4]) # 62.7%

#extract the principal coordinates from the pca object
#prepare plottable dataframe
tmp <- as.data.frame(pca2$CA$u)
tmp$Sample <- rownames(tmp)
#add sample metadata
propc.plot2 <- merge(tmp, metdf2, by.x = "Sample", by.y = "Sample.ID")

#looking at PC1 and PC2 ... what % of the variance do they represent?
one <- round(sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[1])*100, digits = 2)
two <- round(sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[2])*100, digits = 2)

#run through variables and plot
for (i in 1:length(vars)){
  print(vars[i])
  ggplot(data = propc.plot2, 
         aes(x = PC1, y = PC2, colour = as.factor(propc.plot2[,vars[i]]))) +
    geom_point(size = 2, alpha = 0.75) + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("PC1 (", one, "%)", sep = ""),
         y = paste("PC2 (", two, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    guides(alpha = "none") + 
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/PCA_PC1PC2_", title, ".pdf", sep = ""))
}

#looking at pc3 and pc4
three <- round(sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[3])*100, digits = 2)
four <- round(sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[4])*100, digits = 2)

for (i in 1:length(vars)){
  ggplot(data = propc.plot2, 
         aes(x = PC3, y = PC4, colour = as.factor(propc.plot2[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("PC3 (", three, "%)", sep = ""),
         y = paste("PC4 (", four, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/PCA_PC3PC4_", title, ".pdf", sep = ""))
}

#Looking at three axes at once
sum((as.vector(pca2$CA$eig)/sum(pca2$CA$eig))[1:3]) # 54.1%
#may be worth doing some 3D plots
for (i in 1:length(vars)){
  ggplot(data = propc.plot2, 
         aes(x = PC1, y = PC2, z = PC3,
             colour = as.factor(propc.plot2[,vars[i]]))) +
    scale_colour_manual(values = myPal) +
    labs(colour = paste(varslabs[i])) +
    theme_void() +
    axes_3D() +
    stat_3D() +
    labs_3D(labs = c(paste("PC1 (", one, "%)", sep = ""),
                     paste("PC2 (", two, "%)", sep = ""), 
                     paste("PC3 (", three, "%)", sep = "")),
            hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90))
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/PCA_3D_", title, ".pdf", sep = ""))
}


##Step Three: BiPlot ####
###1####
biplot(pca, choices = c(1,2), type = c("text", "points"))

pdf("output/Prokaryote/Composition_Analysis/AllTribes/BiPlot_PC1PC2.pdf")
  biplot(pca, choices = c(1,2), type = c("text", "points"))
dev.off()
#without text
pdf("output/Prokaryote/Composition_Analysis/AllTribes/BiPlot_PC1PC2_Plain.pdf")
  biplot(pca, choices = c(1,2), type = c("points"))
dev.off()

#there's three clear taxa that affect these ordinations (plus a possible fourth)
vip <- c("TRH2", "TRH7", "TRH46", "TRH29")
taxkey[taxkey$taxID %in% vip,]
#Gilliamella (yay), Escherichia, Ralstonia and Chryseobacterium
#Ralstonia and Chryseobacterium downwards, Gilliamella up, and 
#Escherichia to the right
###2####
biplot(pca2, choices = c(1,2), type = c("text", "points"))

pdf("output/Prokaryote/Composition_Analysis/TribeReduced/BiPlot_PC1PC2.pdf")
  biplot(pca2, choices = c(1,2), type = c("text", "points"))
dev.off()
#without text
pdf("output/Prokaryote/Composition_Analysis/TribeReduced/BiPlot_PC1PC2_Plain.pdf")
  biplot(pca2, choices = c(1,2), type = c("points"))
dev.off()

#the same big four are present in the second version, though the biplot is tidier.

##Step Four: Principal Coordinates Analysis ####
###1####
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

#the lowest # reads is 153
str(low_read)
# Using the Bray-Curtis distance metric
#vegan requires microbial hit as columns and samples as rows
cnt.t <- t(cnt)
dist <- avgdist(cnt.t, dmethod="bray", sample=low_read, iterations = 10000)
#The following sampling units were removed because they were below sampling depth: SRR11440494, SRR12053415, SRR2001621

# calculate the principal coordinates
pco <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(pco$values$Relative_eig[1:10])
#again, first eigenvalue explains just over 20% of the variance

# Some distance measures may result in negative eigenvalues. 
# negative eiganvalues need to be corrected
any(pco$values[,"Eigenvalues"] < 0)
# TRUE
# administer correction
pco <- pcoa(dist, correction = "cailliez")
#a "good" pcoa explains ~ 50% of the variance/distance within 3 axes
sum((pco$values$Corr_eig / sum(pco$values$Corr_eig))[1:3])*100
#I'm getting 33 % but ah well

#make the above object more amenable for playing with
tmp <- as.data.frame(pco$vectors)
tmp$Sample <- rownames(tmp)
#add sample metadata
pco.plot <- merge(tmp, metdf, by.x = "Sample", by.y = "Sample.ID")

#looking at axis1 and axis 2 ... 
one <- round(pco$values$Corr_eig[1] / sum(pco$values$Corr_eig)*100, digits = 2)
two <- round(pco$values$Corr_eig[2] / sum(pco$values$Corr_eig)*100, digits = 2)

#run through variables and plot
for (i in 1:length(vars)){
  ggplot(data = pco.plot, 
         aes(x = Axis.1, y = Axis.2, colour = as.factor(pco.plot[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("Axis1 (", one, "%)", sep = ""),
         y = paste("Axis2 (", two, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/PCoA_Ax1Ax2_", title, ".pdf", sep = ""))
}

###2####
table(metdf2$Sociality)
#solitary animals = smallest sociality sample group
#what is the lowest # reads in solitary animals ? 
readCnts2 <- cnt2 %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  group_by(Sample) %>%
  pivot_longer(-Sample) %>%
  mutate(Tot = sum(value)) %>%
  select(Sample, Tot) %>%
  unique() 

low_read2 <- inner_join(metdf, readCnts2, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Solitary") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

# Using the Bray-Curtis distance metric
#vegan requires microbial hit as columns and samples as rows
cnt.t2 <- t(cnt2)
dist2 <- avgdist(cnt.t2, dmethod="bray", sample=low_read2, iterations = 10000)
#note that many samples were dropped
#The following sampling units were removed because they were below sampling depth: CA_N_012, SRR11440494, SRR12053407, SRR12053410, SRR12053413, SRR12053415, SRR12053416, SRR12053423, SRR12053424, SRR12053425, SRR12053426, SRR12053427, SRR12053428, SRR12053429, SRR12053435, SRR12053437, SRR12527935, SRR12527963, SRR2001621, SRR2396653, SRR3383868, SRR3948521, SRR3948548

# calculate the principal coordinates
pco2 <- pcoa(dist2)

# plot the eigenvalues and interpret
barplot(pco2$values$Relative_eig[1:10])
#again, first eigenvalue explains just over 25% of the variance

# Some distance measures may result in negative eigenvalues. 
# negative eiganvalues need to be corrected
any(pco2$values[,"Eigenvalues"] < 0)
# TRUE
# administer correction
pco2 <- pcoa(dist2, correction = "cailliez")
#a "good" pcoa explains ~ 50% of the variance/distance within 3 axes
sum((pco2$values$Corr_eig / sum(pco2$values$Corr_eig))[1:3])*100
#I'm getting 36 % but ah well

#make the above object more amenable for playing with
tmp <- as.data.frame(pco2$vectors)
tmp$Sample <- rownames(tmp)
#add sample metadata
pco.plot2 <- merge(tmp, metdf2, by.x = "Sample", by.y = "Sample.ID")

#looking at axis1 and axis 2 ... 
one <- round(pco2$values$Corr_eig[1] / sum(pco2$values$Corr_eig)*100, digits = 2)
two <- round(pco2$values$Corr_eig[2] / sum(pco2$values$Corr_eig)*100, digits = 2)

#run through variables and plot
for (i in 1:length(vars)){
  ggplot(data = pco.plot2, 
         aes(x = Axis.1, y = Axis.2, colour = as.factor(pco.plot2[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = paste("Axis1 (", one, "%)", sep = ""),
         y = paste("Axis2 (", two, "%)", sep = ""),
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/PCoA_Ax1Ax2_", title, ".pdf", sep = ""))
}

##Step Five: NonMetric Multidimensional Scaling ####
###1####
pro.scree <- scree.NMDS(dist, "Prokaryote", "All", 1)
pro.scree

#below 0.05 is excellent so going for 7 dimensions though I could reduce it to 4 as well 
#as this is still good
nmds <- metaMDS(dist, k = 7, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(nmds)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#save
pdf("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_StressPlot.pdf")
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
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/AllTribes/NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#run through variables and plot
for (i in 1:length(vars)){
  ggplot(data = nmds.plot, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_", title, ".pdf", sep = ""))
}

#run through variables and plot: with ellipses
for (i in 1:length(vars)){
  ggplot(data = nmds.plot, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot[,vars[i]]))) +
    geom_point() + 
    stat_ellipse() +
    scale_colour_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_Ellipses_", title, ".pdf", sep = ""))
}

#run through variables and plot: with polygons
for (i in 1:length(vars)){
  ggplot(data = nmds.plot, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot[,vars[i]]))) +
    geom_point() + 
    stat_ellipse(geom = "polygon",
                 aes(fill = as.factor(nmds.plot[,vars[i]])),
                 alpha = 0.1) +
    scale_colour_manual(values = myPal) +
    scale_fill_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    guides(fill = "none")
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_Polygons_", title, ".pdf", sep = ""))
}

#from looking at the plots so far it looks like the stand out clusters
#that make sense are family and sociality
#getting centroid plots for these for further exploration
#to calculate centroids I need more than / equal to 4 samples per group
#the rarefaction step in producing the matrix removed samples and so I need to 
#check what I have left for each factor (post distance)
metdf.pd <- metdf[metdf$Sample.ID %in% labels(dist),]

#sociality#1
any(table(metdf.pd$Sociality) < 4)

cen.soc <- nmds.plot %>%
  group_by(Sociality) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot, 
       aes(x = NMDS1, y = NMDS2, colour = Sociality)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.soc, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Sociality),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_Sociality_EllipseCentroid.pdf")

#sociality#2
any(table(metdf.pd$Sociality2) < 4)

cen.soc2 <- nmds.plot %>%
  group_by(Sociality2) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot, 
       aes(x = NMDS1, y = NMDS2, colour = Sociality2)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.soc2, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Sociality2),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  labs(colour = "Sociality") +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_Sociality2_EllipseCentroid.pdf")

#family
any(table(metdf.pd$Family) < 4)
#can't be done... unfortunately

#3D Plot
#A quick exploration of what happens if I bring in another NMDS axis
#only for the more interesting factors discussed above
foi <- c("Sociality", "Sociality2", "Family")
foilabs <- c("Sociality", "Sociality ", "Family")

for (i in 1:length(foi)){
  ggplot(data = nmds.plot, 
         aes(x = NMDS1, y = NMDS2, z = NMDS3,
             colour = as.factor(nmds.plot[,foi[i]]))) +
    scale_colour_manual(values = myPal) +
    labs(colour = paste(foilabs[i])) +
    theme_void() +
    axes_3D() +
    stat_3D() +
    labs_3D(labs = c("NMDS1",
                     "NMDS2", 
                     "NMDS3"),
            hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90))
  title <- gsub(" ", "_", foilabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/AllTribes/NMDS_3D_", title, ".pdf", sep = ""))
}

###2####
pro.scree2 <- scree.NMDS(dist2, "Prokaryote", "All", 2)
pro.scree2

#below 0.05 is excellent so going for 6 dimensions though I could reduce it to 4 as well 
#as this is still good
nmds2 <- metaMDS(dist2, k = 6, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(nmds2)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#save
pdf("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_StressPlot.pdf")
  stressplot(nmds2)
dev.off()

#extract dimension axis
#"scores" pulls out the NMDS coordinates
ext <- as.data.frame(scores(nmds2))
#can also use as.data.frame(pro.nmds$points)
#make plottable
ext$Sample <- rownames(ext)
rownames(ext) <- 1:nrow(ext)
nmds.plot2 <- merge(ext, metdf2, by.x = "Sample", by.y = "Sample.ID")

#plot
ggplot(nmds.plot2, aes(x = NMDS1, y = NMDS2))+
  geom_point() 
#check its the same as the vegan version
plot(nmds2, display = "sites")

#save
write.table(nmds.plot2, "output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#run through variables and plot
for (i in 1:length(vars)){
  ggplot(data = nmds.plot2, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot2[,vars[i]]))) +
    geom_point() + 
    scale_colour_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_", title, ".pdf", sep = ""))
}

#run through variables and plot: with ellipses
for (i in 1:length(vars)){
  ggplot(data = nmds.plot2, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot2[,vars[i]]))) +
    geom_point() + 
    stat_ellipse() +
    scale_colour_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Ellipses_", title, ".pdf", sep = ""))
}

#run through variables and plot: with polygons
for (i in 1:length(vars)){
  ggplot(data = nmds.plot2, 
         aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot2[,vars[i]]))) +
    geom_point() + 
    stat_ellipse(geom = "polygon",
                 aes(fill = as.factor(nmds.plot2[,vars[i]])),
                 alpha = 0.1) +
    scale_colour_manual(values = myPal) +
    scale_fill_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(varslabs[i])) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    guides(fill = "none")
  title <- gsub(" ", "_", varslabs[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Polygons_", title, ".pdf", sep = ""))
}

#from looking at the plots so far it looks like the stand out clusters
#that make sense are tribe/family, sociality and possibly platform ...
#getting centroid plots for these for further exploration

#to calculate centroids I need more than / equal to 4 samples per group
#the rarefaction step in producing the matrix removed samples and so I need to 
#check what I have left for each factor (post distance)
metdf.pd2 <- metdf2[metdf2$Sample.ID %in% labels(dist2),]

#sociality#1
any(table(metdf.pd2$Sociality) < 4)

cen.soc2 <- nmds.plot2 %>%
  group_by(Sociality) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot2, 
       aes(x = NMDS1, y = NMDS2, colour = Sociality)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.soc2, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Sociality),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Sociality_EllipseCentroid.pdf")

#sociality#2
any(table(metdf.pd2$Sociality2) < 4)

cen.soc2.2 <- nmds.plot2 %>%
  group_by(Sociality2) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot2, 
       aes(x = NMDS1, y = NMDS2, colour = Sociality2)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.soc2.2, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Sociality2),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  labs(colour = "Sociality") +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Sociality2_EllipseCentroid.pdf")

#family
any(table(metdf.pd2$Family) < 4)

cen.fam2 <- nmds.plot2 %>%
  group_by(Family) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot2, 
       aes(x = NMDS1, y = NMDS2, colour = Family)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.fam2, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Family),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Family_EllipseCentroid.pdf")

#tribe
any(table(metdf.pd2$Tribe) < 4)

cen.tri2 <- nmds.plot2 %>%
  group_by(Tribe) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(data = nmds.plot2, 
       aes(x = NMDS1, y = NMDS2, colour = Tribe)) +
  geom_point() + 
  scale_colour_manual(values = myPal) +
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cen.tri2, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = Tribe),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_Tribe_EllipseCentroid.pdf")

#Exploratory look at 3D plots with 3 NDMS dimensions
#only considering the factors already that look interesting
foi2 <- c("Sociality", "Sociality2", "Family", "Tribe")
foilabs2 <- c("Sociality", "Sociality ", "Family", "Tribe")

for (i in 1:length(foi2)){
  ggplot(data = nmds.plot2, 
         aes(x = NMDS1, y = NMDS2, z = NMDS3,
             colour = as.factor(nmds.plot2[,foi2[i]]))) +
    scale_colour_manual(values = myPal) +
    labs(colour = paste(foilabs2[i])) +
    theme_void() +
    axes_3D() +
    stat_3D() +
    labs_3D(labs = c("NMDS1",
                     "NMDS2", 
                     "NMDS3"),
            hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90))
  title <- gsub(" ", "_", foilabs2[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/TribeReduced/NMDS_3D_", title, ".pdf", sep = ""))
}

##Step Six: Checking Dispersion ####
#if the dispersion of different groups are very variable I'm likely to have false positive
#adonis / anosim results, whatever I decide to pursue
###1####
#first need to make sure the # rows and samples in the metadata match those in the 
#distance matrix 
if (nrow(metdf) != length(labels(dist))){
  metdf <- metdf[metdf$Sample.ID %in% labels(dist),]
}

#run through all variables, run betadisper to compute dispersions from median per
#levels within that variable, run anova to assess whether differences are significant
#control for multiple testing
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

#write up
write.table(disp.df, "output/Prokaryote/Composition_Analysis/AllTribes/VariableDispersion_ANOVA.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

#Sociality (1, with three levels) does not have significant differences in group dispersion
#and nor does Continent.

###2####
#first need to make sure the # rows and samples in the metadata match those in the 
#distance matrix 
if (nrow(metdf2) != length(labels(dist2))){
  metdf2 <- metdf2[metdf2$Sample.ID %in% labels(dist2),]
}

#run through all variables, run betadisper to compute dispersions from median per
#levels within that variable, run anova to assess whether differences are significant
#control for multiple testing
disp.df2 <- data.frame()
for (i in 1:length(vars)){
  print(varslabs[i])
  disp.df2[i,1] <- paste(varslabs[i])
  disp2 <- betadisper(dist2, metdf2[,vars[i]])
  disp.df2[i,2] <- paste(names(disp2$group.distances), collapse = ",")
  disp.df2[i,3] <- paste(as.numeric(disp2$group.distances), collapse = ",")
  a.disp2 <- anova(disp2)
  disp.df2[i,4] <- a.disp2$Df[1]
  disp.df2[i,5] <- a.disp2$`F value`[1]
  disp.df2[i,6] <- a.disp2$`Pr(>F)`[1]
}
names(disp.df2) <- c("Variable", "Factor_Levels", "AvgDistanceToMedian", "ANOVA_DF",
                    "ANOVA_Fvalue", "ANOVA_pvalue")
disp.df2$adj_pvalue <- p.adjust(disp.df2$ANOVA_pvalue, method = "fdr")
disp.df2$Sig <- ifelse(disp.df2$adj_pvalue < 0.05, "Sig", "Non-Sig")

#what variable have significant differences in dispersions ? 
disp.df2 %>%
  filter(adj_pvalue < 0.05) %>%
  select(Variable)

#are there any that don't ? 
disp.df2 %>%
  filter(adj_pvalue > 0.05) %>%
  select(Variable)

#Sociality (1, with three levels) does not have significant diffences in group dispersion
#and nor does Continent, Month Collected, Library Selection or Host Family .

#write up
write.table(disp.df2, "output/Prokaryote/Composition_Analysis/TribeReduced/VariableDispersion_ANOVA.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

##Step Seven: ANOSIM ####
#Analysis of similarities
###1####
#I'm just going to look over all factors for now individually and have a look-see
#ANOSIM won't run if any of the factor levels are NA, so these variables either need
#to be removed or else the NAs removed
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

#make a dataframe ready to populate with test statistic results
multiSIM <- data.frame()
for (i in 1:length(vars.sim$Var)){
  print(vars.sim$Var[i])
  test <- anosim(dist, metdf[,paste(vars.sim$Var[i])], permutations = 9999)
  multiSIM[i,1] <- vars.sim$VarLab[i]
  multiSIM[i,2] <- test$statistic
  multiSIM[i,3] <- test$signif
}
names(multiSIM) <- c("Variable", "R_Statistic", "pvalue")
multiSIM$adjP <-p.adjust(multiSIM$pvalue, method = "fdr")

#consider dispersal
multiSIM <- inner_join(multiSIM, disp.df, by = "Variable") %>%
  select(Variable, R_Statistic, pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

#assess
#what factors are significant and are unlikely to be false positives (no sig diff in
#dispersal)?
multiSIM %>%
  subset(Significance == "*" & Dispersion == "NonSig Dispersed") 
#sociality (simple levels) and continent are both considered significant as explanatory
#factors and neither have significantly different dispersals amongst factor levels
#Sociality has a higher R statistic, meaning there is a higher level of separation
#between the levels of the factor than there is when considering continent as a variable

#write up
write.table(multiSIM, "output/Prokaryote/Composition_Analysis/AllTribes/Multi_ANOSIM_testResults.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

#ANOSIM is noted to be quite sensitive to differences in dispersions so I'm reluctant
#to discuss the other significant variables just yet

#pairwise: sociality
cbn <- combn(x = unique(metdf$Sociality), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Sociality == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  sim <- anosim(subdist, submet$Sociality)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.soc <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.soc

#there's only a sig difference between Eusocial and Polymorphic groups
#no sig difference between sol and eu
#nearly a diff between pol and sol, but non-significant even before adjusting

write.table(pairSIM.soc, "output/Prokaryote/Composition_Analysis/AllTribes/PairwiseANOSIM_Sociality.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

#pairwise: Continent
cbn <- combn(x = unique(metdf$Continent), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Continent == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  sim <- anosim(subdist, submet$Continent)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.con <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.con

#there's only one significant difference after adjustment: Europe and Oceania
#which makes sense

write.table(pairSIM.con, "output/Prokaryote/Composition_Analysis/AllTribes/PairwiseANOSIM_Continent.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

###2####
#I'm just going to look over all factors for now individually and have a look-see
#ANOSIM won't run if any of the factor levels are NA, so these variables either need
#to be removed or else the NAs removed
torem <- vector()
j <- 1
for(i in 1:ncol(metdf2)){
  if(any(is.na(metdf2[,i])) == T){
    torem[j] <- print(names(metdf2[i]))
    j <- j + 1
  }
}
var.key <- data.frame(Var = vars, VarLab = varslabs)
vars.sim <- var.key[!var.key$Var %in% torem,]

#make a dataframe ready to populate with test statistic results
multiSIM2 <- data.frame()
for (i in 1:length(vars.sim$Var)){
  print(vars.sim$Var[i])
  test <- anosim(dist2, metdf2[,paste(vars.sim$Var[i])], permutations = 9999)
  multiSIM2[i,1] <- vars.sim$VarLab[i]
  multiSIM2[i,2] <- test$statistic
  multiSIM2[i,3] <- test$signif
}
names(multiSIM2) <- c("Variable", "R_Statistic", "pvalue")
multiSIM2$adjP <-p.adjust(multiSIM2$pvalue, method = "fdr")

#consider dispersal
multiSIM2 <- inner_join(multiSIM2, disp.df2, by = "Variable") %>%
  select(Variable, R_Statistic, pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

#assess
#what factors are significant and are unlikely to be false positives (no sig diff in
#dispersal)?
multiSIM2 %>%
  subset(Significance == "*" & Dispersion == "NonSig Dispersed") 
#host family and continent are significant and aren't potentially biased by 
#very unequal dispersions
#I've lost sociality which was strong in the previous run ... likely because
#I lost some Polymporphic and Eusocial samples whilst optimising for tribe
#membership, and that was where the big difference was .....
#because it's what I'm considering, I'm still going to do pairwise ANOSIM
#with sociality 

#write up
write.table(multiSIM2, "output/Prokaryote/Composition_Analysis/TribeReduced/Multi_ANOSIM_testResults.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

#ANOSIM is noted to be quite sensitive to differences in dispersions so I'm reluctant
#to discuss the other significant variables just yet

#pairwise: sociality
cbn <- combn(x = unique(metdf2$Sociality), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Sociality == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  sim <- anosim(subdist, submet$Sociality)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.soc2 <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.soc2

#there is  still a significant difference between eusocial and polymorphic,
#just less so.

write.table(pairSIM.soc2, "output/Prokaryote/Composition_Analysis/TribeReduced/PairwiseANOSIM_Sociality.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

#pairwise: Continent
cbn <- combn(x = unique(metdf2$Continent), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Continent == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  sim <- anosim(subdist, submet$Continent)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.con2 <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.con2

#again, there's only one significant difference after adjustment: Europe and Oceania
#which makes sense

write.table(pairSIM.con2, "output/Prokaryote/Composition_Analysis/TribeReduced/PairwiseANOSIM_Continent.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

#pairwise: Family
cbn <- combn(x = unique(metdf2$Family), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Family == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  sim <- anosim(subdist, submet$Family)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.fam2 <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.fam2

#there's a significant difference between Halictidae and everything else, though
#all the others appear to overlap

write.table(pairSIM.fam2, "output/Prokaryote/Composition_Analysis/TribeReduced/PairwiseANOSIM_Family.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

#and finally, as I conserved the tribes for this, I'm going to look at them too
#pairwise: Tribe
cbn <- combn(x = unique(metdf2$Tribe), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Tribe == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  sim <- anosim(subdist, submet$Tribe)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.tri2 <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.tri2 %>%
  subset(padj < 0.05)

#there's a significant difference Augochlorini and Bombini/Osmiini/Andrenini,
#and also between Apini/Euglossini

write.table(pairSIM.tri2, "output/Prokaryote/Composition_Analysis/TribeReduced/PairwiseANOSIM_Tribe.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")


##Step Eight: PERMANOVA / Adonis ####
#with permanovas/ adonis2 (vegan function that is just a permutated anova) I can begin
#to control what I consider fixed / random variables, nestedness, strata control etc etc.
#(I can do some of these with anosim but it doesn't look as simple / well documented)
###1####
#To start with, I'm going to run through every factor considered in isolation and see
#what is giving a significant signal, before looking more closely at sociality/location/
#family(or)tribe
#like anosim, it won't work with NA values, so I'll use the vars.sim object from before
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
#according to this, everything is a significant variable except for Sequencing platform....
#but R2 varies lots
#considering just the most significant (alpha = 0.001) and the dispersals...
multiAOV <- inner_join(multiAOV, disp.df, by = "Variable") %>%
  select(Variable, DF, R2, 'F', pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

multiAOV %>%
  filter(adjP < 0.001) %>%
  arrange(-R2)

#tribe has the least residual error by quite a bit, followed by sequencing platform. 
#Host family follows, However these are nested variables (family:Tribe). 
#Continent comes next with 4 level sociality. Sex and Year collected are also 
#still significant. Continent and Sociality are the only factors that 
#would not be sensitive to the dispersion problem.

#write up
write.table(multiAOV, "output/Prokaryote/Composition_Analysis/AllTribes/Multi_adonis2_results.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

#my main concern continues to be sociality, but will see if I can factor these in
#note: to nest, the factor that is nested within another factor goes last, ie if 
#species is nested in site then it is written Site/Species. For permanova, 
#you then need to designate strata as Site (https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata)
#set permutations to 9999
p <- 9999

#sociality and tribe
soc.tri <- adonis2(dist ~ Sociality + Tribe, 
                   data = metdf, 
                   permutations = p)
soc.tri

#interaction likely won't work as every tribe = 1 sociality.
#nesting also seems to not work, potentially for the same reason

#family?
soc.fam <- adonis2(dist ~ Sociality + Family,
                   data = metdf,
                   permutations = p)
soc.fam

soc.fam.nest <- adonis2(dist ~ Family / Sociality,
                        data = metdf,
                        permutations = p,
                        strata = metdf$Family)
soc.fam.nest

#throwing all phylo in
soc.phy <- adonis2(dist ~ Sociality + Family/Tribe,
                   data = metdf,
                   permutations = p,
                   strata = metdf$Family)
soc.phy
#when controlling for phylogeny by stratifying by family, Tribe remains a significantly 
#differing factor but sociality does not. Though the levels of variance explained still 
#remain very low
soc.phy.int <- adonis2(dist ~ Sociality * Family/Tribe,
                  data = metdf,
                  permutations = p,
                  strata = metdf$Family)
soc.phy.int
#CONSIDER, each tribe has only 1 sociality. Each family is majoritively one sociality
#are these therefore all nested ? 
#regardless, it appears tribe is the most explanatory variable in the model HOWEVER
#it has to be remembered that there is a significant difference in dispersion in the 
#Tribe levels in this version and so this could be a false positive (I don't think it is)

#continent is the only other variable (other than sociality) that does not have 
#significantly difference dispersions
#what about continent?
soc.con <- adonis2(dist ~ Sociality + Continent,
                   data = metdf,
                   permutations = p)
soc.con

soc.con.nest <- adonis2(dist ~ Continent/Sociality,
                        data = metdf,
                        permutations = p, 
                        strata = metdf$Continent)
soc.con.nest

soc.con.int <- adonis2(dist ~ Sociality*Continent,
                       data = metdf,
                       permutations = p)
soc.con.int

#will also consider year of collection
soc.yea <- adonis2(dist ~ Sociality + YearCollected,
                   data = metdf,
                   permutations = p)
soc.yea

adonis2(dist ~ Sociality*YearCollected,
        data = metdf,
        permutations = p)

#considering all 
soc.all <- adonis2(dist ~ Sociality + Family/Tribe + YearCollected + Continent,
                   data = metdf,
                   permutations = p,
                   strata = metdf$Family)
soc.all
#so everything is still significant ... (though perhaps not sociality after correction)
#and this model has the least residual errors. It makes senses really, microbial 
#community should be first and foremost affected by location (especially in non-eusocial
#species), and location/phylogeny are interchangeable in non-global species. Seasonality
#would also affect the types of microbe about (year collected) plus as too would whether
#these were lab reared, caught and put in a lab or else field caught (which may be
#a factor of the project, not actually the year). 
#What would happen if I assessed each sociality individually? Would I see stronger/weaker
#input from locations/phylogeny ? 

####Eusocial Considerations#####
eu.cnt <- cnt[,names(cnt) %in% metdf$Sample.ID[metdf$Sociality == "Eusocial"]]
eu.low_read <- inner_join(metdf, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Eusocial") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

eu.dist <- avgdist(t(eu.cnt), dmethod="bray", sample=eu.low_read, iterations = 10000)
eu.stress <- scree.NMDS(eu.dist, "Prokaryote", "Eusocial", 1)
eu.stress %>%
  filter(Dec == "Excellent")

eu.nmds <- metaMDS(dist, k = 6, trymax = 100, trace = F, distances = "bray")

stressplot(eu.nmds)

pdf("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_StressPlot.pdf")
  stressplot(eu.nmds)
dev.off()

ext <- scores(eu.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
meteu <- metdf %>%
  filter(Sociality == "Eusocial")
eu.plot <- merge(ext, meteu, by.x = "Sample", by.y = "Sample.ID")
eu.plot

#save
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Adonis
#remove factors that now have one level (ie family and sociality)
rownames(vars.sim) <- 1:nrow(vars.sim)
eu.vars <- vars.sim[c(-1,-2,-7),]
eu.multiAOV <- data.frame()
for (i in 1:length(eu.vars$Var)){
    print(eu.vars$Var[i])
    test <- adonis2(eu.dist ~ meteu[,paste(eu.vars$Var[i])], 
                    data = meteu, permutations = 9999)
    eu.multiAOV[i,1] <- eu.vars$VarLab[i]
    eu.multiAOV[i,2] <- test$Df[1]
    eu.multiAOV[i,3] <- test$R2[1] 
    eu.multiAOV[i,4] <- test$F[1]
    eu.multiAOV[i,5] <- test$`Pr(>F)`[1]

}
names(eu.multiAOV) <- c("Variable", "DF", "R2", "F", "pvalue")
eu.multiAOV$adjP <-p.adjust(eu.multiAOV$pvalue, method = "fdr")

eu.multiAOV %>%
  filter(adjP < 0.05) %>%
  arrange(-R2)

#for eusocial animals, the most explanatory variables appear to be Continent,
#platform and sex. Tribe is is significant but doesn't tell much of the story.

write.table(eu.multiAOV, "output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_multiANOVA.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#double check dispersion
eu.disp.df <- data.frame()
for (i in 1:length(eu.vars$Var)){
  print(eu.vars$Var[i])
  eu.disp.df[i,1] <- paste(eu.vars$VarLab[i])
  disp <- betadisper(eu.dist, meteu[,eu.vars$Var[i]])
  eu.disp.df[i,2] <- paste(names(disp$group.distances), collapse = ",")
  eu.disp.df[i,3] <- paste(as.numeric(disp$group.distances), collapse = ",")
  a.disp <- anova(disp)
  eu.disp.df[i,4] <- a.disp$Df[1]
  eu.disp.df[i,5] <- a.disp$`F value`[1]
  eu.disp.df[i,6] <- a.disp$`Pr(>F)`[1]
}
names(eu.disp.df) <- c("Variable", "Factor_Levels", "AvgDistanceToMedian", "ANOVA_DF",
                    "ANOVA_Fvalue", "ANOVA_pvalue")
eu.disp.df$adj_pvalue <- p.adjust(eu.disp.df$ANOVA_pvalue, method = "fdr")
eu.disp.df$Sig <- ifelse(eu.disp.df$adj_pvalue < 0.05, "Sig", "Non-Sig")

eu.disp.df
#consider these
ggplot(data = eu.plot,
       aes(x = NMDS1, y = NMDS2, colour = Continent)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Continent),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Continent") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_Continent.pdf")

ggplot(data = eu.plot,
       aes(x = NMDS1, y = NMDS2, colour = Platform_Spec)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Platform_Spec),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Platform") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_Platform.pdf")

ggplot(data = eu.plot,
       aes(x = NMDS1, y = NMDS2, colour = Sex)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Sex),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Sex") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_Sex.pdf")

ggplot(data = eu.plot,
       aes(x = NMDS1, y = NMDS2, colour = Tribe)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Tribe),
               alpha = 0.1) +
  scale_fill_manual(values = myPal) +
  scale_colour_manual(values = myPal) + 
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Tribe") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_Tribe.pdf")


ggplot(data = eu.plot,
       aes(x = NMDS1, y = NMDS2, colour = as.factor(YearCollected))) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(YearCollected)),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Year Collected") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Eusocial_NMDS_YearCollected.pdf")


####Polymorphic Considerations#####
po.cnt <- cnt[,names(cnt) %in% metdf$Sample.ID[metdf$Sociality == "Polymorphic"]]
po.low_read <- inner_join(metdf, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Polymorphic") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

po.dist <- avgdist(t(po.cnt), dmethod="bray", sample=po.low_read, iterations = 10000)
po.stress <- scree.NMDS(po.dist, "Prokaryote", "Polymorphic", 1)
po.stress %>%
  filter(Dec == "Excellent")

po.nmds <- metaMDS(dist, k = 3, trymax = 100, trace = F, distances = "bray")

stressplot(po.nmds)

pdf("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_StressPlot.pdf")
stressplot(po.nmds)
dev.off()

ext <- scores(po.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
metpo <- metdf %>%
  filter(Sociality == "Polymorphic")
po.plot <- merge(ext, metpo, by.x = "Sample", by.y = "Sample.ID")
po.plot

#save
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Adonis
po.vars <- vars.sim[c(-1,-2),]
po.multiAOV <- data.frame()
for (i in 1:length(po.vars$Var)){
  print(po.vars$Var[i])
  test <- adonis2(po.dist ~ metpo[,paste(po.vars$Var[i])], 
                  data = metpo, permutations = 9999)
  po.multiAOV[i,1] <- po.vars$VarLab[i]
  po.multiAOV[i,2] <- test$Df[1]
  po.multiAOV[i,3] <- test$R2[1] 
  po.multiAOV[i,4] <- test$F[1]
  po.multiAOV[i,5] <- test$`Pr(>F)`[1]
  
}
names(po.multiAOV) <- c("Variable", "DF", "R2", "F", "pvalue")
po.multiAOV$adjP <-p.adjust(po.multiAOV$pvalue, method = "fdr")

po.multiAOV %>%
  filter(adjP < 0.05) %>%
  arrange(-R2)

#for polymorphic animals, the most explanatory variables appear to be Tribe,
#Family, Year Collected and platform. 
#Phylogenetic signal explains over 50% of model

write.table(po.multiAOV, "output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_multiANOVA.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#double check dispersion
po.disp.df <- data.frame()
for (i in 1:length(po.vars$Var)){
  print(po.vars$Var[i])
  po.disp.df[i,1] <- paste(po.vars$VarLab[i])
  disp <- betadisper(po.dist, metpo[,po.vars$Var[i]])
  po.disp.df[i,2] <- paste(names(disp$group.distances), collapse = ",")
  po.disp.df[i,3] <- paste(as.numeric(disp$group.distances), collapse = ",")
  a.disp <- anova(disp)
  po.disp.df[i,4] <- a.disp$Df[1]
  po.disp.df[i,5] <- a.disp$`F value`[1]
  po.disp.df[i,6] <- a.disp$`Pr(>F)`[1]
}
names(po.disp.df) <- c("Variable", "Factor_Levels", "AvgDistanceToMedian", "ANOVA_DF",
                       "ANOVA_Fvalue", "ANOVA_pvalue")
po.disp.df$adj_pvalue <- p.adjust(po.disp.df$ANOVA_pvalue, method = "fdr")
po.disp.df$Sig <- ifelse(po.disp.df$adj_pvalue < 0.05, "Sig", "Non-Sig")

po.disp.df %>%
  filter(Sig == "Non-Sig") %>%
  select(Variable)

#consider these
ggplot(data = po.plot,
       aes(x = NMDS1, y = NMDS2, colour = Tribe)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Tribe),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Tribe") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_Tribe.pdf")

ggplot(data = po.plot,
       aes(x = NMDS1, y = NMDS2, colour = Family)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Family),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Family") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_Family.pdf")

ggplot(data = po.plot,
       aes(x = NMDS1, y = NMDS2, colour = as.factor(YearCollected))) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(YearCollected)),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Year Collected") +
  guides(fill = "none") + 
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_YearCollected.pdf")
#likely to be an artefact of project here.

ggplot(data = po.plot,
       aes(x = NMDS1, y = NMDS2, colour = Platform_Spec)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Platform_Spec),
               alpha = 0.1) +
  scale_fill_manual(values = myPal) +
  scale_colour_manual(values = myPal) + 
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Sequencing Platform") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_Platform_Spec.pdf")

#and for completeness, continent
#Continent
ggplot(data = po.plot,
       aes(x = NMDS1, y = NMDS2, colour = Continent)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Continent),
               alpha = 0.1) +
  scale_fill_manual(values = myPal) +
  scale_colour_manual(values = myPal) + 
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Continent") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_Continent.pdf")

#to compare to others: sex
ggplot(data = po.plot,
      aes(x = NMDS1, y = NMDS2, colour = Sex)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Sex),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Sex") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Polymorphic_NMDS_Sex.pdf")

####Solitary Considerations#####
so.cnt <- cnt[,names(cnt) %in% metdf$Sample.ID[metdf$Sociality == "Solitary"]]
so.low_read <- inner_join(metdf, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Solitary") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

so.dist <- avgdist(t(so.cnt), dmethod="bray", sample=so.low_read, iterations = 10000)
so.stress <- scree.NMDS(so.dist, "Prokaryote", "Solitary", 1)
so.stress %>%
  filter(Dec == "Excellent")
#Warning about the number of samples ... it may not be enought to be informative (stress
# is very close to 0)

so.nmds <- metaMDS(dist, k = 2, trymax = 100, trace = F, distances = "bray")

stressplot(so.nmds)
#still ok r plot ... though it's not the best

pdf("output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_StressPlot.pdf")
  stressplot(so.nmds)
dev.off()

ext <- scores(so.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
metso <- metdf %>%
  filter(Sociality == "Solitary")
so.plot <- merge(ext, metso, by.x = "Sample", by.y = "Sample.ID")

#save
write.table(so.plot, "output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Adonis
torem <- c()
for (i in 1:nrow(vars.sim)){
  noFact <- length(names(table(metso[,vars.sim$Var[i]])))
  if(noFact == 1){
    print(vars.sim$Var[i])
    torem <- c(torem, vars.sim$Var[i])
  }
}
so.vars <- vars.sim[!vars.sim$Var %in% torem,]
so.multiAOV <- data.frame()
for (i in 1:length(so.vars$Var)){
  print(so.vars$Var[i])
  test <- adonis2(so.dist ~ metso[,paste(so.vars$Var[i])], 
                  data = metso, permutations = 9999)
  so.multiAOV[i,1] <- so.vars$VarLab[i]
  so.multiAOV[i,2] <- test$Df[1]
  so.multiAOV[i,3] <- test$R2[1] 
  so.multiAOV[i,4] <- test$F[1]
  so.multiAOV[i,5] <- test$`Pr(>F)`[1]
  
}
names(so.multiAOV) <- c("Variable", "DF", "R2", "F", "pvalue")
so.multiAOV$adjP <-p.adjust(so.multiAOV$pvalue, method = "fdr")

so.multiAOV %>%
  filter(adjP < 0.05) %>%
  arrange(-R2)

#for solitary, the only statistically significant factor is sex.
#will look at this and also those of interest elsewhere:
#tribe, family, continent


write.table(so.multiAOV, "output/Prokaryote/Composition_Analysis/AllTribes/Solitary_multiANOVA.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#double check dispersion
so.disp.df <- data.frame()
for (i in 1:length(so.vars$Var)){
  print(so.vars$Var[i])
  so.disp.df[i,1] <- paste(so.vars$VarLab[i])
  disp <- betadisper(so.dist, metso[,so.vars$Var[i]])
  so.disp.df[i,2] <- paste(names(disp$group.distances), collapse = ",")
  so.disp.df[i,3] <- paste(as.numeric(disp$group.distances), collapse = ",")
  a.disp <- anova(disp)
  so.disp.df[i,4] <- a.disp$Df[1]
  so.disp.df[i,5] <- a.disp$`F value`[1]
  so.disp.df[i,6] <- a.disp$`Pr(>F)`[1]
}
names(so.disp.df) <- c("Variable", "Factor_Levels", "AvgDistanceToMedian", "ANOVA_DF",
                       "ANOVA_Fvalue", "ANOVA_pvalue")
so.disp.df$adj_pvalue <- p.adjust(so.disp.df$ANOVA_pvalue, method = "fdr")
so.disp.df$Sig <- ifelse(so.disp.df$adj_pvalue < 0.05, "Sig", "Non-Sig")

so.disp.df %>%
  filter(Sig == "Non-Sig") %>%
  select(Variable)
ncol(so.cnt)
so.low_read
#sex
ggplot(data = so.plot,
       aes(x = NMDS1, y = NMDS2, colour = Sex)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Sex),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Sex") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_Sex.pdf")

#Tribe
ggplot(data = so.plot,
       aes(x = NMDS1, y = NMDS2, colour = Tribe)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Tribe),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Tribe") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_Tribe.pdf")

#Family
ggplot(data = so.plot,
       aes(x = NMDS1, y = NMDS2, colour = Family)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Family),
               alpha = 0.1) +
  scale_colour_manual(values = myPal) +
  scale_fill_manual(values = myPal) +
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Family") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_Family.pdf")

#Continent
ggplot(data = so.plot,
       aes(x = NMDS1, y = NMDS2, colour = Continent)) +
  geom_point() + 
  stat_ellipse(geom = "polygon",
               aes(fill = Continent),
               alpha = 0.1) +
  scale_fill_manual(values = myPal) +
  scale_colour_manual(values = myPal) + 
  labs(x = "NMDS1",
       y = "NMDS2",
       colour = "Continent") +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave("output/Prokaryote/Composition_Analysis/AllTribes/Solitary_NMDS_Continent_Spec.pdf")


####Pairwise ####
#sociality
cbn <- combn(x=unique(metdf$Sociality), m = 2)
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
#like with the anosim, the only significant difference is between eusocial and 
#polymorphic

write.table(pairwiseAOV.soc, "output/Prokaryote/Composition_Analysis/AllTribes/PairwiseAOV_Sociality.tsv",
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
#like with the anosim, the only significant difference is between eusocial and 
#polymorphic

write.table(pairwiseAOV.fam, "output/Prokaryote/Composition_Analysis/AllTribes/PairwiseAOV_Family.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

###2###
multiAOV2 <- data.frame()
for (i in 1:length(vars.sim$Var)){
  print(vars.sim$Var[i])
  test <- adonis2(dist2 ~ metdf2[,paste(vars.sim$Var[i])], 
                  data = metdf2, permutations = 9999)
  multiAOV2[i,1] <- vars.sim$VarLab[i]
  multiAOV2[i,2] <- test$Df[1]
  multiAOV2[i,3] <- test$R2[1] 
  multiAOV2[i,4] <- test$F[1]
  multiAOV2[i,5] <- test$`Pr(>F)`[1]
}
names(multiAOV2) <- c("Variable", "DF", "R2", "F", "pvalue")
multiAOV2$adjP <-p.adjust(multiAOV2$pvalue, method = "fdr")

multiAOV2
#according to this, everything is a significant variable except for Sequencing platform....
#but R2 varies lots
#considering just the most significant (alpha = 0.001) and the dispersals...
multiAOV2 <- inner_join(multiAOV2, disp.df2, by = "Variable") %>%
  select(Variable, DF, R2, 'F', pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

multiAOV2 %>%
  filter(adjP < 0.001) %>%
  arrange(-R2)

#It appears platform explains the most of the dissmilarities in the matrix when
#tribe numbers are controlled
#otherwise there are 2 nonsignificantly dispersed factors that are significant:
#Host family and continent
#sociality (3 level) is still significant just not at the level of < 0.001
#(its 0.000103...)
#write up
write.table(multiAOV2, "output/Prokaryote/Composition_Analysis/TribeReduced/Multi_adonis2_results.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

#my main concern continues to be sociality, but will see if I can factor these in
#note: to nest, the factor that is nested within another factor goes last, ie if 
#species is nested in site then it is written Site/Species. For permanova, 
#you then need to designate strata as Site (https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata)
#set permutations to 9999
p <- 9999

#sociality and tribe
soc.tri2 <- adonis2(dist2 ~ Sociality + Tribe, 
                   data = metdf2, 
                   permutations = p)
soc.tri2

#interaction likely won't work as every tribe = 1 sociality.

#family?
soc.fam2 <- adonis2(dist2 ~ Sociality + Family,
                   data = metdf2,
                   permutations = p)
soc.fam2

soc.fam.nest2 <- adonis2(dist2 ~ Family / Sociality,
                        data = metdf2,
                        permutations = p,
                        strata = metdf2$Family)
soc.fam.nest2

#throwing all phylo in
soc.phy2 <- adonis2(dist2 ~ Sociality + Family/Tribe,
                   data = metdf2,
                   permutations = p,
                   strata = metdf2$Family)
soc.phy2
#when controlling for phylogeny by stratifying by family, Tribe remains a significantly 
#differing factor but sociality does not. Though the levels of variance explained still 
#remain very low
soc.phy.int2 <- adonis2(dist2 ~ Sociality * Family/Tribe,
                       data = metdf2,
                       permutations = p,
                       strata = metdf2$Family)
soc.phy.int2
#CONSIDER, each tribe has only 1 sociality. Each family is majoritively one sociality
#are these therefore all nested ? 
#here, family appears to be the most explanatory variable and it is not significantly
#dispersed
#continent is the only other variable (other than sociality) that does not have 
#significantly difference dispersions
#what about continent?
soc.con2 <- adonis2(dist2 ~ Sociality + Continent,
                   data = metdf2,
                   permutations = p)
soc.con2

soc.con.nest2 <- adonis2(dist2 ~ Continent/Sociality,
                        data = metdf2,
                        permutations = p, 
                        strata = metdf2$Continent)
soc.con.nest2

soc.con.int2 <- adonis2(dist2 ~ Sociality*Continent,
                       data = metdf2,
                       permutations = p)
soc.con.int2

#will also consider year of collection
soc.yea2 <- adonis2(dist2 ~ Sociality + YearCollected,
                   data = metdf2,
                   permutations = p)
soc.yea2


#considering all 
soc.all2 <- adonis2(dist2 ~ Sociality + Family/Tribe + YearCollected + Continent,
                   data = metdf2,
                   permutations = p,
                   strata = metdf2$Family)
soc.all2
#Sociality is no longer significant when these are all considered together
#and this model has the least residual errors. It makes senses really, microbial 
#community should be first and foremost affected by location (especially in non-eusocial
#species), and location/phylogeny are interchangeable in non-global species. Seasonality
#would also affect the types of microbe about (year collected) plus as too would whether
#these were lab reared, caught and put in a lab or else field caught (which may be
#a factor of the project, not actually the year). 
####Pairwise ####
#sociality
cbn <- combn(x=unique(metdf2$Sociality), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Sociality == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Sociality, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.soc2 <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.soc2 %>%
  filter(padj < 0.05)
#like with the anosim, and version 1,
#the only significant difference is between eusocial and 
#polymorphic

write.table(pairwiseAOV.soc2, 
            "output/Prokaryote/Composition_Analysis/TribeReduced/PairwisedANOVA_Sociality.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

#in version 2, there are other non-significantly dispersed factors to consider
#continent
cbn <- combn(x=unique(metdf2$Continent), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Continent == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Continent, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.con2 <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.con2 %>%
  filter(padj < 0.05)

#the only significant difference between continents was between oceania and europe
write.table(pairwiseAOV.con2, 
            "output/Prokaryote/Composition_Analysis/TribeReduced/PairwisedANOVA_Continent.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")


#family
cbn <- combn(x=unique(metdf2$Family), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf2$Sample.ID[metdf2$Family == cbn[,i]]
  subcnt <- cnt2[names(cnt2) %in% sub]
  submet <- metdf2[metdf2$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read2, iterations = 1e4)
  test <- adonis2(subdist ~ submet$Family, 
                  data = submet, permutations = 9999)
  pvalue <- c(pvalue, test$`Pr(>F)`[1])
}

p.adj <- p.adjust(pvalue, method = "BH")
pairwiseAOV.fam2 <- cbind.data.frame(t(cbn), pvalue=pvalue, padj=p.adj)
pairwiseAOV.fam2 %>%
  filter(padj < 0.05)

#Halictidae is significantly different from everything and Apidae and
#megachilidae are different from each other
#note: all of halictidae is polymorphic, and all of Megachilidae is solitary
#so what is the real factor driving these differences?
table(metdf2$Sociality, metdf2$Family)


write.table(pairwiseAOV.fam2, 
            "output/Prokaryote/Composition_Analysis/TribeReduced/PairwisedANOVA_Family.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")
