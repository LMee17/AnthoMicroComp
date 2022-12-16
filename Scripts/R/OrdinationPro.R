#8th December 2022
#Making Ordination Plots; Simple; Prokaryote

#Load Libraries
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)

set.seed(1119)

dir.create("output/Plots/Ordinations/")
dir.create("output/Plots/Ordinations/Prokaryote/")

##Load Metadata####
#sample metadata
met <- read.table("input/Metadata/MetaData_Edit_Nov22.tsv", 
                  sep = "\t", header = T)
#microbial metadata (species level)
microkey <- read.table("input/Metadata/MicrobialSpeciesKey.tsv", 
                       sep = "\t", header = T)
#microbial metadata (genus and family level)
genkey <- read.table("input/Metadata/GenusFamilySpeciesKey_Nov22.tsv", 
                       sep = "\t", header = T)


##Step One: Convert raw reads to relative abundances####
#read in raw read count tables
#for now, just genus-level for Euk, Pro and Viral tables
rawpath <- "input/Counts/RawReads/"
cnt.files <- list.files(path = rawpath)[c(4,6,8)]

data.raw <- list(length = length(cnt.files))

#prepare to remove dna samples from the data
dna <- met$Sample.ID[met$NucleotideType == "DNA"]

#read in files and remove DNA samples
#remove any samples that have less than 100 reads
for (i in 1:length(cnt.files)){
  data.raw[[i]] <- read.table(paste(rawpath, cnt.files[i], sep = ""))
  data.raw[[i]] <- data.raw[[i]][, !names(data.raw[[i]]) %in% dna]
  data.raw[[i]] <- data.raw[[i]][, ! colSums(data.raw[[i]]) < 100]
}

#prepare to populate with relative abundance matrices
data.rel <- vector(mode = "list", length = length(data.raw))

for (l in 1:length(data.rel)){
  data.rel[[l]] <- matrix(nrow = nrow(data.raw[[l]]), ncol = ncol(data.raw[[l]]))
  samples <- names(data.raw[[l]])
  for (i in 1:length(samples)){
    x <- vector(length = nrow(data.rel[[l]]))
    tot <- sum(data.raw[[l]][,i])
    if (tot == 0){
      x[1:length(x)] <- 0
    } else {
      for(j in 1:length(x)){
        x[j] <- (data.raw[[l]][j,i] / tot) * 100
      }
    }
    data.rel[[l]][,i] <- x
  }
  data.rel[[l]] <- as.data.frame(data.rel[[l]])
  names(data.rel[[l]]) <- names(data.raw[[l]])
  rownames(data.rel[[l]]) <- rownames(data.raw[[l]])
}

#prepare to split list into dataframes in the global environment
names(data.rel) <- c("euk", "pro", "vir")
list2env(data.rel, envir = .GlobalEnv)

##Step Two: Filtering Bacterial Data #####
#starting with just bacteria
#keep only microbes that appear in at least 5% of the samples
#!keep only samples that have 3 or more microbial genera 
#remove rows of 0 counts
pro <- pro[!rowSums(pro) == 0,]

#make an incidence matrix for ease of summing rows and columns
pro.inc <- pro
pro.inc[pro.inc > 0] <- 1

#get microbial IDs and sample IDs of rows and columns that are to be removed
#store what 5% of samples are
five <- ncol(pro) * 0.05
pro <- pro[rownames(pro.inc)[rowSums(pro.inc) > 13.25],]
#final check for empty columns 
pro <- pro[,!colSums(pro) == 0]

#store all sample IDs that have less than 3 microbials present
lt3 <- apply(pro, 2, function(x) sum(x > 0)) < 3
pro <- pro[, !names(pro) %in% names(lt3[lt3 == TRUE])]

#write up
dir.create("input/Counts/RelativeAbundance/")
write.table(pro, "input/Counts/RelativeAbundance/Prokaryote_RA_filtered_Dec22.tsv")

#sanity check: these filters should have gotten rid of any clear outliers
#i.e. columns or rows that are mostly zero
#to check, do a quick and dirty NMDS plot using Bray-Curtis
pro.t <- t(pro)

pro.t %>%
  metaMDS(trace = F) %>%
  ordiplot(type = "none") %>%
  text("sites")

#plot are in the same space but looking at them in this way is not super informative
#but at least there's no outliers messing up the picture

##Step Three: Determining Ordination with Bacterial Results ####
###Principal Components Analysis####
#not using scales as the variables are not on different scales
pro.pca <- rda(pro.t, scale = FALSE)

# Now plot a bar plot of relative eigenvalues. 
#This is the percentage variance explained by each axis
barplot(as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig)) 
# over 25% of variance explained by first PC (would probably expect a good pca to be up
# to 60 but there is so much noise in this dataset)

# Calculate the percent of variance explained by first two axes
sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[1:2]) # 43%
sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[1:4]) # 66%

#this plots a pca with points per microbial genus taxa, and points per
#samples
plot(pro.pca)

#looking at just species
sam.pro.pca <- pro.pca$CA$u
#hmm.. if only I could see what this means.

###Playing with PCA Plots ####
#make the above object more amenable for playing with
tmp <- as.data.frame(sam.pro.pca)
tmp$Sample <- rownames(tmp)
#add sample metadata
propc.plot <- merge(tmp, met, by.x = "Sample", by.y = "Sample.ID")

#looking at PC1 and PC2 ... what % of the variance do they represent?
one <- round(sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[1])*100, digits = 2)
two <- round(sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[2])*100, digits = 2)

myPal <- c("#ff8c00","#008080", "#00bfff", "#7300e6")

#look at socialities
ggplot(data = propc.plot, aes(x = PC1, y = PC2, colour = Sociality))+
  geom_point(aes(fill = Sociality)) + 
  scale_colour_manual(values = myPal) +
  labs(title = "Prokaryotic PCA", 
       x = paste("PC1 (", one, "%)", sep = ""),
       y = paste("PC2 (", two, "%)", sep = ""))
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sociality_PC1PC2.pdf")

#tissue type
ggplot(data = propc.plot, aes(x = PC1, y = PC2, colour = Tissue3))+
  geom_point() + 
  #scale_colour_manual(values = myPal5) +
  scale_colour_brewer(palette = "Dark2")+
  labs(title = "Prokaryotic PCA", 
       x = paste("PC1 (", one, "%)", sep = ""),
       y = paste("PC2 (", two, "%)", sep = ""),
       colour = "Tissue Type")
#save
ggsave("output/Plots/Ordinations/Prokaryote/TissueType_PC1PC2.pdf")


#Continental Location
ggplot(data = propc.plot, aes(x = PC1, y = PC2, colour = Continent))+
  geom_point() + 
  scale_colour_brewer(palette = "Dark2")+
  labs(title = "Prokaryotic PCA", 
       x = paste("PC1 (", one, "%)", sep = ""),
       y = paste("PC2 (", two, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Continent_PC1PC2.pdf")

#Sex
ggplot(data = propc.plot, aes(x = PC1, y = PC2, colour = Sex))+
  geom_point() + 
  scale_colour_manual(values = myPal)+
  labs(title = "Prokaryotic PCA", 
       x = paste("PC1 (", one, "%)", sep = ""),
       y = paste("PC2 (", two, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sex_PC1PC2.pdf")

#looking at pc3 and pc4
three <- round(sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[3])*100, digits = 2)
four <- round(sum((as.vector(pro.pca$CA$eig)/sum(pro.pca$CA$eig))[4])*100, digits = 2)

#look at socialities
ggplot(data = propc.plot, aes(x = PC3, y = PC4, colour = Sociality))+
  geom_point(aes(fill = Sociality)) + 
  scale_colour_manual(values = myPal)+
  labs(title = "Prokaryotic PCA", 
       x = paste("PC3 (", three, "%)", sep = ""),
       y = paste("PC4 (", four, "%)", sep = ""))
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sociality_PC3PC4.pdf")

#tissue type
ggplot(data = propc.plot, aes(x = PC3, y = PC4, colour = Tissue3))+
  geom_point() + 
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Prokaryotic PCA", 
       x = paste("PC3 (", three, "%)", sep = ""),
       y = paste("PC4 (", four, "%)", sep = ""),
       colour = "Tissue Type")
#save
ggsave("output/Plots/Ordinations/Prokaryote/TissueType_PC3PC4.pdf")

#Continental Location
ggplot(data = propc.plot, aes(x = PC3, y = PC4, colour = Continent))+
  geom_point() + 
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Prokaryotic PCA", 
       x = paste("PC3 (", three, "%)", sep = ""),
       y = paste("PC4 (", four, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Continent_PC3PC4.pdf")

#Sex / Queen 
ggplot(data = propc.plot, aes(x = PC3, y = PC4, colour = Sex))+
  geom_point() + 
  scale_fill_manual(values = myPal) +
  labs(title = "Prokaryotic PCA", 
       x = paste("PC3 (", three, "%)", sep = ""),
       y = paste("PC4 (", four, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sex_PC3PC4.pdf")

###Bi Plots#####
#see if there's any sense in trying a biplot with this many samples
#biplot of axis 1 verus axis 2
biplot(pro.pca, choices = c(1,2), type = c("text", "points"),
       xlim = c(-5,10))

#hmm.
pdf("output/Plots/Ordinations/Prokaryote/BiPlot_BacterialGenus_PC1PC2.pdf")
  biplot(pro.pca, choices = c(1,2), type = c("text", "points"))
dev.off()

#Looks to be 3 microbial genera that are causing much of the separation
big3 <- c("GRH1", "GRH5", "GRH2")
genkey[genkey$GenHitID %in% big3,]
#GRH1 = Pseudomonas; GRH2 = Acinetobacter; GRH5 = Escherichia

##Principal Coordinates Analysis ####
# First step is to calculate a distance matrix. 
# Here we use Bray-Curtis distance metric
dist <- vegdist(pro.t,  method = "bray")

# calculate the princ
pro.pcoa <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(pro.pcoa$values$Relative_eig[1:10])

# Some distance measures may result in negative eigenvalues. 
#In that case, add a correction:
pro.pcoa$values[pro.pcoa$values[,"Eigenvalues"] < 0,]
#there are negative eigenvalues
pro.pcoa <- pcoa(dist, correction = "cailliez")
#a "good" pcoa explains ~ 50% of the variance/distance within 3 axes
sum((pro.pcoa$values$Corr_eig / sum(pro.pcoa$values$Corr_eig))[1:3])*100
#I'm getting 25 % but ah well

# Plot your results
biplot.pcoa(pro.pcoa, pro.t)

#Ok, let's make this more legible
###Playing with PCoA Plots####
#make the above object more amenable for playing with
tmp <- as.data.frame(pro.pcoa$vectors)
tmp$Sample <- rownames(tmp)

#add sample metadata
propco.plot <- merge(tmp, met, by.x = "Sample", by.y = "Sample.ID")

#looking at axis1 and axis 2 ... 
one <- round(pro.pcoa$values$Corr_eig[1] / sum(pro.pcoa$values$Corr_eig)*100, digits = 2)
two <- round(pro.pcoa$values$Corr_eig[2] / sum(pro.pcoa$values$Corr_eig)*100, digits = 2)

#look at socialities
ggplot(data = propco.plot, aes(x = Axis.1, y = Axis.2, colour = Sociality))+
  geom_point(aes(fill = Sociality)) + 
  scale_colour_manual(values = myPal) +
  labs(title = "Prokaryotic PCoA", 
       x = paste("Axis 1 (", one, "%)", sep = ""),
       y = paste("Axis 2 (", two, "%)", sep = ""))
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sociality_PCoA_Axis12.pdf")

#tissue type
ggplot(data = propco.plot, aes(x = Axis.1, y = Axis.2, colour = Tissue3))+
  geom_point() + 
  #scale_colour_manual(values = myPal5) +
  scale_colour_brewer(palette = "Dark2")+
  labs(title = "Prokaryotic PCA", 
       x = paste("Axis 1 (", one, "%)", sep = ""),
       y = paste("Axis 2 (", two, "%)", sep = ""),
       colour = "Tissue Type")
#save
ggsave("output/Plots/Ordinations/Prokaryote/TissueType_PCoA_Axis12.pdf")


#Continental Location
ggplot(data = propco.plot, aes(x = Axis.1, y = Axis.2, colour = Continent))+
  geom_point() + 
  scale_colour_brewer(palette = "Dark2")+
  labs(title = "Prokaryotic PCA", 
       x = paste("Axis 1 (", one, "%)", sep = ""),
       y = paste("Axis 2 (", two, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Continent_PCoA_Axis12.pdf")

#Sex
ggplot(data = propco.plot, aes(x = Axis.1, y = Axis.2, colour = Sex))+
  geom_point() + 
  scale_colour_manual(values = myPal)+
  labs(title = "Prokaryotic PCA", 
       x = paste("Axis 1 (", one, "%)", sep = ""),
       y = paste("Axis 2 (", two, "%)", sep = ""),
       colour = "Continent of Collection")
#save
ggsave("output/Plots/Ordinations/Prokaryote/Sex_PCoA_Axes12.pdf")

#as the most explanatory looks to be maybe sociality and sex here, let's combine
ggplot(data = propco.plot, 
       aes(x = Axis.1, y = Axis.2, colour = Sociality, shape = Sex))+
  geom_point(aes(fill = Sociality)) + 
  scale_colour_manual(values = myPal) +
  labs(title = "Prokaryotic PCoA", 
       x = paste("Axis 1 (", one, "%)", sep = ""),
       y = paste("Axis 2 (", two, "%)", sep = ""))
#save
ggsave("output/Plots/Ordinations/Prokaryote/SocialityNSex_PCoA_Axis12.pdf")

##NonMetric MultiDimensional Scaling####
#using the bray-curtis dissimilarity matrix from before
#make a function that automatically performs an NMDS for 1-10 dimensions
#and plots the n of dimensions versus the stress of the run
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#now run
#NMDS.scree(dist)
# ok so that's very cute but not very helpful in terms of keeping a good track of what's
# going on. I would like the plot saved and a table of stress versus number of dimensions
# for my records.

#distmat = distance/dissimilarity matrix, title = file name to be saved
scree.NMDS <- function(distmat, title){
  strs <- vector(length = 10)
  for(i in 1:10){
    mds <- metaMDS(dist, distance = "bray", autotransform = F, k = i)
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
  ggsave(paste("output/Plots/Ordinations/Prokaryote/",
               title, "_StressVDimensionScreePlot.pdf", sep = ""))
  return(out)
}

pro.scree <- scree.NMDS(dist, "Prokaryote")
pro.scree

#below 0.1 is good, 0.05 is preferred
#there is a trade-off between number of dimensions and stress, and so I'm going to 
#go with 6, as its the lowest number of dimensions that have stress below 0.1

pro.nmds <- metaMDS(dist, k = 6, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(pro.nmds)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#plot results
#requires original input to run through in one step to keep metadata
pro.nmds2 <- metaMDS(pro.t, k = 6, 
                     trymax = 100, trace = F, autotransform = F, 
                     distance = "bray")

plot(pro.nmds2, display = "sites")

#not very informative
#to apply sample metadata, I need a dataframe with only the factors I want to look at
#and sample ids as rownames
#not sure if this will work with factors other than numbers ....
met2 <- met[,c(2:6,9,10,13,15,16)]
rownames(met2) <- met$Sample.ID
met2 <- met2[rownames(met2) %in% rownames(pro.t),]

#trying to plot metadata ontop of the simple plot that this tutorial does, doesn't work
#perhaps cos most of my metadata is not numeric
#can try with polygons though
met <- met[met$Sample.ID %in% rownames(pro.t),]
met.grp <- data.frame(rownames = met$Sample.ID,
                      Sociality = met$Sociality,
                      Tissue = met$Tissue3,
                      Location = met$Continent)

ordihull(pro.nmds, draw = "polygon",
         groups = met.grp$Sociality, col = c("orange", "green", "blue"),
         label = F)

#god I hate this ugly, ugly graphs
#though there is complete overlap between sociality
#tissue type and continent doesn't work (too many polygons ?)


