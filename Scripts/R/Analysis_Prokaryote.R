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

dir.create("output/Prokaryote/")
dir.create("output/Prokaryote/Composition_Analysis/")

## Libraries ####
library(tidyverse)
library(vegan)
library(ape)
library(gg3D)

## Functions ####

#distmat = distance/dissimilarity matrix, 
#title = analysis level for saving (prokaryote/eukaryote/virus)
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

#reduce the counts by tribes = ie keep only tribes that have at least 4 members
tri.tab <- table(met$Tribe)
metdf <- met[!met$Tribe %in% names(tri.tab)[tri.tab < 4],]

cnt <- cnt[,names(cnt) %in% metdf$Sample.ID]
cnt.rel <- as.data.frame(cnt.rel)
cnt.rel <- cnt.rel[,names(cnt.rel) %in% metdf$Sample.ID]

write.table(cnt, "input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

#sample metadata
#keep only the samples that are in the table
metdf <- met[met$Sample.ID %in% colnames(cnt),]

#add a second sociality field (complex versus primitive eusocial)
#this is just to perhaps balance the sample sets more fairly
metdf$Sociality2 <- metdf$Sociality
metdf$Sociality2[metdf$Tribe == "Apini"] <- "Complex Eusocial"
metdf$Sociality2[metdf$Tribe == "Meliponini"] <- "Complex Eusocial"
metdf$Sociality2[metdf$Tribe == "Bombini"] <- "Primitive Eusocial"

##Step Two: Principal Components ####
#most of these steps assume that rows are samples / sites, not columns
#using relative abundance data
rel.t <- t(cnt.rel)

#not using scales as the variables are not on different scales
pca <- rda(rel.t, scale = FALSE)

# Now plot a bar plot of relative eigenvalues. 
#This is the percentage variance explained by each axis
barplot(as.vector(pca$CA$eig)/sum(pca$CA$eig)) 
# over 30% of variance explained by first PC (would probably expect a good pca to be up
# to 60 but there is so much noise in this dataset)

# Calculate the percent of variance explained by first two to four axes
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:2]) # 43.5%
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:4]) # 64.0%

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
  ggsave(paste("output/Prokaryote/Composition_Analysis/PCA_PC1PC2_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Prokaryote/Composition_Analysis/PCA_PC3PC4_", title, ".pdf", sep = ""))
}

#Looking at three axes at once
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:3]) # 55.3%
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
  ggsave(paste("output/Prokaryote/Composition_Analysis/PCA_3D_", title, ".pdf", sep = ""))
}

##Step Three: BiPlot ####
biplot(pca, choices = c(1,2), type = c("text", "points"))

pdf("output/Prokaryote/Composition_Analysis/BiPlot_PC1PC2.pdf")
  biplot(pca, choices = c(1,2), type = c("text", "points"))
dev.off()
#without text
pdf("output/Prokaryote/Composition_Analysis/BiPlot_PC1PC2_Plain.pdf")
  biplot(pca, choices = c(1,2), type = c("points"))
dev.off()

#there's three clear taxa that affect these ordinations (plus a possible fourth)
vip <- c("TRH2", "TRH7", "TRH46", "TRH29")
taxkey[taxkey$taxID %in% vip,]
#Gilliamella (yay), Escherichia, Ralstonia and Chryseobacterium
#Ralstonia and Chryseobacterium downwards, Gilliamella up, and 
#Escherichia to the right

##Step Four: Principal Coordinates Analysis ####
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
#The following sampling units were removed because they were below sampling depth: CA_N_012, SRR11440494, SRR12053407, SRR12053410, SRR12053413, SRR12053415, SRR12053416, SRR12053423, SRR12053424, SRR12053425, SRR12053426, SRR12053427, SRR12053428, SRR12053429, SRR12053435, SRR12053437, SRR12527935, SRR12527963, SRR2001621, SRR2396653, SRR3383868, SRR3948521, SRR3948548

#write up
d <- as.matrix(dist)
write.table(d, "output/Prokaryote/Composition_Analysis/BC_DistanceMatrix.tsv",
            sep = "\t", row.names = T, col.names = T, quote = F)

# calculate the principal coordinates
pco <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(pco$values$Relative_eig[1:10])
#again, first eigenvalue explains over 25% of the variance

# Some distance measures may result in negative eigenvalues. 
# negative eiganvalues need to be corrected
if (any(pco$values[,"Eigenvalues"] < 0) == TRUE){
  print("Negative Eigenvalues present")
  #administer correction
  pco <-pcoa(dist, correction = "cailliez")
}
#a "good" pcoa explains ~ 50% of the variance/distance within 3 axes
sum((pco$values$Corr_eig / sum(pco$values$Corr_eig))[1:3])*100
#I'm getting 37 % but ah well

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
  ggsave(paste("output/Prokaryote/Composition_Analysis/PCoA_Ax1Ax2_", title, ".pdf", sep = ""))
}

##Step Five: NonMetric Multidimensional Scaling ####
pro.scree <- scree.NMDS(dist, "Prokaryote", "Prokaryote")
pro.scree

#below 0.05 is excellent so going for 7 dimensions though I could reduce it to 4 as well 
#as this is still good
nmds <- metaMDS(dist, k = 6, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(nmds)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#save
pdf("output/Prokaryote/Composition_Analysis/NMDS_StressPlot.pdf")
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
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/NMDS_meta.tsv",
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
  ggsave(paste("output/Prokaryote/Composition_Analysis/NMDS_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Prokaryote/Composition_Analysis/NMDS_Ellipses_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Prokaryote/Composition_Analysis/NMDS_Polygons_", title, ".pdf", sep = ""))
}

#getting centroid plots for these for further exploration
#to calculate centroids I need more than / equal to 4 samples per group
#the rarefaction step in producing the matrix removed samples and so I need to 
#check what I have left for each factor (post distance)
metdf.pd <- metdf[metdf$Sample.ID %in% labels(dist),]

#sociality#1
for (i in 1:length(vars)){
  print(vars[i])
  if(any(table(metdf.pd[,vars[i]]) < 4)){
    print("Cannot compute")
  } else{
    v <- unname(unlist(vars[i]))
    print(v)
    cen <- nmds.plot %>%
      group_by(across(all_of(v))) %>%
      summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
      as.data.frame()
    print(cen)
    ggplot(data = nmds.plot, 
           aes(x = NMDS1, y = NMDS2, colour = as.factor(nmds.plot[,vars[i]]))) +
      geom_point() + 
      scale_colour_manual(values = myPal) +
      stat_ellipse(show.legend = FALSE) +
      geom_point(data = cen, size = 3, 
                 shape = 21, colour = "black", 
                 aes(fill = cen[,1]),
                 show.legend = FALSE) + 
      scale_fill_manual(values = myPal) +
      labs(colour = paste(varslabs[i])) +
      theme_classic() +
      theme(strip.background = element_blank())
    title <- gsub(" ", "_", varslabs[i])
    ggsave(paste("output/Prokaryote/Composition_Analysis/NMDS_", title, 
                 "_CentroidEllipsis.pdf", sep = ""))
  }
}

#3D Plot
#A quick exploration of what happens if I bring in another NMDS axis
#only for the more interesting factors discussed above
foi <- c("Sociality", "Sociality2", "Family", "Tribe", "Platform_Spec")
foilabs <- c("Sociality", "Sociality ", "Family", "Tribe", "Sequencing Platform")

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
  ggsave(paste("output/Prokaryote/Composition_Analysis/NMDS_3D_", 
               title, ".pdf", sep = ""))
}

##Step Six: Checking Dispersion ####
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
write.table(disp.df, "output/Prokaryote/Composition_Analysis/VariableDispersion_ANOVA.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

#Sociality (1, with three levels) does not have significant differences in group dispersion
#and nor does Continent, Host Family or Month Collected or Library Selection.

##Step Seven: ANOSIM ####
#Analysis of similarities
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
#host family and continent are significant contributors

#write up
write.table(multiSIM, "output/Prokaryote/Composition_Analysis/Multi_ANOSIM_testResults.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

#ANOSIM is noted to be quite sensitive to differences in dispersions so I'm reluctant
#to discuss the other significant variables just yet

#pairwise: sociality
#I'm still the most interesting in sociality so will still consider it
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

write.table(pairSIM.soc, 
            "output/Prokaryote/Composition_Analysis/PairwiseANOSIM_Sociality.tsv",
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

write.table(pairSIM.con, "output/Prokaryote/Composition_Analysis/PairwiseANOSIM_Continent.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

#pairwise: Family
cbn <- combn(x = unique(metdf$Family), m = 2)
pvalue <- c()

for (i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Family == cbn[,i]]
  subcnt <- cnt[names(cnt) %in% sub]
  submet <- metdf[metdf$Sample.ID %in% names(subcnt),]
  subdist <- avgdist(t(subcnt), dmethod = "bray", sample = low_read, iterations = 1e4)
  sim <- anosim(subdist, submet$Family)
  pvalue <- c(pvalue, sim$signif[1])
}
padj <- p.adjust(pvalue, method = "BH")

pairSIM.fam <- cbind.data.frame(t(cbn), pvalue = pvalue, padj = padj)
pairSIM.fam

#there's a sig difference between Halictidae and everything else

write.table(pairSIM.fam, 
            "output/Prokaryote/Composition_Analysis/PairwiseANOSIM_Continent.tsv",
            col.names = T, row.names = F, quote = F,
            sep = "\t")

##Step Eight: PERMANOVA / Adonis ####
#with permanovas/ adonis2 (vegan function that is just a permutated anova) I can begin
#to control what I consider fixed / random variables, nestedness, strata control etc etc.
#(I can do some of these with anosim but it doesn't look as simple / well documented)
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
#according to this, everything is a significant variable 
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
#Host family follows, However these are basically just phylogeny (perhaps only consider family)
#Continent comes next with 4 level sociality. Sex and Year collected are also 
#still significant. Continent, Family and Sociality are the only factors that 
#would not be sensitive to the dispersion problem.

#write up
write.table(multiAOV, "output/Prokaryote/Composition_Analysis/Multi_adonis2_results.tsv",
            row.names = F, col.names = T, quote = F,
            sep = "\t")

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

write.table(pairwiseAOV.soc, 
            "output/Prokaryote/Composition_Analysis/PairwiseAOV_Sociality.tsv",
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
#differences between Apidae and Megachilidae, and Halicidae and everything

write.table(pairwiseAOV.fam, 
            "output/Prokaryote/Composition_Analysis/PairwiseAOV_Family.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

###Modellling #####

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
#using the tribe reduced metadata, not the further reduced set that occurs with 
#the dissimilarity matrix
met <- met[met$Sample.ID %in% names(cnt),]
eu.cnt <- cnt[,names(cnt) %in% met$Sample.ID[met$Sociality == "Eusocial"]]
eu.low_read <- inner_join(met, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Eusocial") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

eu.dist <- avgdist(t(eu.cnt), dmethod="bray", sample=eu.low_read, iterations = 10000)
eu.stress <- scree.NMDS(eu.dist, "Prokaryote", "Eusocial")
eu.stress %>%
  filter(Dec == "Excellent")

eu.nmds <- metaMDS(dist, k = 6, trymax = 100, trace = F, distances = "bray")

stressplot(eu.nmds)

pdf("output/Prokaryote/Composition_Analysis/Eusocial_NMDS_StressPlot.pdf")
  stressplot(eu.nmds)
dev.off()

ext <- scores(eu.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
meteu <- met %>%
  filter(Sociality == "Eusocial")
eu.plot <- merge(ext, meteu, by.x = "Sample", by.y = "Sample.ID")

#save
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/Eusocial_NMDS_meta.tsv",
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
#platform and sex

write.table(eu.multiAOV, 
            "output/Prokaryote/Composition_Analysis/Eusocial_multiANOVA.tsv",
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
#everything is significantly dispersed, except for Year Collected

#considered continent, platform, tribe, sex, year collected
eu.vars2 <- eu.vars[c(1:2,4:6),]
for (i in 1:nrow(eu.vars2)){
  ggplot(data = eu.plot,
         aes(x = NMDS1, y = NMDS2, colour = as.factor(eu.plot[,eu.vars$Var[i]]))) +
    geom_point() + 
    stat_ellipse(geom = "polygon",
                 aes(fill = as.factor(eu.plot[,eu.vars$Var[i]])),
                 alpha = 0.1) +
    scale_colour_manual(values = myPal) +
    scale_fill_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(eu.vars2$VarLab[i])) +
    guides(fill = "none") +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", eu.vars2$VarLab[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/Eusocial_NMDS_",
               title, ".pdf", sep = ""))
}

####Polymorphic Considerations#####
po.cnt <- cnt[,names(cnt) %in% met$Sample.ID[met$Sociality == "Polymorphic"]]
po.low_read <- inner_join(met, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Polymorphic") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

po.dist <- avgdist(t(po.cnt), dmethod="bray", sample=po.low_read, iterations = 10000)
po.stress <- scree.NMDS(po.dist, "Prokaryote", "Polymorphic")
po.stress %>%
  filter(Dec == "Excellent")

po.nmds <- metaMDS(dist, k = 3, trymax = 100, trace = F, distances = "bray")

stressplot(po.nmds)

pdf("output/Prokaryote/Composition_Analysis/Polymorphic_NMDS_StressPlot.pdf")
  stressplot(po.nmds)
dev.off()

ext <- scores(po.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
metpo <- met %>%
  filter(Sociality == "Polymorphic")
po.plot <- merge(ext, metpo, by.x = "Sample", by.y = "Sample.ID")

#save
write.table(nmds.plot, "output/Prokaryote/Composition_Analysis/Polymorphic_NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Adonis
vars.sim2 <- vars.sim[c(-1,-2),]
torem <- c()
for (i in 1:nrow(vars.sim2)){
  noFact <- length(names(table(metpo[,vars.sim2$Var[i]])))
  if(noFact == 1){
    print(vars.sim2$Var[i])
    torem <- c(torem, vars.sim2$Var[i])
  }
}
po.vars <- vars.sim2[!vars.sim2$Var %in% torem,]
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
#Family (both phylogeny) and Year collected
#Phylogenetic signal explains over 50% of model

write.table(po.multiAOV, "output/Prokaryote/Composition_Analysis/Polymorphic_multiANOVA.tsv",
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
#continent, year collected, host tribe/family are all not significantly dispersed
for (i in 1:nrow(po.vars)){
  ggplot(data = po.plot,
         aes(x = NMDS1, y = NMDS2, colour = as.factor(po.plot[,po.vars$Var[i]]))) +
    geom_point() + 
    stat_ellipse(geom = "polygon",
                 aes(fill = as.factor(po.plot[,po.vars$Var[i]])),
                 alpha = 0.1) +
    scale_colour_manual(values = myPal) +
    scale_fill_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(po.vars$VarLab[i])) +
    guides(fill = "none") +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", po.vars$VarLab[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/Polymorphic_NMDS_",
               title, ".pdf", sep = ""))
}


####Solitary Considerations#####
so.cnt <- cnt[,names(cnt) %in% met$Sample.ID[met$Sociality == "Solitary"]]
so.low_read <- inner_join(met, readCnts, by = c("Sample.ID" = "Sample")) %>%
  filter(Sociality == "Solitary") %>%
  arrange(Tot) %>%
  head(n = 1) %>%
  select(Tot) %>%
  as.numeric()

so.dist <- avgdist(t(so.cnt), dmethod="bray", sample=so.low_read, iterations = 10000)
so.stress <- scree.NMDS(so.dist, "Prokaryote", "Solitary")
so.stress %>%
  filter(Dec == "Excellent")
#Warning about the number of samples ... it may not be enought to be informative (stress
# is very close to 0)
#1 seems too low to use so I'm going for 2
so.nmds <- metaMDS(dist, k = 2, trymax = 100, trace = F, distances = "bray")

stressplot(so.nmds)
#still ok r plot ... though it's not the best

pdf("output/Prokaryote/Composition_Analysis/Solitary_NMDS_StressPlot.pdf")
  stressplot(so.nmds)
dev.off()

ext <- scores(so.nmds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
metso <- met %>%
  filter(Sociality == "Solitary")
so.plot <- merge(ext, metso, by.x = "Sample", by.y = "Sample.ID")

#save
write.table(so.plot, "output/Prokaryote/Composition_Analysis/Solitary_NMDS_meta.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Adonis
torem <- c()
for (i in 1:nrow(vars.sim2)){
  noFact <- length(names(table(metso[,vars.sim2$Var[i]])))
  if(noFact == 1){
    print(vars.sim2$Var[i])
    torem <- c(torem, vars.sim2$Var[i])
  }
}
so.vars <- vars.sim2[!vars.sim2$Var %in% torem,]
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


write.table(so.multiAOV, "output/Prokaryote/Composition_Analysis/Solitary_multiANOVA.tsv",
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
#pretty much everything is non-significanty dispersed, but that may be due to low
#sample sizes

#just to compare to others, I want to look at sex, family, tribe, continent
#and platform

for (i in 1:nrow(so.vars)){
  ggplot(data = so.plot,
         aes(x = NMDS1, y = NMDS2, colour = as.factor(so.plot[,so.vars$Var[i]]))) +
    geom_point() + 
    stat_ellipse(geom = "polygon",
                 aes(fill = as.factor(so.plot[,so.vars$Var[i]])),
                 alpha = 0.1) +
    scale_colour_manual(values = myPal) +
    scale_fill_manual(values = myPal) +
    labs(x = "NMDS1",
         y = "NMDS2",
         colour = paste(so.vars$VarLab[i])) +
    guides(fill = "none") +
    theme_classic() +
    theme(strip.background = element_blank())
  title <- gsub(" ", "_", so.vars$VarLab[i])
  ggsave(paste("output/Prokaryote/Composition_Analysis/Solitary_NMDS_",
               title, ".pdf", sep = ""))
}

##SessionLogs####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Analysis_Prokaryote_Jan23.txt")
