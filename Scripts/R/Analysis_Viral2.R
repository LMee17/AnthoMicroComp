#13th January 2023
#Analysis : Viral Count Data
#See Prokaryote Script for Resource details
set.seed(1436)

dir.create("output/Viral/")
dir.create("output/Viral/Composition_Analysis/")

## Libraries ####
library(tidyverse)
library(vegan)
library(ape)
library(gg3D)

## Functions ####

#distmat = distance/dissimilarity matrix, 
#title = analysis level for saving (Viral/Viral/virus)
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
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Feb23.tsv",
                  sep = "\t", header = T)
#microbial metadata
taxkey <- read.table("input/Metadata/Taxa_HitKey_Dec22.tsv",
                     sep ="\t", header = T)
#taxomy
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t", quote = "")

##Step One: Readying Counts and MetaData ####
#count table
cnt <- read.table("input/Counts/Viral_RawCounts.tsv")

#remove any empty taxa rows that appear in less than 5% samples
lt5 <- apply(cnt, 1, function(x) sum(x > 0)) <= ncol(cnt)*0.05
cnt <- cnt[!rownames(cnt) %in% names(lt5)[lt5 == T],]
#remove any samples that have < 100 reads
cnt <- cnt[colSums(cnt)>100]

#write up for ease of later iterations
write.table(cnt, "input/Counts/Viral_Filtered_Raw.tsv",
            col.names = T, row.names = T, quote = F)

#for PCA / negative binomial regression, I will need a cnt matrix of relative abundances
cnt.rel <- apply(cnt,2, FUN=function(x){ x / sum(x)})

#write up
write.table(cnt.rel, "input/Counts/Viral_RelativeAbundance.tsv",
            col.names = T, row.names = T, quote = F)

#this dataset is already VERY small. So will try to do my best to keep socialities
#over anything else
fam.tab <- table(met$Family)
metdf <- met[!met$Family %in% names(fam.tab)[fam.tab < 4],]
#I have more solitary samples than I did with prokaryote analysis!

cnt <- cnt[,names(cnt) %in% metdf$Sample.ID]
cnt.rel <- as.data.frame(cnt.rel)
cnt.rel <- cnt.rel[,names(cnt.rel) %in% metdf$Sample.ID]

write.table(cnt, "input/Counts/Viral_Filtered_FamilyReduced_Raw.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

#sample metadata
#keep only the samples that are in the table
metdf <- met[met$Sample.ID %in% colnames(cnt),]
table(metdf$Sociality)
#there's no point keeping 1 polymorphic
metdf <- metdf[!metdf$Sociality == "Polymorphic",]
cnt <- cnt[, names(cnt) %in% metdf$Sample.ID]

##Step Two: Principal Components ####
#most of these steps assume that rows are samples / sites, not columns
#using relative abundance data
rel.t <- t(cnt.rel)

#not using scales as the variables are not on different scales
pca <- rda(rel.t, scale = FALSE)

# Now plot a bar plot of relative eigenvalues. 
#This is the percentage variance explained by each axis
barplot(as.vector(pca$CA$eig)/sum(pca$CA$eig)) 
# over 25% of variance explained by first PC (would probably expect a good pca to be up
# to 60 but there is so much noise in this dataset)

# Calculate the percent of variance explained by first two to four axes
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:2]) # 44.7%
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:4]) # 67.8%

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
vars <- c("Sociality", "Family", "Tribe", "Continent", "Tissue2")
#and better labels to use
varslabs <- c("Sociality", "Host Family", "Host Tribe", "Continent", "Tissue Type")

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
  ggsave(paste("output/Viral/Composition_Analysis/PCA_PC1PC2_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Viral/Composition_Analysis/PCA_PC3PC4_", title, ".pdf", sep = ""))
}

#Looking at three axes at once
sum((as.vector(pca$CA$eig)/sum(pca$CA$eig))[1:3]) # 58.1%
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
  ggsave(paste("output/Viral/Composition_Analysis/PCA_3D_", 
               title, ".pdf", sep = ""))
}

##Step Three: BiPlot ####
biplot(pca, choices = c(1,2), type = c("text", "points"))

pdf("output/Viral/Composition_Analysis/BiPlot_PC1PC2.pdf")
  biplot(pca, choices = c(1,2), type = c("text", "points"))
dev.off()
#without text
pdf("output/Viral/Composition_Analysis/BiPlot_PC1PC2_Plain.pdf")
  biplot(pca, choices = c(1,2), type = c("points"))
dev.off()

#there's three clear taxa that affect these ordinations (plus a possible fourth)
vip <- c("TRH788", "TRH777", "TRH778")
taxkey[taxkey$taxID %in% vip,]
#Parvoviridae, Iflavirus, Triatovirus

##Step Four: Principal Coordinates Analysis ####
# First step is to calculate a distance matrix. 
# using avgdist from vegan to rarefy counts to the lowest number of the smallest sample
# set by sociality (to try and conserve sample number) with 10,000 iterations

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

#write up
d <- as.matrix(dist)
write.table(d, "output/Viral/Composition_Analysis/BC_DistanceMatrix.tsv",
            sep = "\t", row.names = T, col.names = T, quote = F)

# The following sampling units were removed because they were below sampling depth: CA_N_007, CA_N_030, SRR11426433, SRR12053407, SRR12053408, SRR12053426, SRR13442822, SRR15496107, SRR3948521, SRR3948533, SRR3948536, SRR3948544, SRR3948562, SRR3948566, SRR3948570, SRR3948575, SRR8754044, SRR8867394
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
#I'm getting 50 % .... well well

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
  ggsave(paste("output/Viral/Composition_Analysis/PCoA_Ax1Ax2_", title, ".pdf", sep = ""))
}

##Step Five: NonMetric Multidimensional Scaling ####
#there's a problem further down the line
#in the current iteration there are three singletons that fuck with 
#computing dispersions
#1 f.eusocial, 1 Halictidae and Andrenidae and 1 South American
torem <- c()
torem <- c(metdf$Sample.ID[metdf$Sociality == "F. Eusocial"], torem) #this is also the halictid
torem <- c(metdf$Sample.ID[metdf$Family == "Andrenidae"], torem) #this is also the halictid
torem <- c(metdf$Sample.ID[metdf$Continent == "South America"], torem) #this is also the halictid

metdf <- metdf[!metdf$Sample.ID %in% torem,]
cnt <- cnt[, ! names(cnt) %in% torem]

cnt.t <- t(cnt)
dist <- avgdist(cnt.t, dmethod="bray", sample=low_read, iterations = 10000)

pro.scree <- scree.NMDS(dist, "Viral", "Viral")
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


nmds <- metaMDS(dist, k = kval, trymax = 100, trace = F, distances = "bray")

#check the results of this NMDS with a stress plot
#here I want the scatter to not fall too far from the plotted line
stressplot(nmds)
#both R2 are above .9, same as examples in tutorials, so I'm going to guess that 
#this is good

#save
pdf("output/Viral/Composition_Analysis/NMDS_StressPlot.pdf")
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
write.table(nmds.plot, "output/Viral/Composition_Analysis/NMDS_meta.tsv",
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
  ggsave(paste("output/Viral/Composition_Analysis/NMDS_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Viral/Composition_Analysis/NMDS_Ellipses_", title, ".pdf", sep = ""))
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
  ggsave(paste("output/Viral/Composition_Analysis/NMDS_Polygons_", title, ".pdf", sep = ""))
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
                 aes(fill = as.factor(cen[,1])),
                 show.legend = FALSE) + 
      scale_fill_manual(values = myPal) +
      labs(colour = paste(varslabs[i])) +
      theme_classic() +
      theme(strip.background = element_blank())
    title <- gsub(" ", "_", varslabs[i])
    ggsave(paste("output/Viral/Composition_Analysis/NMDS_", title, 
                 "_CentroidEllipsis.pdf", sep = ""))
  }
}


##Step Six: Checking Dispersion ####
#if the dispersion of different groups are very variable I'm likely to have false positive
#adonis / anosim results, whatever I decide to pursue
#first need to make sure the # rows and samples in the metadata match those in the 
#distance matrix 
if (nrow(metdf) != length(labels(dist))){
  metdf <- metdf[metdf$Sample.ID %in% labels(dist),]
}

table(metdf$Continent)
table(metdf$Sociality)
table(metdf$Family)

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
#sociality (both versions), year collected, month collected and library selection

#write up
write.table(disp.df, 
            "output/Viral/Composition_Analysis/VariableDispersion_ANOVA.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

##Step Seven: PERMANOVA / Adonis ####
#it won't work with NA values, so I'll need to remove
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

multiAOV <- inner_join(multiAOV, disp.df, by = "Variable") %>%
  select(Variable, DF, R2, 'F', pvalue, adjP, Sig) %>%
  mutate(Significance = ifelse(adjP < 0.05, "*", "")) %>%
  mutate(Dispersion = ifelse(Sig == "Non-Sig", "NonSig Dispersed", "Sig Dispersed")) %>%
  select(-Sig)

multiAOV %>%
  filter(adjP < 0.05) %>%
  arrange(-R2)

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
            "output/Viral/Composition_Analysis/PairwiseAOV_Sociality.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

#family
cbn <- combn(x=unique(metdf$Family[!metdf$Family == "Andrenidae"]), m = 2)
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
            "output/Viral/Composition_Analysis/PairwiseAOV_Family.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

#continent
cbn <- combn(x=unique(metdf$Continent), m = 2)
pvalue <- c()

for(i in 1:ncol(cbn)){
  sub <- metdf$Sample.ID[metdf$Continent == cbn[,i]]
  print(cbn[,i])
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
            "output/Viral/Composition_Analysis/PairwiseAOV_Continent.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")

###Modellling #####
#set permutations to 9999
p <- 9999

#sociality and tribe
soc <- adonis2(dist ~ Sociality, 
                   data = metdf, 
                   permutations = p)
soc
disp <- betadisper(dist, metdf$Sociality)
anova(disp)
#family?
fam <- adonis2(dist ~ Family,
                   data = metdf,
                   permutations = p)
fam
disp <- betadisper(dist, metdf$Family)
anova(disp)
#what about continent?
con <- adonis2(dist ~ Continent,
                   data = metdf,
                   permutations = p)
con

disp <- betadisper(dist, metdf$Continent)
anova(disp)

##SessionLogs####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Analysis_Viral_Jan23.txt")
