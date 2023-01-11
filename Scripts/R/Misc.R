#9th January 2023
#Looking at some exploratory plots (prokaryote)

#Libaries
library(tidyverse)
library(ggplot2)
library(reshape)

dir.create("output/Prokaryote/Miscellaneous/")

##Counts and Metadata ####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_RNAOnly_Dec22.tsv",
                  sep = "\t", header = T)
#microbial taxa metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv", 
                  header = T, sep = "\t")

#prokaryote count table (filtered)
pro <- read.table("input/Counts/Prokaryote_Filtered_Raw.tsv")
#prokaryote count table (tribe reduced)
pro2 <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")
#prokaryote relative abundance
pro.rel <- read.table("input/Counts/Prokaryote_RelativeAbundance.tsv")
pro.rel2 <- pro.rel[,names(pro.rel) %in% names(pro2)]

#manually set palette
#colour-blind friendly palette sourced: 
#https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
myPal <- c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
           "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
           "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
           "#9acd32")

##Species Richness vs Sequencing Platform ####
#species richness = total number of species
#per sample, calculate how many different species are present
#probably be easier working with incidence data for a second
pro.inc <- pro
pro.inc[pro.inc > 0] <- 1

specCnt <- data.frame(Sample = names(pro.inc),
                      SpeciesTot = colSums(pro.inc))
pro.Spec.Cnt <- pro.inc %>%
  t() %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") 
names(pro.Spec.Cnt)[2] <- "Total"
pro.Spec.Plot <- inner_join(pro.Spec.Cnt, met, by = c("Sample" = "Sample.ID"))

readCnts <- pro %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  group_by(Sample) %>%
  pivot_longer(-Sample) %>%
  mutate(Tot = sum(value)) %>%
  select(Sample, Tot) %>%
  unique() 
pro.Spec.Plot <- inner_join(pro.Spec.Plot, readCnts, by = "Sample")

#tribe reduced version
pro.Spec.Plot2 <- pro.Spec.Plot[pro.Spec.Plot$Sample %in% names(pro2),]

#platformspec
pro.plots <- list(pro.Spec.Plot, pro.Spec.Plot2)

for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = Platform_Spec, y = Total)) +
    geom_boxplot(aes(colour = Platform_Spec)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Sequencing Platform",
         y = "Number of Unique Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Platform.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Platform_TriRed.pdf")
  }
}

#what about by sociality / genus / family ? 
for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = Sociality, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Sociality",
         y = "Number of Unique Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Sociality.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Sociality_TriRed.pdf")
  }
}

#family
for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = Family, y = Total)) +
    geom_boxplot(aes(colour = Family)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Family",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Family.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Family_TriRed.pdf")
  }
}

#for genus, order genera by sociality
pro.plots[[1]]$Genus <- factor(pro.plots[[1]]$Genus, 
                              levels = c("Xylocopa", "Osmia", "Habropoda",
                                         "Epeolus", "Dufourea", "Colletes",
                                         "Anthophora", "Andrena",
                                         "Megalopta", "Lasioglossum", "Halictus",
                                         "Exoneura", "Euglossa", "Eufriesea",
                                         "Ceratina",
                                         "Tetragonula", "Tetragonisca", "Bombus", 
                                         "Apis"))
pro.plots[[2]]$Genus <- factor(pro.plots[[2]]$Genus,
                               levels = c("Osmia", "Andrena",
                                          "Megalopta",
                                          "Euglossa", "Eufriesea",
                                          "Ceratina",
                                          "Tetragonula", "Tetragonisca", "Bombus", 
                                          "Apis"))

#genus
for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = Genus, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Genus",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus_TriRed.pdf")
  }
}

#tribe
#again, order by sociality
pro.plots[[1]]$Tribe <- factor(pro.plots[[1]]$Tribe,
                               levels = c("Xylocopini", "Rophitini", "Osmiini", 
                                          "Epeolini", "Colletini", "Anthophorini",
                                          "Andrenini",
                                          "Halictini", "Euglossini", "Ceratini", 
                                          "Augochlorini", "Allodapini",
                                          "Meliponini", "Bombini", "Apini"))

pro.plots[[2]]$Tribe <- factor(pro.plots[[2]]$Tribe,
                               levels = c("Osmiini", 
                                          "Andrenini",
                                          "Euglossini", "Ceratini", 
                                          "Augochlorini", "Allodapini",
                                          "Meliponini", "Bombini", "Apini"))

for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = Tribe, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Tribe",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus_TriRed.pdf")
  }
}

#library prep
for (i in 1:length(pro.plots)){
  print(i)
  ggplot(data = pro.plots[[i]], aes(x = LibrarySelection, y = Total)) +
    geom_boxplot(aes(colour = LibrarySelection)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Method of Library Selection",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_LibrarySelection.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_LibrarySelection_TriRed.pdf")
  }
}

#total number of reads versus number of species
for (i in 1:length(pro.plots)){
  ggplot(data = pro.plots[[i]], aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = Sociality, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = Sociality), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_BySoc.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_BySoc_TriRed.pdf")
  }
}

#total number of reads versus number of species (by library selection)
for (i in 1:length(pro.plots)){
  ggplot(data = pro.plots[[i]], aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = LibrarySelection, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = LibrarySelection), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none") +
    labs(colour = "Method of Library Selection")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByLibSel.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByLibSel_TriRed.pdf")
  }
}

#total number of reads versus number of species (by family)
for (i in 1:length(pro.plots)){
  ggplot(data = pro.plots[[i]], aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = Family, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = Family), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none") +
    labs(colour = "Host Family")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByFam.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByFam_TriRed.pdf")
  }
}

#interesting.

##Bacterial Family / Order Facet Plot ####
#first, check taxa included here all have family entries
taxkey <- data.frame(TaxID = rownames(pro))
taxkey <- inner_join(taxkey, tax, by = c("TaxID" = "ID"))
taxkey %>%
  filter(family == "")

taxkey$family[taxkey$TaxID == "TRH117"] <- "Enterobacteriaceae"
taxkey$family[taxkey$TaxID == "TRH119"] <- "Acetobacteraceae"
taxkey$family[taxkey$TaxID == "TRH395"] <- "Methylothermaceae"

#begin to create plot table
pro.micro.plot <- pro.rel %>%
  rownames_to_column() %>%
  melt() 
names(pro.micro.plot) <- c("TaxID", "Sample", "RelAbundance")
pro.micro.plot <- inner_join(pro.micro.plot, taxkey, by = c("TaxID"))
#don't want clashes with the host/micro metadata. Add micro prefix before continuing
for(i in 4:length(names(pro.micro.plot))){
  print(names(pro.micro.plot[i]))
  names(pro.micro.plot)[i] <- sub("^", "Micro_", names(pro.micro.plot[i]))
}
pro.micro.plot <- inner_join(pro.micro.plot, met, by = c("Sample" = "Sample.ID"))

#differentiate between host and micro - adding host prefix
names(pro.micro.plot)[c(13,27:30)] <- sub("^", "Host_", names(pro.micro.plot[c(13,27:30)]))

pro.micro.plot2 <- pro.micro.plot[pro.micro.plot$Sample %in% names(pro2),]

pro.micro.plots <- list(pro.micro.plot, pro.micro.plot2)
#plot
#x = bacterial family / order
#y - abundance
#facet = sociality
for(i in 1:length(pro.micro.plots)){
  ggplot(data = pro.micro.plots[[i]], 
         aes(x = fct_rev(factor(Micro_family)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    #ylim(0,5)+
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialFam_RelAbu_SocialFacet.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialFam_RelAbu_SocialFacet_TriRed.pdf")
  }
}

#by order
for(i in 1:length(pro.micro.plots)){
  ggplot(data = pro.micro.plots[[i]], 
         aes(x = fct_rev(factor(Micro_order)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialOrd_RelAbu_SocialFacet.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialOrd_RelAbu_SocialFacet_TriRed.pdf")
  }
}

#fuck it let's do genus
for(i in 1:length(pro.micro.plots)){
  ggplot(data = pro.micro.plots[[i]], 
         aes(x = fct_rev(factor(Micro_genus)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    #ylim(0,5)+
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialGen_RelAbu_SocialFacet.pdf")
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialGen_RelAbu_SocialFacet_TriRed.pdf")
  }
}

##Prokaryote Heatmap####
heat.plot.inc <- pro.inc %>%
  rownames_to_column() %>%
  melt() 
names(heat.plot.inc) <- c("MicroTax", "Sample", "Inc")
heat.plot.inc <- inner_join(heat.plot.inc, tax, by = c("MicroTax" = "ID")) %>%
        select(MicroTax, Sample, Inc, family)
heat.plot.inc <- inner_join(heat.plot.inc, met, by = c("Sample" = "Sample.ID")) %>%
  select(MicroTax, Sample, Inc, family, Sociality)

eu.samps <- unique(heat.plot.inc$Sample[heat.plot.inc$Sociality == "Eusocial"])
po.samps <- unique(heat.plot.inc$Sample[heat.plot.inc$Sociality == "Polymorphic"])
so.samps <- unique(heat.plot.inc$Sample[heat.plot.inc$Sociality == "Solitary"])
samps <- c(eu.samps,po.samps,so.samps)

heat.plot.inc$Sample <- factor(heat.plot.inc$Sample,
                               levels = samps)

ggplot(data = heat.plot.inc, aes(x = fct_rev(family), y = Sample, fill = factor(Inc))) +
  geom_tile() + 
  coord_flip()

#ok, this is getting me nowhere. Let's try collapsing by species
#for each species per taxa family I want average relative abundance and prevalence

#start with prevalence
#make a table
genera <- unique(met$Genus[met$Sample.ID %in% names(pro)])
micros <- unique(heat.plot.inc$MicroTax)
prevmat <- matrix(nrow = length(micros),
                  ncol = length(genera))

samplabs <- c()
for (i in 1:length(genera)){
  print(genera[i])
  samps <- unique(pro.plots[[1]]$Sample[pro.plots[[1]]$Genus == genera[i]])
  tot <- length(samps)
  label <- paste(genera[i], " (n = ", tot, ")", sep = "")
  samplabs <- c(samplabs, label)
  #prepare to populate with prevalence per microtaxa
  x <- c()
  for (j in 1:length(micros)){
    m.gen <- tax$genus[tax$ID == micros[j]]
    prev <- (sum(pro.inc[micros[j], names(pro.inc) %in% samps]) / tot)*100
    x <- c(x, prev)
  }
  prevmat[,i] <- x
}
prevmat <- as.data.frame(prevmat)
rownames(prevmat) <- micros
names(prevmat) <- genera

hpi <- prevmat %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) 
names(hpi) <- c("TaxID", "Host_Genus", "Prevalence")
hpi <- inner_join(hpi, taxkey, by = c("TaxID")) %>%
  mutate(Taxa = genus) %>%
  select(TaxID, Host_Genus, Prevalence, Taxa) 
hpi <- inner_join(hpi, met, by = c("Host_Genus" = "Genus")) %>%
  select(TaxID, Host_Genus, Prevalence, Taxa, Sociality) %>%
  unique()
#add labels
hpi <- hpi %>%
  mutate(Labels = rep(samplabs, 65)) %>%
  as.data.frame()
#order bee genera by sociality
levelz <- hpi %>%
  select(Host_Genus, Sociality, Labels) %>%
  arrange(Sociality) %>%
  unique() %>%
  select(Labels) %>%
  unlist()

hpi$Labels <- factor(hpi$Labels, levels = c(levelz))

#order microbial taxa so core are together
core <- factor(c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
            "Bartonella", "Apibacter"))
noncore <- factor(unique(hpi$Taxa[!hpi$Taxa %in% core]))

hpi$Taxa <- factor(hpi$Taxa, levels = c(core, noncore))

#factorise prevalence
hpi$Prev2[hpi$Prevalence > 75] <- "> 75%"
hpi$Prev2[hpi$Prevalence <= 75 &
            hpi$Prevalence > 50] <- "51 - 75%"
hpi$Prev2[hpi$Prevalence <= 50 &
            hpi$Prevalence > 30] <- "31 - 50%"
hpi$Prev2[hpi$Prevalence <= 30 &
            hpi$Prevalence > 10] <- "11 - 30%"
hpi$Prev2[hpi$Prevalence < 10 &
            hpi$Prevalence > 0] <- "< 11%"
hpi$Prev2[hpi$Prevalence == 0] <- " "
  
hpi$Prev2 <- factor(hpi$Prev2, levels =c("> 75%", "51 - 75%", "31 - 50%",
                                         "11 - 30%", "< 11%", " "))

ggplot(data = hpi, aes(x = fct_rev(Taxa), y = Labels, fill = Prev2)) + 
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 80, hjust = -.1,
                                   face = "italic"),
        axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("black",
                               "#004545",
                               "#006b6b",
                               "#00cece",
                               "#00f6f6",
                               "white")) +
  labs(x = "Prokaryote Genus",
       y = "Host Genus",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 


ggsave("output/Prokaryote/Miscellaneous/Pro_PrevalenceHeatmap_Factored_HostGenus.pdf",
      height = 28, units = "cm")

#now to do that with tribe reduced :(

genera <- unique(met$Genus[met$Sample.ID %in% names(pro2)])
micros <- unique(heat.plot.inc$MicroTax)
prevmat <- matrix(nrow = length(micros),
                  ncol = length(genera))

samplabs <- c()
for (i in 1:length(genera)){
  print(genera[i])
  samps <- unique(pro.plots[[2]]$Sample[pro.plots[[2]]$Genus == genera[i]])
  tot <- length(samps)
  label <- paste(genera[i], " (n = ", tot, ")", sep = "")
  samplabs <- c(samplabs, label)
  #prepare to populate with prevalence per microtaxa
  x <- c()
  for (j in 1:length(micros)){
    m.gen <- tax$genus[tax$ID == micros[j]]
    prev <- (sum(pro.inc[micros[j], names(pro.inc) %in% samps]) / tot)*100
    x <- c(x, prev)
  }
  prevmat[,i] <- x
}

prevmat <- as.data.frame(prevmat)
rownames(prevmat) <- micros
names(prevmat) <- genera

hpi <- prevmat %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) 
names(hpi) <- c("TaxID", "Host_Genus", "Prevalence")
hpi <- inner_join(hpi, taxkey, by = c("TaxID")) %>%
  mutate(Taxa = genus) %>%
  select(TaxID, Host_Genus, Prevalence, Taxa) 
hpi <- inner_join(hpi, met, by = c("Host_Genus" = "Genus")) %>%
  select(TaxID, Host_Genus, Prevalence, Taxa, Sociality) %>%
  unique()
#add labels
hpi <- hpi %>%
  mutate(Labels = rep(samplabs, 65)) %>%
  as.data.frame()
#order bee genera by sociality
levelz <- hpi %>%
  select(Host_Genus, Sociality, Labels) %>%
  arrange(Sociality) %>%
  unique() %>%
  select(Labels) %>%
  unlist()

hpi$Labels <- factor(hpi$Labels, levels = c(levelz))

#order microbial taxa so core are together
core <- factor(c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
                 "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
                 "Bartonella", "Apibacter"))
noncore <- factor(unique(hpi$Taxa[!hpi$Taxa %in% core]))

hpi$Taxa <- factor(hpi$Taxa, levels = c(core, noncore))

#factorise prevalence
hpi$Prev2[hpi$Prevalence > 75] <- "> 75%"
hpi$Prev2[hpi$Prevalence <= 75 &
            hpi$Prevalence > 50] <- "51 - 75%"
hpi$Prev2[hpi$Prevalence <= 50 &
            hpi$Prevalence > 30] <- "31 - 50%"
hpi$Prev2[hpi$Prevalence <= 30 &
            hpi$Prevalence > 10] <- "11 - 30%"
hpi$Prev2[hpi$Prevalence < 10 &
            hpi$Prevalence > 0] <- "< 11%"
hpi$Prev2[hpi$Prevalence == 0] <- " "

hpi$Prev2 <- factor(hpi$Prev2, levels =c("> 75%", "51 - 75%", "31 - 50%",
                                         "11 - 30%", "< 11%", " "))

ggplot(data = hpi, aes(x = fct_rev(Taxa), y = Labels, fill = Prev2)) + 
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 80, hjust = -.1,
                                   face = "italic"),
        axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("black",
                               "#004545",
                               "#006b6b",
                               "#00cece",
                               "#00f6f6",
                               "white")) +
  labs(x = "Prokaryote Genus",
       y = "Host Genus",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Prokaryote/Miscellaneous/Pro_PrevalenceHeatmap_Factored_HostGenus_TriRed.pdf",
       height = 28, units = "cm")

##Sociality vs SharedMicrobial Species #####
#this will be difficult to code.....
cntz <- list(pro, pro2)

for (i in 1:length(cntz)){
  soc.inc <- cntz[[i]] %>%
    rownames_to_column(var = "MicroTax") %>%
    pivot_longer(-MicroTax) %>%
    mutate(inc = ifelse(value > 0, 1, 0)) %>%
    mutate(Sample = name) %>%
    select(MicroTax, Sample, inc)
  soc.inc2 <- inner_join(soc.inc, met, by = c("Sample" = "Sample.ID")) %>%
    select(MicroTax, Sample, inc, Sociality) %>%
    mutate(SocID = substr(Sociality, 1, 1)) %>%
    select(-Sample) %>%
    filter(inc > 0) %>%
    unique() %>%
    select(-inc) %>%
    arrange(SocID) %>%
    group_by(MicroTax) %>%
    mutate(Combo = paste0(SocID, collapse = "")) %>%
    group_by(Sociality, Combo) %>%
    count(name = "Count")
  #order factors
  soc.inc2$Combo <- factor(soc.inc2$Combo,
                           levels = c("EPS", 
                                      "PS", "ES", "EP",
                                      "E", "P", "S"))
  ggplot(data = soc.inc2,
         aes(x = Sociality, y = Count, fill = (Combo), label = Count))+
    geom_bar(stat = "identity", colour = "grey") + 
    geom_text(size = 3, 
              position = position_stack(vjust = 0.5)) +
    scale_fill_manual(labels = c("All Socialities", "Eusocial and Solitary",
                                 "Eusocial and Polymorphic", "Eusocial only",
                                 "Polymorphic Only"),
                      values = c("#6893b4",
                                 "#42282e",
                                 "#a2ce47",
                                 "#db6d00",
                                 "#009292")) +
    theme(panel.grid =  element_blank()) +
    labs(y = "Number of Prokaryote Genera",
         x = "",
         fill = "")
  if (i == 1){
    ggsave("output/Prokaryote/Miscellaneous/Prok_MicroGenera_bySocCombo.pdf")  
  } else {
    ggsave("output/Prokaryote/Miscellaneous/Prok_MicroGenera_bySocCombo_TriRed.pdf")  
  }
  #get list of microbial taxa for each combination
  soc.inc3 <- inner_join(soc.inc, met, by = c("Sample" = "Sample.ID")) %>%
    select(MicroTax, Sample, inc, Sociality) %>%
    select(-Sample) %>%
    filter(inc > 0) %>%
    unique() %>%
    select(-inc) %>%
    arrange(Sociality) %>%
    group_by(MicroTax) %>%
    mutate(Combo = paste0(Sociality, collapse = ", ")) %>%
    select(MicroTax, Combo)
  names(soc.inc3) <- c("Tax_ID", "Sociality")
  soc.inc3 <- inner_join(soc.inc3, taxkey, by = c("Tax_ID" = "TaxID")) %>%
    mutate(Microbial_Genus = genus) %>%
    select(Tax_ID, Sociality, Microbial_Genus) %>%
    relocate(Microbial_Genus, .before = "Sociality") %>%
    as.data.frame() 
  if (i == 1){
    write.table(soc.inc3,
                "output/Prokaryote/Miscellaneous/Pro_MicrobialTaxa_by_SocialityPresent.tsv",
                  col.names = T, row.names = F, quote = F, sep = "\t")
  } else {
    write.table(soc.inc3, 
                "output/Prokaryote/Miscellaneous/Pro_MicrobialTaxa_by_SocialityPresent_TriRed.tsv",
                col.names = T, row.names = F, quote = F, sep = "\t")
  }
}

##Microbial Taxa by Species #####
spec.inc <- pro2 %>%
  rownames_to_column(var = "TaxID") %>%
  pivot_longer(-TaxID) %>%
  filter(value > 0)
names(spec.inc)[2:3] <- c("Sample", "ReadCnt")
spec.inc <- inner_join(spec.inc, met, by = c("Sample" = "Sample.ID")) %>%
  select(TaxID, Genus) %>%
  group_by(TaxID) %>%
  unique() %>%
  mutate(Host_Genera = paste0(Genus, collapse = ", "))
spec.inc <- inner_join(spec.inc, tax, by = c("TaxID" = "ID")) %>%
  mutate(Microbe_Genus = genus) %>%
  select(TaxID, Host_Genera, Microbe_Genus) %>%
  relocate(Microbe_Genus, .before = Host_Genera) %>%
  unique()

write.table(spec.inc,
            "output/Prokaryote/Miscellaneous/Pro_BacterialGen_by_HostGen.tsv",
            sep = "\t",
            row.names = F, col.names = T, quote = F)
