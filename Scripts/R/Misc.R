#9th January 2023
#Looking at some exploratory plots (prokaryote)

#Libaries
library(tidyverse)
library(ggplot2)
library(reshape)

dir.create("output/Prokaryote/Miscellaneous/")

##Counts and Metadata ####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  sep = "\t", header = T)
#microbial taxa metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv", 
                  header = T, sep = "\t")

#prokaryote count table (tribe reduced)
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")
#prokaryote relative abundance
pro.rel <- read.table("input/Counts/Prokaryote_RelativeAbundance.tsv")
pro.rel <- pro.rel[,names(pro.rel) %in% names(pro)]

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

#plot
ggplot(data = pro.Spec.Plot, aes(x = Platform_Spec, y = Total)) +
    geom_boxplot(aes(colour = Platform_Spec)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Sequencing Platform",
         y = "Number of Unique Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Platform.pdf")

#what about by sociality / genus / family ? 
ggplot(data = pro.Spec.Plot, aes(x = Sociality, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Sociality",
         y = "Number of Unique Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Sociality.pdf")
    
#family
ggplot(data = pro.Spec.Plot, aes(x = Family, y = Total)) +
    geom_boxplot(aes(colour = Family)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Family",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Family.pdf")


#for genus, order genera by sociality
pro.Spec.Plot$Genus <- factor(pro.Spec.Plot$Genus,
                               levels = c("Osmia", "Andrena",
                                          "Megalopta",
                                          "Euglossa", "Eufriesea",
                                          "Ceratina",
                                          "Tetragonula", "Tetragonisca", "Bombus", 
                                          "Apis"))

#genus
ggplot(data = pro.Spec.Plot, aes(x = Genus, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Genus",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus.pdf")

#tribe
#again, order by sociality
pro.Spec.Plot$Tribe <- factor(pro.Spec.Plot$Tribe,
                               levels = c("Osmiini", 
                                          "Andrenini",
                                          "Euglossini", "Ceratini", 
                                          "Augochlorini", "Allodapini",
                                          "Meliponini", "Bombini", "Apini"))

ggplot(data = pro.Spec.Plot, aes(x = Tribe, y = Total)) +
    geom_boxplot(aes(colour = Sociality)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Host Tribe",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_Genus.pdf")

#library prep
ggplot(data = pro.Spec.Plot, aes(x = LibrarySelection, y = Total)) +
    geom_boxplot(aes(colour = LibrarySelection)) +
    coord_flip() +
    scale_colour_manual(values = myPal) +
    labs(x = "Method of Library Selection",
         y = "Number of Unique Prokaryote Species") +
    guides(colour = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_SpeciesRichness_Vs_LibrarySelection.pdf")

#total number of reads versus number of species
ggplot(data = pro.Spec.Plot, aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = Sociality, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = Sociality), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none")
ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_BySoc.pdf")

#total number of reads versus number of species (by library selection)
ggplot(data = pro.Spec.Plot, aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = LibrarySelection, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = LibrarySelection), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none") +
    labs(colour = "Method of Library Selection")
ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByLibSel.pdf")

#total number of reads versus number of species (by family)
ggplot(data = pro.Spec.Plot, aes(x = log(Tot), y = Total)) +
    geom_point(aes(colour = Family, alpha = 0.85)) +
    labs(x = "log(Total Number of Reads)",
         y = "Number of Unique Species") +
    geom_smooth(aes(colour = Family), 
                method = lm, se = FALSE, fullrange = TRUE,
                linetype = "twodash") +
    guides(alpha = "none") +
    labs(colour = "Host Family")
ggsave("output/Prokaryote/Miscellaneous/Pro_TotSpec_Vs_TotalReads_ByFam.pdf")

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

#plot
#x = bacterial family / order
#y - abundance
#facet = sociality
ggplot(data = pro.micro.plot, 
         aes(x = fct_rev(factor(Micro_family)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    #ylim(0,5)+
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialFam_RelAbu_SocialFacet.pdf")

#by order
ggplot(data = pro.micro.plot, 
         aes(x = fct_rev(factor(Micro_order)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialOrd_RelAbu_SocialFacet.pdf")

#fuck it let's do genus
ggplot(data = pro.micro.plot, 
         aes(x = fct_rev(factor(Micro_genus)), y = RelAbundance)) +
    geom_bar(aes(fill = Sociality), stat = "identity") +
    labs(y = "Cumulative Relative Abundance",
         x = "Prokaryote Family") +
    #ylim(0,5)+
    scale_fill_manual(values = myPal) +
    facet_grid(~Sociality) +
    coord_flip()
ggsave("output/Prokaryote/Miscellaneous/Pro_BacterialGen_RelAbu_SocialFacet.pdf",
       height = 28, units = "cm")

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
  samps <- unique(pro.Spec.Plot$Sample[pro.Spec.Plot$Genus == genera[i]])
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
#repeat the labels x number of microbial taxa
notax <- length(unique(hpi$TaxID))
hpi <- hpi %>%
  mutate(Labels = rep(samplabs, notax)) %>%
  as.data.frame()

#add line breaks into labels
hpi$Labels <- str_replace(hpi$Labels, " ", "\n")

#order samples by sociality
hpi$Labels <- as.character(hpi$Labels)
levelz <- hpi %>%
  select(Host_Genus, Sociality, Labels) %>%
  arrange(Labels) %>%
  arrange(Sociality) %>%
  unique() %>%
  select(Labels) %>%
  unlist()
#have euglossini near corbiculates
levelz <- levelz[c(1:4,6:7,5,8:10)]

hpi$Labels <- factor(hpi$Labels, levels = c(levelz))

#order microbial taxa so core are together ? 
core <- factor(c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
                 "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
                 "Bartonella", "Apibacter", 
                 "Bombella", "Parasaccharibacter", "Commensalibacter"))
#everything else by bacterial order
pro.tax <- tax[tax$ID %in% rownames(pro),]

bigones <- c("Actinobacteria",  "Bacteroidota", "Firmicutes", "Proteobacteria")
baclist <- vector(mode = "list", length = length(bigones))
for (i in 1:length(bigones)){
  baclist[[i]] <- pro.tax %>%
    subset(phylum == bigones[i]) %>%
    select(genus) %>%
    unlist() %>%
    as.vector() %>%
    sort()
}
other <- pro.tax %>%
  subset(!phylum %in% bigones) %>%
  select(genus) %>%
  unlist() %>%
  as.vector() %>%
  sort()

hpi$Taxa <- factor(hpi$Taxa, levels = c(baclist[[1]], baclist[[2]],
                                        baclist[[3]], bacllist[[4]],
                                        noncore))

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

ggplot(data = hpi, aes(x = fct_rev(Taxa), y = Labels, fill = Prev2)) + 
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 80, hjust = -.1,
                                   face = "italic"),
        axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("#440088",
                               "#6a00d4",
                               "#8307ff",
                               "#b66eff",
                               "#d0a1ff",
                               "white")) +
  labs(x = "",
       y = "",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Prokaryote/Miscellaneous/Pro_PrevalenceHeatmap_Factored_HostGenus_Purp.pdf",
       height = 28, units = "cm")

##Sociality vs SharedMicrobial Species #####
#this will be difficult to code.....
soc.inc <- pro %>%
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
ggsave("output/Prokaryote/Miscellaneous/Pro_MicroGenera_bySocCombo.pdf")  
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
write.table(soc.inc3,
      "output/Prokaryote/Miscellaneous/Pro_MicrobialTaxa_by_SocialityPresent.tsv",
                  col.names = T, row.names = F, quote = F, sep = "\t")

##Microbial Taxa by Species #####
spec.inc <- pro %>%
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

##SessionLog####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Misc_Pro_Jan23.txt")
