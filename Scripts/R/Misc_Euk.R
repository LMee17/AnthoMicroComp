#13th January 2023
#Looking at some exploratory plots (eukaryote)

#Libaries
library(tidyverse)
library(ggplot2)
library(reshape)

dir.create("output/Eukaryote/Miscellaneous/")

##Counts and Metadata ####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Feb23.tsv",
                  sep = "\t", header = T)
#microbial taxa metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv", 
                  header = T, sep = "\t")

#Eukaryote count table (family reduced)
euk <- read.table("input/Counts/Eukaryote_Filtered_FamilyReduced_Raw.tsv")
#Eukaryote relative abundance
euk.rel <- read.table("input/Counts/Eukaryote_RelativeAbundance.tsv")
euk.rel <- euk.rel[,names(euk.rel) %in% names(euk)]

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
#eukbably be easier working with incidence data for a second
euk.inc <- euk
euk.inc[euk.inc > 0] <- 1

specCnt <- data.frame(Sample = names(euk.inc),
                      SpeciesTot = colSums(euk.inc))
euk.Spec.Cnt <- euk.inc %>%
  t() %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") 
names(euk.Spec.Cnt)[2] <- "Total"
euk.Spec.Plot <- inner_join(euk.Spec.Cnt, met, by = c("Sample" = "Sample.ID"))

readCnts <- euk %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  group_by(Sample) %>%
  pivot_longer(-Sample) %>%
  mutate(Tot = sum(value)) %>%
  select(Sample, Tot) %>%
  unique() 
euk.Spec.Plot <- inner_join(euk.Spec.Plot, readCnts, by = "Sample")

#plot
ggplot(data = euk.Spec.Plot, aes(x = Platform_Spec, y = Total)) +
  geom_boxplot(aes(colour = Platform_Spec)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Sequencing Platform",
       y = "Number of Unique Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_Platform.pdf")

#what about by sociality / genus / family ? 
euk.Spec.Plot$Sociality <- factor(euk.Spec.Plot$Sociality,
                                  levels = c("O. Eusocial",
                                             "F. Eusocial",
                                             "Solitary"))
ggplot(data = euk.Spec.Plot, aes(x = Sociality, y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Sociality",
       y = "Number of Unique Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_Sociality.pdf")

#family
ggplot(data = euk.Spec.Plot, aes(x = Family, y = Total)) +
  geom_boxplot(aes(colour = Family)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Family",
       y = "Number of Unique Eukaryote Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_Family.pdf")

#for genus, order genera by sociality
euk.Spec.Plot$Genus <- factor(euk.Spec.Plot$Genus,
                              levels = c("Apis", "Bombus", "Tetragonisca", 
                                         "Tetragonula",
                                         "Ceratina", "Euglossa", "Exoneura", 
                                         "Halictus", "Lasioglossum", "Megalopta",
                                         "Anthophora", "Andrena", "Dufourea", 
                                         "Epeolus", "Habropoda",
                                         "Nomada", "Osmia"))
#genus
ggplot(data = euk.Spec.Plot, aes(x = fct_rev(Genus), y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Genus",
       y = "Number of Unique Eukaryote Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_Genus.pdf")

table(euk.Spec.Plot$Sociality, euk.Spec.Plot$Tribe)

#tribe
#again, order by sociality
euk.Spec.Plot$Tribe <- factor(euk.Spec.Plot$Tribe,
                              levels = c("Apini", "Bombini", "Meliponini",
                                         "Allodapini", "Augochlorini", "Ceratini",
                                         "Euglossini", "Halictini",
                                         "Andrenini", "Anthophorini", "Epeolini",
                                         "Nomadini", "Osmiini", "Rophitini"))

ggplot(data = euk.Spec.Plot, aes(x = fct_rev(Tribe), y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Tribe",
       y = "Number of Unique Eukaryote Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_Genus.pdf")

#library prep
ggplot(data = euk.Spec.Plot, aes(x = LibrarySelection, y = Total)) +
  geom_boxplot(aes(colour = LibrarySelection)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Method of Library Selection",
       y = "Number of Unique Eukaryote Species") +
  guides(colour = "none")
ggsave("output/Eukaryote/Miscellaneous/euk_SpeciesRichness_Vs_LibrarySelection.pdf")

#total number of reads versus number of species
ggplot(data = euk.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = Sociality, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = Sociality), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  scale_colour_manual(values = myPal)
ggsave("output/Eukaryote/Miscellaneous/euk_TotSpec_Vs_TotalReads_BySoc.pdf")

#total number of reads versus number of species (by library selection)
ggplot(data = euk.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = LibrarySelection, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = LibrarySelection), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  labs(colour = "Method of Library Selection")
ggsave("output/Eukaryote/Miscellaneous/euk_TotSpec_Vs_TotalReads_ByLibSel.pdf")

#total number of reads versus number of species (by family)
ggplot(data = euk.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = Family, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = Family), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  labs(colour = "Host Family")
ggsave("output/Eukaryote/Miscellaneous/euk_TotSpec_Vs_TotalReads_ByFam.pdf")

#interesting.

##eukaryotic Family / Order Facet Plot ####
#first, check taxa included here all have family entries
taxkey <- data.frame(TaxID = rownames(euk))
taxkey <- inner_join(taxkey, tax, by = c("TaxID" = "ID"))
taxkey %>%
  filter(family == "")

#Starmerella has not been assigned a family phylogenetically
#below is taken from ncbi taxoniy 13th Januar 2023
taxkey$family[taxkey$TaxID == "TRH1"] <- "Saccharomycetales incertae sedis"

#check order
taxkey %>%
  filter(order == "")

#nosema has not been fully classified (see notes for references)
#so will refer to them as Nosemitadae [Order Unclassified]

taxkey$order[taxkey$TaxID == "TRH462"] <- "Nosematidae \n(Order Unclassified)"

#begin to create plot table
euk.micro.plot <- euk.rel %>%
  rownames_to_column() %>%
  melt() 
names(euk.micro.plot) <- c("TaxID", "Sample", "RelAbundance")
euk.micro.plot <- inner_join(euk.micro.plot, taxkey, by = c("TaxID"))
#don't want clashes with the host/micro metadata. Add micro prefix before continuing
for(i in 4:length(names(euk.micro.plot))){
  print(names(euk.micro.plot[i]))
  names(euk.micro.plot)[i] <- sub("^", "Micro_", names(euk.micro.plot[i]))
}
euk.micro.plot <- inner_join(euk.micro.plot, met, by = c("Sample" = "Sample.ID"))

#differentiate between host and micro - adding host prefix
names(euk.micro.plot)[c(13,27:30)] <- sub("^", "Host_", names(euk.micro.plot[c(13,27:30)]))

#plot
#x = bacterial family / order
#y - abundance
#facet = sociality
#I would like to order the eukaryotes by what they are
#
fun.fam <- sort(unique(taxkey$family[taxkey$kingdom == "Fungi"]))
fun.ord <- sort(unique(taxkey$order[taxkey$kingdom == "Fungi"]))
other.fam <- sort(unique(taxkey$family[!taxkey$kingdom == "Fungi"]))
other.ord <- sort(unique(taxkey$order[!taxkey$kingdom == "Fungi"]))

##family
euk.micro.plot$Micro_family <- factor(euk.micro.plot$Micro_family, 
                                      levels = c(fun.fam, other.fam))

#sociality
euk.micro.plot$Sociality <- factor(euk.micro.plot$Sociality,
                                levels = c("O. Eusocial",
                                           "F. Eusocial",
                                           "Solitary"))

ggplot(data = euk.micro.plot, 
       aes(x = fct_rev(factor(Micro_family)), y = RelAbundance)) +
  geom_bar(aes(fill = Sociality), stat = "identity") +
  labs(y = "Cumulative Relative Abundance",
       x = "Eukaryote Family") +
  #ylim(0,5)+
  scale_fill_manual(values = myPal) +
  facet_grid(~Sociality) +
  coord_flip()
ggsave("output/Eukaryote/Miscellaneous/euk_BacterialFam_RelAbu_SocialFacet.pdf")

#by order
euk.micro.plot$Micro_order <- factor(euk.micro.plot$Micro_order, 
                                      levels = c(fun.ord, other.ord))

ggplot(data = euk.micro.plot, 
       aes(x = fct_rev(factor(Micro_order)), y = RelAbundance)) +
  geom_bar(aes(fill = Sociality), stat = "identity") +
  labs(y = "Cumulative Relative Abundance",
       x = "Eukaryote Family") +
  scale_fill_manual(values = myPal) +
  facet_grid(~Sociality) +
  coord_flip()
ggsave("output/Eukaryote/Miscellaneous/euk_BacterialOrd_RelAbu_SocialFacet.pdf")

#fuck it let's do genus
fun.gen <- sort(unique(taxkey$genus[taxkey$kingdom == "Fungi"]))
other.gen <- sort(unique(taxkey$genus[!taxkey$kingdom == "Fungi"]))

euk.micro.plot$Micro_genus <- factor(euk.micro.plot$Micro_genus,
                                     levels = c(fun.gen, other.gen))

ggplot(data = euk.micro.plot, 
       aes(x = fct_rev(factor(Micro_genus)), y = RelAbundance)) +
  geom_bar(aes(fill = Sociality), stat = "identity") +
  labs(y = "Cumulative Relative Abundance",
       x = "Eukaryote Family") +
  #ylim(0,5)+
  scale_fill_manual(values = myPal) +
  facet_grid(~Sociality) +
  coord_flip()
ggsave("output/Eukaryote/Miscellaneous/euk_BacterialGen_RelAbu_SocialFacet.pdf",
       height = 28, units = "cm")

##Eukaryote Heatmap####
heat.plot.inc <- euk.inc %>%
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
genera <- unique(met$Genus[met$Sample.ID %in% names(euk)])
micros <- unique(heat.plot.inc$MicroTax)
prevmat <- matrix(nrow = length(micros),
                  ncol = length(genera))

samplabs <- c()
for (i in 1:length(genera)){
  print(genera[i])
  samps <- unique(euk.Spec.Plot$Sample[euk.Spec.Plot$Genus == genera[i]])
  tot <- length(samps)
  label <- paste(genera[i], " (n = ", tot, ")", sep = "")
  samplabs <- c(samplabs, label)
  #prepare to populate with prevalence per microtaxa
  x <- c()
  for (j in 1:length(micros)){
    m.gen <- tax$genus[tax$ID == micros[j]]
    prev <- (sum(euk.inc[micros[j], names(euk.inc) %in% samps]) / tot)*100
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
hpi$Arrange[hpi$Sociality == "O. Eusocial"] <- 1
hpi$Arrange[hpi$Sociality == "F. Eusocial"] <- 2
hpi$Arrange[hpi$Sociality == "Solitary"] <- 3
levelz <- hpi %>%
  select(Host_Genus, Labels, Arrange) %>%
  arrange(Labels) %>%
  arrange(Arrange) %>%
  unique() %>%
  select(Labels) %>%
  unlist()
#have euglossini near corbiculates
levelz <- levelz[c(1:4,6,5,7:17)]

hpi$Labels <- factor(hpi$Labels, levels = c(levelz))

#arrange micobial taxa by order, ie all fungus, all tryps, etc
eu.tax <- tax[tax$ID %in% rownames(euk),]

fungi <- eu.tax %>%
  subset(kingdom == "Fungi") %>%
  select(genus) %>%
  unlist() %>%
  as.vector() %>%
  sort()
tryps <- eu.tax %>%
  subset(family == "Trypanosomatidae") %>%
  select(genus) %>%
  unlist() %>%
  as.vector() %>%
  sort()
other <-  eu.tax %>%
  subset(family != "Trypanosomatidae") %>%
  subset(kingdom != "Fungi") %>%
  select(genus) %>%
  unlist() %>%
  as.vector() %>%
  sort()

hpi$Taxa <- factor(hpi$Taxa, levels = c(fungi, tryps, other))

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
  labs(x = "Eukaryote Genus",
       y = "Host Genus",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Eukaryote/Miscellaneous/euk_PrevalenceHeatmap_Factored_HostGenus.pdf",
       height = 28, units = "cm")

ggplot(data = hpi, aes(x = fct_rev(Taxa), y = Labels, fill = Prev2)) + 
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 80, hjust = -.1,
                                   face = "italic"),
        axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("#130000",
                               "#460000",
                               "#790000",
                               "#e00000",
                               "#ff6161",
                               "white")) +
  labs(x = " ",
       y = " ",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Eukaryote/Miscellaneous/euk_PrevalenceHeatmap_Factored_HostGenus_Red.pdf",
       height = 28, units = "cm")

##Sociality vs SharedMicrobial Species #####
#this will be difficult to code.....
soc.inc <- euk %>%
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
                         levels = c("FOS", 
                                    "FS", "OS", "FO",
                                    "O", "F", "S"))
soc.inc2$Sociality <- factor(soc.inc2$Sociality,
                             levels = c("O. Eusocial",
                                        "F. Eusocial", 
                                        "Solitary"))
ggplot(data = soc.inc2,
       aes(x = Sociality, y = Count, fill = (Combo), label = Count))+
  geom_bar(stat = "identity", colour = "grey") + 
  geom_text(size = 3, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(labels = c("All Socialities", "O. Eusocial and Solitary",
                               "O. Eusocial and F. Eusocial", "O. Eusocial only",
                               "F. Eusocial Only"),
                    values = c("#6893b4",
                               "#42282e",
                               "#a2ce47",
                               "#db6d00",
                               "#009292")) +
  theme(panel.grid =  element_blank()) +
  labs(y = "Number of Eukaryote Genera",
       x = "",
       fill = "")
ggsave("output/Eukaryote/Miscellaneous/euk_MicroGenera_bySocCombo.pdf")  
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
            "output/Eukaryote/Miscellaneous/euk_MicrobialTaxa_by_SocialityPresent.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

##Microbial Taxa by Species #####
spec.inc <- euk %>%
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
            "output/Eukaryote/Miscellaneous/euk_BacterialGen_by_HostGen.tsv",
            sep = "\t",
            row.names = F, col.names = T, quote = F)

##SessionLog####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Misc_Euk_Jan23.txt")
 