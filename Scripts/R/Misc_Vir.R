#13th January 2023
#Looking at some exploratory plots (Viral)

#Libaries
library(tidyverse)
library(ggplot2)
library(reshape)

dir.create("output/Viral/Miscellaneous/")

##Counts and Metadata ####
#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Feb23.tsv",
                  sep = "\t", header = T)
#microbial taxa metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv", 
                  header = T, sep = "\t")

#Viral count table (family reduced)
vir <- read.table("input/Counts/Viral_Filtered_FamilyReduced_Raw.tsv")
#Viral relative abundance
vir.rel <- read.table("input/Counts/Viral_RelativeAbundance.tsv")
vir.rel <- vir.rel[,names(vir.rel) %in% names(vir)]

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
#virbably be easier working with incidence data for a second
vir.inc <- vir
vir.inc[vir.inc > 0] <- 1

specCnt <- data.frame(Sample = names(vir.inc),
                      SpeciesTot = colSums(vir.inc))
vir.Spec.Cnt <- vir.inc %>%
  t() %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") 
names(vir.Spec.Cnt)[2] <- "Total"
vir.Spec.Plot <- inner_join(vir.Spec.Cnt, met, by = c("Sample" = "Sample.ID"))

readCnts <- vir %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  group_by(Sample) %>%
  pivot_longer(-Sample) %>%
  mutate(Tot = sum(value)) %>%
  select(Sample, Tot) %>%
  unique() 
vir.Spec.Plot <- inner_join(vir.Spec.Plot, readCnts, by = "Sample")

#plot
ggplot(data = vir.Spec.Plot, aes(x = Platform_Spec, y = Total)) +
  geom_boxplot(aes(colour = Platform_Spec)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Sequencing Platform",
       y = "Number of Unique Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_Platform.pdf")

#what about by sociality / genus / family ? 
vir.Spec.Plot$Sociality <- factor(vir.Spec.Plot$Sociality,
                                  levels = c("O. Eusocial",
                                             "F. Eusocial",
                                             "Solitary"))
ggplot(data = vir.Spec.Plot, aes(x = Sociality, y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Sociality",
       y = "Number of Unique Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_Sociality.pdf")

#family
ggplot(data = vir.Spec.Plot, aes(x = Family, y = Total)) +
  geom_boxplot(aes(colour = Family)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Family",
       y = "Number of Unique Viral Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_Family.pdf")

#for genus, order genera by sociality
table(vir.Spec.Plot$Genus, vir.Spec.Plot$Sociality)

vir.Spec.Plot$Genus <- factor(vir.Spec.Plot$Genus,
                              levels = c("Apis", "Bombus", "Tetragonisca", 
                                         "Lasioglossum",
                                         "Andrena", "Osmia"))
#genus
ggplot(data = vir.Spec.Plot, aes(x = fct_rev(Genus), y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Genus",
       y = "Number of Unique Viral Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_Genus.pdf")

#tribe
#again, order by sociality
table(vir.Spec.Plot$Sociality, vir.Spec.Plot$Tribe)

vir.Spec.Plot$Tribe <- factor(vir.Spec.Plot$Tribe,
                              levels = c("Apini", "Bombini", "Meliponini",
                                         "Halictini",
                                         "Andrenini", "Osmiini"))

ggplot(data = vir.Spec.Plot, aes(x = fct_rev(Tribe), y = Total)) +
  geom_boxplot(aes(colour = Sociality)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Host Tribe",
       y = "Number of Unique Viral Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_Genus.pdf")

#library prep
ggplot(data = vir.Spec.Plot, aes(x = LibrarySelection, y = Total)) +
  geom_boxplot(aes(colour = LibrarySelection)) +
  coord_flip() +
  scale_colour_manual(values = myPal) +
  labs(x = "Method of Library Selection",
       y = "Number of Unique Viral Species") +
  guides(colour = "none")
ggsave("output/Viral/Miscellaneous/vir_SpeciesRichness_Vs_LibrarySelection.pdf")

#total number of reads versus number of species
ggplot(data = vir.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = Sociality, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = Sociality), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  scale_colour_manual(values = myPal)
ggsave("output/Viral/Miscellaneous/vir_TotSpec_Vs_TotalReads_BySoc.pdf")

#total number of reads versus number of species (by library selection)
ggplot(data = vir.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = LibrarySelection, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = LibrarySelection), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  labs(colour = "Method of Library Selection")
ggsave("output/Viral/Miscellaneous/vir_TotSpec_Vs_TotalReads_ByLibSel.pdf")

#total number of reads versus number of species (by family)
ggplot(data = vir.Spec.Plot, aes(x = log(Tot), y = Total)) +
  geom_point(aes(colour = Family, alpha = 0.85)) +
  labs(x = "log(Total Number of Reads)",
       y = "Number of Unique Species") +
  geom_smooth(aes(colour = Family), 
              method = lm, se = FALSE, fullrange = TRUE,
              linetype = "twodash") +
  guides(alpha = "none") +
  labs(colour = "Host Family")
ggsave("output/Viral/Miscellaneous/vir_TotSpec_Vs_TotalReads_ByFam.pdf")
#so few datapoints, a lot of these are irrelevant

##viraryotic Family / Order Facet Plot ####
#first, check taxa included here all have family entries
taxkey <- data.frame(TaxID = rownames(vir))
taxkey <- inner_join(taxkey, tax, by = c("TaxID" = "ID"))
taxkey %>%
  filter(family == "")

#check order
taxkey %>%
  filter(order == "")

#begin to create plot table
vir.micro.plot <- vir.rel %>%
  rownames_to_column() %>%
  melt() 
names(vir.micro.plot) <- c("TaxID", "Sample", "RelAbundance")
vir.micro.plot <- inner_join(vir.micro.plot, taxkey, by = c("TaxID"))
#don't want clashes with the host/micro metadata. Add micro prefix before continuing
for(i in 4:length(names(vir.micro.plot))){
  print(names(vir.micro.plot[i]))
  names(vir.micro.plot)[i] <- sub("^", "Micro_", names(vir.micro.plot[i]))
}
vir.micro.plot <- inner_join(vir.micro.plot, met, by = c("Sample" = "Sample.ID"))

#differentiate between host and micro - adding host prefix
names(vir.micro.plot)[c(13,27:30)] <- sub("^", "Host_", names(vir.micro.plot[c(13,27:30)]))

#plot
#x = bacterial family / order
#y - abundance
#facet = sociality
vir.micro.plot$Sociality <- factor(vir.micro.plot$Sociality,
                                  levels = c("O. Eusocial",
                                             "F. Eusocial",
                                             "Solitary"))
ggplot(data = vir.micro.plot, 
       aes(x = fct_rev(factor(Micro_family)), y = RelAbundance)) +
  geom_bar(aes(fill = Sociality), stat = "identity") +
  labs(y = "Cumulative Relative Abundance",
       x = "Viral Family") +
  #ylim(0,5)+
  scale_fill_manual(values = myPal) +
  facet_grid(~Sociality) +
  coord_flip()
ggsave("output/Viral/Miscellaneous/vir_BacterialFam_RelAbu_SocialFacet.pdf")

#by order
ggplot(data = vir.micro.plot, 
       aes(x = fct_rev(factor(Micro_order)), y = RelAbundance)) +
  geom_bar(aes(fill = Sociality), stat = "identity") +
  labs(y = "Cumulative Relative Abundance",
       x = "Viral Family") +
  scale_fill_manual(values = myPal) +
  facet_grid(~Sociality) +
  coord_flip()
ggsave("output/Viral/Miscellaneous/vir_BacterialOrd_RelAbu_SocialFacet.pdf")

#can't do genus for viruses

##Viral Heatmap####
#start with prevalence
#make a table
genera <- unique(met$Genus[met$Sample.ID %in% names(vir)])
micros <- unique(heat.plot.inc$MicroTax)
prevmat <- matrix(nrow = length(micros),
                  ncol = length(genera))

samplabs <- c()
for (i in 1:length(genera)){
  print(genera[i])
  samps <- unique(vir.Spec.Plot$Sample[vir.Spec.Plot$Genus == genera[i]])
  tot <- length(samps)
  label <- paste(genera[i], " (n = ", tot, ")", sep = "")
  samplabs <- c(samplabs, label)
  #prepare to populate with prevalence per microtaxa
  x <- c()
  for (j in 1:length(micros)){
    m.gen <- tax$genus[tax$ID == micros[j]]
    prev <- (sum(vir.inc[micros[j], names(vir.inc) %in% samps]) / tot)*100
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
  mutate(Taxa = family) %>%
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

hpi$Labels <- factor(hpi$Labels, levels = c(levelz))

#order viruses by order
vir.tax <- tax[tax$ID %in% rownames(vir),]

pisu <- vir.tax %>%
  subset(phylum == "Pisuviricota") %>%
  select(family) %>%
  unlist() %>%
  as.vector() %>%
  sort()

other <- vir.tax %>%
  subset(phylum != "Pisuviricota") %>%
  select(family) %>%
  unlist() %>%
  as.vector() %>%
  sort()

hpi$Taxa <- factor(hpi$Taxa, levels = c(pisu, other))

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
  labs(x = "Viral Genus",
       y = "Host Family",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Viral/Miscellaneous/vir_PrevalenceHeatmap_Factored_HostGenus.pdf",
       height = 28, units = "cm")

ggplot(data = hpi, aes(x = fct_rev(Taxa), y = Labels, fill = Prev2)) + 
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 80, hjust = -.1,
                                   face = "italic"),
        axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("#004488",
                               "#006ad4",
                               "#0884ff",
                               "#6fb6ff",
                               "#bcddff",
                               "white")) +
  labs(x = " ",
       y = " ",
       fill = " ") +
  theme(panel.background = element_blank()) +
  coord_flip() 

ggsave("output/Viral/Miscellaneous/vir_PrevalenceHeatmap_Factored_HostGenus_Red.pdf",
       units = "cm")

##Sociality vs SharedMicrobial Species #####
#this will be difficult to code.....
soc.inc <- vir %>%
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
  scale_fill_manual(labels = c("Eusocial and Solitary",
                               "Eusocial and Polymorphic", "Eusocial only"),
                    values = c("#42282e",
                               "#a2ce47",
                               "#db6d00",
                               "#009292")) +
  theme(panel.grid =  element_blank()) +
  labs(y = "Number of Viral Families",
       x = "",
       fill = "")
ggsave("output/Viral/Miscellaneous/virk_MicroGenera_bySocCombo.pdf")  
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
            "output/Viral/Miscellaneous/vir_MicrobialTaxa_by_SocialityPresent.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

##Microbial Taxa by Species #####
spec.inc <- vir %>%
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
            "output/Viral/Miscellaneous/vir_BacterialGen_by_HostGen.tsv",
            sep = "\t",
            row.names = F, col.names = T, quote = F)

##SessionLog####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/Misc_vir_Jan23.txt")
