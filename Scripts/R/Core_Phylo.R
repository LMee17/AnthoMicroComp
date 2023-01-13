#14th November 2022
#Updated 12th Jan 2023
#Assessing the Core Phylotypes in Data

library(ggplot2)
library(reshape)
library(tidyverse)
#library(RColorBrewer)

dir.create("output/")
dir.create("output/Prokaryote/CorePhylo/")

#####Load Necessary Counts and Metadata####
#Sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_RNAOnly_Dec22.tsv",
                  sep = "\t", header = T)

#Microbial metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t")

#prokaryote count table
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")

#####Quick Incidence Overview: Core Phylo + All Samples####
#extract the core microbiota phylotype
beecore <- c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
             "Bartonella", "Bombiscardovia", "Schmidhempelia", "Apibacter")
coreID <- tax %>%
  subset(genus %in% beecore) %>%
  select(ID) %>%
  unlist %>% as.vector()

#convert count table to 3 columns: Microbial Taxa, Sample, NoReads
#subset so that only the taxa from the core phylotypes are kept
core.df <- pro %>%
  rownames_to_column(var = "TaxID") %>%
  pivot_longer(-TaxID) %>%
  subset(TaxID %in% coreID)
#rename variables
names(core.df)[2:3] <- c("Sample", "NoReads")
#add sample metada (just tribe and sociality)
core.df <- inner_join(core.df, met, by = c("Sample" = "Sample.ID")) %>%
  select(TaxID, Sample, NoReads, Tribe, Sociality) %>%
  mutate(Category = Tribe)
#add categories (tribes for eusocial, sociality for other socialities)
core.df$Category[core.df$Sociality == "Polymorphic"] <- "Polymorphic Species"
core.df$Category[core.df$Sociality == "Solitary"] <- "Solitary Species"
#get n samples per category
catno <- c()
cats <- unique(core.df$Category)
for (i in 1:length(unique(core.df$Category))){
  print(cats[i])
  n <- as.numeric(length(unique(core.df$Sample[core.df$Category == cats[i]])))
  catno <- c(catno, n)
}
catkey <- as.data.frame(cbind(catno, cats))
catkey$catno <- as.numeric(catkey$catno)
for (i in 1:nrow(core.df)){
  core.df$CatN[i] <- as.numeric(catkey$catno[catkey$cats == core.df$Category[i]])
}
#for some reason, despite many directions to the contrary, CatN is seen as a
#character vector
core.df$CatN <- as.numeric(core.df$CatN)
core.df <- core.df %>%
  mutate(Label = paste(Category, " (n = ", CatN, ")", sep ="")) %>%
  group_by(Sample) %>%
  mutate(TotReads = sum(NoReads), .after = NoReads) %>%
  group_by(TaxID) %>%
  mutate(RelAbundance = NoReads / TotReads, .after = TotReads) %>%
  #remove NA values
  mutate_all(~replace(., is.nan(.), 0)) %>%
  #remove unnecessary columns
  select(-Tribe, -Sociality) %>%
  #get average abundance by Category (when it's present, and across all samples)
  group_by(TaxID, Category) %>%
  mutate(AvgRelAbundancePres = sum(RelAbundance) / length(unique(Sample[RelAbundance!=0])),
         .after = RelAbundance) %>%
  mutate(AvgRelAbundanceAll = sum(RelAbundance) / CatN,
         .after = AvgRelAbundancePres) %>%
  #again replace NA
  mutate_all(~replace(., is.nan(.), 0)) %>%
  #now to compute prevalence
  ungroup() %>% mutate(inc = ifelse(NoReads > 0, 1, 0)) %>%
  group_by(TaxID, Category) %>%
  mutate(TotInc = sum(inc)) %>%
  ungroup() %>%
  mutate(Prevalence = TotInc / CatN, .after = AvgRelAbundanceAll) %>%
  #remove unnecessary columns
  select(TaxID, AvgRelAbundanceAll, AvgRelAbundancePres, Prevalence, Category, Label)
#add microbiota metadata
core.df <- inner_join(core.df, tax, by = c("TaxID" = "ID")) %>%
  select(-species, -tax_name, -kingdom, -phylum, -superkingdom) %>%
  unique() %>%
  as.data.frame()
#add new lines in labels
core.df$Label <- str_replace(core.df$Label, " ", "\n")

#arrange microbials
core.df$genus <- factor(core.df$genus, levels = c("Gilliamella","Snodgrassella", 
                                                  "Lactobacillus: Firm-5","Bifidobacterium", 
                                                  "Frischella", "Bartonella", "Apibacter", 
                                                  "Apilactobacillus","Bombilactobacillus", 
                                                  "Bombiscardovia"))
#factorise prevalence
core.df$Prev2[core.df$Prevalence > 0.8] <- "81 - 100%"
core.df$Prev2[core.df$Prevalence <= 0.8 &
                core.df$Prevalence > 0.6] <- "61 - 80%"
core.df$Prev2[core.df$Prevalence <= 0.6 &
                core.df$Prevalence > 0.4] <- "41 - 60%"
core.df$Prev2[core.df$Prevalence <= 0.4 &
                core.df$Prevalence > 0.2] <- "21 - 40%"
core.df$Prev2[core.df$Prevalence <= 0.2 &
                core.df$Prevalence >= 0.05 ] <- "5 - 20%"
core.df$Prev2[core.df$Prevalence < 0.05] <- "< 5%"
core.df$Prev2[core.df$Prevalence == 0] <- " "

#order
core.df$Prev2 <- factor(core.df$Prev2,
                        levels = c("81 - 100%", "61 - 80%", "41 - 60%", 
                                   "21 - 40%", "5 - 20%", "< 5%", " "))
#prepare average abundances to printed over tiles (need to end up character vectors)
for (i in 2:3){}
  core.df[,i] <- core.df[,i]*100
  core.df[,i] <- round(core.df[,i], digits = 1)
  core.df[,i] <- as.character(core.df[,i])
  core.df[,i] <- paste0(core.df[,i], "%")
}

#plot (all )
ggplot(data = core.df, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  geom_text(aes(label = AvgRelAbundanceAll), color = "white", size = 4) +
  scale_fill_manual(values = c("#004d4d",
                                 "#006767", 
                                 "#008080",
                                 "#009a9a",
                                 "#00b3b3",
                                 "#00cdcd",
                                 "white")) +
  labs(fill = "Prevalence",
       x = "",
       y = "Core Phylotypes") + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                                 face = "italic"))
ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_PrevAndAvgAbu.pdf")


#plot (without average relative abundance)
ggplot(data = core.df, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  scale_fill_manual(values = c("#004d4d",
                               "#006767", 
                               "#008080",
                               "#009a9a",
                               "#00b3b3",
                               "#00cdcd",
                               "white")) +
  labs(fill = "Prevalence",
       x = "",
       y = "Core Phylotypes") + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   face = "italic"))
ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_Prev.pdf")

