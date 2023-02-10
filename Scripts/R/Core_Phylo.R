#14th November 2022
#Updated 12th Jan 2023
#Updated 10th Feb 2023: changed names of social categories
#Assessing the Core Phylotypes in Data

library(ggplot2)
library(reshape)
library(tidyverse)
#library(RColorBrewer)

dir.create("output/")
dir.create("output/Prokaryote/CorePhylo/")

#####Load Necessary Counts and Metadata####
#Sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Feb23.tsv",
                  sep = "\t", header = T)

#Microbial metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t")

#prokaryote count table
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")

pro.rel <- read.table("input/Counts/Prokaryote_RelativeAbundance.tsv")
pro.rel <- pro.rel[, names(pro.rel) %in% names(pro)]

#####Quick Incidence Overview: Core Phylo + All Samples####
#extract the core microbiota phylotype
beecore <- c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
             "Bartonella", "Bombiscardovia", "Schmidhempelia", "Apibacter", 
             "Bombella", "Parasaccharibacter", "Commensalibacter")
coreID <- tax %>%
  subset(tax$genus %in% beecore) %>%
  select(ID) %>%
  unlist %>% as.vector()

#convert count table to 3 columns: Microbial Taxa, Sample, NoReads
core.df <- pro %>%
  rownames_to_column(var = "TaxID") %>%
  pivot_longer(-TaxID)
#rename variables
names(core.df)[2:3] <- c("Sample", "NoReads")
#add sample metada (just tribe and sociality)
core.df <- inner_join(core.df, met, by = c("Sample" = "Sample.ID")) %>%
  select(TaxID, Sample, NoReads, Tribe, Sociality) %>%
  mutate(Category = Tribe)
#add categories (tribes for eusocial, sociality for other socialities)
core.df$Category[core.df$Sociality == "F. Eusocial"] <- "F.Eusocial Species"
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
  subset(genus %in% beecore) %>%
  unique() %>%
  as.data.frame()
#add new lines in labels
core.df$Label <- str_replace(core.df$Label, " ", "\n")

#arrange microbials
core.df$genus <- factor(core.df$genus, levels = c("Bifidobacterium",
                                                  "Gilliamella","Snodgrassella", 
                                                  "Lactobacillus: Firm-5", 
                                                  "Bombilactobacillus", 
                                                  "Apilactobacillus",
                                                  "Apibacter", 
                                                  "Frischella", "Bartonella",
                                                  "Commensalibacter", 
                                                  "Bombella", "Parasaccharibacter",
                                                  "Bombiscardovia"))
#arrange tribe/social categories
core.df$Label <- factor(core.df$Label, 
                           levels = c("Apini\n(n = 86)", "Bombini\n(n = 57)", 
                                      "Meliponini\n(n = 7)",
                                      "F.Eusocial\nSpecies (n = 53)",
                                      "Solitary\nSpecies (n = 24)"))

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
for (i in 2:3){
  core.df[,i] <- core.df[,i]*100
  core.df[,i] <- round(core.df[,i], digits = 2)
  core.df[,i] <- as.character(core.df[,i])
  core.df[,i] <- paste0(core.df[,i], "%")
}
#check if there's any that need fixing
fixrow <- core.df %>%
  subset(AvgRelAbundanceAll == "0%" & Prevalence > 0) %>%
  rownames()
core.df[c(fixrow), 2] <- "< 0.01%"

#there's a straight up "4%" in bifido/meliponini I wanna fix too
core.df %>%
  subset(genus == "Bifidobacterium" & Category == "Meliponini")
core.df[24,2] <- "4.00%"

#also Frischella and Solitary
core.df %>%
  subset(genus == "Frischella" & Category == 'Solitary Species')
core.df[30,2] <- "0.90%"

#plot (all )
ggplot(data = core.df, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  theme_classic() +
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
                                                 face = "italic"),
        panel.background = element_blank())

ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_PrevAndAvgAbu.pdf")


#plot (without average relative abundance)
ggplot(data = core.df, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  theme_classic() +
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
                                   face = "italic"),
        panel.background = element_blank()) 
ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_Prev.pdf")

##The same as above ... just with Euglossini in it####
#convert count table to 3 columns: Microbial Taxa, Sample, NoReads
#subset so that only the taxa from the core phylotypes are kept
core.df2 <- pro %>%
  rownames_to_column(var = "TaxID") %>%
  pivot_longer(-TaxID)
#rename variables
names(core.df2)[2:3] <- c("Sample", "NoReads")
#add sample metada (just tribe and sociality)
core.df2 <- inner_join(core.df2, met, by = c("Sample" = "Sample.ID")) %>%
  select(TaxID, Sample, NoReads, Tribe, Sociality) %>%
  mutate(Category = Tribe)
#add categories (tribes for eusocial, sociality for other socialities)
core.df2$Category[core.df2$Sociality == "F. Eusocial"] <- "F.Eusocial Non-Corbiculates"
core.df2$Category[core.df2$Sociality == "Solitary"] <- "Solitary Non-Corbiculates"
core.df2$Category[core.df2$Tribe == "Euglossini"] <- "Euglossini"
#get n samples per category
catno <- c()
cats <- unique(core.df2$Category)
for (i in 1:length(unique(core.df2$Category))){
  print(cats[i])
  n <- as.numeric(length(unique(core.df2$Sample[core.df2$Category == cats[i]])))
  catno <- c(catno, n)
}
catkey <- as.data.frame(cbind(catno, cats))
catkey$catno <- as.numeric(catkey$catno)
for (i in 1:nrow(core.df2)){
  core.df2$CatN[i] <- as.numeric(catkey$catno[catkey$cats == core.df2$Category[i]])
}
#for some reason, despite many directions to the contrary, CatN is seen as a
#character vector
core.df2$CatN <- as.numeric(core.df2$CatN)
core.df2 <- core.df2 %>%
  mutate(Label = paste(Category, " (n=", CatN, ")", sep ="")) %>%
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
core.df2 <- inner_join(core.df2, tax, by = c("TaxID" = "ID")) %>%
  filter(TaxID %in% coreID) %>%
  select(-species, -tax_name, -kingdom, -phylum, -superkingdom) %>%
  unique() %>%
  as.data.frame()

#add new lines in labels
core.df2$Label <- str_replace_all(core.df2$Label, " ", "\n")

#arrange microbials
core.df2$genus <- factor(core.df2$genus, levels = c("Bifidobacterium",
                                                  "Gilliamella","Snodgrassella", 
                                                  "Lactobacillus: Firm-5", 
                                                  "Bombilactobacillus", 
                                                  "Apilactobacillus",
                                                  "Apibacter", 
                                                  "Frischella", "Bartonella",
                                                  "Commensalibacter", 
                                                  "Bombella", "Parasaccharibacter",
                                                  "Bombiscardovia"))#factorise prevalence
core.df2$Prev2[core.df2$Prevalence > 0.8] <- "81 - 100%"
core.df2$Prev2[core.df2$Prevalence <= 0.8 &
                 core.df2$Prevalence > 0.6] <- "61 - 80%"
core.df2$Prev2[core.df2$Prevalence <= 0.6 &
                 core.df2$Prevalence > 0.4] <- "41 - 60%"
core.df2$Prev2[core.df2$Prevalence <= 0.4 &
                 core.df2$Prevalence > 0.2] <- "21 - 40%"
core.df2$Prev2[core.df2$Prevalence <= 0.2 &
                 core.df2$Prevalence >= 0.05 ] <- "5 - 20%"
core.df2$Prev2[core.df2$Prevalence < 0.05] <- "< 5%"
core.df2$Prev2[core.df2$Prevalence == 0] <- " "

#order prevalence factors
core.df2$Prev2 <- factor(core.df2$Prev2,
                         levels = c("81 - 100%", "61 - 80%", "41 - 60%", 
                                    "21 - 40%", "5 - 20%", "< 5%", " "))

#order categories
labs <- unique(core.df2$Label)
labs <- labs[c(1:2,4,5,3,6)]

core.df2$Label <- factor(core.df2$Label,
                         levels = labs)

#prepare average abundances to printed over tiles (need to end up character vectors)
for (i in 2:3){
  core.df2[,i] <- core.df2[,i]*100
  core.df2[,i] <- round(core.df2[,i], digits = 2)
  core.df2[,i] <- as.character(core.df2[,i])
  core.df2[,i] <- paste0(core.df2[,i], "%")
}
#add point with too many decimals
fixrow <- core.df2 %>%
  subset(AvgRelAbundanceAll == "0%" & Prevalence > 0) %>%
  rownames()
core.df2[c(fixrow), 2] <- "< 0.01%"

#plot (all )
ggplot(data = core.df2, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  theme_classic() +
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
                                   face = "italic"),
        panel.background = element_blank())
  
ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_PrevAndAvgAbu_Eug.pdf",
       height = 15, unit = "cm")


#plot (without average relative abundance)
ggplot(data = core.df2, aes(x = Label, y = fct_rev(genus), fill = Prev2)) +
  geom_tile() +
  theme_classic() +
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
                                   face = "italic"),
        panel.background = element_blank())
ggsave("output/Prokaryote/CorePhylo/CorePhylos_vs_Cat_Prev_Eug.pdf")


#writeup
#need to remove \n
core.df2$Label <- str_replace_all(core.df2$Label, "\n", " ")
write.table(core.df2, "output/Prokaryote/CorePhylo/CorePhylotypes_dataframe.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

