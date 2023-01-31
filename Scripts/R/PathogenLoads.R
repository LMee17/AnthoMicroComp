#29th January 2023
#Pathogen versus potentially beneficial microbes

library(tidyverse)

dir.create("output/Pathogen/")

#load tables
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")
pro.rel <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_RelAb_MicroTaxaLabelled.tsv",
                    sep = "\t", header = T)
euk <- read.table("input/Counts/Eukaryote_Filtered_FamilyReduced_Raw_MicroTaxaLabelled.tsv")
euk.rel <- apply(euk, 2, FUN=function(x){ x / sum(x)})
euk.inc <- euk
euk.inc[euk.inc > 0] <- 1

vir <- read.table("input/Counts/Viral_Filtered_FamilyReduced_Raw.tsv",
                   header = T, sep = "\t")
vir.rel <- apply(vir, 2, FUN=function(x){x / sum(x)})
vir.inc <- vir
vir.inc[vir.inc > 0] <- 1

met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  sep ="\t", header = T)

tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  sep = "\t", header = T)

#list pathogens
bact <- c("Melissococcus", "Paenibacillus", "Brevibacillus", "Enteroccocus",
          "Achromobacter", "Spiroplasma", "Serratia")
badEu <- c("Crithidia", "Leishmania", "Leptomonas", "Lotmaria", "Trypanosoma", "Apicystis", "Nosema")

#all viruses will be considered as a form of pathogen load
corePop <- c("Gilliamella", "Apilactobacillus", "Bombilactobacillus",
             "Lactobacillus: Firm-5", "Bifidobacterium", "Snodgrassella")

#starting with pro: get list of all apis and bombus samples
#get relative abudances
#get 25% and 75% percentiles to grade low/mid/high levels of corbiculates relative
#to the average r.a. for each bee genus.
#also get a presence/absence in case that ends up easiest.
path.df <- pro.rel %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  subset(Genus == "Bombus" | Genus == "Apis") %>%
  select(id, name, value, Genus) %>%
  group_by(id, Genus) %>%
  mutate(avgRel = mean(value)) %>%
  mutate(lowBound = quantile(value, probs = c( .25)),
         highBound = quantile(value, (probs = c(.75))),
         inc = ifelse(value > 0, "Present", "Absent"),
         IncStatus = paste(id, inc)) %>%
          select(-avgRel) %>%
  filter(id %in% corePop)
#asses low/mid/high
for(i in 1:nrow(path.df)){
  if(path.df$value[i] < path.df$lowBound[i]){
    callit <- "Low"
  } else {
    if(path.df$value[i] > path.df$highBound[i]){
      callit <- "High"
    } else {
      callit <- "Mid"
    }
  }
  path.df$Level[i] <- callit
}  

#add bacterial load
bact.df <- pro.rel %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  subset(Genus == "Bombus" | Genus == "Apis") %>%
  filter(id %in% bact) %>%
  select(id, name, value, Species) %>%
  mutate(inc = ifelse(value > 0, 1, 0)) %>%
  group_by(name) %>%
  mutate(Bact_totA = sum(value),
         Bact_totInc = sum(inc)) %>%
  select(name, Bact_totA, Bact_totInc) %>%
  unique() %>%
  inner_join(., path.df, by = "name")

bact.df$Level <- factor(bact.df$Level,
                        levels = c("Low", "Mid", "High"))

#test plot
ggplot(data = bact.df, 
       aes(x = id, y = log(Bact_totA), fill = Level)) +
  geom_boxplot() 
#doesn't look like levels x relative abundnace of bacterial load is all that
#informative so sticking with incidence x incidence

ggplot(data = bact.df, 
       aes(x = inc, y = Bact_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#800080", "#008080"))
ggsave("output/Pathogen/BacterialOnly.pdf")

ggplot(data = bact.df, 
       aes(x = inc, y = Bact_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#800080", "#ff6dff"))
ggsave("output/Pathogen/BacterialOnly_purp.pdf")
#it looks like only Apilactobacillus and Bombilactobacillus may have any effect
#maybe bacterial ones should be taken with a pinch of salt - more reads for
#bacteria in general may mean both counts of "good" and "bad" taxa increase together

#eukaryote load
euk.df <- euk.rel %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  subset(Genus == "Bombus" | Genus == "Apis") %>%
  filter(id %in% badEu) %>%
  select(id, name, value, Species) %>%
  mutate(inc = ifelse(value > 0, 1, 0)) %>%
  group_by(name) %>%
  mutate(Euk_totA = sum(value),
         Euk_totInc = sum(inc)) %>%
  select(name, Euk_totA, Euk_totInc) %>%
  unique() %>%
  inner_join(., path.df, by = "name")

#plot
ggplot(data = euk.df, 
       aes(x = inc, y = Euk_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#800080", "#008080"))
ggsave("output/Pathogen/EukaryoteOnly.pdf")
#only affect seems to be with Apilactobacillus
#as this is eukaryote rna it shouldn't have been too affected by poly(A) enrichment
#and thus there is more robustness to any conclusions taken from this

ggplot(data = euk.df, 
       aes(x = inc, y = Euk_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#790000", "#ef0000"))
ggsave("output/Pathogen/EukaryoteOnly_red.pdf")

#viral load
vir.df <- vir.rel %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  subset(Genus == "Bombus" | Genus == "Apis") %>%
  select(id, name, value, Species) %>%
  mutate(inc = ifelse(value > 0, 1, 0)) %>%
  group_by(name) %>%
  mutate(Vir_totA = sum(value),
         Vir_totInc = sum(inc)) %>%
  select(name, Vir_totA, Vir_totInc) %>%
  unique() %>%
  inner_join(., path.df, by = "name")

#plot
ggplot(data = vir.df, 
       aes(x = inc, y = Vir_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#800080", "#008080"))
ggsave("output/Pathogen/VirusOnly.pdf")
#again, only affect seems to be with Apilactobacillus
#as this is eukaryote rna it shouldn't have been too affected by poly(A) enrichment
#and thus there is more robustness to any conclusions taken from this

ggplot(data = vir.df, 
       aes(x = inc, y = Vir_totInc, fill = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(~id,
             labeller = label_wrap_gen(width = 16)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#0884ff", "#7ebfff"))
ggsave("output/Pathogen/VirusOnly_blue.pdf")

#combine?
names(bact.df)[2:3] <- c("totA", "totInc")
names(euk.df)[2:3] <- c("totA", "totInc")
names(vir.df)[2:3] <- c("totA", "totInc")

bact.df$king <- paste("Bacteria")
euk.df$king <- paste("Eukaryote")
vir.df$king <- paste("Virus")

all.df <- rbind(bact.df, euk.df, vir.df)

ggplot(data = all.df, 
       aes(x = inc, y = totInc, fill = inc, colour = inc)) +
  theme_grey() +  
  geom_bar(stat = "identity") +
  facet_grid(king~id,
             labeller = label_wrap_gen(width = 16),
             scale = "free_y") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Incidence of Detected Pathogens",
       x = "",
       fill = "Corbiculate Core\nTaxa Status") +
  scale_fill_manual(values = c("#800080", "#008080")) +
  scale_colour_manual(values = c("#800080", "#008080")) +
  guides(colour = "none")
ggsave("output/Pathogen/All_Incidence.pdf")

#write up: nead more readable format
all.df <- all.df %>%
  ungroup() %>%
  select(totInc, id, inc, king) %>%
  group_by(id, inc, king) %>%
  mutate(sumTot = sum(totInc)) %>%
  select(-totInc) %>%
  unique()
  

names(all.df) <- c("CoreMicrobe_Taxa", "CoreMicrobe_Presence", 
                   "Pathogen_Kingdom", "Pathogen_No_Detected")

write.table(all.df, 
            "output/Pathogen/Pathogen_df.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#Apilactobacillus scatterplots
api.df <- all.df %>%
  subset(id == "Apilactobacillus")

ggplot(api.df, aes(x = inc, y = (totInc), fill = inc)) +
  geom_bar(stat = "identity") +
  facet_grid(~king) +
  labs(fill = "Apilactobacillus\nStatus",
       y =  "Number of Detected Pathogens",
       x = "") +
  theme(strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("#800080", "#008080"))
ggsave("output/Pathogen/Apilactobacillus_byInc.pdf")

#chisq ? 
chi.api <- api.df %>%
  ungroup() %>%
  select(totInc, inc, king) %>%
  group_by(king, inc) %>%
  mutate(sumPath = sum(totInc)) %>%
  select(-totInc) %>%
  unique() %>%
  pivot_wider(names_from = inc, values_from = sumPath) %>%
  as.data.frame() %>%
  column_to_rownames(var = "king") %>%
  as.matrix() %>%
  t()
all.test <- chisq.test(chi.api)
bac.test <- chisq.test(chi.api[,1])
euk.test <- chisq.test(chi.api[,2])
vir.test <- chisq.test(chi.api[,3])

all.test
#no sig when considered altogether
bac.test
#sig when considering just bacteria
euk.test
#sig when considering just eukaryota
vir.test
#sig when considering just viruses .... 
#this makes no sense to me

#fisher ? 
