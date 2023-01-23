#20th January 2023
#Overlaps

dir.create("output/Overlaps/")

#Libraries

library(tidyverse)
library(UpSetR)
library(ggplot2)
library(ggVennDiagram)


#####Load Necessary Counts and Metadata####
#Sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_RNAOnly_Dec22.tsv",
                  sep = "\t", header = T)

#Microbial metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  header = T, sep = "\t")

#prokaryote count table
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")

##All Microbial Taxa, by Tribe####
tribe.df <- pro %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(-ID)
names(tribe.df)[2:3] <- c("Sample", "Count")
tribe.mat <- inner_join(tribe.df, met, by = c("Sample" = "Sample.ID")) %>%
  select(ID, Sample, Count, Tribe) %>%
  group_by(Tribe, ID) %>%
  mutate(tri.Count = sum(Count)) %>%
  select(ID, Tribe, tri.Count) %>%
  unique() %>%
  mutate(tri.Inc = ifelse(tri.Count > 0, 1, 0)) %>%
  ungroup() %>%
  select(-tri.Count) %>%
  pivot_wider(names_from = Tribe,
              values_from = tri.Inc) %>%
  column_to_rownames("ID")
tribe.mat <- as.data.frame(tribe.mat)

upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
              "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      order.by = "freq", keep.order = T)

pdf("output/Overlaps/Tribe_v_ProCom_upsetr.pdf")
  upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      order.by = "freq", keep.order = T)
dev.off()

upset.met <- data.frame(Tribe = names(tribe.mat))
#upset.met$Sociality <- c("Eusocial", "Eusocial", "Polymorphic", "Eusocial",
 #                        "Polymorphic", "Solitary", "Solitary", "Polymorphic")
for (i in 1:nrow(upset.met)){
  upset.met$Family[i] <- unique(met$Family[met$Tribe == upset.met$Tribe[i]])
}
str(upset.met)

beecore <- c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
             "Bartonella", "Bombiscardovia", "Schmidhempelia", "Apibacter")
coreID <- tax %>%
  subset(genus %in% beecore) %>%
  select(ID) %>%
  unlist %>% as.vector()

tribe.mat$Core <- ifelse(rownames(tribe.mat) %in% coreID, "yes", "no")
head(tribe.mat)

pdf("output/Overlaps/Tribe_v_ProCom_ColorCodedAll_upsetr.pdf")
  upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      nintersects = 15,
      set.metadata = list(data = upset.met, 
                          plots = list(list(type = "matrix_rows", 
                                            column = "Family", assign = 10, 
                                            colors = c(Apidae = "#cfe2f3",
                                                              Megachilidae = "#d9ead3",
                                                              Halictidae = "#fff2cc",
                                                              Andrenidae = "#ea9999")))),
      shade.alpha = 1,
      sets.bar.color = c("#3b3bc4", "#009292",
                         "#3b3bc4",
                         "#009292", "#009292",
                         "#db6d00","#db6d00","#db6d00"),
      sets.x.label = "No Detected \nProkaryote Species",
      queries = list(list(query = elements,
                              params = list("Core", "yes"), color = "maroon", 
                              active = T,
                              query.name = "Core Phylotypes"))) 
dev.off()

pdf("output/Overlaps/Tribe_v_ProCom_ColorCoded_SocOnly_upsetr.pdf")
upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      nintersects = 15,
      shade.alpha = 1,
      sets.bar.color = c("#3b3bc4", "#009292",
                         "#3b3bc4",
                         "#009292", "#009292",
                         "#db6d00","#db6d00","#db6d00"),
      sets.x.label = "No Detected \nProkaryote Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

##Sociality####
#venn diagrams baby

soc.df <- pro %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(-ID)
names(soc.df)[2:3] <- c("Sample", "Count")
soc.mat <- inner_join(tribe.df, met, by = c("Sample" = "Sample.ID")) %>%
  select(ID, Sample, Count, Sociality) %>%
  group_by(Sociality, ID) %>%
  mutate(soc.Count = sum(Count)) %>%
  select(ID, Sociality, soc.Count) %>%
  unique() %>%
  mutate(soc.Inc = ifelse(soc.Count > 0, 1, 0)) %>%
  ungroup() %>%
  select(-soc.Count) %>%
  pivot_wider(names_from = Sociality,
              values_from = soc.Inc) %>%
  as.data.frame()
eu <- soc.mat %>%
  filter(Eusocial == 1) %>%
  select(ID) %>%
  unlist %>%
  as.vector()
po <- soc.mat %>%
  filter(Polymorphic == 1) %>%
  select(ID) %>%
  unlist %>%
  as.vector()
so <- soc.mat %>%
  filter(Solitary == 1) %>%
  select(ID) %>%
  unlist %>%
  as.vector()

soc.plot <- list(Solitary = so,
                 Eusocial = eu,
                 Polymorphic = po)

ggVennDiagram(soc.plot) +
  scale_x_continuous(expand = expansion(mult = .25))

#attempt #2
soc.venn <- Venn(soc.plot)
soc.venn <- process_data(soc.venn)

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(soc.venn), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(soc.venn), show.legend = FALSE) + 
  geom_sf_text(aes(label = name), size = 5, hjust = .4, vjust = -.1,
               fontface = "bold",
               data = venn_setlabel(soc.venn)) +
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(soc.venn)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("#3b3bc4", "#8b5462", "#5d6972", "#1e67ab",
                               "#db6d00", "#6e8049", "#009292"))

ggsave("output/Overlaps/Venn_Socs_AllPro_Labelled.pdf")

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(soc.venn), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(soc.venn), show.legend = FALSE) + 
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(soc.venn)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("#3b3bc4", "#8b5462", "#5d6972", "#1e67ab",
                               "#db6d00", "#6e8049", "#009292"))

ggsave("output/Overlaps/Venn_Socs_AllPro_UnLabelled.pdf")

##Using Top Tribe Hits Only ####
top <- read.table("output/Prokaryote/Tribal_Descriptives/TopPrevMicrobes_byTribe.tsv",
                  header = T, sep = "\t")

top.mat <- top %>%
  select(genus, Tribe) %>%
  group_by(Tribe, genus) %>%
  mutate(Inc = 1) %>%
  pivot_wider(names_from = Tribe, 
              values_from = Inc,
              values_fill = 0) %>%
  ungroup() %>%
 column_to_rownames("genus")

upset(top.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      order.by = "freq", keep.order = T)

pdf("output/Overlaps/Tribe_v_TopPro_upsetr.pdf")
upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      order.by = "freq", keep.order = T)
dev.off()

top.mat$Core <- ifelse(rownames(top.mat) %in% beecore, "yes", "no")
head(top.mat)

pdf("output/Overlaps/Tribe_v_Top_ColorCodedAll_upsetr.pdf")
upset(top.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      set.metadata = list(data = upset.met, 
                          plots = list(list(type = "matrix_rows", 
                                            column = "Family", assign = 10, 
                                            colors = c(Apidae = "#cfe2f3",
                                                       Megachilidae = "#d9ead3",
                                                       Halictidae = "#fff2cc",
                                                       Andrenidae = "#ea9999")))),
      shade.alpha = 1,
      sets.bar.color = c("#3b3bc4", "#009292",
                         "#3b3bc4",
                         "#009292", "#009292",
                         "#db6d00","#db6d00","#db6d00"),
      sets.x.label = "No Potential \nTribe-specific Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

pdf("output/Overlaps/Tribe_v_Top_ColorCoded_SocOnly_upsetr.pdf")
upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      nintersects = 15,
      shade.alpha = 1,
      sets.bar.color = c("#3b3bc4", "#009292",
                         "#3b3bc4",
                         "#009292", "#009292",
                         "#db6d00","#db6d00","#db6d00"),
      sets.x.label = "No Detected \nProkaryote Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

##Using Top Tribe Hits Only ####
for (i in 1:nrow(top)){
  top$Sociality[i] <- unique(met$Sociality[met$Tribe == top$Tribe[i]])
}

socs <- unique(top$Sociality)

topz <- vector(mode = "list", length = length(socs))
for(i in 1:length(socs)){
  topz[[i]] <- subset(top, Sociality == socs[i]) %>% select(genus) %>% unlist %>% as.vector()
}
names(topz) <- socs
topz <- topz[c("Solitary", "Eusocial", "Polymorphic")]
top.venn <- Venn(topz)
top.venn <- process_data(top.venn)

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(top.venn), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(top.venn), show.legend = FALSE) + 
  geom_sf_text(aes(label = name), size = 5, hjust = .4, vjust = -.1,
               fontface = "bold",
               data = venn_setlabel(top.venn)) +
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(top.venn)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("#3b3bc4", "#8b5462", "#5d6972", "#1e67ab",
                               "#db6d00", "#6e8049", "#009292"))

ggsave("output/Overlaps/Venn_Socs_Top_Labelled.pdf")

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(top.venn), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(top.venn), show.legend = FALSE) + 
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(top.venn)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("#3b3bc4", "#8b5462", "#5d6972", "#1e67ab",
                               "#db6d00", "#6e8049", "#009292")) +
  labs(title = "Associated Species\nby Sociality")

ggsave("output/Overlaps/Venn_Socs_Top_UnLabelled.pdf")

##Everything ... but by Family #####
fam.df <- pro %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(-ID)
names(fam.df)[2:3] <- c("Sample", "Count")
fam.mat <- inner_join(fam.df, met, by = c("Sample" = "Sample.ID")) %>%
  select(ID, Sample, Count, Family) %>%
  group_by(Family, ID) %>%
  mutate(fam.Count = sum(Count)) %>%
  select(ID, Family, fam.Count) %>%
  unique() %>%
  mutate(fam.Inc = ifelse(fam.Count > 0, 1, 0)) %>%
  ungroup() %>%
  select(-fam.Count) %>%
  pivot_wider(names_from = Family,
              values_from = fam.Inc) %>%
  column_to_rownames("ID")
fam.mat <- as.data.frame(fam.mat)

upset(fam.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Megachilidae",
                "Halictidae",
               "Apidae",
               "Andrenidae"),
      order.by = "freq", keep.order = T)

pdf("output/Overlaps/Family_v_ProCom_upsetr.pdf")
upset(fam.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Megachilidae",
               "Halictidae",
               "Apidae",
               "Andrenidae"),
      order.by = "freq", keep.order = T)
dev.off()

fam.mat$Core <- ifelse(rownames(fam.mat) %in% coreID, "yes", "no")
head(fam.mat)

pdf("output/Overlaps/Fam_v_ProCom_ColorCoded_FamOnly_upsetr.pdf")
upset(fam.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Megachilidae",
               "Halictidae",
               "Apidae",
               "Andrenidae"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      nintersects = 15,
      shade.alpha = 1,
      sets.bar.color = c("#7fffbf", 
                         "#ffff00",
                         "blue", 
                         "red"),
      sets.x.label = "No Detected \nProkaryote Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

fam.top.mat <- inner_join(top, met, by = "Tribe") %>%
  select(genus, Family) %>%
  unique() %>%
  group_by(genus, Family) %>%
  mutate(Inc = 1) %>%
  pivot_wider(names_from = Family, 
              values_from = Inc,
              values_fill = 0) %>%
  ungroup() %>%
  column_to_rownames("genus")

pdf("output/Overlaps/Fam_v_TopPro_upsetr.pdf")
upset(fam.top.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Megachilidae",
               "Halictidae",
               "Apidae",
               "Andrenidae"),
      order.by = "freq", keep.order = T)
dev.off()

fam.top.mat$Core <- ifelse(rownames(fam.top.mat) %in% beecore, "yes", "no")
head(fam.top.mat)

pdf("output/Overlaps/Fam_v_Top_ColorCodedAll_upsetr.pdf")
upset(fam.top.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Megachilidae",
               "Halictidae",
               "Apidae",
               "Andrenidae"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      shade.alpha = 1,
      sets.bar.color = c("#7fffbf", 
                         "#ffff00",
                         "blue", 
                         "red"),
      sets.x.label = "No Potential \nTribe-specific Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

pdf("output/Overlaps/Tribe_v_Top_ColorCoded_SocOnly_upsetr.pdf")
upset(tribe.mat, mb.ratio = c(0.55, 0.45),
      sets = c("Osmiini",
               "Augochlorini",
               "Andrenini",
               "Ceratini", "Euglossini","Meliponini","Bombini","Apini"),
      query.legend = "top",
      order.by = "freq", keep.order = T,
      nintersects = 15,
      shade.alpha = 1,
      sets.bar.color = c("#3b3bc4", "#009292",
                         "#3b3bc4",
                         "#009292", "#009292",
                         "#db6d00","#db6d00","#db6d00"),
      sets.x.label = "No Detected \nProkaryote Species",
      queries = list(list(query = elements,
                          params = list("Core", "yes"), color = "maroon", 
                          active = T,
                          query.name = "Core Phylotypes"))) 
dev.off()

##As a Venn #####

for (i in 1:nrow(top)){
  top$Family[i] <- unique(met$Family[met$Tribe == top$Tribe[i]])
}
top
fams <- unique(top$Family)

topz2 <- vector(mode = "list", length = length(fams))
for(i in 1:length(fams)){
  topz2[[i]] <- subset(top, Family == fams[i]) %>% select(genus) %>% unlist %>% as.vector()
}
names(topz2) <- fams

top.venn2 <- Venn(topz2)
top.venn2 <- process_data(top.venn2)

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(top.venn2), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(top.venn2), show.legend = FALSE) + 
  geom_sf_text(aes(label = name), size = 5, hjust = .4, vjust = -.1,
               fontface = "bold",
               data = venn_setlabel(top.venn2)) +
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(top.venn2)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("blue","#4080df","#7f5595","#9f8070","#7faa95",
                               "#800080","#aa5555","#808080","#7FFFBF","#bf8060",
                               "#d4aa40","#bfff60","red","#ff8000","yellow"))

ggsave("output/Overlaps/Venn_Fams_Top_Labelled.pdf")

ggplot() +
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(top.venn2), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(top.venn2), show.legend = FALSE) + 
  geom_sf_label(aes(label = count), fontface = "bold", size = 6,
                label.size = 0.1, data = venn_region(top.venn2)) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_fill_manual(values = c("blue","#4080df","#7f5595","#9f8070","#7faa95",
                               "#800080","#aa5555","#808080","#7FFFBF","#bf8060",
                               "#d4aa40","#bfff60","red","#ff8000","yellow"))

ggsave("output/Overlaps/Venn_Fams_Top_UnLabelled.pdf")
