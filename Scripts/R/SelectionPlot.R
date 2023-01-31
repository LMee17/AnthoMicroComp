#31st January 2023
#Assess library selection method vs number of reads per kingdom plot

library(tidyverse)

#load tables
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")
euk <- read.table("input/Counts/Eukaryote_Filtered_FamilyReduced_Raw_MicroTaxaLabelled.tsv")
vir <- read.table("input/Counts/Viral_Filtered_FamilyReduced_Raw.tsv",
                  header = T, sep = "\t")

met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  sep ="\t", header = T)

tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  sep = "\t", header = T)

p.plot <- pro %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value)) %>%
  ggplot(aes(x = LibrarySelection, y = log(tot))) +
  geom_boxplot(aes(fill = LibrarySelection)) +
  scale_fill_manual(values = myPal) +
  labs(y = "log(Total Reads)",
       x = "Method of Library Selection") +
  guides(fill = "none")
ggsave("output/SampleStats/LibrarySelectionvsTotalRead_Pro.pdf")

e.plot <- euk %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value)) %>%
  ggplot(aes(x = LibrarySelection, y = log(tot))) +
  geom_boxplot(aes(fill = LibrarySelection)) +
  scale_fill_manual(values = myPal) +
  labs(y = "log(Total Reads)",
       x = "Method of Library Selection") +
  guides(fill = "none")
ggsave("output/SampleStats/LibrarySelectionvsTotalRead_Euk.pdf")

v.plot <- vir %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value)) %>%
  ggplot(aes(x = LibrarySelection, y = log(tot))) +
  geom_boxplot(aes(fill = LibrarySelection)) +
  scale_fill_manual(values = myPal[c(1,3:10)]) +
  labs(y = "log(Total Reads)",
       x = "Method of Library Selection") +
  guides(fill = "none")
ggsave("output/SampleStats/LibrarySelectionvsTotalRead_Vir.pdf")

library(ggpubr)

ggarrange(p.plot, e.plot, v.plot, ncol = 1, label.x = "Library Selection Method")

myPal <- 	c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
            "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
            "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
            "#9acd32")


v <- vir %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value)) %>%
  mutate(kingdom = paste("Virus"))
e <- euk %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value))%>%
  mutate(kingdom = paste("Eukaryote"))

p <- pro %>%
  rownames_to_column(var = "id") %>%
  pivot_longer(-id) %>%
  inner_join(., met, by = c("name" = "Sample.ID")) %>%
  select(name, value, LibrarySelection) %>%
  group_by(name) %>%
  mutate(tot = sum(value)) %>%
  mutate(kingdom = paste("Bacteria"))
all <- rbind(e,p,v)
all %>%
  ggplot(aes(x = LibrarySelection, y = log(tot))) +
  geom_boxplot(aes(fill = LibrarySelection)) +
  scale_fill_manual(values = myPal[c(1,3:10)]) +
  labs(y = "log(Total Reads)",
       x = "Method of Library Selection") +
  guides(fill = "none") +
  facet_grid(kingdom~.)
ggsave("output/SampleStats/LibraryPrep_All.pdf")

