#20th January 2023
#Overlaps

#Libraries

library(tidyverse)
library(UpSetR)

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

