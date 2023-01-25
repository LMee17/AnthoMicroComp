#11th January 2023
#Shallow look at Tribal microbial (prokaryote) complements

##Libraries####
library(tidyverse)

dir.create("output/Prokaryote/Tribal_Descriptives/")

##Counts and Metadata####
#count
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")

#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  sep = "\t",
                  header = T)
met <- met[met$Sample.ID %in% names(pro),]

#microbial metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  sep = "\t",
                  header = T)
tax <- tax[tax$ID %in% rownames(pro),]

##Tribe Descriptive Loop ####ß
tribes <- unique(met$Tribe)

triPrev <- data.frame()
for (i in 1:length(tribes)){
  samps <- met %>%
    subset(Tribe == tribes[i]) %>%
    select(Sample.ID) %>%
    unlist() %>%
    as.vector()
  ex.df <- pro %>%
    rownames_to_column(var = "TaxID") %>%
    pivot_longer(-TaxID) %>%
    filter(name %in% samps) %>%
    group_by(name) %>%
    mutate(TotReads = sum(value)) %>%
    mutate(RelAbundance = value / TotReads) %>%
    arrange(-RelAbundance) %>%
    group_by(TaxID) %>%
    mutate(AvgRelAbundance = sum(RelAbundance) / length(unique(name))) %>%
    mutate(Inc = ifelse(value > 0, 1, 0)) %>%
    mutate(Prevalence = sum(Inc) / length(unique(name))) %>%
    select(-Inc) %>%
    filter(Prevalence > 0)
  names(ex.df)[2:3] <- c("Sample", "NoReads")
  #add sample metadata
  #remove columns with only 1 level 
  ex.df <- inner_join(ex.df, met, by = c("Sample" = "Sample.ID")) %>%
    select(-Family, -SubFamily, -Tissue2, -Tissue3)
  #add microbial metadata
  ex.df <- inner_join(ex.df, tax, by = c("TaxID"= "ID")) %>%
    as.data.frame() %>%
    select(-species, -kingdom, -superkingdom)
  over50 <- ex.df %>%
    filter(Prevalence >= 0.5) %>%
    filter(AvgRelAbundance >= 0.01) %>%
    select(genus, AvgRelAbundance ,Prevalence) %>%
    mutate(Tribe = paste(tribes[i])) %>%
    mutate(n = length(samps)) %>%
    unique()
  triPrev <- rbind(triPrev, over50)
  #write up
  write.table(ex.df,
              paste("output/Prokaryote/Tribal_Descriptives/TribeBreakdown_",
              tribes[i], ".tsv", sep = ""),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

write.table(triPrev,
            "output/Prokaryote/Tribal_Descriptives/TopPrevMicrobes_byTribe.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

##SessionLogs####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/DescribeTribe_Jan23.txt")
