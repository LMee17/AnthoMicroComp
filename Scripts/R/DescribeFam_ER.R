#26th January 2023
#Shallow look at Tribal microbial (eukaryote and virus) complements

##Libraries####
library(tidyverse)

dir.create("output/Eukaryote/Associates/")
dir.create("output/Viral/Associates/")

##Counts and Metadata####
#count
euk <- read.table("input/Counts/Eukaryote_Filtered_FamilyReduced_Raw.tsv")
vir <- read.table("input/Counts/Viral_Filtered_FamilyReduced_Raw.tsv")

#sample metadata
met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  sep = "\t",
                  header = T)

#microbial metadata
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  sep = "\t",
                  header = T)

##Tribe Descriptive Loop: eukaryote ####
e.tribes <- unique(met$Tribe[met$Sample.ID %in% names(euk)])
e.met <- met[met$Sample.ID %in% names(euk),]

triPrev <- data.frame()
for (i in 1:length(e.tribes)){
  samps <- e.met %>%
    subset(Tribe == e.tribes[i]) %>%
    select(Sample.ID) %>%
    unlist() %>%
    as.vector()
  ex.df <- euk %>%
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
  ex.df <- inner_join(ex.df, e.met, by = c("Sample" = "Sample.ID")) %>%
    select(-Family, -SubFamily, -Tissue2, -Tissue3)
  #add microbial metadata
  ex.df <- inner_join(ex.df, tax, by = c("TaxID"= "ID")) %>%
    as.data.frame() %>%
    select(-species, -kingdom, -superkingdom)
  over50 <- ex.df %>%
    filter(Prevalence >= 0.5) %>%
    filter(AvgRelAbundance >= 0.01) %>%
    select(genus, AvgRelAbundance ,Prevalence) %>%
    mutate(Tribe = paste(e.tribes[i])) %>%
    mutate(n = length(samps)) %>%
    unique()
  triPrev <- rbind(triPrev, over50)
  #write up
  write.table(ex.df,
              paste("output/Eukaryote/Associates/TribeBreakdown_",
                    e.tribes[i], ".tsv", sep = ""),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

write.table(triPrev,
            "output/Eukaryote/Associates/TopPrevMicrobes_byTribe.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#by family
e.fam <- unique(met$Family[met$Sample.ID %in% names(euk)])

famPrev <- data.frame()
for (i in 1:length(e.fam)){
  samps <- e.met %>%
    subset(Family == e.fam[i]) %>%
    select(Sample.ID) %>%
    unlist() %>%
    as.vector()
  ex.df <- euk %>%
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
  ex.df <- inner_join(ex.df, e.met, by = c("Sample" = "Sample.ID")) %>%
    select(-Tribe, -SubFamily, -Tissue2, -Tissue3)
  #add microbial metadata
  ex.df <- inner_join(ex.df, tax, by = c("TaxID"= "ID")) %>%
    as.data.frame() %>%
    select(-species, -kingdom, -superkingdom)
  over50 <- ex.df %>%
    filter(Prevalence >= 0.5) %>%
    filter(AvgRelAbundance >= 0.01) %>%
    select(genus, AvgRelAbundance ,Prevalence) %>%
    mutate(Family = paste(e.fam[i])) %>%
    mutate(n = length(samps)) %>%
    unique()
  famPrev <- rbind(famPrev, over50)
  #write up
  write.table(ex.df,
              paste("output/Eukaryote/Associates/FamilyBreakdown_",
                    e.fam[i], ".tsv", sep = ""),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

write.table(famPrev,
            "output/Eukaryote/Associates/TopPrevMicrobes_byFamily.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")


##Tribe Descriptive Loop: Virus ####
v.met <- met[met$Sample.ID %in% names(vir),]
v.tribes <- unique(met$Tribe[met$Sample.ID %in% names(euk)])

triPrev <- data.frame()
for (i in 1:length(v.tribes)){
  samps <- v.met %>%
    subset(Tribe == v.tribes[i]) %>%
    select(Sample.ID) %>%
    unlist() %>%
    as.vector()
  ex.df <- vir %>%
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
  ex.df <- inner_join(ex.df, v.met, by = c("Sample" = "Sample.ID")) %>%
    select(-Family, -SubFamily, -Tissue2, -Tissue3)
  #add microbial metadata
  ex.df <- inner_join(ex.df, tax, by = c("TaxID"= "ID")) %>%
    as.data.frame() %>%
    select(-species, -kingdom, -superkingdom)
  over50 <- ex.df %>%
    filter(Prevalence >= 0.5) %>%
    filter(AvgRelAbundance >= 0.01) %>%
    select(family, AvgRelAbundance ,Prevalence) %>%
    mutate(Tribe = paste(v.tribes[i])) %>%
    mutate(n = length(samps)) %>%
    unique()
  triPrev <- rbind(triPrev, over50)
  #write up
  write.table(ex.df,
              paste("output/Viral/Associates/TribeBreakdown_",
                    v.tribes[i], ".tsv", sep = ""),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

write.table(triPrev,
            "output/Viral/Associates/TopPrevMicrobes_byTribe.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#by family
v.fam <- unique(met$Family[met$Sample.ID %in% names(euk)])

famPrev <- data.frame()
for (i in 1:length(e.fam)){
  samps <- v.met %>%
    subset(Family == e.fam[i]) %>%
    select(Sample.ID) %>%
    unlist() %>%
    as.vector()
  ex.df <- vir %>%
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
  ex.df <- inner_join(ex.df, v.met, by = c("Sample" = "Sample.ID")) %>%
    select(-Tribe, -SubFamily, -Tissue2, -Tissue3)
  #add microbial metadata
  ex.df <- inner_join(ex.df, tax, by = c("TaxID"= "ID")) %>%
    as.data.frame() %>%
    select(-species, -kingdom, -superkingdom)
  over50 <- ex.df %>%
    filter(Prevalence >= 0.5) %>%
    filter(AvgRelAbundance >= 0.01) %>%
    select(family, AvgRelAbundance ,Prevalence) %>%
    mutate(Family = paste(e.fam[i])) %>%
    mutate(n = length(samps)) %>%
    unique()
  famPrev <- rbind(famPrev, over50)
  #write up
  write.table(ex.df,
              paste("output/Viral/Associates/FamilyBreakdown_",
                    e.fam[i], ".tsv", sep = ""),
              sep = "\t", col.names = T, row.names = F, quote = F)
}
write.table(famPrev,
            "output/Viral/Associates/TopPrevMicrobes_byFamily.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")
