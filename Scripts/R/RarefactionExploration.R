#3rd January 2023
#Do I need to rarefy ?

#Prepared in preparation for discussion with Seth or justification for viva
#Following tutorial : https://www.youtube.com/watch?v=ht3AX5uZTTQ&list=PLmNrK_nkqBpJuhS93PYC-Xr5oqur7IIWf
#Original script found: https://github.com/riffomonas/distances/blob/0f64c79deabf30603205ffc983abc6674cd121b4/code/distances.R

library(tidyverse)
library(vegan)


#make a null community model from my data (mixing up all samples and microbial
#taxa so there is no effect re host / microbe)
cnt.null <- read.table("input/Counts/Prokaryote_RawCounts.tsv") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var= "Sample") %>%
  pivot_longer(-Sample) %>%
  #group by sample and remove any that have less than 100 total reads
  group_by(Sample) %>%
  mutate(total = sum(value)) %>%
  filter(total > 100) %>%
  #group by microbe and remove any empty rows
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  #remove zero counts
  filter(value != 0)

rand <- cnt.null %>%
  #a function that produces rows per number of of a sample
  #ie if there are 35 counts of TRH2 in AM_N_072, then there were will be
  #72 rows of AM_N_072 / TRH2
  uncount(value) %>%
  #randomise the microbial hit
  mutate(rand_name = sample(name)) %>%
  #remove the actual hit
  select(-name) %>%
  #produce the new count
  count(Sample, rand_name)
  
#produce a tibble with the number of reads per sample
shared_samp_count <- cnt.null %>% 
  group_by(Sample) %>%
  summarize(n = sum(value))

#the same as above but with the randomised object
rand_samp_count <- rand %>%
  group_by(Sample) %>%
  summarize(n = sum(n))

#check that we still have the same # of reads per sample
inner_join(shared_samp_count, rand_samp_count,  by="Sample")

#check the same with microbial taxa
shared_trh_count <- cnt.null %>%
  group_by(name) %>%
  summarize(n = sum(value))

rand_trh_count <- rand %>%
  group_by(rand_name) %>%
  summarize(n = sum(n))

inner_join(shared_trh_count, rand_trh_count,  by=c("name" = "rand_name" ))

#plot the distribution fo reads
rand_samp_count %>%
  ggplot(aes(x=n)) + geom_histogram()

#convert the null sample set into a dataframe
rand_df <- rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill = 0) %>%
  as.data.frame()
#remove the group variable and convert into a matrix for use in vegan
rownames(rand_df) <- rand_df$Sample
rand_df <- rand_df[, -1]
rand_matrix <- as.matrix(rand_df)

#run the distance on the randomised set with vegdist (no rarefaction)
norare_dist_matrix <- vegdist(rand_matrix, method="bray")

#run the distance function on the randomised set with rarefaction (153, 
#lowest number of reads in solitary samples), with 10,000 iterations (avgdist takes
#the average distance of each iteration)
rare_dist_matrix <- avgdist(rand_matrix, dmethod="bray", sample=153, iterations = 10000)
#The following sampling units were removed because they were below sampling depth: SRR12053415, SRR12053434, SRR12527932, SRR13442814

norare_dist_tibble <- norare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

rare_dist_tibble <- rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

comparison <- inner_join(norare_dist_tibble, rare_dist_tibble, by=c("sample", "name")) %>%
  select(sample, name, norare=value.x, rare=value.y) %>%
  inner_join(., rand_samp_count, by=c("sample" = "Sample")) %>%
  inner_join(., rand_samp_count, by=c("name" = "Sample")) %>%
  mutate(n_diff = abs(n.x-n.y)) %>%
  select(-n.x, -n.y)

comparison %>%
  ggplot(aes(x=norare, y=rare, color=n_diff)) +
  geom_point(size=0.25, alpha=0.25) +
  geom_smooth()

comparison %>%
  pivot_longer(cols=c("norare", "rare"), names_to="type", values_to="dist") %>%
  ggplot(aes(x=n_diff,  y=dist)) +
  geom_point(size=0.25, alpha=0.25) +
  facet_wrap(~type, nrow=2)
#shows that the difference in number of sequences drives the distances within the 
#matrices

