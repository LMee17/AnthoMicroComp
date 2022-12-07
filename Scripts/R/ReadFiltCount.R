#7th December 2022
#A script to read in report files, filter according to parameters (NR/NT percent identity
#> 90%, NR/NT rPM > 1 for non-viruses, the same parameters but just for NR for viruses)
#and produce count tables per blast hit taxon (to start) and then collapse to
#genera (family?)

#load necessary libraries
library(stringr)
library(dplyr)

###load necessary metadata and taxonomic files####
#sample metadata 
met <- read.table("input/Metadata/MetaData_Edit_Nov22.tsv",
                  sep = "\t", header = T)
#eukaryotic tax classifications
euk <- read.csv("input/Phylo_Misc/Eukaryota_Classifications_Nov22.csv")
#taxonomy file
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22_v2.tsv",
                  sep = "\t", header = T)
#updated taxonomic terms (those that have changed since CZID last updated their databases)
tax.up <- read.table("input/Phylo_Misc/TaxoUpdate.txt",
                     header = T, sep = "\t", quote = "")
#microbial metadata
microkey <- read.table("input/Metadata/MicrobialSpeciesKey.tsv",
                       sep = "\t", header = T)
#higher level phylogeny microbial metadata
genkey <- read.table("input/Metadata/GenusFamilySpeciesKey_Nov22.tsv",
                     sep = "\t", header = T)

###read in czid report files####
#list taxon report files
csvs <- list.files(path = "input/CZID_TaxonReports/", pattern = "*.csv")

#prepare a list ready to be populated with report files
data.raw <- vector(mode = "list", length = length(csvs))

#iterate through list of files, reading in each into the data list
for (i in 1:length(csvs)){
  data.raw[[i]] <- read.csv(paste("input/CZID_TaxonReports/", csvs[i], sep = ""),
                            header = T)
}

###starting with the non-viral filter#####
#prepare a list space to populate
data.filt.nv <- vector(mode = "list", length = length(data.raw))

for (i in 1:length(data.filt.nv)){
  #only keep necessary columns
  data.filt.nv[[i]] <- data.raw[[i]][,c(4,6,8,11,12,15,23,24,27)]
  #get sample ID from the filename
  samp <- paste(strsplit(csvs[i], "_2[0-9]")[[1]][1])
  samp <- str_remove(samp, "_[1-2]")
  #add sample ID to dataframe for ease of identification when these are combined.
  data.filt.nv[[i]]$SampleID <- paste(samp)
  #get species  and genus information from the metadata object
  data.filt.nv[[i]]$Species <- unique(met$Species[met$Sample.ID == samp])
  data.filt.nv[[i]]$Genus <- unique(met$Genus[met$Sample.ID == samp])
  #filter the results by the parameters outlined above
  data.filt.nv[[i]] <- subset(data.filt.nv[[i]], nt_rpm > 1 &
                                nt_percent_identity > 90 &
                                nr_rpm > 1 &
                                nr_percent_identity > 90)
}

#combine into a single dataframe
filt.df.nv <- bind_rows(data.filt.nv)
#remove viruses
filt.df.nv <- subset(filt.df.nv, category != "viruses")

#add genus information as a separate column
for(i in 1:nrow(filt.df.nv)){
  gen <- strsplit(filt.df.nv$name[i], " ")
  filt.df.nv$NonHostPhylo[i] <- paste(gen[[1]][1])
}

#for some reason the dataframe has amassed rows of all NAs. 
#These are causing issues and need to be removed
filt.df.nv <- filt.df.nv[grepl("^NA", rownames(filt.df.nv))==F,]

#how to deal with the parts that don't have genera?
#the genus will instead be the Family name (these have been added manually to the euk
#filtering file)
for (i in 1:nrow(filt.df.nv)){
  if (filt.df.nv$NonHostPhylo[i] == "non-genus-specific"){
    filt.df.nv$NonHostPhylo[i] <- strsplit(filt.df.nv$name[i], " ")[[1]][5]
  }
}

#there's also this guy in there
filt.df.nv$NonHostPhylo[filt.df.nv$name == "parasitid 'Pas'"] <- "Parasitoidea"

#make lists of wanted eukaryota
wantedClass <- c("Fungus", "Trypanosome", "Unicellular")
tokeep <- euk$Genus[euk$Classification %in% wantedClass]

#keep only wanted taxa
filt.df.nv <- filt.df.nv[filt.df.nv$NonHostPhylo %in% tokeep |
                           filt.df.nv$category == "bacteria" |
                           filt.df.nv$category == "archaea",]

#to allow for easier recall of fungal hits
fungi <- euk$Genus[euk$Classification == "Fungus"]

#add fungus to classifications
filt.df.nv$category[filt.df.nv$NonHostPhylo %in% fungi] <- "fungus"

#change eukaryota classification to make more sense
filt.df.nv$category[filt.df.nv$category == "eukaryota"] <- "other_eukaryota"

#Now that that filtering is done, I don't need the NonHostPhylo column
filt.df.nv$NonHostPhylo <- NULL



###viral filter####
data.filt.v <- vector(mode = "list", length = length(data.raw))

for (i in 1:length(data.filt.v)){
  #only keep necessary columns
  data.filt.v[[i]] <- data.raw[[i]][,c(4,6,8,11,12,15,23,24,27)]
  #get sample ID from the filename
  samp <- paste(strsplit(csvs[i], "_2[0-9]")[[1]][1])
  samp <- str_remove(samp, "_[1-2]")
  #add sample ID to dataframe for ease of identification when these are combined.
  data.filt.v[[i]]$SampleID <- paste(samp)
  #get species  and genus information from the metadata object
  data.filt.v[[i]]$Species <- unique(met$Species[met$Sample.ID == samp])
  data.filt.v[[i]]$Genus <- unique(met$Genus[met$Sample.ID == samp])
  #filter the results by the parameters outlined above
  data.filt.v[[i]] <- subset(data.filt.v[[i]], nr_rpm > 1 &
                               nr_percent_identity > 90)
}

filt.df.v <- bind_rows(data.filt.v)
filt.df.v <- subset(filt.df.v, category == "viruses")

#combine with non-viral filtered data
filt.df <- rbind(filt.df.nv, filt.df.v)

#And lastly, remove an unnecessary entry 
filt.df <- subset(filt.df, 
                  name != "all taxa with neither family nor genus classification")

###Raw Read Hits Count Tables: By Species ####
tax.hit <- unique(microkey$TaxHit)
hitids <- unique(microkey$HitID)

#compile the samples
samples <- unique(filt.df$SampleID)

#prepare an empty matrix to populate
cnt <- matrix(nrow = length(tax.hit), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(tax.hit))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(microkey$MicroHit[microkey$TaxHit == tax.hit[j]])
    t <- filt.df$nr_count[filt.df$name %in% y &
                            filt.df$SampleID == samples[i]]
    if (length(t)==0){
      t <- 0
    }
    #check for duplicates within the sample
    if (length(t) == 2){
      print("Here!")
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
    #add sample's counts to count matrix
  }
  cnt[,i] <- x
}
#add sample and microbe IDs
rownames(cnt) <- unique(microkey$HitID)
colnames(cnt) <- samples

#convert to dataframe for further manipulation
cnt <- as.data.frame(cnt)

#prepare to subset: add class to microkey
for (i in 1:nrow(microkey)){
  microkey$Class[i] <- unique(filt.df$category[filt.df$name == microkey$MicroHit[i]])
}

#extract the MRH ids that correspond to bacteria/archaea
ext <- unique(microkey$HitID[microkey$Class == "bacteria" | microkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(microkey$HitID[microkey$Class == "fungus" | 
                               microkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(microkey$HitID[microkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

#create a directory ready to fill 
dir.create("input/Counts/RawReads")
path <- "input/Counts/RawReads/"

#combine subsets into a listed object for ease of write up
cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)
#store subset file names
filez <- c("All", "Prokaryote", "Eukaryote", "Viral")
#write up by iterating through the list of subsetted dataframes and file dinstinguishers
for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "rawCounts", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}

#save the microbial metadata / key 
dir.create("input/Metadata/")
write.table(microkey, "input/Metadata/MicrobialSpeciesKey.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)


###Raw Read Hit Counts. Collapsed by Genus/Family####
#prepare an empty matrix to populate
phyloz <- unique(genkey$PhyloHit)
genhitids <- paste0("GRH", 1:length(phyloz), sep = "")

cnt <- matrix(nrow = length(phyloz), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(phyloz))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(genkey$MicroHit[genkey$PhyloHit == phyloz[j]])
    t <- filt.df$nr_count[filt.df$name %in% y &
                            filt.df$SampleID == samples[i]]
    if (length(t)==0){
      t <- 0
    }
    #check for duplicates within the sample
    if (length(t) > 1){
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
    #add sample's counts to count matrix
    cnt[,i] <- x
  }
}

rownames(cnt) <- unique(genkey$GenHitID)
colnames(cnt) <- samples

cnt <- as.data.frame(cnt)

#extract the GRH ids that correspond to bacteria/archaea
ext <- unique(genkey$GenHitID[genkey$Class == "bacteria" | genkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(genkey$GenHitID[genkey$Class == "fungus" | 
                                genkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(genkey$GenHitID[genkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

#set path
path <- "input/Counts/RawReads/"

#make list of the count dataframes
cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)

filez <- c("All", "Prokaryote", "Eukaryote", "Viral")

for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "rawCounts_PhyloCollapsed", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}


###Reads per Million: By Species####
#prepare an empty matrix to populate
cnt <- matrix(nrow = length(tax.hit), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(tax.hit))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(microkey$MicroHit[microkey$TaxHit == tax.hit[j]])
    t <- filt.df$nr_rpm[filt.df$name %in% y &
                          filt.df$SampleID == samples[i]]
    if (length(t)==0){
      t <- 0
    }
    #check for duplicates within the sample
    if (length(t) == 2){
      print("Here!")
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
  }
  #add sample's counts to count matrix
  cnt[,i] <- x
}

rownames(cnt) <- unique(microkey[microkey$HitID])
colnames(cnt) <- samples

cnt <- as.data.frame(cnt)

#extract the MRH ids that correspond to bacteria/archaea
ext <- unique(microkey$HitID[microkey$Class == "bacteria" |
                               microkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(microkey$HitID[microkey$Class == "fungus" | 
                               microkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(microkey$HitID[microkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

dir.create("input/Counts/Raw_rPM")
path <- "input/Counts/Raw_rPM/"

cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)

filez <- c("All", "Prokaryote", "Eukaryote", "Viral")

for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "raw_rPM", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}
###Reads per Million: Collapsed by Genus/Family####
#prepare an empty matrix to populate
cnt <- matrix(nrow = length(phyloz), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(phyloz))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(genkey$MicroHit[genkey$PhyloHit == phyloz[j]])
    t <- filt.df$nr_rpm[filt.df$name %in% y &
                          filt.df$SampleID == samples[i]]
    if (length(t)==0){
      t <- 0
    }
    #check for duplicates within the sample
    if (length(t) > 1){
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
    #add sample's counts to count matrix
    cnt[,i] <- x
  }
}

rownames(cnt) <- unique(genkey$GenHitID)
colnames(cnt) <- samples

cnt <- as.data.frame(cnt)

#extract the GRH ids that correspond to bacteria/archaea
ext <- unique(genkey$GenHitID[genkey$Class == "bacteria" | genkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(genkey$GenHitID[genkey$Class == "fungus" | 
                                genkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(genkey$GenHitID[genkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

#set path
path <- "input/Counts/Raw_rPM/"


#make list of the count dataframes
cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)

filez <- c("All", "Prokaryote", "Eukaryote", "Viral")

for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "raw_rPM_PhyloCollapsed", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}
###Incidence: Per Species#####
#prepare an empty matrix to populate
cnt <- matrix(nrow = length(tax.hit), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(tax.hit))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(microkey$MicroHit[microkey$TaxHit == tax.hit[j]])
    t <- filt.df$nr_count[filt.df$name %in% y &
                            filt.df$SampleID == samples[i]]
    if (length(t)==0){
      t <- 0
    }
    if (length(t) == 2){
      t <- sum(t)
    }
    #override counts that aren't 0 with a 1
    if (t >= 1){
      t <- 1
    }
    #populate sample's vector with counts
    x[j] <- t
  }
  #add sample's incidences to count matrix
  cnt[,i] <- x
}

rownames(cnt) <- hitids
colnames(cnt) <- samples

cnt <- as.data.frame(cnt)

#extract the MRH ids that correspond to bacteria/archaea
ext <- unique(microkey$HitID[microkey$Class == "bacteria" | 
                               microkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(microkey$HitID[microkey$Class == "fungus" | 
                               microkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(microkey$HitID[microkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

#create directory to save reads into
dir.create("input/Counts/Incidence")
path <- "input/Counts/Incidence/"

#compile the count tables
cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)

filez <- c("All", "Prokaryote", "Eukaryote", "Viral")

for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "Incidence", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}

###Incidence: Collapsed by Genus/Family ####
#prepare an empty matrix to populate
cnt <- matrix(nrow = length(phyloz), ncol = length(samples))

for (i in 1:length(samples)){
  print(paste(i, "/317: ", samples[i], sep = ""))
  x <- vector(length = length(phyloz))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    y <- unique(genkey$MicroHit[genkey$PhyloHit == phyloz[j]])
    t <- filt.df$nr_count[filt.df$name %in% y &
                            filt.df$SampleID == samples[i]]
    #in instances were there was nothing to extract, mark the absence with a 0
    if (length(t)==0){
      t <- 0
      #however, if there was data to extract, assess
    } else {
      #run through the reads extracted and replace values with the presence marker 1
      for (k in 1:length(t)){
        if (t[k] >= 1){
          t[k] <- 1
        }
      }
      #sum the number of presences
      t <- sum(t)
    }
    #populate the count vector
    x[j] <- t
  }
  #add this sample counts to the matrix
  cnt[,i] <- x
}  

rownames(cnt) <- genhitids
colnames(cnt) <- samples

cnt <- as.data.frame(cnt)

#extract the GRH ids that correspond to bacteria/archaea
ext <- unique(genkey$GenHitID[genkey$Class == "bacteria" | genkey$Class == "archaea"])
#extract these reads from the main count table
cnt.ab <- cnt[rownames(cnt) %in% ext,]

#repeat for fungi and other eukaryotes
ext <- unique(genkey$GenHitID[genkey$Class == "fungus" | 
                                genkey$Class == "other_eukaryota"])
cnt.fe <- cnt[rownames(cnt) %in% ext,]  

#and finally, viruses
ext <- unique(genkey$GenHitID[genkey$Class == "viruses"])
cnt.v <- cnt[rownames(cnt) %in% ext,]

#set path
path <- "input/Counts/Incidence/"

#make list of the count dataframes
cnt.list <- list(cnt, cnt.ab, cnt.fe, cnt.v)

filez <- c("All", "Prokaryote", "Eukaryote", "Viral")

for (i in 1:length(cnt.list)){
  nombre <- paste(filez[i], "Incidence_PhyloCollapsed", "Nov22.tsv", sep = "_")
  write.table(cnt.list[i], paste(path, nombre, sep = ""),
              sep = "\t", row.names = T, col.names = T, quote = F)
}

###SessionLog####
dir.create("SessionLogs/")
writeLines(capture.output(sessionInfo()),
           "SessionLogs/ReadFiltCount_Dec22.txt")
