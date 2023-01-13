#13th December
#Read in czid.org taxon reports, filter, produce count tables

##Libraries####
library(dplyr)
library(stringr)

dir.create("SessionLogs/")

##Functions####
#get maximum value from a vector of numbers that may or may not contain NAs (or only NAs)
getMax <- function(counts){
  if (length(unique(counts)) > 1){
    out <- max(counts, na.rm = T)
  } else {
    out <- max(counts)
  } 
  return(out)
  print(out)
}

#get minimum value from a vector of numbers that may or may not contain NAs (or only NAs)
getMin <- function(counts){
  if (length(unique(counts)) > 1){
    out <- min(counts, na.rm = T)
  } else {
    out <- min(counts)
  } 
  return(out)
  print(out)
}


##load necessary metadata and taxonomic files####
#sample metadata 
met <- read.table("input/Metadata/SampleMetaData_Edit_RNAOnly_Dec22.tsv",
                  sep = "\t", header = T)
#eukaryotic tax classifications
euk <- read.csv("input/Phylo_Misc/Eukaryota_Classifications_Dec22.csv")

#taxonomy
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22_v2.tsv",
                  sep ="\t", header = T)
#outdated and updated taxonomy terms
#only including genera names or species I'm interested in, 
#i.e Apilactobacillus bombintestini
taxup <- read.table("input/Phylo_Misc/TaxoUpdate_2_Dec22.txt", sep = "\t", header = T)

##Step One: Read in Taxon Reports####
#list taxon report files
csvs <- list.files(path = "input/CZ_TaxonReports/", pattern = "*.csv")

#prepare a list ready to be populated with report files
data.raw <- vector(mode = "list", length = length(csvs))

#iterate through list of files, reading in each into the data list
#updating outdated taxa with the current version (or at least the version
#that matches my NCBI taxa file, downloaded Oct 2022)
#replace empty cells with NA
for (i in 1:length(csvs)){
  data.raw[[i]] <- read.csv(paste("input/CZ_TaxonReports/", csvs[i], sep = ""),
                            header = T)
  for (j in 1:nrow(data.raw[[i]])){
    if (data.raw[[i]]$name[j] %in% taxup$Old){
      data.raw[[i]]$name[j] <- taxup$Current[taxup$Old == data.raw[[i]]$name[j]]
    }
  }
  data.raw[[i]][data.raw[[i]] == ""] <- NA
}

#fix evalues
#nt evalues are in the format 10^-n, but has decimals in which seem to have
#stopped r from recognising them. I'm just going to remove the decimals
#and everything afterwards and that should fix NT evalue
#nr evalues on the other hand appear to be in the format 1e-183 etc,
#but are prefixed by 10^ which I think is making it a character string to r.
#Also need to replace all 10^0 or 10^0.0 to 1
for (i in 1:length(data.raw)){
  #nt evalues
  #replace 10^0/10^0.0 with 1
  data.raw[[i]]$nt_e_value[data.raw[[i]]$nt_e_value == "10^0"] <- 1
  data.raw[[i]]$nt_e_value[data.raw[[i]]$nt_e_value == "10^0.0"] <- 1
  #split the entries by "."
  et <- strsplit(data.raw[[i]]$nt_e_value, ".", fixed = T)
  #overwrite with just the first part of the split values
  data.raw[[i]]$nt_e_value <- sapply(et, "[", 1)
  #replace 0^ with an e, converting it to e-notation
  data.raw[[i]]$nt_e_value <- gsub("0\\^", "e", data.raw[[i]]$nt_e_value)
  #make numeric
  data.raw[[i]]$nt_e_value <- as.numeric(data.raw[[i]]$nt_e_value)
  #nr evalues
  #replace 10^0/10^0.0 with 1
  data.raw[[i]]$nr_e_value[data.raw[[i]]$nr_e_value == "10^0"] <- 1
  data.raw[[i]]$nr_e_value[data.raw[[i]]$nr_e_value == "10^0.0"] <- 1
  #remove 10^
  data.raw[[i]]$nr_e_value <- gsub("10\\^", "", data.raw[[i]]$nr_e_value)
  #convert to numeric
  data.raw[[i]]$nr_e_value <- as.numeric(data.raw[[i]]$nr_e_value)
}

##Step Two: The Lactobacillus Consideration####
#two things need to be checked for corrected
#one: if Apilactobacillus bombintestini is present, the Apilactobacillus totals need to be recalculated
#two: Lactobacillus needs to split into Lactobacillus and Lactobacillus Firm-5
#this would require an entire new entry for Firm-5 with appropriate entries and
#either replacement or editing of the previous Lactobacillus entry

#sort through Apilactobacillus hits first
for (i in 1:length(data.raw)){
  #see if apilactobacillus bombestini was found in this sample
  if (any(data.raw[[i]] == "Apilactobacillus bombintestini", na.rm = T) == T){
    #change Tax_id to 2767877
    data.raw[[i]]$tax_id[data.raw[[i]]$name == "Apilactobacillus bombintestini"] <- 2767877
    #check to see if there's already a Apilactobacillus entry
    if(any(data.raw[[i]] == "Apilactobacillus", na.rm = T) == T){
      #log the row (location) of the Apilactobacillus entry
      loc <- rownames(data.raw[[i]][data.raw[[i]]$name == "Apilactobacillus",])
      #extract all Apilactobacillus species-level entries
      api <- data.raw[[i]][grepl("Apilactobacillus", data.raw[[i]]$name) 
                           & data.raw[[i]]$tax_level == 1,]
      #agg score is the highest agg score
      data.raw[[i]][loc, 8] <- max(api$agg_score)
      #for evalues, take the minimum evalues
      data.raw[[i]][[loc,17]] <- getMin(api$nt_e_value)
      data.raw[[i]][[loc,29]] <- getMin(api$nr_e_value)
      #for nt/nr_rpm, nt/nr_counts, contigs, r_contig
      #the genus value is the sum of all the species
      #within that genus thus
      data.raw[[i]][loc, 11] <- sum(api$nt_rpm, na.rm = T)
      data.raw[[i]][loc, 12] <- sum(api$nt_count, na.rm = T)
      data.raw[[i]][loc, 13] <- sum(api$nt_contigs, na.rm = T)
      data.raw[[i]][loc, 14] <- sum(api$nt_contig_r, na.rm = T)
      data.raw[[i]][loc, 23] <- sum(api$nr_rpm, na.rm = T)
      data.raw[[i]][loc, 24] <- sum(api$nr_count, na.rm = T)
      data.raw[[i]][loc, 25] <- sum(api$nr_contigs, na.rm = T)
      data.raw[[i]][loc, 26] <- sum(api$nr_contig_r, na.rm = T)
    } else {
      #when A. bombintestini is present, but Apilactobacillus is not, 
      #I can just duplicate the A. bombintestini metrics as a new genus entry
      #changing the name and taxon level
      lastrow <- nrow(data.raw[[i]]) + 1
      data.raw[[i]][lastrow,] <- data.raw[[i]][data.raw[[i]]$name == "Apilactobacillus bombintestini",]
      #upgrade taxa level to genus
      #change name to Apilactobacillus
      data.raw[[i]][lastrow,2] <- 2
      data.raw[[i]][lastrow,4] <- "Apilactobacillus"
    }
  }
}

#split lactobacillus genus
#list the firm-5 species
firm5 <- c("Lactobacillus apis", "Lactobacillus melliventris", "Lactobacillus kimbladii",
           "Lactobacillus kullabergensis", "Lactobacillus panisapium", 
           "Lactobacillus bombicola", "Lactobacillus helsingborgensis")

#start with the firm-5 : making a new genus grouping
#then move on to fixing the lactobacillus that are left
for (i in 1:length(data.raw)){
  #extract firm5, if present
  five <- data.raw[[i]][data.raw[[i]]$name %in% firm5,]
  if(nrow(five) > 0){
    ##genus tax id will be a placeholder 5 (check on NCBI and appears to not be in use)
    for (r in 1:nrow(five)){
      #log the rownames of where these microbes appear in the data
      loc <- rownames(five[r,])
      #replace the taxid there with my placeholder
      data.raw[[i]][loc,3] <- 5
    }
    #prepare to make new row for firm 5 "genus" entry
    lastrow <- nrow(data.raw[[i]]) + 1
    #tax level = 2, tax/Taxa_tax_level = placeholder 5
    #name = Lactobacillus: Firm-5, category = bacteria, is_phage = FALSE
    data.raw[[i]][lastrow,1] <- 5
    data.raw[[i]][lastrow,2] <- 2
    data.raw[[i]][lastrow,3] <- 5
    data.raw[[i]][lastrow,4] <- "Lactobacillus: Firm-5"
    data.raw[[i]][lastrow,5] <- NA #common name
    data.raw[[i]][lastrow,6] <- "bacteria"
    data.raw[[i]][lastrow,7] <- "false"
    #agg score is max agg score of the lactobacilli
    data.raw[[i]][lastrow,8] <- max(five$agg_score)
    #same for z, for my purposes ( I can't compute these so this is a quick workaround)
    data.raw[[i]][lastrow,9] <- getMax(five$max_z_score)
    data.raw[[i]][lastrow,10] <- getMax(five$nt_z_score)
    data.raw[[i]][lastrow,22] <- getMax(five$nr_z_score)
    #nt counts, rpm, contigs, contig_r are all sums of each of these
    data.raw[[i]][lastrow,11] <- sum(five$nt_rpm, na.rm = T)
    data.raw[[i]][lastrow,12] <- sum(five$nt_count, na.rm = T)
    data.raw[[i]][lastrow,13] <- sum(five$nt_contigs, na.rm = T)
    data.raw[[i]][lastrow,14] <- sum(five$nt_contig_r, na.rm = T)
    data.raw[[i]][lastrow,23] <- sum(five$nr_rpm, na.rm = T)
    data.raw[[i]][lastrow,24] <- sum(five$nr_count, na.rm = T)
    data.raw[[i]][lastrow,25] <- sum(five$nr_contigs, na.rm = T)
    data.raw[[i]][lastrow,26] <- sum(five$nr_contig_r, na.rm = T)
    #I'll just be taking max percent ID and alignment length for ease
    data.raw[[i]][lastrow,15] <- getMax(five$nt_percent_identity)
    data.raw[[i]][lastrow,16] <- getMax(five$nt_alignment_length)
    data.raw[[i]][lastrow,27] <- getMax(five$nr_percent_identity)
    data.raw[[i]][lastrow,28] <- getMax(five$nr_alignment_length)
    #take the minimum for evalue
    data.raw[[i]][lastrow,17] <- getMin(five$nt_e_value)
    data.raw[[i]][lastrow,29] <- getMin(five$nr_e_value)
  }
}

#now for the remaining lactobacilli
for (i in 1:length(data.raw)){
  #check if this frame has been affected by firm5 changes
  if(any(data.raw[[i]] == "Lactobacillus: Firm-5", na.rm = T)){
    #extract remaining lactobacilli
    rem <- data.raw[[i]][data.raw[[i]]$Taxa_tax_id == 1578 &
                           !data.raw[[i]]$name %in% firm5 &
                           data.raw[[i]]$tax_level == 1,]
    #is there any remaininng lactobacillus to consider? 
    if (nrow(rem) == 0){
      #remove the lactobacillus genus entry
      data.raw[[i]] <- subset(data.raw[[i]], name != "Lactobacillus")
    } else {
      #get the location of the lactobacillus genus entry
      loc <- rownames(data.raw[[i]][data.raw[[i]]$name == "Lactobacillus",])
      #agg score is the highest agg score
      data.raw[[i]][loc, 8] <- max(rem$agg_score)
      #for evalues, take the minimum evalues
      data.raw[[i]][[loc,17]] <- getMin(rem$nt_e_value)
      data.raw[[i]][[loc,29]] <- getMin(rem$nr_e_value)
      #again, max alignment length and percent ID
      data.raw[[i]][loc,15] <- getMax(rem$nt_percent_identity)
      data.raw[[i]][loc,16] <- getMax(rem$nt_alignment_length)
      data.raw[[i]][loc,27] <- getMax(rem$nr_percent_identity)
      data.raw[[i]][loc,28] <- getMax(rem$nr_alignment_length)
      #I'm going to leave z scores as they are
      #nt counts, rpm, contigs, contig_r are all sums of each of these
      data.raw[[i]][loc,11] <- sum(rem$nt_rpm, na.rm = T)
      data.raw[[i]][loc,12] <- sum(rem$nt_count, na.rm = T)
      data.raw[[i]][loc,13] <- sum(rem$nt_contigs, na.rm = T)
      data.raw[[i]][loc,14] <- sum(rem$nt_contig_r, na.rm = T)
      data.raw[[i]][loc,23] <- sum(rem$nr_rpm, na.rm = T)
      data.raw[[i]][loc,24] <- sum(rem$nr_count, na.rm = T)
      data.raw[[i]][loc,25] <- sum(rem$nr_contigs, na.rm = T)
      data.raw[[i]][loc,26] <- sum(rem$nr_contig_r, na.rm = T)
    }
  }
}

##Step Three: Logging Microbes of Interest####
#read in list of microbes previously identified and characterised in various bee species
#(essentially corbiculates)
moi <- read.table("input/Metadata/MicrobesOfInterest.txt", sep = "\t")

#prepare dataframe to be populated
moi.df <- data.frame()

#loop through each taxon report, tracking each time these m.o.i appear
#set counter for rows to populate
tic <- 1
for (i in 1:length(data.raw)){
  if (any(data.raw[[i]]$name %in% moi$V1)){
    #extract microbes of interest from dataframe
    tar <- unique(data.raw[[i]]$name[data.raw[[i]]$name %in% moi$V1])
    #get sample name 
    samp <- strsplit(csvs[i], "_2", fixed = T)[[1]][1]
    #get sample species
    spec <- met$Species[met$Sample.ID == samp]
    print(paste(samp,spec))
    for (j in 1:length(tar)){
      #get rowname of this hit
      loc <- rownames(data.raw[[i]][data.raw[[i]]$name == tar[j],])
      #start to put information into dataframe
      moi.df[tic, 1:3] <- c(tar[j], samp, spec)
      #agg score
      moi.df[tic, 4] <- data.raw[[i]]$agg_score[data.raw[[i]]$name == tar[j]]
      #count and contig count (nt)
      moi.df[tic, 5:6] <- data.raw[[i]][loc, 12:13]
      #nt_percent and alignment length and evalue
      moi.df[tic, 7:9] <- data.raw[[i]][loc, 15:17]
      #count and contig (nr)
      moi.df[tic, 10:11] <- data.raw[[i]][loc, 24:25]
      #nr_percent and alignment length and evalue
      moi.df[tic, 12:14] <- data.raw[[i]][loc, 27:29]
      tic <- tic + 1
    }
  }
}
#fix dataframe column names
names(moi.df)[c(1:4)] <- c("Microbe", 
                           "Sample",
                           "Sample_Species",
                           "AggregateScore")

#write up
write.table(moi.df, "output/MicrobesOfInterest.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")


##Step Four: Filtering#####
###starting with the non-viral filter#####
#filter parameters for nonvirals:
#nt_count > 1, nt_count >= 10, NT L >= 50, evalue > 1e-6, agg > 0, 
#max_zscore > 0, %id >= 90%, and tax level 2

#check the number of reads per sample
#Any with total < 1000 will be removed
for(i in 1:length(data.raw)){
  s <- sum(data.raw[[i]]$nt_count, na.rm = T)
  if (s < 1000){
    print(paste(csvs[i], i))
  }
}

#dataframes 11 and 12 which correspond to samples SRR072891 and SRR072892 have less
#than 1000 (in fact less than 50) reads each and will be removed before going further

#prepare list ready to populated
data.filt.nv <- data.raw

#run through dataframes, filtering  and adding information as required
for (i in 1:length(data.filt.nv)){
  #get sample ID from the filename
  samp <- paste(strsplit(csvs[i], "_2[0-9]")[[1]][1])
  samp <- str_remove(samp, "_[1-2]")
  #add sample ID to dataframe for ease of identification when these are combined.
  data.filt.nv[[i]]$SampleID <- paste(samp)
  #get species and genus information from the metadata object
  data.filt.nv[[i]]$Species <- unique(met$Species[met$Sample.ID == samp])
  data.filt.nv[[i]]$Genus <- unique(met$Genus[met$Sample.ID == samp])
  #filter the results by the parameters outlined above
  data.filt.nv[[i]] <- subset(data.filt.nv[[i]], nt_count > 1 &
                                nt_rpm >= 5 &
                                nt_alignment_length >= 50 &
                                nt_e_value <= 1e-6 &
                                agg_score > 0 &
                                max_z_score > 0 &
                                nt_percent_identity >= 90 &
                                tax_level == 2)
  #only keep necessary columns
  data.filt.nv[[i]] <- data.filt.nv[[i]][,c(4,6,8:12,15:17,23:24,27:29,35:37)]
}

#combine into a single dataframe
filt.df.nv <- bind_rows(data.filt.nv)
#remove viruses
filt.df.nv <- subset(filt.df.nv, category != "viruses")

#for some reason the dataframe has amassed rows of all NAs. 
#These are causing issues and need to be removed
filt.df.nv <- filt.df.nv[grepl("^NA", rownames(filt.df.nv))==F,]

#add genus information as a separate column
filt.df.nv$MicroPhylo <- filt.df.nv$name
#for those that don't have a specific genus
for (i in 1:nrow(filt.df.nv)){
  if (grepl("non-genus-specific", filt.df.nv$name[i]) == T){
    filt.df.nv$MicroPhylo[i] <- strsplit(filt.df.nv$name[i], " ")[[1]][5]
  }
}

#make lists of wanted eukaryota
wantedClass <- c("Fungus", "Trypanosome", "Unicellular")
tokeep <- euk$Genus[euk$Classification %in% wantedClass]

#keep only wanted taxa
filt.df.nv <- filt.df.nv[filt.df.nv$MicroPhylo %in% tokeep |
                           filt.df.nv$category == "bacteria" |
                           filt.df.nv$category == "archaea",]

###viral filter####
#prepare list ready to populated
data.filt.v <- data.raw

#run through dataframes, filtering  and adding information as required
#for viruses use the NR side of the taxon report tables
for (i in 1:length(data.filt.v)){
  #get sample ID from the filename
  samp <- paste(strsplit(csvs[i], "_2[0-9]")[[1]][1])
  samp <- str_remove(samp, "_[1-2]")
  #add sample ID to dataframe for ease of identification when these are combined.
  data.filt.v[[i]]$SampleID <- paste(samp)
  #get species and genus information from the metadata object
  data.filt.v[[i]]$Species <- unique(met$Species[met$Sample.ID == samp])
  data.filt.v[[i]]$Genus <- unique(met$Genus[met$Sample.ID == samp])
  #filter the results by the parameters outlined above
  data.filt.v[[i]] <- subset(data.filt.v[[i]], nr_count > 1 &
                               nr_rpm >= 5 &
                               nr_alignment_length >= 50 &
                               nr_e_value <= 1e-6 &
                               agg_score > 0 &
                               max_z_score > 0 &
                               nr_percent_identity >= 90 &
                               tax_level == 2)
  #only keep necessary columns
  data.filt.v[[i]] <- data.filt.v[[i]][,c(4,6,8:12,15:17,23:24,27:29,35:37)]
}

filt.df.v <- bind_rows(data.filt.v)
filt.df.v <- subset(filt.df.v, category == "viruses")

#add a microphylo column to line up with the nonviral table
filt.df.v$MicroPhylo <- filt.df.v$name
for (i in 1:nrow(filt.df.v)){
  if (grepl("non-genus-specific", filt.df.v$name[i]) == T){
    filt.df.v$MicroPhylo[i] <- strsplit(filt.df.v$name[i], " ")[[1]][5]
  }
}

#combine with non-viral filtered data
filt.df <- rbind(filt.df.nv, filt.df.v)

#And lastly, remove an unnecessary entry 
filt.df <- subset(filt.df, 
                  name != "all taxa with neither family nor genus classification")

#again, NAs have amassed
filt.df <- filt.df[grepl("^NA", rownames(filt.df.nv))==F,]

##Step Five: Making the Count Table ####
#this will require an assignment of arbitrary taxa ids in place of genus information
### Producing the Genus Hit Key ####
taxkey <- data.frame(Taxa = unique(filt.df$name))

#add classifications
for (i in 1:nrow(taxkey)){
  taxkey$Class[i] <- unique(filt.df$category[filt.df$name == taxkey$Taxa[i]])
}

#remove "non-genus-specific"
for (i in 1:nrow(taxkey)){
  if (grepl("non-genus-specific", taxkey$Taxa[i]) == T){
    taxkey$Taxa[i] <- strsplit(taxkey$Taxa[i], " ")[[1]][5]
  }
}
#are there now duplicates ? 
taxkey[duplicated(taxkey$Taxa),]

#no, so we can continue (and this makes it much simpler)
#populate taxa read hit ids (TRH)
taxkey$taxID <- paste0("TRH", 1:length(taxkey$Taxa), sep = "")

#save
write.table(taxkey, "input/Metadata/Taxa_HitKey_Dec22.tsv",
            sep = "\t", row.names = F, quote = F)

### Producing the Count Table ###
#check the number of reads per sample
#Any with total < 1000 will be removed
for(i in 1:length(data.raw)){
  s <- sum(data.raw[[i]]$nt_count, na.rm = T)
  if (s < 1000){
    print(paste(csvs[i], i))
  }
}

#dataframes 11 and 12 which correspond to samples SRR072891 and SRR072892 have less
#than 1000 (in fact less than 50) reads each and will be removed before going further
lt100 <- c("SRR072891", "SRR072892")
filt.df <- filt.df[!filt.df$SampleID %in% lt100,]

####Prokaryotes#####
#starting with the prokaryote read counts (nt)
pro <- filt.df[filt.df$category == "bacteria" |
                 filt.df$category == "archaea",]

proID <- taxkey$taxID[taxkey$Class == "bacteria" | 
                        taxkey$Class == "archaea"]

#get column headers ready
#keeping only samples that have microbes after filtering
samples <- unique(pro$SampleID)

#record length of samples
totSam <- length(samples)

#prepare matrix to be populated with counts
cnt <- matrix(ncol = length(samples), nrow = length(proID))

for (i in 1:length(samples)){
  print(paste(i, "/", totSam, ": ", samples[i], sep = ""))
  #iterate through the genera found across all samples, pulling out counts where they're
  #present, populating with a 0 where absent
  x <- vector(length = length(proID))
  y <- subset(filt.df, SampleID == samples[i])
  for (g in 1:length(proID)){
    micro <- taxkey$Taxa[taxkey$taxID == proID[g]]
    if (micro %in% y$MicroPhylo){
      x[g] <- y$nt_count[y$SampleID == samples[i] & y$MicroPhylo == micro]
    } else {
      x[g] <- 0
    }
  }
  #add to count matrix
  cnt[,i] <- x
}

#convert to dataframe for more amenable editing
cnt <- as.data.frame(cnt)
names(cnt) <- samples
rownames(cnt) <- proID

#write up
write.table(cnt, "input/Counts/Prokaryote_RawCounts.tsv",
            sep = "\t", quote = F,
            col.names = T, row.names = T)

####Eukaryotes#####
#I'll be doing this twice: once with all eukaryotes, once with fungus + other eukaryotes
#split
#there is an instance where a sample has two Candida hits with different taxa
#ids that lead back to the same taxonomy on NCBI that causes issues later on 
#so first I'll remove the extra candida entry, choosing the entry with the lowest counts
#and highest evalue
for (i in 1:length(samples)){
  y <- filt.df[filt.df$SampleID == samples[i] & filt.df$MicroPhylo == "Candida",]
  if (nrow(y) > 1){
    print(samples[i])
  }
}
#fix rownames so they can be used to pull rows out
rownames(filt.df) <- paste(1:nrow(filt.df))
#determine location of duplicates
filt.df[filt.df$name == "Candida" & filt.df$SampleID == "SRR14567222",]
#remove
filt.df <- filt.df[-3161,]

#continuing with the table
euks <- filt.df[filt.df$category == "eukaryota",]

euID <- taxkey$taxID[taxkey$Class == "eukaryota" &
                       taxkey$Taxa %in% filt.df$MicroPhylo]

#get column headers ready
#keeping only samples that have microbes after filtering
samples <- unique(euks$SampleID)

#record length of samples
totSam <- length(samples)

#prepare matrix to be populated with counts
cnt <- matrix(ncol = length(samples), nrow = length(euID))

for (i in 1:length(samples)){
  print(paste(i, "/", totSam, ": ", samples[i], sep = ""))
  #iterate through the genera found across all samples, pulling out counts where they're
  #present, populating with a 0 where absent
  x <- vector(length = length(euID))
  y <- subset(filt.df, SampleID == samples[i])
  for (g in 1:length(euID)){
    micro <- taxkey$Taxa[taxkey$taxID == euID[g]]
    if (micro %in% y$MicroPhylo){
      x[g] <- y$nt_count[y$SampleID == samples[i] & y$MicroPhylo == micro]
    } else {
      x[g] <- 0
    }
  }
  #add to count matrix
  cnt[,i] <- x
}


#convert to dataframe for more amenable editing
cnt <- as.data.frame(cnt)
names(cnt) <- samples
rownames(cnt) <- euID

#write up
write.table(cnt, "input/Counts/Eukaryote_RawCounts.tsv",
            sep = "\t", quote = F,
            col.names = T, row.names = T)

#Fungus only
fungus <- euk$Genus[euk$Classification == "Fungus"]
#continuing with the table
fun <- filt.df[filt.df$MicroPhylo %in% fungus,]

funID <- taxkey$taxID[taxkey$Class == "eukaryota" &
                        taxkey$Taxa %in% fungus]

#get column headers ready
#keeping only samples that have microbes after filtering
samples <- unique(fun$SampleID)

#record length of samples
totSam <- length(samples)

#prepare matrix to be populated with counts
cnt <- matrix(ncol = length(samples), nrow = length(funID))

for (i in 1:length(samples)){
  print(paste(i, "/", totSam, ": ", samples[i], sep = ""))
  #iterate through the genera found across all samples, pulling out counts where they're
  #present, populating with a 0 where absent
  x <- vector(length = length(funID))
  y <- subset(filt.df, SampleID == samples[i])
  for (g in 1:length(funID)){
    micro <- taxkey$Taxa[taxkey$taxID == funID[g]]
    if (micro %in% y$MicroPhylo){
      x[g] <- y$nt_count[y$SampleID == samples[i] & y$MicroPhylo == micro]
    } else {
      x[g] <- 0
    }
  }
  #add to count matrix
  cnt[,i] <- x
}


#convert to dataframe for more amenable editing
cnt <- as.data.frame(cnt)
names(cnt) <- samples
rownames(cnt) <- funID

#write up
write.table(cnt, "input/Counts/Fungus_RawCounts.tsv",
            sep = "\t", quote = F,
            col.names = T, row.names = T)

#Other eukaryotes
other <- unique(filt.df$MicroPhylo[filt.df$category == "eukaryota" &
                                     !filt.df$MicroPhylo %in% fungus])
#continuing with the table
other.df <- filt.df[filt.df$MicroPhylo %in% other,]

othID <- taxkey$taxID[taxkey$Taxa %in% other]

#get column headers ready
#keeping only samples that have microbes after filtering
samples <- unique(other.df$SampleID)

#record length of samples
totSam <- length(samples)

#prepare matrix to be populated with counts
cnt <- matrix(ncol = length(samples), nrow = length(othID))

for (i in 1:length(samples)){
  print(paste(i, "/", totSam, ": ", samples[i], sep = ""))
  #iterate through the genera found across all samples, pulling out counts where they're
  #present, populating with a 0 where absent
  x <- vector(length = length(othID))
  y <- subset(filt.df, SampleID == samples[i])
  for (g in 1:length(othID)){
    micro <- taxkey$Taxa[taxkey$taxID == othID[g]]
    if (micro %in% y$MicroPhylo){
      x[g] <- y$nt_count[y$SampleID == samples[i] & y$MicroPhylo == micro]
    } else {
      x[g] <- 0
    }
  }
  #add to count matrix
  cnt[,i] <- x
}


#convert to dataframe for more amenable editing
cnt <- as.data.frame(cnt)
names(cnt) <- samples
rownames(cnt) <- othID

#write up
write.table(cnt, "input/Counts/OtherEukaryote_RawCounts.tsv",
            sep = "\t", quote = F,
            col.names = T, row.names = T)

####Viruses#####
#virus read counts are from the nr database
vir <- filt.df[filt.df$category == "viruses",]

virID <- taxkey$taxID[taxkey$Class == "viruses"]

#get column headers ready
#keeping only samples that have microbes after filtering
samples <- unique(vir$SampleID)

#record length of samples
totSam <- length(samples)

#prepare matrix to be populated with counts
cnt <- matrix(ncol = length(samples), nrow = length(virID))

for (i in 1:length(samples)){
  print(paste(i, "/", totSam, ": ", samples[i], sep = ""))
  #iterate through the genera found across all samples, pulling out counts where they're
  #present, populating with a 0 where absent
  x <- vector(length = length(virID))
  y <- subset(filt.df, SampleID == samples[i])
  for (g in 1:length(virID)){
    micro <- taxkey$Taxa[taxkey$taxID == virID[g]]
    if (micro %in% y$MicroPhylo){
      x[g] <- y$nr_count[y$SampleID == samples[i] & y$MicroPhylo == micro]
    } else {
      x[g] <- 0
    }
  }
  #add to count matrix
  cnt[,i] <- x
}

#convert to dataframe for more amenable editing
cnt <- as.data.frame(cnt)
names(cnt) <- samples
rownames(cnt) <- virID

#write up
write.table(cnt, "input/Counts/Viral_RawCounts.tsv",
            sep = "\t", quote = F,
            col.names = T, row.names = T)

###SessionLog ####
writeLines(capture.output(sessionInfo()),
           "SessionLogs/ReadFiltCount_Dec22.txt")
