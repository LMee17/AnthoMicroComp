#13th December
#Read in czid.org taxon reports, filter, produce count tables

##Libraries####
library(dplyr)
library(stringr)

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
euk <- read.csv("input/Phylo_Misc/Eukaryota_Classifications_Nov22.csv")
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
  data.raw[[i]]$nt_e_value <- gsub("10^0", 1, data.raw[[i]]$nt_e_value)
  data.raw[[i]]$nt_e_value <- gsub("10^0.0", 1, data.raw[[i]]$nt_e_value)
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
  data.raw[[i]]$nr_e_value <- gsub("10^0", 1, data.raw[[i]]$nr_e_value)
  data.raw[[i]]$nr_e_value <- gsub("10^0.0", 1, data.raw[[i]]$nr_e_value)
}




data <- bind_rows(data.raw)
testR <- data$nr_e_value[1:20]

testR <- gsub("10\\^", "", testR)
testR2 <- as.numeric(testR)

testT <- data$nt_e_value[690:700]
pows <- strsplit(testT, ".", fixed = T)
testT <- sapply(pows, "[", 1)
testT <- gsub("0\\^", "e", testT)
str(testT)
min(testT, na.rm = T)

as.numeric(testT)
testone <- "10^-27"
str_replace(testone, "0[^]", "e")

for(i in 1:length(data.raw)){
  et <- unique(grepl("e", data.raw[[i]]$nt_e_value))
  if(length(et) >1){
    print(paste(i, "e present in NT"))
  }
  er <- unique(grepl("e", data.raw[[i]]$nr_e_value))
  if(length(er) >1){
    print(paste(i, "e present in NR"))
  } else {
    print(paste(i, "not present in NR"))
  }
}

csvs[11:12]

met[met$Sample.ID == "SRR072892" | met$Sample.ID == "SRR072891",]

data.raw[[1]][8,]
str(data$nt_e_value)
10^189
10^0
data <- bind_rows(data.raw)
data[grepl("-", data$nr_e_value),]
data[(grepl("e", data$nr_e_value)),]
tail(data[!grepl("-", data$nt_e_value),])
head(data[grepl("-", data$nt_e_value),])
csvs[1]
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
    #change genus_taxon_id to 2767877
    data.raw[[i]]$genus_tax_id[data.raw[[i]]$name == "Apilactobacillus bombintestini"] <- 2767877
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

check <- data.raw
data.raw <- check

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
    #tax level = 2, tax/genus_tax_level = placeholder 5
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
    #nt counts, rpm, contigs, contig_r tare all sums of each of these
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

check2 <- data.raw

#now for the remaining lactobacilli
for (i in 1:10){#length(data.raw)
  print(i)
  #check if this frame has been affected by firm5 changes
  if(any(data.raw[[i]] == "Lactobacillus: Firm-5", na.rm = T)){
    #extract remaining lactobacilli
    rem <- data.raw[[i]][data.raw[[i]]$genus_tax_id == 1578 &
                           !data.raw[[i]]$name %in% firm5 &
                           data.raw[[i]]$tax_level == 1,]
    #get the location of the lactobacillus genus entry
    loc <- rownames(data.raw[[i]][data.raw[[i]]$name == "Lactobacillus",])
    #agg score is the highest agg score
    data.raw[[i]][loc, 8] <- max(rem$agg_score)
    #for evalues, take the minimum evalues
    print(data.raw[[i]][loc,c(17,29)])
    data.raw[[i]][[loc,17]] <- getMin(rem$nt_e_value)
    data.raw[[i]][[loc,29]] <- getMin(rem$nr_e_value)
    print(data.raw[[i]][loc,c(17,29)])
  }
}

pre <- check2[[6]]
t <- pre[pre$genus_tax_id == 1578,]
min(t$nr_e_value, rm.na = T)
str(t$nr_e_value)
str(check[[6]]$nr_e_value)
str(pre)
#for evalues, take the minimum evalues
data.raw[[i]][[loc,17]] <- getMin(rem$nt_e_value)
data.raw[[i]][[loc,29]] <- getMin(rem$nr_e_value)
#for nt/nr_rpm, nt/nr_counts, contigs, r_contig
#the genus value is the sum of all the species
#within that genus thus
data.raw[[i]][loc, 11] <- sum(rem$nt_rpm, na.rm = T)
data.raw[[i]][loc, 12] <- sum(rem$nt_count, na.rm = T)
data.raw[[i]][loc, 13] <- sum(rem$nt_contigs, na.rm = T)
data.raw[[i]][loc, 14] <- sum(rem$nt_contig_r, na.rm = T)
data.raw[[i]][loc, 23] <- sum(rem$nr_rpm, na.rm = T)
data.raw[[i]][loc, 24] <- sum(rem$nr_count, na.rm = T)
data.raw[[i]][loc, 25] <- sum(rem$nr_contigs, na.rm = T)
data.raw[[i]][loc, 26] <- sum(rem$nr_contig_r, na.rm = T)
print(data.raw[[i]][loc,])


prework <- check[[6]]
postwork <- data.raw[[6]]
tail(postwork)
rem1 <- prework[prework$genus_tax_id == 1578 &
                  !prework$name %in% firm5 &
                  prework$tax_level == 1,]

rem2 <- postwork[postwork$genus_tax_id == 1578 &
                   !postwork$name %in% firm5 &
                   postwork$tax_level == 1,]
prework[prework$name == "Lactobacillus: Firm-5",]
postwork[postwork$name == "Lactobacillus: Firm-5",]

#keep only hits at genus level (tax_level = 2)
#data.raw[[i]] <- subset(data.raw[[i]], tax_level == 2)
data <- bind_rows(data.raw)
data[order(data$nt_e_value, decreasing = FALSE),]

data %>%
  #select(nt_e_value) %>%
  arrange(desc(nt_e_value)) %>%
  head()

###starting with the non-viral filter#####
#prepare list ready to populated
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
#maybe remove negative aggregates ???
head(filt.df.nv)
filt.df.nv2 <- filt.df.nv[filt.df.nv$agg_score > 0,]
nrow(filt.df.nv[filt.df.nv$name=="Apilactobacillus",])
nrow(filt.df.nv2[filt.df.nv2$name=="Apilactobacillus",])


