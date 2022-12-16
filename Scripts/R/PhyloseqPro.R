#12th December 2022
#PhyloSeq Analysis; Prokaryote

library(phyloseq)
library(ggplot2)

set.seed(1549)

dir.create("output/Plots/PhyloSeq")

names(cnt)[colSums(cnt) > 7500000]
met[met$Sample.ID == "SRR5117446",]

##Functions####
#To filter data so that: 
# keep microbes that appear in at least 5% of samples
# keep samples that have 3 or more microbes in the community
# remove rows and columns that summed = 0
filter.cnt <- function(cntmat){
  #remove DNA samples
  dna <- met$Sample.ID[met$NucleotideType == "DNA"]
  filt <- cnt[, ! names(cnt) %in% dna]
  #remove any rows that equal 0 (first round)
  filt <- filt[!rowSums(filt) == 0,]
  #convert cntmat to incidence for ease of summing etc
  inc <- filt
  inc[inc > 0] <- 1
  #get microbial IDs and sample IDs of rows and columns that are to be removed
  #store what 5% of samples are
  five <- ncol(filt) * 0.05
  filt <- filt[rownames(inc)[rowSums(inc) > five],]
  #store all sample IDs that have less than 3 microbials present
  lt3 <- apply(filt, 2, function(x) sum(x > 0)) < 3
  filt <- filt[, !names(filt) %in% names(lt3[lt3 == TRUE])]
  #remove any columns that == 0, then double check for any empty rows
  filt <- filt[, ! colSums(filt) == 0]
  filt <- filt[! rowSums(filt) == 0, ]
  return(filt)
}


##Load Metadata####
#sample metadata
met <- read.table("input/Metadata/MetaData_Edit_Nov22.tsv", 
                  sep = "\t", header = T)
#microbial metadata (species level)
microkey <- read.table("input/Metadata/MicrobialSpeciesKey.tsv", 
                       sep = "\t", header = T)
#microbial metadata (genus and family level)
genkey <- read.table("input/Metadata/GenusFamilySpeciesKey_Nov22.tsv", 
                     sep = "\t", header = T)

#taxonomy file
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22_v2.tsv",
                  sep = "\t", header = T)

##Make PhyloSeq Object####
###Counts####
#load relative abundance count table for prokaryotes
cnt <- read.table("input/Counts/RawReads/Prokaryote_rawCounts_PhyloCollapsed_Nov22.tsv")

names(cnt)[colSums(cnt) > 7500000]

#apply filtering
cnt.filt <- filter.cnt(cnt)

#make matrix
cnt.mat <- as.matrix(cnt.filt)

###Taxonomy####
#I only need genus information
tax <- tax[,3:9]

#extract the genera that is present in the count table
gentab <- data.frame(Genus = unique(genkey$PhyloHit[genkey$GenHitID %in% rownames(cnt.mat)]),
                     GenID = rownames(cnt.mat))

#extract these from tax and see if there are any mismatches
tax2 <- unique(tax[tax$genus %in% gentab$Genus,])

#merge and see what's missing
taxmat <- merge(gentab, tax2, by.x = "Genus", by.y = "genus", all.x = T)
taxmat[is.na(taxmat$family),]
#Massilia, Ochrobactrum and Orbus were in the tax file but as genus only entries
#meaning the genus was in the species column.
#massilia
taxmat[64,3:8] <- c("Oxalobacteraceae", "Burkholderiales", "Betaproteobacteria", "Proteobacteria", NA, "Bacteria")
#ochrobactrum
taxmat[78,3:8] <- c("Brucellaceae", "Hyphomicrobiales", "Alphaproteobacteria", "Proteobacteria", NA, "Bacteria")
#orbus
taxmat[79,3:8] <- c("Orbaceae", "Orbales", "Gammaproteobacteria", "Proteobacteria", NA, "Bacteria")
#finally, the unspecified bacteria
taxmat[109,8] <- "Bacteria"

#make ids rownames and convert to matrix
rownames(taxmat) <- taxmat$GenID
taxmat$GenID <- NULL
taxmat <- as.matrix(taxmat)

###Make Simple PhyoSeq Object####
cnts <- otu_table(cnt.mat, taxa_are_rows = T)
tax <- tax_table(taxmat)

phy <- phyloseq(cnts, tax)

#plot absolute abundances ?
plot_bar(phy, fill = "phylum")
colnames(taxmat)


cnt2 <- read.table("input/Counts/RelativeAbundance/Prokaryote_RA_filtered_Dec22.tsv")