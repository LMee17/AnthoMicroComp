#Lactobacillus Firm-5: Compilation
#11 November 2022

library(edgeR)
library(ggplot)
library(reshape)

#Sample metadata
met <- read.table("input/SRA/MetaData_Edit_Oct22.tsv",
                  header = T, sep ="\t", quote = "")

#Microbial metadata
genkey <- read.table("output/Keys/GenusFamilySpeciesKey.tsv",
                     sep = "\t", header = T)
#Taxonomy
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22.tsv",
                  header = T, sep = "\t", quote = "")

#core phylogroup lactobacillus firm -5
firm5 <- c("Lactobacillus apis", "Lactobacillus melliventris", "Lactobacillus kimbladii",
           "Lactobacillus kullabergensis", "Lactobacillus panisapium", 
           "Lactobacillus bombicola", "Lactobacillus helsingborgensis")

#load in species level count table
spec.cnt <- read.table("output/Counts/RawReads/All_rawCounts_Nov22.tsv")

#check the length of the IDs to make a new one
length(unique(genkey$GenHitID))
#apply this to the firm5 group
genkey$GenHitID[genkey$MicroHit %in% firm5] <- "GRH761"
genkey[genkey$MicroHit %in% firm5,]
#change phylogeny
genkey$PhyloHit[genkey$MicroHit %in% firm5] <- "Lactobacillus: Firm-5"

#add to taxonomy
tax[tax$genus == "Lactobacillus",]
n <- nrow(tax)+1
tax[n,] <- c("Lactobacillus: Firm-5", "", "Lactobacillus: Firm-5", "Lactobacillaceae",
             "Lactobacillales", "Bacilli", "Firmicutes", NA, "Bacteria")

tail(tax)

#Add the microhit ids to the genkey so I can recount the single species count table
microkey <- read.table("output/Keys/MicrobialSpeciesKey.tsv",
                       sep = "\t", header = T)
for (i in 1:nrow(genkey)){
  genkey$MicroHitID[i] <- microkey$HitID[microkey$MicroHit == genkey$MicroHit[i]]
}
#write up
write.table(genkey, "output/Keys/GenusFamilySpeciesKey_Nov22.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(tax, "input/Phylo_Misc/rankedlin_Edit_Nov22_v2.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)


#Now make a new, phylo collapsed count table ....
#compile samples
samples <- names(spec.cnt)
#compile microbial genus hit ids
phyloz <- unique(genkey$GenHitID)
#make a matrix of sufficient size
cnt <- matrix(ncol = length(samples), nrow = length(phyloz))
#run through each of the samples
for (i in 1:length(samples)){
  #print to screen so progress can be checked
  print(paste(i, "/317: ", samples[i], sep = ""))
  #make a vector for the sample read to populate with counts per microbial genus group
  x <- vector(length = length(phyloz))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    #extract the microbial IDs within this microbial genus grouping
    y <- unique(genkey$MicroHitID[genkey$GenHitID == phyloz[j]])
    #pull out the rows of the count table that match these microbial ids
    t <- spec.cnt[rownames(spec.cnt) %in% y, i]
    #combine multiple species from a phylogenetic classification into one
    if (length(t) > 1){
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
    #add sample's counts to count matrix
    cnt[,i] <- x
  }
}
#covert to a dataframe
cnt <- as.data.frame(cnt)
#add row and colomn names
rownames(cnt) <- phyloz
names(cnt) <- samples
#write up
write.table(cnt, "output/Counts/RawReads/All_rawCounts_CorePhylo_PhyloCollapsed_Nov22.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

cnt
cnt2 <- cnt
cnt2[cnt2 > 0] <- 1
cnt2
write.table(cnt2, "output/Counts/Incidence/All_Incidence_CorePhylo_PhyloCollapsed_Nov22.tsv")