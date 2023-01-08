#20th December
#Make taxonomy file using NCBI DUMP file and genera detected post filtering of CZID
#taxon reports

#outside of this script the bash script SplitTax.sh was ran on the file 
#rankedlineage_181022.tsv to split it into smaller, manageable files that can
#be processed through r without breaking it (this is a very large file)

#list split files
taxfiles <- list.files(path = "input/Phylo_Misc/", pattern = "*.tsv")
#remove the original file / previous taxonomy versions
#double check before running
taxfiles <- taxfiles[c(-1,-2)]

#loop through the files, editing them as I go
#the edited files will be written back into the PhyloSeq directory 
for (i in 1:length(taxfiles)){
  x <- read.table(paste("input/Phylo_Misc/", taxfiles[i], sep = ""), 
                  comment.char = "", header = F, sep = "\t", quote = "")
  #remove header from the first file
  if(i == 1){
    x <- x[-1,]
  }
  #only keep lines that have either the species / common name entry that matches
  #one of the microbial taxa hits I already have
  x2 <- x[x$V2 %in% genkey$Genus ,]
  write.table(x2, paste("input/Phylo_Misc/", taxfiles[i], sep = ""),
              col.names = F, row.names = F, sep = "\t", quote = F)
}

#back in the directory space I catted these files together into a new taxonomy file
#splitfiles were then removed
#read this table in the final touches
tax <- read.table("input/Phylo_Misc/rankedlin_Dec22.tsv",
                  sep = "\t", header = F)
#remove the taxid column
tax <- tax[,-1]
#add the taxonomy ranks
names(tax) <- c("tax_name",	"species", "genus", 
                "family", "order", "class", "phylum",	
                "kingdom", "superkingdom")

#add the lactobacillus firm-5 entry
lastRow <-  nrow(tax) +1
tax[lastRow,] <- c("Lactobacillus: Firm-5", NA, "", "Lactobacillaceae",
                   "Lactobacillales", "Bacilli", "Firmicutes", NA, "Bacteria")

#there will be duplicated entries I need to address
tax[duplicated(tax$tax_name),]

#there should be no plants or metazoans in there
#before subsetting, I need to replace NA's with empty spaces as otherwise
#dplyr filters it out
tax[is.na(tax)] <- ""
tax <- tax %>%
  subset(kingdom != "Metazoa") %>%
  subset(class != "Magnoliopsida")

#remove duplicated entries (Leishmania has two entries, one with a genus entry)
rownames(tax) <- paste0(1:nrow(tax))
tax[tax$tax_name == "Leishmania",]
tax <- tax[-459,]

#remove duplicated Graphium entries
#the only fungal Graphium in my data is from the Graphiaceae family
rownames(tax) <- paste0(1:nrow(tax))
tax[tax$tax_name == "Graphium",]
tax <- tax[c(-557,-580),]

#there are two Morganella entries
#only the bacterial hit remains in my data after filtering
rownames(tax) <- paste0(1:nrow(tax))
tax[tax$tax_name == "Morganella",]
tax <- tax[c(-646),]

#write up
write.table(tax, "input/Phylo_Misc/rankedlin_Dec22.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")
