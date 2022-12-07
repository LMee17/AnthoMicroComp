#October 2022
#Iterate through the split files of large taxonomy dump file and process each, 
#pulling out only wanted microbial genera

#First: Viruses
taxfiles <- list.files(path = "input/Phylo_Misc/", pattern = "*.tsv")
#remove the original file
taxfiles <- taxfiles[c(-1,-2)]
taxfiles

#make space to store viral entries
#no idea how many there will be, so let's guess maybe 10? Safe side.
vir.list <- vector(mode = "list", length = 10)

j <- 1
for (i in 1:length(taxfiles)){
  x <- read.table(paste("input/Phylo_Misc/", taxfiles[i], sep = ""), 
                  comment.char = "", header = F, sep = "\t", quote = "")
  #take the header off the first file
  if(i == 1){
    x <- x[-1,]
  }
  x <- subset(x, V10 == "Viruses")
  v.check <- nrow(x)
  if (v.check > 1){
    print(j)
    vir.list[[j]] <- x
    j <- j + 1
  }
}

vir.df <- bind_rows(vir.list)

#replace blanks with NA and then remove any lines with NA
#then remove unnecessary columns
vir.df <- vir.df %>% 
  mutate_all(na_if, "")

vir.df <- vir.df[complete.cases(vir.df),]

vir.df <- vir.df[,c(2:10)]

names(vir.df) <-  c("Common name", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "SuperKingdom")

write.table(vir.df, "input/Phylo_Misc/ViralTaxo_Oct22.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")

#And now for the rest of the taxonomy
#define what's getting kept: microbial (bacteria/archaea/viruses/single celled)
wantedSuper <- c("Bacteria", "Archaea")
wantedKing <- c("Fungi", "Ciliophora", "Choanoflagellata", "Heterolobosea",
                "Evosea", "Apicomplexa", "Oomycota", "Bigyra")
wantedOrder <- c("Trypanosomatida", "Trombidiformes")
wantedPhylum <- c("Nematoda")
wantedGenus <- c("Apocephalus", "Aethina")

for (i in 1:length(taxfiles)){
  x <- read.table(paste("input/Phylo_Misc/", taxfiles[i], sep = ""), 
                  comment.char = "", header = F, sep = "\t", quote = "")
  if(i == 1){
    x <- x[-1,]
  }
  x <- x[,c(4:10)]
  #don't include lines that don't have a genus entry
  x <- x[!x$V4 == "",]
  x2 <- x[x$V10 %in% wantedSuper |
            x$V9 %in% wantedKing |
            x$V6 %in% wantedOrder |
            x$V4 %in% wantedGenus |
            x$V8 %in% wantedPhylum,]
  write.table(x2, paste("input/Phylo_Misc/", taxfiles[i], sep = ""),
              col.names = F, row.names = F, sep = "\t", quote = F)
}

#Cat the .tsv files together in commandline