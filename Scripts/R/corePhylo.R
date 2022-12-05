#CorePhylogroups #2
#14th November 2022

library(ggplot2)
library(reshape)
library(tidyverse)
library(RColorBrewer)
library(dplyr)

dir.create("output/Plots/")
dir.create("output/Plots/CorePhylo/")

#####Functions#####
#make plottable df, removing samples that have no counts at all
plotGen <- function(microbes,samples,cntmat){
  ex.cnt <- cntmat[rownames(cntmat) %in% microbes,
                   names(cntmat) %in% samples]
  ex.cnt <- ex.cnt[,!colSums(ex.cnt) == 0]
  ex.cnt$MicroID <- rownames(ex.cnt)
  tmp <- melt(ex.cnt)
  ex.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
  names(ex.melt)[c(1,3)] <- c("Sample", "Count")
  return(ex.melt)
}

#make plottable df, but with 0s kept in
plotGen0 <- function(microbes,samples,cntmat){
  ex.cnt <- cntmat[rownames(cntmat) %in% microbes,
                   names(cntmat) %in% samples]
  ex.cnt$MicroID <- rownames(ex.cnt)
  tmp <- melt(ex.cnt)
  ex.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
  names(ex.melt)[c(1,3)] <- c("Sample", "Count")
  return(ex.melt)
}

#make plot showing the prevalence of each microbial species per host taxonomic specification
#make sure that host taxonomic specification is input in quotes
#requires melted dataframe such as those made above
#CAUTION: this function doesn't work if the number of sample species is the same
#as the microbe species
plotPrev <- function(df.melt,microbes,HostTax, key){
  stat <- data.frame(Host = rep(unique(df.melt[,HostTax]), length(microbes)),
                     MicroID = rep(unique(df.melt$MicroID), length(unique(df.melt[,HostTax]))))
  for (i in 1:nrow(stat)){
    stat$Total[i] <- nrow(df.melt[df.melt[,HostTax] == stat$Host[i],])/length(microbes) 
    x <- df.melt[df.melt[,HostTax] == stat$Host[i] & df.melt$MicroID == stat$MicroID[i],]
    stat$Present[i] <- sum(x$Count)
    stat$PresProp[i] <- (sum(x$Count)/stat$Total[i])
    stat$Absent[i] <- stat$Total[i] - stat$Present[i]
    stat$AbsProp[i] <- 1 - stat$PresProp[i]
  }
  plot.df <- data.frame(Host = rep(unique(df.melt[,HostTax]),length(microbes)),
                        MicroID = rep(unique(df.melt$MicroID),length(unique(df.melt[,HostTax]))),
                        Prop = stat$PresProp,
                        Inc = paste("Present"))
  tmp <- data.frame(Host = rep(unique(df.melt[,HostTax]),length(microbes)),
                    MicroID = rep(unique(df.melt$MicroID),length(unique(df.melt[,HostTax]))),
                    Prop = stat$AbsProp,
                    Inc = paste("Absent"))
  plot.df <- rbind(plot.df, tmp)
  if (key == "genus"){
    for (i in 1:nrow(plot.df)){
      plot.df$MicroPhy[i] <- unique(genkey$PhyloHit[genkey$GenHitID == plot.df$MicroID[i]])  
    }
  }
  if (key == "micro"){
    for (i in 1:nrow(plot.df)){
      plot.df$MicroPhy[i] <- microkey$TaxHit[microkey$HitID == plot.df$MicroID[i]]
    }
  }
  ggplot(data = plot.df, aes(x = Host, y = Prop, fill = Inc)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.text = element_text(face = "italic"),
          strip.text = element_text(face = "italic")) + 
    labs(x = paste("Host", HostTax),
         fill = "",
         y = "Prevalence") +
    facet_wrap(~MicroPhy, labeller = label_wrap_gen(width = 2, multi_line = TRUE))
}

#####Load Necessary Counts and Metadata####
#Sample metadata
met <- read.table("input/SRA/MetaData_Edit_Oct22.tsv",
                  header = T, sep ="\t", quote = "")

#fix tissue field in metadata
met$Tissue2 <- met$Tissue
met$Tissue2[grepl("pooled whole bodies ", met$Tissue)] <- "whole body, pooled"
met$Tissue2[grepl("whole body ", met$Tissue)] <- "whole body, pooled"
met$Tissue2[grepl("pooled whole bodies", met$Tissue)] <- "whole body, pooled"
met$Tissue2[grepl("head and abdomen ", met$Tissue)] <- "abdomen and head, pooled"
met$Tissue2[grepl("pool of heads and abdomen", met$Tissue)] <- 
  "abdomen and head, pooled"
met$Tissue2[grepl("pooled gut", met$Tissue)] <- "gut, pooled"
met$Tissue2[grepl("midgut", met$Tissue)] <- "gut"
met$Tissue2[grepl("pool of abdomens and fat bodies", met$Tissue)] <- 
  "abdomen and fatbody, pooled"
met$Tissue2[grepl("abdomen with fat bodies", met$Tissue)] <- 
  "abdomen and fatbody"

#and finally
met$Tissue3 <- met$Tissue
unique(met$Tissue3)
met$Tissue3[grepl("whole bod", met$Tissue)] <- "whole body"
met$Tissue3[grepl("gut", met$Tissue)] <- "gut"
met$Tissue3[grepl("abdomen with fat", met$Tissue)] <- "abdomen and fat body"
met$Tissue3[grepl("head and abdomen", met$Tissue)] <- "abdomen and head"
met$Tissue3[grepl("abdomens and fat", met$Tissue)] <- "abdomen and fat body"
met$Tissue3[grepl("heads and abdomen", met$Tissue)] <- "abdomen and head"

#add continent values to metadata
#fix whitespace issue
met$Location <- trimws(met$Location)

met$Continent[met$Location == "Australia" | met$Location == "New Zealand"] <- "Oceania"
met$Continent[met$Location == "Belgium" | met$Location == "France" |
                met$Location == "Germany" | met$Location == "Italy" |
                met$Location == "Scotland" |
                met$Location == "United Kingdom"  |
                met$Location == "Switzerland" |
                met$Location == "Netherlands"  ] <- "Europe"
met$Continent[met$Location == "Brazil"  ] <- "South America"
met$Continent[met$Location == "California" | met$Location == "Canada" |
                met$Location == "Florida" | met$Location == "Hawaii" |
                met$Location == "California" | met$Location == "Canada" | 
                met$Location == "Panama" | met$Location == "Mexico" |
                met$Location == "North Dakota" | met$Location == "Utah"] <- "North America"
met$Continent[met$Location == "China" | met$Location == "Israel" |
                met$Location == "South Korea" ] <- "Asia"
met$Continent[met$Location == "South Africa" ] <- "Africa"

#save
write.table(met, "input/SRA/MetaData_Edit_Nov22.tsv",
            sep = "\t", quote = F, col.names = T, row.names = F)

#Microbial metadata
genkey <- read.table("input/Keys/GenusFamilySpeciesKey_Nov22.tsv",
                     sep = "\t", header = T)
#Taxonomy
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22_v2.tsv",
                  header = T, sep = "\t", quote = "")

#read in incidence and rPM tables
cnt.inc <- read.table("input/Counts/Incidence/All_Incidence_CorePhylo_PhyloCollapsed_Nov22.tsv")
cnt.rpm <- read.table("input/Counts/Raw_rPM/All_raw_rPM_CorePhylo_PhyloCollapsed_Nov22.tsv")


#####Quick Incidence Overview: Core Phylo + All Samples####
#remove DNA samples
dna <- met$Sample.ID[met$NucleotideType =="DNA"]
cnt.inc <- cnt.inc[,! names(cnt.inc) %in% dna]
cnt.rpm <- cnt.rpm[,! names(cnt.rpm) %in% dna]

#extract the core microbiota phylotypes
beecore <- c("Lactobacillus: Firm-5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
             "Bartonella", "Bombiscardovia", "Schmidhempelia", "Apibacter")
bcoreID <- unique(genkey$GenHitID[genkey$PhyloHit %in% beecore])
bc.inc <- cnt.inc[rownames(cnt.inc) %in% bcoreID,]

#make plottable table
bc.inc$MicrobialHit <- rownames(bc.inc)
tmp <- melt(bc.inc)
#add metadata
#for host
bc.melt <- merge(tmp, met, by.x= "variable", by.y = "Sample.ID")
#for microbes
for (i in 1:nrow(bc.melt)){
  bc.melt$PhyloHit[i] <- unique(genkey$PhyloHit[genkey$GenHitID == bc.melt$MicrobialHit[i]])
}

#rename variable to sample
names(bc.melt)[1] <- "Sample"

#factorise presence/absence
bc.melt$Inc <- ifelse(bc.melt$value == 1, "Present", "Absent")
bc.melt$Inc <- factor(bc.melt$Inc, levels = c("Present", "Absent"))

#arrange microbials
bc.melt$PhyloHit <- factor(bc.melt$PhyloHit, levels = c("Bombiscardovia",
                                      "Apibacter", "Bartonella", "Frischella",
                                      "Bombilactobacillus", "Apilactobacillus",
                                      "Bifidobacterium", "Lactobacillus: Firm-5",
                                      "Snodgrassella", "Gilliamella"))

#order samples ... this will be big
#most social to least social
bc.melt$Sample <- factor(bc.melt$Sample, levels = c(
  #apis
  "SRR072891","SRR072892","SRR10034516","SRR10034517","SRR10034518",
  "SRR10034523","SRR12097444","SRR12097445","SRR12097446","SRR12097447",
  "SRR13189057","SRR13189058","SRR13189059","SRR13189063","SRR13189066",
  "SRR13189067","SRR13404632","SRR13404634","SRR13404636","SRR13404641",
  "SRR13442815","SRR13442819","SRR13442827","SRR13442828","SRR13442829",
  "SRR13442831","SRR1380984","SRR15496106","SRR15496107","SRR15496108",
  "SRR15496109","SRR15496114","SRR15496115","SRR17475010","SRR17475011",
  "SRR17475012","SRR17475013","SRR17475014","SRR1790685","SRR1790686",
  "SRR1790687","SRR1790690","SRR1790691","SRR18500042","SRR18500043",
  "SRR18500044","SRR18500045","SRR18500050","SRR18500061","SRR18500072",
  "SRR18500076","SRR18500077","SRR18500078","SRR18500079","SRR18500080",
  "SRR18500081","SRR18500082","SRR18500083","SRR18500084","SRR18500085",
  "SRR18500086","SRR18500087","SRR18500088","SRR18500089","SRR18500090",
  "SRR18500091","SRR5109820","SRR5109826","SRR5109829","SRR5109831","SRR5109832",
  "SRR5109833","SRR5117442","SRR5117443","SRR5117444","SRR5117445","SRR5117446",
  "SRR5117447","SRR5117448","SRR5117449","SRR5117450","SRR5387735","SRR5387737",
  "SRR5387738","SRR8335255","SRR8335256","SRR8335257","SRR8335258","SRR8754036",
  "SRR8754037","SRR8754038","SRR8754039","SRR8754040","SRR8754041","SRR8754042",
  "SRR8754043","SRR8754044","SRR8867392","SRR8867393","SRR8867394","SRR8867395",
  "AM_N_056","AM_N_068","AM_N_072",
  #tetragonula
  "SRR13442814","SRR13442818","SRR13442830",
  #tetragonisca
  "SRR11426431","SRR11426432","SRR11426433","SRR11440494","SRR11440495","SRR11440496",
  #bombus
  "SRR12245258","SRR12245523","SRR12527928","SRR12527929","SRR12527930","SRR12527931",
  "SRR12527932","SRR12527933","SRR12527934","SRR12527935","SRR12527936","SRR12527937",
  "SRR12527938","SRR12527939","SRR12527940","SRR12527941","SRR12527942","SRR12527943",
  "SRR12527944","SRR12527945","SRR12527946","SRR12527947","SRR12527948","SRR12527949",
  "SRR12527950","SRR12527951","SRR12527952","SRR12527953","SRR12527954","SRR12527955",
  "SRR12527956","SRR12527957","SRR12527958","SRR12527959","SRR12527961","SRR12527962",
  "SRR12527963","SRR12527964","SRR12527965","SRR12527966","SRR12527967","SRR12527968",
  "SRR12527969","SRR12527970","SRR12527971","SRR12527972","SRR12527974","SRR12527975",
  "SRR12527977","SRR13769013","SRR13769014","SRR13769015","SRR13769016","SRR13769017",
  "SRR13769018","SRR13769019","SRR14567209","SRR14567210","SRR14567211","SRR14567212",
  "SRR14567213","SRR14567214","SRR14567217","SRR14567218","SRR14567219","SRR14567220",
  "SRR14567221","SRR1503088","SRR2396633","SRR2396634","SRR2396640","SRR2396641",
  "SRR2396644","SRR2396653","SRR2396654","SRR2396655","SRR2396656","SRR2396657",
  "SRR6148366","SRR6148368","SRR6148369","SRR6148370","SRR6148372","SRR6148374",
  "SRR6148376","SRR6148378","SRR12235348","SRR12244814","ERR10123638","SRR3383868",
  "SRR3383869","SRR11445089","SRR11445090","SRR11448239","SRR11448240","SRR11448241",
  "SRR12233573","SRR12233820","SRR6370988","SRR14567222","BT_N_001","BT_N_025",
  "BT_N_029","SRR12527960","SRR11445091","SRR12527973","SRR12527976",
  #euglossa
  "SRR12053406","SRR12053407","SRR12053408","SRR12053409","SRR12053410","SRR12053411",
  "SRR12053412","SRR12053413","SRR12053414","SRR12053415","SRR12053416","SRR12053417",
  "SRR12053418","SRR12053419","SRR12053420","SRR12053421","SRR12053422","SRR12053423",
  "SRR12053424","SRR12053425","SRR12053426","SRR12053427","SRR12053428","SRR12053429",
  "SRR12053430","SRR12053431","SRR12053432","SRR12053433","SRR12053434","SRR12053435",
  "SRR12053436","SRR12053437","SRR1503113",
  #lasioglossum
  "SRR13442813","SRR13442825",
  #exoneura
  "SRR13442817","SRR13442822","SRR13442824",
  #ceratina
  "SRR2915407","SRR2954866","CA_N_007","CA_N_012","CA_N_030",
  #Eufriesea
  "SRR2001621","SRR1945063","SRR1945064","SRR1945065",
  #Megalopta
  "SRR3948521","SRR3948528","SRR3948531","SRR3948533","SRR3948536","SRR3948537",
  "SRR3948540","SRR3948542","SRR3948544","SRR3948548","SRR3948550","SRR3948554",
  "SRR3948556","SRR3948560","SRR3948562","SRR3948563","SRR3948564","SRR3948566",
  "SRR3948568","SRR3948570","SRR3948572","SRR3948575","SRR3948577",
  #Halictus
  "SRR5253659",
  #Epeolus
  "SRR1503066",
  #Thyreus
  "SRR1503122",
  #Anthophora
  "SRR1503129",
  #Nomada
  "SRR1503133",
  #xylocopa
  "SRR1503136",
  #Osmia
  "SRR6148364","SRR6148371","SRR6148375","SRR6148379","SRR3383870","SRR2895245",
  "SRR2895246","SRR2895247","SRR2895248","SRR2895249","SRR2895250","SRR2895251",
  "SRR9676054","SRR9676058","SRR9676059",
  #Andrena
  "SRR6148365","SRR6148367","SRR6148373","SRR6148377","SRR3383871","SRR8335251",
  "SRR8335252","SRR8335253","SRR8335254","SRR13404631","SRR13404633","SRR13404635",
  "SRR13404640",
  #habropoda
  "SRR2001612","SRR1943074","SRR1943075","SRR1943076",
  #Colletes
  "SRR10765021",
  #Dufourea
  "SRR2001766", "SRR1945093", "SRR1945094", "SRR1945095"
))
  
#plot 
#x axis begins with most social and moves over to solitary (left to right)
ggplot(data = bc.melt, aes(x = Sample, y = PhyloHit, fill = Inc)) +
  geom_tile() +
  scale_fill_manual(values = c("#601a4a","#63abce"))+
  labs(fill = "",
       x = "Host Samples",
       y = "Microbial Taxa") + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_blank())


#####Collapsing by Host Genera / Sociality####
#collapse samples into genera
genera <- unique(met$Genus)

genmat <- matrix(nrow = nrow(bc.inc), ncol = length(genera))
spepro <- vector(length = length(genera))
for (i in 1:length(genera)){
 x <- met$Sample.ID[met$Genus == genera[i]]
 print(x)
 for (j in 1:length(bcoreID))
   genmat[j,i] <- sum(bc.inc[bcoreID[j], names(bc.inc) %in% x])
}

gen.df <- as.data.frame(genmat)
names(gen.df) <- genera
rownames(gen.df) <- bcoreID

#from looking at this, it's pretty pointless to have it broken down by genus 
#keep genera for the social species, then collapse solitary / polymorphic together

#store the solitary species' genera as one group
sol <- met$Genus[met$Sociality == "Solitary"]
#store the polymorphic species' genera as one group
pol <- met$Genus[met$Sociality == "Polymorphic"]
#store Meliponini species as one group
mel <- c("Tetragonula", "Tetragonisca")
#combine these into a list to run through
wantedGen <- list("Apis", "Bombus", mel, pol, sol)

#create a matrix ready to populate
genmat2 <- matrix(nrow = nrow(bc.inc), ncol = length(wantedGen))
#iterate through each item of genera in the wantedGen list
for (i in 1:length(wantedGen)){
  #extract the sample IDs of samples within this genus/group of genera
  x <- met$Sample.ID[met$Genus %in% wantedGen[[i]]]
  #for each of the core phylogroup stored in bcoreID, populate the column corresponding
  #to this species (i) with each incidence value per phylogroup (j)
  for (j in 1:length(bcoreID)){
    #fi tehre is more that one incidence, sum them
    genmat2[j,i] <- sum(bc.inc[bcoreID[j], names(bc.inc) %in% x])
  }
}
#change into a dataframe for ease of data manipulation
gen.df2 <- as.data.frame(genmat2)
names(gen.df2) <- c("Apini", "Bombini", "Meliponini", 
                    "Polymorphic\nSpecies", "Solitary\nSpecies")
rownames(gen.df2) <- bcoreID

#make plottable table (dplyr melt)
#make the rownames of the count table a separate variable
gen.df2$MicrobialHit <- rownames(gen.df2)
#melt the count table into a dataframe broken down by variable values per incidence
#(column headers become variables)
gen.melt <- melt(gen.df2)

#add metadata
for (i in 1:nrow(gen.melt)){
  #pull out any socialities associated with "variable" = either a tribe name for 
  #the corbiculates, else a manual annotation of poly or solitary species
  x <- unique(met$Sociality[met$Tribe == gen.melt$variable[i]])
  #if there is a hit, ie, if the "variable" is a tribe, then add the extracted
  #sociality to the dataframe
  if(length(x) > 0){
    gen.melt$Sociality[i] <- x
  }
  #if there isn't a hit, meaning the "variable" is not a tribe, then leave the 
  #value blank
  if(length(x) ==0){
    gen.melt$Sociality[i] <- ""
  }
}
#now fix these manual entries
gen.melt$Sociality[gen.melt$variable == "Polymorphic\nSpecies"] <- "Polymorphic"
gen.melt$Sociality[gen.melt$variable == "Solitary\nSpecies"] <- "Solitary"

#convert to proportions
#get the number of samples per genera class ready
#pull out the column names of the genus dataframe, but remove the Microbial Hit column
genmet <- data.frame(GenClass = names(gen.df2)[-6])
#empty vector ready to be filled - one per classification outlined above (5, 3 tribes
#plus all polymorphic and all solitary)
totsamp <- vector(length = nrow(genmet))
#use the wantedGen list object as a tool to count how many samples there are per 
#classification
for (i in 1:length(wantedGen)){
  #extract the samples that belong to this classification
  x <- met$Sample.ID[met$Genus %in% wantedGen[[i]]]
  #keep those that are actually in the analysis (the metadata has samples that I removed above)
  y <- x[x %in% names(bc.inc)]
  #count the number of samples there are per classification
  totsamp[i] <- length(y)
}
#add to the genus metadata object
genmet$TotalSamples <- totsamp

ncol(cnt.inc[,names(cnt.inc) %in% sol2])

#get proportion into the gen plot dataframe
for (i in 1:nrow(gen.melt)){
  #extract the total number of samples from the gen metadata object
  tot <- genmet$TotalSamples[genmet$GenClass == gen.melt$variable[i]]
  #add another column with the proportion of incidences per classification
  gen.melt$PercPresent[i] <- (gen.melt$value[i] / tot)*100
}

#make some factors for ease of colouring the heatmap
gen.melt$PercPresentFact[gen.melt$PercPresent > 60] <- ">60%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 60 &
                           gen.melt$PercPresent > 50] <- "51-60%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 50 &
                           gen.melt$PercPresent > 40] <- "41-50%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 40 &
                           gen.melt$PercPresent > 30] <- "31-40%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 30 &
                           gen.melt$PercPresent > 20] <- "21-30%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 20 &
                           gen.melt$PercPresent > 10] <- "10-20%"
gen.melt$PercPresentFact[gen.melt$PercPresent < 10] <- "<10%"
gen.melt$PercPresentFact[gen.melt$PercPresent == 0.000000] <- " "

#order
gen.melt$PercPresentFact <- factor(gen.melt$PercPresentFact,
                                   levels = c(" ", "<10%", "10-20%", "21-30%",
                                              "31-40%", "41-50%", "51-60%", ">60%"))

#bring the core phylo names back in
for (i in 1:nrow(gen.melt)){
  gen.melt$MicroPhyLabel[i] <- unique(genkey$PhyloHit[genkey$GenHitID == gen.melt$MicrobialHit[i]])
}

#order microbes
gen.melt$MicroPhyLabel <- factor(gen.melt$MicroPhyLabel, 
                                 levels = c("Bombiscardovia",
                                            "Bartonella", "Frischella", 
                                            "Apibacter",
                                            "Bombilactobacillus", "Apilactobacillus",
                                            "Bifidobacterium", "Lactobacillus: Firm-5",
                                            "Snodgrassella", "Gilliamella"))

#plot
#heatmap with darker colours for the higher the prevalence of each microbial taxa
ggplot(data = gen.melt, aes(x = variable, y = MicroPhyLabel)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        panel.background = element_blank()) +
  geom_tile(aes(fill = (PercPresentFact))) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("white",
                               "#cfffff",
                               "#00e2e2",
                               "#00b8b8",
                               "#006b6b",
                               "#005151", 
                               "#001e1e",
                               "black")) +
  labs(fill = "Prevalence")+
  theme(panel.border=element_rect(fill = NA, colour=alpha('black', 
                                                              .5),size=1)) 

ggsave("output/Plots/CorePhylo/CorePhylo_AcrossBees_Heatmap.pdf")

#Looking at just the corbiculate species
ggplot(data = gen.melt[gen.melt$Sociality == "Eusocial",], 
       aes(x = variable, y = MicroPhyLabel)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        panel.background = element_blank()) +
  geom_tile(aes(fill = (PercPresentFact))) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("white",
                               "#cfffff",
                               "#00e2e2",
                               "#00b8b8",
                               "#006b6b",
                               "#005151", 
                               "#001e1e",
                               "black")) +
  labs(fill = "Prevalence")+
  theme(panel.border=element_rect(fill = NA, colour=alpha('black', 
                                                          .5),size=1))

ggsave("output/Plots/CorePhylo/CorePhylo_AcrossBees_CorbOnly_Heatmap.pdf")


#####Frischella and Bombus#####
#Frischella is found quite prevalently in Bombus, despite it apparently not being 
#found previously in literature (from a very shallow reading)
#Where do I find it?
bom.fri <- bc.melt %>%
  subset(PhyloHit == "Frischella" & Genus == "Bombus" & Inc == "Present") %>%
  select(Sample)

#look at species level
spec.rpm <- read.table("input/Counts/Raw_rPM/All_raw_rPM_Nov22.tsv")
#load the required key for species-level IDs
microkey <- read.table("input/Keys/MicrobialSpeciesKey.tsv", 
                       sep = "\t", header = T)
#extract the frischella species IDs
friID <- microkey$HitID[grepl("Frischella", microkey$TaxHit)]

#pull out rpm from counts tables
bom.fri.rpm <- spec.rpm[rownames(spec.rpm) %in% friID,
                        names(spec.rpm) %in% bom.fri$Sample]

#make plottable
bom.fri.rpm$MicroID <- rownames(bom.fri.rpm)
bom.fri.melt <- melt(bom.fri.rpm)
#add species info for microbes
bom.fri.melt$MicrobialSpecies[bom.fri.melt$MicroID == "MRH295"] <- "Other Frischella"
bom.fri.melt$MicrobialSpecies[bom.fri.melt$MicroID == "MRH1388"] <- "Frischella perrara"
#add species info for samples
for (i in 1:nrow(bom.fri.melt)){
  bom.fri.melt$HostSpecies[i] <- met$Species[met$Sample.ID == bom.fri.melt$variable[i]]
}
#rename "variable" and "value"
names(bom.fri.melt)[2:3] <- c("Sample", "rPM")

#plot - barplot
ggplot(data = bom.fri.melt, aes(x = Sample, y = rPM, fill = HostSpecies))+
  geom_bar(stat = "identity") +
  facet_grid(MicrobialSpecies~.) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic")) +
  labs(fill = "Host Species")

ggsave("output/Plots/CorePhylo/Bombus_Frisch_rPM.pdf")

#whats the average Frischella rpm in Apis ? (where it's present)
#extract apis samples and frischella rows from species-level rpm table
apis <- met$Sample.ID[met$Genus == "Apis"]
apis.fri.rpm <- spec.rpm[rownames(spec.rpm) %in% friID, 
                        names(spec.rpm) %in% apis]
#melt into something more useable
apis.fri.rpm$MicroID <- rownames(apis.fri.rpm)
apis.fri.melt <- melt(apis.fri.rpm)

#remove samples that have no frischella 
apis.fri.melt <- filter(apis.fri.melt, value != 0)
#assess average rpm
apis.fri.melt %>%
  group_by(MicrobialSpecies) %>%
  summarise(average = mean(value), max(value), min(value))

#compared to the bombus ? 
#add species info for microbes
apis.fri.melt$MicrobialSpecies[apis.fri.melt$MicroID == "MRH295"] <- "Other Frischella"
apis.fri.melt$MicrobialSpecies[apis.fri.melt$MicroID == "MRH1388"] <- "Frischella perrara"
#add species info for samples
for (i in 1:nrow(apis.fri.melt)){
  apis.fri.melt$HostSpecies[i] <- met$Species[met$Sample.ID == apis.fri.melt$variable[i]]
}
#rename "variable" and "value"
names(apis.fri.melt)[2:3] <- c("Sample", "rPM")

#combine
api.melt <- rbind(bom.fri.melt, apis.fri.melt)
api.melt$Species2[grepl("Apis", api.melt$HostSpecies)] <- "Apis"
api.melt$Species2[grepl("Bombus", api.melt$HostSpecies)] <- "Bombus"
#plot
ggplot(data = api.melt, aes(x = Species2, y = log(rPM), colour = MicrobialSpecies)) +
  geom_boxplot(size = 1) +
  geom_jitter() +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_colour_manual(values = c("#7ccc00", "#812fff")) +
  labs(fill = "Microbial Species",
       x = "Host Genus")

ggsave("output/Plots/CorePhylo/ApisVBombus_Frischella.pdf")

api.melt %>%
  group_by(Species2) %>%
  summarise(average = mean(rPM), max(rPM), min(rPM))

#####Gilliamella and Snodgrassella #####
#both of these are found to be more at higher abundance in bumbles than honeys
#do we replicate that ? 
#will also add meliponini samples in for completeness
#snod = GRH8, gill = GRH7
#looking at reads per million

#extract bombus samples IDs
bombus <- met$Sample.ID[met$Genus == "Bombus"]
meli <- met$Sample.ID[met$Tribe == "Meliponini"]

sg.rpm <- cnt.rpm[rownames(cnt.rpm) == "GRH8" | rownames(cnt.rpm) == "GRH7",
                  names(cnt.rpm) %in% bombus | names(cnt.rpm) %in% apis |
                    names(cnt.rpm) %in% meli]

#remove samples that have no counts
sg.rpm <- sg.rpm[,!colSums(sg.rpm) == 0]

#make plottable and add metadata
sg.rpm$MicrobeID <- rownames(sg.rpm)
tmp <- melt(sg.rpm)
sg.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
#add microbiota phylo labels
sg.melt$MicroPhylo[sg.melt$MicrobeID == "GRH8"] <- "Snodgrassella"
sg.melt$MicroPhylo[sg.melt$MicrobeID == "GRH7"] <- "Gilliamella"
#make more informative column headers
names(sg.melt)[c(1,3)] <- c("Sample", "rPM")

#plot (boxplot)
ggplot(data = sg.melt, aes(x = Genus, y = log(rPM), fill = MicroPhylo)) +
  geom_boxplot() +
  theme(axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#e6e600", "#0000e6")) +
  labs(fill = "Microbial Genus", 
       x = "Host Genus")

#woop. We do recapture it.
ggsave("output/Plots/CorePhylo/ApisVBombus_GillnSnod.pdf")



#assess average rPM between the two genera (where these microbes are found)


#####Assessing potential species-level differences in Bombilactobacillus####
#make species-level incidence table
spec.inc <- spec.rpm
spec.inc[spec.inc > 0] <- 1
#remove DNA samples
spec.inc <- spec.inc[, !names(spec.inc) %in% dna]

#extract species of bactobacillus
blID <- genkey$MicroHitID[genkey$PhyloHit == "Bombilactobacillus"]
#combine previous lists of corb samples
corb <- c(apis, bombus, meli)

#make amenable for plotting
#going to make this a function as I'm sick of redoing it
plotGen <- function(microbes,samples,cntmat){
  ex.cnt <- cntmat[rownames(cntmat) %in% microbes,
                   names(cntmat) %in% samples]
  ex.cnt <- ex.cnt[,!colSums(ex.cnt) == 0]
  ex.cnt$MicroID <- rownames(ex.cnt)
  tmp <- melt(ex.cnt)
  ex.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
  names(ex.melt)[c(1,3)] <- c("Sample", "Count")
  return(ex.melt)
}

bl.melt <- plotGen(blID, corb, spec.inc)

#add species information
microkey[microkey$HitID %in% bl.melt$MicroID,]
for (i in 1:nrow(bl.melt)){
  bl.melt$MicroPhylo[i] <- unique(microkey$TaxHit[microkey$HitID == bl.melt$MicroID[i]])
}

#order by genus
bl.melt$Sample <- factor(bl.melt$Sample, 
                          levels = c(#Apis cerana
                            "SRR10034516", "SRR10034523",
                            #Apis mellifera
                            "SRR13189057", "SRR13189058", "SRR13189059", "SRR13404632", 
                            "SRR13404636", "SRR13404641", "SRR13442819", "SRR18500042", 
                            "SRR18500043", "SRR18500044", "SRR18500045", "SRR18500050", 
                            "SRR18500061", "SRR18500072", "SRR18500076", "SRR18500077", 
                            "SRR18500078", "SRR18500079", "SRR18500080", "SRR18500081", 
                            "SRR18500082", "SRR18500083", "SRR18500084", "SRR18500085",
                            "SRR18500086", "SRR18500087", "SRR18500088", "SRR18500090",
                            "SRR5109820",  "SRR5109826",  "SRR5109829",  "SRR5109831", 
                            "SRR5109832",  "SRR5109833",
                            #bombus lucorum
                            "SRR13769016",
                            #Bombus pascurorum
                            "SRR13769014",
                            #Bombus terrestris
                            "SRR11445089", "SRR11445090", "SRR11445091", "SRR11448239", 
                            "SRR11448240", "SRR11448241", "SRR13769018","SRR2396655"))


#plot
ggplot(data = bl.melt, aes(x = Sample, y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(MicroPhylo~.) +
  theme(legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Host Species",
       y = "") +
  scale_fill_manual(values = c("red",
                               "orange",
                               "skyblue",
                               "navy",
                               "blue"))

ggsave("output/Plots/CorePhylo/Bombilactobacillus_ByMicrobeSpecies.pdf")

#####Lactobacillus Firm-5#####
#extract firm-5 ids
f5 <- genkey$MicroHitID[genkey$PhyloHit == "Lactobacillus: Firm-5"]

#make plottable df, but with 0s kept in
plotGen0 <- function(microbes,samples,cntmat){
  ex.cnt <- cntmat[rownames(cntmat) %in% microbes,
                   names(cntmat) %in% samples]
  ex.cnt$MicroID <- rownames(ex.cnt)
  tmp <- melt(ex.cnt)
  ex.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
  names(ex.melt)[c(1,3)] <- c("Sample", "Count")
  return(ex.melt)
}


#make plottable df
f5.melt <- plotGen0(f5, corb, spec.inc)

#make stacked barchart with presence / absence (prevalence) per species
#start a new dataframe based off the unique species included in the extracted count table
#I need a proportion per host species per microbe (ie 24 * 7)
f5.stat <- data.frame(Species = rep(unique(f5.melt$Species),length(f5)),
                      MicroID = rep(unique(f5.melt$MicroID),24))
#add total count of species from the melt dataframe
for (i in 1:nrow(f5.stat)){
  f5.stat$Total[i] <- nrow(f5.melt[f5.melt$Species == f5.stat$Species[i],])/7  
  x <- subset(f5.melt, Species == f5.stat$Species[i] & MicroID == f5.stat$MicroID[i])
  f5.stat$Present[i] <- sum(x$Count)
  f5.stat$PresProp[i] <- sum(x$Count)/f5.stat$Total[i]
  f5.stat$Absent[i] <- f5.stat$Total[i] - f5.stat$Present[i]
  f5.stat$AbsProp[i] <- 1 - f5.stat$PresProp[i]
}

f5.df2 <- data.frame(Species = rep(unique(f5.melt$Species),length(f5)),
                    MicroID = rep(unique(f5.melt$MicroID),24),
                    Prop = f5.stat$PresProp,
                    Inc = paste("Present"))
tmp <- data.frame(Species = rep(unique(f5.melt$Species),length(f5)),
                   MicroID = rep(unique(f5.melt$MicroID),24),
                   Prop = f5.stat$AbsProp,
                   Inc = paste("Absent"))
f5.df2 <- rbind(f5.df2, tmp)

#add microbe phylo labels
for (i in 1:nrow(f5.df2)){
  f5.df2$MicroPhy[i] <- microkey$TaxHit[microkey$HitID == f5.df2$MicroID[i]]
}

#remove any species that have no hits anywhere
tmp <- f5.df2 %>%
  filter(Inc == "Present") %>%
  group_by(Species) %>%
  summarise(sum(Prop))
remove <- tmp$Species[tmp$`sum(Prop)` == 0]
f5.df3 <- f5.df2[!f5.df2$Species %in% remove,]

#plot
ggplot(data = f5.df3, aes(x = Species, y = Prop, fill = Inc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(MicroPhy~.) +
  theme(axis.text.y = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic")) + 
  labs(x = "Host Species",
       fill = "",
       y = "Prevalence")

ggsave("output/Plots/CorePhylo/Lactobacillus_speciesBreakdown_Corb.pdf")

#####Bifidobacterium#####
#extract bifo ids
bifo <- unique(genkey$MicroHitID[genkey$PhyloHit == "Bifidobacterium"])

#make plottable df
bifo.melt <- plotGen0(bifo, corb, spec.inc)

#make stacked barchart with presence / absence (prevalence) per tribe
#start a new dataframe based off the unique species included in the extracted count table
#I need a proportion per host species per microbe (ie 24 * 7)
bifo.stat <- data.frame(Tribe = rep(unique(bifo.melt$Tribe),length(bifo)),
                      MicroID = rep(unique(bifo.melt$MicroID),length(unique(bifo.melt$Tribe))))
#add total count of species from the melt dataframe
for (i in 1:nrow(bifo.stat)){
  bifo.stat$Total[i] <- nrow(bifo.melt[bifo.melt$Tribe == bifo.stat$Tribe[i],])/length(bifo) 
  x <- subset(bifo.melt, Tribe == bifo.stat$Tribe[i] & MicroID == bifo.stat$MicroID[i])
  bifo.stat$Present[i] <- sum(x$Count)
  bifo.stat$PresProp[i] <- (sum(x$Count)/bifo.stat$Total[i])
  bifo.stat$Absent[i] <- bifo.stat$Total[i] - bifo.stat$Present[i]
  bifo.stat$AbsProp[i] <- 1 - bifo.stat$PresProp[i]
}

bifo.df2 <- data.frame(Tribe = rep(unique(bifo.melt$Tribe),length(bifo)),
                       MicroID = rep(unique(bifo.melt$MicroID),length(unique(bifo.melt$Tribe))),
                       Prop = bifo.stat$PresProp,
                       Inc = paste("Present"))
tmp <- data.frame(Tribe = rep(unique(bifo.melt$Tribe),length(bifo)),
                  MicroID = rep(unique(bifo.melt$MicroID),length(unique(bifo.melt$Tribe))),
                  Prop = bifo.stat$AbsProp,
                  Inc = paste("Absent"))
bifo.df2 <- rbind(bifo.df2, tmp)

#add microbe phylo labels
for (i in 1:nrow(bifo.df2)){
  bifo.df2$MicroPhy[i] <- microkey$TaxHit[microkey$HitID == bifo.df2$MicroID[i]]
}

#plot
ggplot(data = bifo.df2, aes(x = Tribe, y = Prop, fill = Inc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic")) + 
  labs(x = "Host Species",
       fill = "",
       y = "Prevalence") +
  facet_wrap(~MicroPhy)

ggsave("output/Plots/CorePhylo/Bifido_SpeciesBreakdown_Corbs.pdf")

#####Assessing Bifidobacterium Variability across Apis, Bombus, Meliponini#####
bi.gen <- genkey$GenHitID[genkey$PhyloHit == "Bifidobacterium"]
bi.melt <- plotGen(bi.gen, corb, cnt.rpm)

#plot
ggplot(bi.melt, aes(x = Species, y = log(Count), fill = Tribe))+
  geom_boxplot() +
  labs(y = "Log(rPM)") +
  theme(axis.text.y = element_text(face = "italic")) +
  coord_flip() 

ggsave("output/Plots/CorePhylo/Bifi_byrPM_corbs.pdf")

#####Apibacter: Prevalence Across Honey and BumbleBees#####
#make plottable dataframe
apib <- genkey$MicroHitID[genkey$PhyloHit == "Apibacter"]
#make list of samples for honey and bumble bees only
apinae <- c(apis, bombus)
apib.melt <- plotGen0(microbes = apib,
                      samples = apinae,
                      cntmat = spec.inc)

#plot
plotPrev(apib.melt, apib, "Species")
ggsave("output/Plots/CorePhylo/Apibacter_microSpecvApinaeSpec_Prev.pdf")

#again, but with Apibacter in general, not broken down by species
apib2 <- unique(genkey$GenHitID[genkey$PhyloHit == "Apibacter"])
apib.melt2 <- plotGen0(microbes = apib2,
                      samples = corb,
                      cntmat = cnt.inc)

plotPrev(apib.melt2, apib2, "Genus", "genus")
ggsave("output/Plots/CorePhylo/Apibacter_genusvCorb_Prev.pdf")

#####Apibacter versus Snodgrassella and Crithidia####
#apibacter = GRH128, snodgrassella = GRH8, crithidia = GRH434
#apibacter and snodgrassella
#looking at rpm
avs <- c("GRH128", "GRH8")

#extract rpm
avs.rpm <- cnt.rpm[rownames(cnt.rpm) %in% avs,]
#remove samples with neither present
avs.rpm <- avs.rpm[,! colSums(avs.rpm) == 0]

#make plottable
avs.rpm$MicroID <- rownames(avs.rpm)
tmp <- melt(avs.rpm)
avs.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
names(avs.melt)[c(1,3)] <- c("Sample", "rPM")

#annotate each sample as either Apibacter dominant, snodgrassella dominant or 
#absent in either case
samps <- unique(avs.melt$Sample)
#prepare vector to populate with classifications
class <- vector(length = length(samps))
#for each rpm table, row 1 = snodgrassella, row 2 = apibacter
for (i in 1:length(class)){
  s <- avs.rpm[1, names(avs.rpm) == paste(samps[i])]
  a <- avs.rpm[2, names(avs.rpm) == paste(samps[i])]
  if (s == 0){
    class[i] <- "Snodgrassella\nAbsent"
  }
  if (a == 0){
    class[i] <- "Apibacter\nAbsent"
  }
  if (s > a & a != 0){
    class[i] <- "Snodgrassella\nDominant"
  }
  if (a > s & s != 0){
    class[i] <- "Apibacter\nDominant"
  }
}
classkey <- data.frame(Sample = samps, Composition = class)

#add composition classification and microphylo information
for(i in 1:nrow(avs.melt)){
  avs.melt$Composition[i] <- classkey$Composition[classkey$Sample == avs.melt$Sample[i]]
  avs.melt$MicroPhylo[i] <- genkey$PhyloHit[genkey$GenHitID == avs.melt$MicroID[i]]
}

avs.melt$Composition <- factor(avs.melt$Composition, 
                               levels = c("Apibacter\nAbsent",
                                          "Apibacter\nDominant",
                                          "Snodgrassella\nDominant",
                                          "Snodgrassella\nAbsent"))

myPal <- brewer.pal(4, "Dark2")

ggplot(data = avs.melt, aes(x = Composition, y = log(rPM), fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa"))

ggsave("output/Plots/CorePhylo/ApibactervsSnodgrassella_AllBees.pdf")

ggplot(data = avs.melt, aes(x = Composition, y = log(rPM), fill = MicroPhylo, colour = Genus))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa")) +
  scale_colour_manual(values=myPal)

ggsave("output/Plots/CorePhylo/ApibactervsSnodgrassella_AllBees_IncludeGen.pdf")

#look at just corbiculates
#there only Andrena that's in here with corbiculates anyway
avs.melt2 <- subset(avs.melt, Genus != "Andrena")

ggplot(data = avs.melt2, aes(x = Composition, y = log(rPM), fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa"))

ggsave("output/Plots/CorePhylo/ApibactervsSnodgrassella_CorbBees.pdf")

myPal <- myPal[c(2,1,3)]

ggplot(data = avs.melt2, aes(x = Composition, y = log(rPM), fill = MicroPhylo, colour = Genus))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa")) +
  scale_colour_manual(values=myPal)

ggsave("output/Plots/CorePhylo/ApibactervsSnodgrassella_CorbBees_IncludeGen.pdf")

#apibacter and crithidia
#apibacter = GRH128, snodgrassella = GRH8, crithidia = GRH434
#looking at rpm
avc <- c("GRH128", "GRH434")

#extract rpm
avc.rpm <- cnt.rpm[rownames(cnt.rpm) %in% avc,]
#remove samples with neither present
avc.rpm <- avc.rpm[,! colSums(avc.rpm) == 0]

#make plottable
avc.rpm$MicroID <- rownames(avc.rpm)
tmp <- melt(avc.rpm)
avc.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
names(avc.melt)[c(1,3)] <- c("Sample", "rPM")

#annotate each sample as either Apibacter present or absent
samps <- unique(avc.melt$Sample)
#prepare vector to populate with classifications
class <- vector(length = length(samps))
#for each rpm table, row 1 = apibacter, row 2 = crithidia
for (i in 1:length(class)){
  a <- avc.rpm[1, names(avc.rpm) == paste(samps[i])]
  c <- avc.rpm[2, names(avc.rpm) == paste(samps[i])]
  if (a > 0){
    class[i] <- "Apibacter\nPresent"
  }
  if (a == 0){
    class[i] <- "Apibacter\nAbsent"
  }
}
classkey <- data.frame(Sample = samps, Composition = class)

#add composition classification and microphylo information
for(i in 1:nrow(avc.melt)){
  avc.melt$Composition[i] <- classkey$Composition[classkey$Sample == avc.melt$Sample[i]]
  avc.melt$MicroPhylo[i] <- unique(genkey$PhyloHit[genkey$GenHitID == avc.melt$MicroID[i]])
}

myPal <- brewer.pal(4, "Dark2")

ggplot(data = avc.melt[avc.melt$MicroID == "GRH434",], 
       aes(x = Composition, y = log(rPM)))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "black"))

ggsave("output/Plots/CorePhylo/ApibactervsCrithidia_AllBees.pdf")

myPal <- brewer.pal(4, "Dark2")
myPal <- myPal[c(2,1,3,4)]

ggplot(data = avc.melt[avc.melt$MicroID == "GRH434",], 
       aes(x = Composition, y = log(rPM), colour = Genus, fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("black")) +
  scale_colour_manual(values=myPal)

ggsave("output/Plots/CorePhylo/ApibactervsCrithidia_AllBees_IncludeGen.pdf")

#look at just corbiculates
#there only Osmia that's in here with corbiculates anyway
avc.melt2 <- subset(avc.melt, Genus != "Osmia")

ggplot(data = avc.melt2, aes(x = Composition, y = log(rPM), fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "black"))

ggsave("output/Plots/CorePhylo/ApibactervsCrithidia_CorbBees.pdf")


ggplot(data = avc.melt2, aes(x = Composition, y = log(rPM), fill = MicroPhylo, colour = Genus))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "black")) +
  scale_colour_manual(values=myPal)

ggsave("output/Plots/CorePhylo/ApibactervsCrithidia_CorbBees_IncludeGen.pdf")

#what about just numbers of cases of crithidia when apibacter is present versus when
#its absent ? 
avc.melt2$Inc <- ifelse(avc.melt2$rPM > 0, 1, 0)

avc.melt2 %>%
  filter(MicroID == "GRH434") %>%
  group_by(Composition, Genus) %>%
  summarise(sum(Inc))

#####Apilactobacillus: An overview#####
#make plottable dataframe
apil <- genkey$MicroHitID[genkey$PhyloHit == "Apilactobacillus"]

#make dataframe of apil microbes versus corbiculate bees
apil.melt <- plotGen0(microbes = apil,
                      samples = corb,
                      cntmat = spec.inc)

#plot
plotPrev(apil.melt, apil, "Tribe", "micro")
ggsave("output/Plots/CorePhylo/Apilactobacillus_microSpecvCorbSpec_Prev.pdf")

#again, but with Apilactobacillus in general, not broken down by species
apil2 <- unique(genkey$GenHitID[genkey$PhyloHit == "Apilactobacillus"])
apil.melt2 <- plotGen0(microbes = apil2,
                       samples = corb,
                       cntmat = cnt.inc)

plotPrev(apil.melt2, apil2, "Genus", "genus")
ggsave("output/Plots/CorePhylo/Apilactobacillus_genusvCorb_Prev.pdf")

#####Apilactobacillus and Disease Microbes#####
#looking at rpm
#get disease genus ids
dis <- c("Paenibacillus", "Nosema", "Melissococcus")
disID <- unique(genkey$GenHitID[genkey$PhyloHit %in% dis])
#there is no melissococcus post filtering, and no P. larvae within the Paenibacillus
#breakdown
#apilactobacillus id = GRH465
avd <- c(disID, "GRH465")

#extract rpm
avd.rpm <- cnt.rpm[rownames(cnt.rpm) %in% avd,]
#remove samples with neither present
avd.rpm <- avd.rpm[,! colSums(avd.rpm) == 0]

#make plottable
avd.rpm$MicroID <- rownames(avd.rpm)
tmp <- melt(avd.rpm)
avd.melt <- merge(tmp, met, by.x = "variable", by.y = "Sample.ID")
names(avd.melt)[c(1,3)] <- c("Sample", "rPM")

#annotate each sample as either Apibacter present or absent
samps <- unique(avd.melt$Sample)
#prepare vector to populate with classifications
class <- vector(length = length(samps))
#for each rpm table, row 1 = nosema, row 2 = Paenibacillus, row 3 = apilactobacillus
for (i in 1:length(class)){
  a <- avd.rpm[3, names(avd.rpm) == paste(samps[i])]
  if (a == 0){
    class[i] <- "Apilactobacillus\nAbsent"
  }
  if (a > 0){
    class[i] <- "Apilactobacillus\nPresent"
  }
}
classkey <- data.frame(Sample = samps, Composition = class)

#add composition classification and microphylo information
for(i in 1:nrow(avd.melt)){
  avd.melt$Composition[i] <- classkey$Composition[classkey$Sample == avd.melt$Sample[i]]
  avd.melt$MicroPhylo[i] <- unique(genkey$PhyloHit[genkey$GenHitID == avd.melt$MicroID[i]])
}



ggplot(data = avd.melt[!avd.melt$MicroPhylo == "Apilactobacillus",],
       aes(x = Composition, y = log(rPM), fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa"))

ggsave("output/Plots/CorePhylo/ApilactobacillusvsDiseases_AllBees.pdf")

myPal <- brewer.pal(6, "Dark2")

avd.melt2 <- subset(avd.melt, MicroPhylo != "Apilactobacillus")

ggplot(data = avd.melt2,
       aes(x = Composition, y = log(rPM), fill = MicroPhylo, colour = Genus))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa", "pink")) +
  scale_colour_manual(values=myPal)



ggsave("output/Plots/CorePhylo/AplactobacillusvsDisease_AllBees_IncludeGen.pdf")

#look at just corbiculates
#only apis and bombus have any apilactobacillus out of corbiculate species
avd.melt3 <- subset(avd.melt, Genus == "Apis" | Genus == "Bombus")


ggplot(data = subset(avd.melt3, MicroPhylo != "Apilactobacillus"),
       aes(x = Composition, y = log(rPM), fill = MicroPhylo))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa", "pink"))

ggsave("output/Plots/CorePhylo/ApilactobacillusvsDiseases_CorbBees.pdf")


ggplot(data = subset(avd.melt3, MicroPhylo != "Apilactobacillus"),
       aes(x = Composition, y = log(rPM), fill = MicroPhylo, colour = Genus))+
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle =69, hjust = 1),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ,
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#fab387", "#87cefa", "pink")) +
  scale_colour_manual(values=myPal)


ggsave("output/Plots/CorePhylo/ApilactobacillusvsDisease_CorbBees_IncludeGen.pdf")

#what about just numbers of cases of diseases when apilactobacillus is present 
#versus when its absent ? 
avd.melt3$Inc <- ifelse(avd.melt3$rPM > 0, 1, 0)

avd.df <- as.data.frame(avd.melt2 %>%
  group_by(Composition, Genus, MicroPhylo) %>%
  summarise(sum(Inc)))
#change a column header that has a function in it
names(avd.df)[4] <- "Inc"
#remove Apilactobacillus
avd.df <- subset(avd.df, MicroPhylo != "Apilactobacillus")

myPal <- myPal[c(2,1,3,4)]

#plot
ggplot(data = avd.df,
       aes(x = MicroPhylo, y = Inc, fill = Genus)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Composition) +
  scale_fill_manual(values = myPal) +
  theme(axis.text.x = element_text(face = "italic"),
        strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) + 
  labs(x = "Pathogen",
       y = "Number of Samples\nWhere Microbe was found")

ggsave("output/Plots/CorePhylo/ApilactobacillusvsDisease_Inc_Corbs.pdf")

#repeated but with all bees
avd.melt2$Inc <- ifelse(avd.melt2$rPM > 0, 1, 0)

avd.df <- as.data.frame(avd.melt2 %>%
                          group_by(Composition, Genus, MicroPhylo) %>%
                          summarise(sum(Inc)))
#change a column header that has a function in it
names(avd.df)[4] <- "Inc"
#remove Apilactobacillus
avd.df <- subset(avd.df, MicroPhylo != "Apilactobacillus")
#remove samples that only had Apilactobacillus and no pathogens
avd.df <- subset(avd.df, Genus == "Apis" | Genus == "Bombus" |
                   Genus == "Euglossa" | Genus == "Exoneura")

myPal <- brewer.pal(6, "Dark2")
myPal <- myPal[c(2,1,3,4,6,5)]

#plot
ggplot(data = avd.df,
       aes(x = MicroPhylo, y = Inc, fill = Genus)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Composition) +
  scale_fill_manual(values = myPal) +
  theme(axis.text.x = element_text(face = "italic"),
        strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) + 
  labs(x = "Pathogen",
       y = "Number of Samples\nWhere Microbe was found")

ggsave("output/Plots/CorePhylo/ApilactobacillusvsDisease_Inc_Corbs.pdf")

#####PolyMorphic Species and Apilactobacillus####
#Apilactobacillus occurs in <5% of Polymorphic species. Which species?
#extract sample ids of polymorphic species,
#make melted dataframe, plot prevalence for just apilactobacillus
pol2 <- met$Sample.ID[met$Sociality == "Polymorphic"]
pol.melt <- plotGen0(apil2, pol2, cnt.inc)
plotPrev(pol.melt, apil2, "Species", "genus")

ggsave("output/Plots/CorePhylo/Apilactobacillus_PolymorphicPrev.pdf")

#again but with different apilactobacillus species
pol.melt2 <- plotGen0(apil, pol2, spec.inc)
length(apil)
length(pol2)
plotPrev(pol.melt2, apil, "Genus", "micro")

ggsave("output/Plots/CorePhylo/Apilactobacillus_Polymorphic_PrevByMicroSpec.pdf")

pol.melt[pol.melt$Count > 0,]
met[met$Project.Accession == "PRJNA374528",]

#####Snod, Bombi, Frisc in Solitary Species####
#extract solitary sample IDs
sol2 <- met$Sample.ID[met$Sociality == "Solitary"]
#extract targeted microbe ids
sbf <- unique(genkey$GenHitID[genkey$PhyloHit == "Snodgrassella" |
                         genkey$PhyloHit == "Frischella" |
                         genkey$PhyloHit == "Bombilactobacillus"])
sbf.melt <- plotGen0(sbf, sol2, cnt.inc)
plotPrev(sbf.melt, sbf, "Species", "genus")

sbf.melt$Sample[sbf.melt$Count > 0]

#####Bifidobacterium in Solitary Species ####
bifi <- unique(genkey$GenHitID[genkey$PhyloHit == "Bifidobacterium"])
bifsol.melt <- plotGen0(bifi, sol2, cnt.inc)
plotPrev(bifsol.melt, bifi, "Genus", "genus")
length(sol2)

bifsol.melt %>%
  filter(Count > 0)

#what actual bifidobacterium species ? 
bifsol.melt2 <- plotGen0(bifo, sol2, spec.inc)
plotPrev(bifsol.melt2, bifo, "Genus", "micro")

#####Sessionlog#####
dir.create("SessionLogs/")
writeLines(capture.output(sessionInfo()),
           "SessionLogs/CorePhylo_Dec22.txt")
