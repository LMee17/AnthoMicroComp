#28th January 2023
#Describe tribe plots

library(tidyverse)
library(ggpubr)

met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Jan23.tsv",
                  header = T, sep = "\t")
pro <- read.table("input/Counts/Prokaryote_Filtered_TribeReduced_Raw.tsv")

met2 <- met[met$Sample.ID %in% names(pro),]

tribes <- sort(unique(met2$Tribe))
noAssociates <- c(1, 8, 4, 5, 2, 3, 4, 2)
n <- c(12, 86, 22, 57, 5, 26, 7, 11)

assoc <- data.frame(Tribe = tribes, NoAssoc = noAssociates, Sample_n = n)
assoc$Sociality <- c("Solitary", "Eusocial", "Polymorphic", "Eusocial", 
                     rep("Polymorphic",2), "Eusocial", "Solitary")
assoc$Sociality <- factor(assoc$Sociality, levels = c("Solitary",
                                                      "Polymorphic",
                                                      "Eusocial"))
myPal <- 	c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
            "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
            "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
            "#9acd32") 

one <- ggplot(assoc, aes(x = Sociality, y = NoAssoc)) + 
  theme_minimal() +
  geom_bar(stat = "identity", aes(fill = Tribe)) + 
  scale_fill_manual(values = c(myPal[3], myPal[1], myPal[2],
                               "#ffa852",
                               "#57ffff",
                               "#005757",
                               "#793c00",
                               "#8686db" )) +
  labs(y = "Number of Bacterial Associates") +
  theme(legend.position = "bottom")+
  guides(guide_legend(nrow = 2))
ggsave("output/Prokaryote/Tribal_Descriptives/SocialityvNoAssociates.pdf")

two <- ggplot(assoc, aes(x = Sample_n, y = noAssociates)) +
  theme_minimal() +
  geom_smooth(method = "glm", se = T, colour = "grey", alpha = 0.3) +
  labs(x = "Number of Samples",
       y = "Number of Bacterial Associates") +
  geom_point(aes(colour = Sociality), size = 3)  +
  scale_colour_manual(values = myPal[c(3,2,1)]) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks=seq(0,8, 2)) + 
  guides(guide_legend(nrow = 1))
ggsave("output/Prokaryote/Tribal_Descriptives/NoSamplesvNoAssociates.pdf")  



ggarrange(one, two,  ncol = 2, align = "h")
ggsave("output/Prokaryote/Tribal_Descriptives/AssociateCombined.pdf",
       height = 10, width = 20, units = "cm")  