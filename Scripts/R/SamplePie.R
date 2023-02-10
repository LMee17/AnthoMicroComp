#24th January 2023
#Updated 10th January 2023: New Social Cat Names
#Sample pie charts

library(tidyverse)

dir.create("output/SampleStats/")

met <- read.table("input/Metadata/SampleMetaData_Edit_Final_Feb23.tsv",
                  header = T, sep = "\t")

##Sociality ####

myPal <- 	c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
            "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
            "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
            "#9acd32")

met %>%
  select(Sociality) %>%
  group_by(Sociality) %>%
  count() %>%
  mutate(Sociality = factor(Sociality)) %>%
  mutate(Sociality = fct_relevel(Sociality, c("O. Eusocial", "F. Eusocial", "Solitary"))) %>%
  ggplot(aes(x = "", y = n, fill = Sociality)) +
  geom_bar(stat = "identity", width = 1, colour = "black") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = myPal) +
  theme_void() +
  theme(legend.title = element_text(face = "bold"))
ggsave("output/SampleStats/Sociality_PieChart.pdf")

##Family, Tribe, Genus(?)

met %>%
  select(Family) %>%
  group_by(Family) %>%
  count() %>%
  ggplot(aes(x = "", y = n, fill = Family)) +
  geom_bar(stat = "identity", width = 1, colour = "black") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = myPal) +
  theme_void() +
  theme(legend.title = element_text(face = "bold"))
ggsave("output/SampleStats/Family_PieChart.pdf")

met %>%
  select(Tribe) %>%
  group_by(Tribe) %>%
  count() %>%
  ggplot(aes(x = "", y = n, fill = Tribe)) +
  geom_bar(stat = "identity", width = 1, colour = "black") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = myPal) +
  theme_void() +
  theme(legend.title = element_text(face = "bold"))
ggsave("output/SampleStats/Tribe_PieChart.pdf")

##Location
met %>%
  select(Continent) %>%
  group_by(Continent) %>%
  count() %>%
  ggplot(aes(x = "", y = n, fill = Continent)) +
  geom_bar(stat = "identity", width = 1, colour = "black") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = myPal) +
  theme_void() +
  theme(legend.title = element_text(face = "bold"))
ggsave("output/SampleStats/Continent_PieChart.pdf")

