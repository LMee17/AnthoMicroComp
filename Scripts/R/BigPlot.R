#27th January 2023
#Combine NMDS plots for chapter

library(tidyverse)
library(ggpubr)
library("RColorBrewer")         

#unfortunately, I can't loop these. Not just because some need individual adjustmments
#for colours and centroid etc, but also I can only have every third plot have a legend
#and ggarange won't accept concatenated lists. So expect A LOT of repeated code.

pro.n <- read.table("output/Prokaryote/Composition_Analysis/NMDS_meta.tsv",
                    header = T, sep = "\t")
euk.n <- read.table("output/Eukaryote/Composition_Analysis/NMDS_meta.tsv",
                    header = T, sep = "\t")
vir.n <- read.table("output/Viral/Composition_Analysis/NMDS_meta.tsv",
                    header = T, sep = "\t")

myPal <- 	c("#db6d00", "#009292", "#3b3bc4", "#bb00bb", "#920000",
            "#000000", "#ff6db6", "#6db6ff", "#24ff24", "#00e6e6",
            "#ffb6db", "#b66dff", "#924900", "#b6dbff", "#ffff6d",
            "#9acd32")

famPal <- c("#db7093", "#00bfff", "#e6c000", "green")
conPal <- brewer.pal(6, "Dark2")

factors <- c("Sociality", "Continent", "Family")

all <- vector(mode = "list", length = 9) # for centroid versions
all2 <- all #for no centroid versions
#plot 1. viral sociality
cenvs <- vir.n %>%
  group_by(Sociality) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame()
all[[1]] <- ggplot(data = vir.n, #check data 
       aes(x = NMDS1, y = NMDS2, colour = vir.n$Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal[c(2,3)]) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenvs, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenvs[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal[c(1,3)]) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[1]] <- ggplot(data = vir.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = vir.n$Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal[c(2,3)]) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

#plot 2. prokaryote sociality
cenps <- pro.n %>%
  group_by(Sociality) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame()
all[[2]] <- ggplot(data = pro.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = pro.n$Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenps, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenps[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[2]] <- ggplot(data = pro.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())


#plot 3. eukaruote sociality with guide
cenes <- euk.n %>%
  group_by(Sociality) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame()
all[[3]] <- ggplot(data = euk.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenes, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenes[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = myPal) +
  labs(colour = "Sociality") +
  theme_classic() +
  theme(strip.background = element_blank())+
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))


all2[[3]] <- ggplot(data = euk.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Sociality)) + #check factor
  geom_point() + 
  scale_colour_manual(values = myPal) + #check colour scheme
  labs(colour = "Sociality") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))


#plot 4. viral Family
cenvf <- vir.n %>%
  group_by(Family) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame() %>%
  subset(Family != "Andrenidae")
all[[4]] <- ggplot(data = vir.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal[c(1:2,4)]) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenvf, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenvf[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = famPal[c(2,4)]) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[4]] <- ggplot(data = vir.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal[c(1:2,4)]) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

#plot 5. prokaryote Family
cenpf <- pro.n %>%
  group_by(Family) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame()
all[[5]] <- ggplot(data = pro.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenpf, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenpf[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = famPal) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[5]] <- ggplot(data = pro.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())


#plot 6. eukaruote Family with guide
cenef <- euk.n %>%
  group_by(Family) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame() %>%
  subset(Family != "Andrenidae")
all[[6]] <- ggplot(data = euk.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenef, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenef[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = famPal[c(2:4)]) +
  labs(colour = "Family") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))


all2[[6]] <- ggplot(data = euk.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Family)) + #check factor
  geom_point() + 
  scale_colour_manual(values = famPal) + #check colour scheme
  labs(colour = "Family") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))


#plot 4. viral Continent
cenvc <- vir.n %>%
  group_by(Continent) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame() %>%
  subset(Continent != "South America")
all[[7]] <- ggplot(data = vir.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenvc, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenvc[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = conPal) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[7]] <- ggplot(data = vir.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

#plot 5. prokaryote Continent
cenpc <- pro.n %>%
  group_by(Continent) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame()
all[[8]] <- ggplot(data = pro.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenpc, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenpc[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = conPal) +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())

all2[[8]] <- ggplot(data = pro.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(strip.background = element_blank())


#plot 6. eukaruote Continent with guide
cenec <- euk.n %>%
  group_by(Continent) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2)) %>%
  as.data.frame() %>%
  subset(Continent != "South America")
all[[9]] <- ggplot(data = euk.n, #check data 
                   aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  stat_ellipse(show.legend = FALSE) +
  geom_point(data = cenec, size = 3, 
             shape = 21, colour = "black", 
             aes(fill = cenec[,1]),
             show.legend = FALSE) + 
  scale_fill_manual(values = conPal) +
  labs(colour = "Continent") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))


all2[[9]] <- ggplot(data = euk.n, #check data 
                    aes(x = NMDS1, y = NMDS2, colour = Continent)) + #check factor
  geom_point() + 
  scale_colour_manual(values = conPal) + #check colour scheme
  labs(colour = "Continent") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 3, byrow = T))

all.v2 <- all[c(1,4,7,2,5,8,3,6,9)]
all2.v2 <- all2[c(1,4,7,2,5,8,3,6,9)]
#arrange
ggarrange(plotlist = all.v2, heights = c(1,1,1.75))
ggsave("All_Centroid_Ellipsis_1guide_v2.pdf")
ggarrange(plotlist = all2.v2, heights = c(1,1,1.75))
ggsave("All_Simple_1guide_v2.pdf")

#just sociality
soc <- all[c(1,2,3)]
ggarrange(plotlist = soc, heights = c(1,1,1.75), ncol = 1)

theme_get()
