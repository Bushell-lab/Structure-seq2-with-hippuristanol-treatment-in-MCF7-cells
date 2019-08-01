###Imports
library(gridExtra)
library(grid)
library(tidyverse)

setwd('\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Paper/Figures/R/Retic/panels')

my_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_blank())

#read in data----
data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/JWALDRON/Beatson_translation_assays/test/treated_retic_RNA_titration_results.csv", col_names = T)

n <- 3

#sort data----
data %>%
  group_by(Concentration) %>%
  summarise(mean_Fluc = mean(Fluc),
            SD_Fluc = sd(Fluc),
            SE_Fluc = (sd(Fluc)) / sqrt(n)) %>%
  ungroup() -> Fluc

#plot data----
Fluc %>%
  ggplot(aes(x = Concentration, y = mean_Fluc))+
  geom_line(size=1)+
  geom_errorbar(aes(x = Concentration, ymin = mean_Fluc - SE_Fluc, ymax = mean_Fluc + SE_Fluc), size = 0.5, width = 5)+
  ylab("Fluc")+
  scale_x_continuous(limits = c(20, 360), breaks = c(20,40,80,160,320), trans='log2')+
  scale_y_continuous(labels = scales::comma)+
  xlab("RNA concentration (ng/ul)")+
  my_theme -> CAA_Fluc_plot
pdf(file = "RNA_titration.pdf", width = 6, height = 4)
CAA_Fluc_plot
dev.off()

#harringtonine----
#read in data----
data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Translation_assays/retic/harringtonine/harringtonine_results.csv", col_names = T)

n <- 3

#sort data----
data %>%
  group_by(Harringtonine, Time) %>%
  summarise(mean_Fluc = mean(Fluc),
            SD_Fluc = sd(Fluc),
            SE_Fluc = (sd(Fluc)) / sqrt(n)) %>%
  ungroup() -> Fluc

Fluc %>%
  ggplot(aes(x = Time, y = mean_Fluc, colour = factor(Harringtonine)))+
  geom_line(size = 1)+
  geom_errorbar(aes(x = Time, ymin = mean_Fluc - SE_Fluc, ymax = mean_Fluc + SE_Fluc, colour = factor(Harringtonine), width = 0.1))+
  scale_y_continuous(labels = scales::comma)+
  xlim(c(10,20))+
  ylab("Fluc")+
  xlab("Time (min)")+
  my_theme+
  theme(legend.text = element_text(size = 20)) -> Har_plot

pdf(file = "harringtonine_titration.pdf", width = 8, height = 4)
Har_plot
dev.off()

#hipp----
#read in data----
data <- read_csv(file = '\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Translation_assays/retic/CAA_GAA/CAA_Hipp_results.csv', col_names = T)

n <- 4

#sort data----
data %>%
  filter(Reporter == "CAA") %>%
  group_by(Hipp_conc) %>%
  summarise(mean_Fluc = mean(Fluc),
            SD_Fluc = sd(Fluc),
            SE_Fluc = (sd(Fluc)) / sqrt(n)) %>%
  ungroup() -> Fluc

Fluc %>%
  ggplot(aes(x = Hipp_conc, y = mean_Fluc))+
  geom_line(size=1)+
  geom_errorbar(aes(x = Hipp_conc, ymin = mean_Fluc - SE_Fluc, ymax = mean_Fluc + SE_Fluc), size = 0.5, width = 0.1)+
  ylab("Fluc")+
  scale_y_continuous(labels = scales::comma)+
  xlab("Hippuristanol concentration (uM)")+
  ggtitle("CAA reporters")+
  my_theme -> Hipp_plot

pdf(file = "hipp_titration.pdf", width = 6, height = 4)
Hipp_plot
dev.off()
