###This script was written by Joseph A. Waldron and produces Figure S3 from Waldron et al. (2019) Genome Biology
###Note that the data required for this script is generated following the mapping of the sequencing reads to two additional transcriptomes
###To avoid confusion and difficulties in naming the files appropriately, it was decided not to host the data for this script in this repository
###Please contact us if you require this data as we would be happy to share

#load packages----
library(gridExtra)
library(grid)
library(tidyverse)

#set coverage threshold----
coverage <- 1

#load data----
#binned_rtsc
#the following for loop reads in binned rtsc data across the whole transcript for all transcripts in the 3 different transcriptomes
df_list <- list()
for (transcriptome in c("CAGE", "MCF7_2015", "refseq")) {
  df <- read_csv(file = file.path(transcriptome, "control_minus_DMS_all_binned.csv"), col_names = T)
  df$transcriptome <- as.factor(rep(transcriptome, nrow(df)))
  df$position <- 1:nrow(df)
  df_list[[transcriptome]] <- df
}
binned_data <- do.call("rbind", df_list)

#5' coverage_data
#the following for loop reads in 5' coverage data with increasing lengths for 3 different transcriptomes
df_list <- list()
for (transcriptome in c("CAGE", "MCF7_2015", "refseq")) {
  for (n in c(1:5, 10, 15, 20, 25, 50, 100, 200)) {
    df <- read_csv(file = file.path(transcriptome, paste0('control_minus_DMS_FP_', n, '.csv')), col_names = T)
    df$transcriptome <- as.factor(rep(transcriptome, nrow(df)))
    df$n <- rep(n, nrow(df))
    df_list[[paste(transcriptome, n, sep = "_")]] <- df
  }
}
FP_coverage_data <- do.call("rbind", df_list)

#3' coverage_data
#the following for loop reads in 3' coverage data with increasing lengths for 3 different transcriptomes
df_list <- list()
for (transcriptome in c("MCF7_2015", "refseq")) {
  for (n in c(seq(0,300,25))) {
    df <- read_csv(file = file.path(transcriptome, paste0('control_minus_DMS_TP_50_300_', n, '.csv')), col_names = T)
    df$transcriptome <- as.factor(rep(transcriptome, nrow(df)))
    df$n <- rep(n, nrow(df))
    df_list[[paste(transcriptome, n, sep = "_")]] <- df
  }
}
TP_coverage_data <- do.call("rbind", df_list)

#coverage_data
#the following for loop reads in coverage data for 3 different transcriptomes
df_list <- list()
for (transcriptome in c("CAGE", "MCF7_2015", "refseq")) {
  df <- read_csv(file = file.path(transcriptome, "control_plus_DMS_coverage.csv"), col_names = T)
  df$transcriptome <- as.factor(rep(transcriptome, nrow(df)))
  df_list[[transcriptome]] <- df
}
coverage_data <- do.call("rbind", df_list)

#filter by coverage----
FP_coverage_data %>%
  inner_join(coverage_data, by = c("transcript", "transcriptome")) %>%
  filter(control_plus_DMS_coverage > coverage) -> filtered_FP_coverage_data

TP_coverage_data %>%
  inner_join(coverage_data, by = c("transcript", "transcriptome")) %>%
  filter(control_plus_DMS_coverage > coverage) -> filtered_TP_coverage_data

#calculate medians----
filtered_FP_coverage_data %>%
  group_by(transcriptome, n) %>%
  summarise(FP_coverage_median = median(FP_coverage)) -> summarised_FP_coverage_data

filtered_TP_coverage_data %>%
  group_by(transcriptome, n) %>%
  summarise(TP_coverage_median = median(TP_coverage, na.rm = T)) ->summarised_TP_coverage_data

#plot data----
#binned_data
ylim <- max(binned_data$value)

CAGE_plot <- ggplot(data = binned_data[binned_data$transcriptome == "CAGE",], aes(x = position, y = value))+
  geom_bar(stat='identity', fill = "#619CFF")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none",
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"))+
  scale_y_continuous(limits = c(0,ylim))+
  coord_trans(limx = c(0, 101))

MCF7_2015_plot <- ggplot(data = binned_data[binned_data$transcriptome == "MCF7_2015",], aes(x = position, y = value))+
  geom_bar(stat='identity', fill = "#F8766D")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none",
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"))+
  scale_y_continuous(limits = c(0,ylim))+
  coord_trans(limx = c(0, 101))

refseq_plot <- ggplot(data = binned_data[binned_data$transcriptome == "refseq",], aes(x = position, y = value))+
  geom_bar(stat='identity', fill = "#00BA38")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none",
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"))+
  scale_y_continuous(limits = c(0,ylim))+
  coord_trans(limx = c(0, 101))

pdf(file = 'DMS_minus_binned_stops.pdf', height = 10, width = 20)
grid.arrange(CAGE_plot, MCF7_2015_plot, refseq_plot, nrow = 3)
dev.off()

#5' coverage
CAGE_plot <- ggplot(data = summarised_FP_coverage_data[summarised_FP_coverage_data$transcriptome == "CAGE",], aes(x = n, y = FP_coverage_median))+
  geom_point(size=2, colour = "#619CFF")+ 
  geom_path(size = 1, lineend = "round", colour = "#619CFF")+
  xlab("n")+
  ylab("median 5\' end coverage")+
  theme_bw()+
  theme(axis.text = element_text(colour="grey20", size=14, face="plain"), 
        axis.title = element_text(colour="grey20", size=18, face="plain"), 
        legend.position = 'none')

MCF7_2015_plot <- ggplot(data = summarised_FP_coverage_data[summarised_FP_coverage_data$transcriptome == "MCF7_2015",], aes(x = n, y = FP_coverage_median))+
  geom_point(size=2, colour = "#F8766D")+ 
  geom_path(size = 1, lineend = "round", colour = "#F8766D")+
  xlab("n")+
  ylab("median 5\' end coverage")+
  theme_bw()+
  theme(axis.text = element_text(colour="grey20", size=14, face="plain"), 
        axis.title = element_text(colour="grey20", size=18, face="plain"), 
        legend.position = 'none')

refseq_plot <- ggplot(data = summarised_FP_coverage_data[summarised_FP_coverage_data$transcriptome == "refseq",], aes(x = n, y = FP_coverage_median))+
  geom_point(size=2, colour = "#00BA38")+ 
  geom_path(size = 1, lineend = "round", colour = "#00BA38")+
  xlab("n")+
  ylab("median 5\' end coverage")+
  theme_bw()+
  theme(axis.text = element_text(colour="grey20", size=14, face="plain"), 
        axis.title = element_text(colour="grey20", size=18, face="plain"), 
        legend.position = 'none')

pdf(file = 'FP_coverage.pdf', height = 5, width = 21)
grid.arrange(CAGE_plot, MCF7_2015_plot, refseq_plot, nrow = 1)
dev.off()

#3' coverage
ylim <- max(summarised_TP_coverage_data$TP_coverage_median)
MCF7_2015_plot <- ggplot(data = summarised_TP_coverage_data[summarised_TP_coverage_data$transcriptome == "MCF7_2015",], aes(x = n, y = TP_coverage_median))+
  geom_point(size=2, colour = "#F8766D")+ 
  geom_path(size = 1, lineend = "round", colour = "#F8766D")+
  xlab("trim (nt)")+
  ylim(0, ylim)+
  ylab("median 3\' end coverage")+
  theme_bw()+
  theme(axis.text = element_text(colour="grey20", size=14, face="plain"), 
        axis.title = element_text(colour="grey20", size=18, face="plain"), 
        legend.position = 'none')

refseq_plot <- ggplot(data = summarised_TP_coverage_data[summarised_TP_coverage_data$transcriptome == "refseq",], aes(x = n, y = TP_coverage_median))+
  geom_point(size=2, colour = "#00BA38")+ 
  geom_path(size = 1, lineend = "round", colour = "#00BA38")+
  xlab("trim (nt)")+
  ylim(0, ylim)+
  ylab("median 3\' end coverage")+
  theme_bw()+
  theme(axis.text = element_text(colour="grey20", size=14, face="plain"), 
        axis.title = element_text(colour="grey20", size=18, face="plain"), 
        legend.position = 'none')

pdf(file = 'TP_coverage.pdf', height = 5, width = 14)
grid.arrange(MCF7_2015_plot, refseq_plot, nrow = 1)
dev.off()
