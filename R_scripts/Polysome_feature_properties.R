###Imports
library(grid)
library(gridExtra)
library(tidyverse)

#import functions----
source("N:\\JWALDRON/R_scripts/functions.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/functions.R")

#import variables----
source("N:\\JWALDRON/R_scripts/structure_seq_variables.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/structure_seq_variables.R")

#set home and working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Polysomes/panels')
#setwd('\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Paper/Figures/R/Polysomes/panels')
#home <- '\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/MCF7_2015'

#make theme----
myTheme <- theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#common data
source("N:\\JWALDRON/R_scripts/Structure_seq_common_data.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/Structure_seq_common_data.R")

FASTA_compositions %>%
  mutate(GC_content = GC_content *100,
         GA_content = GA_content *100,
         G_content = G_content *100,
         A_content = A_content *100,
         C_content = C_content *100,
         T_content = T_content *100) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_data, by = "gene") %>%
  filter(translation == "4A-dep" | translation == "4A-indep") %>%
  inner_join(abundance_data, by = "transcript") %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() -> merged_data

fourAdep_transcript_n <- n_distinct(merged_data$transcript[merged_data$translation == "4A-dep"])

merged_data %>%
  group_by(translation, region) %>%
  top_n(wt = -posterior_probability, n = fourAdep_transcript_n) %>%
  ungroup() -> merged_data 

summary(merged_data$translation[merged_data$region=="fpUTR"])
summary(merged_data$translation[merged_data$region=="CDS"])
summary(merged_data$translation[merged_data$region=="tpUTR"])

#plot data----
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- merged_data[merged_data$region == region,]
  
  #plot length
  t <- wilcox.test(data = df, length ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$length[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$length[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = length, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('length')+
    scale_x_log10(breaks=c(10,100,1000,10000),limits=c(10, 10000))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_length_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot GC_content
  t <- wilcox.test(data = df, GC_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$GC_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$GC_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = GC_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('GC content')+
    xlim(c(25,100))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_GC_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot GA_content
  t <- wilcox.test(data = df, GA_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$GA_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$GA_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = GA_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('GA content')+
    xlim(c(25,100))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_GA_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot G_content
  t <- wilcox.test(data = df, G_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$G_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$G_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = G_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('G content')+
    xlim(c(0,50))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_G_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot C_content
  t <- wilcox.test(data = df, C_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$C_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$C_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = C_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('C content')+
    xlim(c(0,50))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_C_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot A_content
  t <- wilcox.test(data = df, A_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$A_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$A_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = A_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('A content')+
    xlim(c(0,50))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_A_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
  
  #plot T_content
  t <- wilcox.test(data = df, T_content ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_mean <- mean(df$T_content[df$translation == "4A-dep"])
  fourAindep_mean <- mean(df$T_content[df$translation == "4A-indep"])
  
  density_plot <- ggplot(data = df, aes(x = T_content, colour = translation, fill = translation))+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    scale_colour_manual(values=c("#74add1", "#fdae61"))+
    xlab('U content')+
    xlim(c(0,50))+
    geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
    geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
    myTheme+
    ggtitle(p_label)
  pdf(file = paste0(region, "_T_content_density.pdf"), width = 4, height = 4)
  print(density_plot)
  dev.off()
}

#export legend----
density_plot <- ggplot(data = merged_data, aes(x = G_content, colour = translation, fill = translation))+
  geom_density()+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  theme(legend.title = element_blank())

density_legend <- myLegend(density_plot)
pdf(file = 'density_legend.pdf', height = 1, width = 1)
grid.arrange(density_legend)
dev.off()

#calculate GC3 content
library(seqinr)

CDS_FASTA <- read.fasta(file = file.path(home, 'fasta/spliced/MCF7_2015_CDSs.fasta'))

fourAdepIDs <- pull(merged_data[merged_data$translation == "4A-dep" & merged_data$region == "CDS",], transcript)
fourAindepIDs <- pull(merged_data[merged_data$translation == "4A-indep" & merged_data$region == "CDS",], transcript)

CDS_fourAdep_fasta <- CDS_FASTA[names(CDS_FASTA) %in% fourAdepIDs]
CDS_fourAindep_fasta <- CDS_FASTA[names(CDS_FASTA) %in% fourAindepIDs]

#fourAdep
fourAdep_GC1 <- lapply(CDS_fourAdep_fasta, GC1)
fourAdep_GC2 <- lapply(CDS_fourAdep_fasta, GC2)
fourAdep_GC3 <- lapply(CDS_fourAdep_fasta, GC3)

fourAdep_GC1 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC1) %>%
  mutate(translation = rep("4A-dep")) -> fourAdep_GC1

fourAdep_GC2 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC2) %>%
  mutate(translation = rep("4A-dep")) -> fourAdep_GC2

fourAdep_GC3 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC3) %>%
  mutate(translation = rep("4A-dep")) -> fourAdep_GC3

#fourAindep
fourAindep_GC1 <- lapply(CDS_fourAindep_fasta, GC1)
fourAindep_GC2 <- lapply(CDS_fourAindep_fasta, GC2)
fourAindep_GC3 <- lapply(CDS_fourAindep_fasta, GC3)

fourAindep_GC1 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC1) %>%
  mutate(translation = rep("4A-indep")) -> fourAindep_GC1

fourAindep_GC2 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC2) %>%
  mutate(translation = rep("4A-indep")) -> fourAindep_GC2

fourAindep_GC3 %>%
  as.data.frame() %>%
  gather(key = transcript, value = GC3) %>%
  mutate(translation = rep("4A-indep")) -> fourAindep_GC3

GC1_data <- bind_rows(fourAdep_GC1, fourAindep_GC1)
GC2_data <- bind_rows(fourAdep_GC2, fourAindep_GC2)
GC3_data <- bind_rows(fourAdep_GC3, fourAindep_GC3)

#plot
#GC1
t <- wilcox.test(data = GC1_data, GC1 ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

density_plot <- ggplot(data = GC1_data, aes(x = GC1, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('GC1 content')+
  xlim(c(0.25, 0.75))+
  myTheme+
  ggtitle(p_label)
pdf(file = "GC1_content_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

#GC2
t <- wilcox.test(data = GC2_data, GC2 ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

density_plot <- ggplot(data = GC2_data, aes(x = GC2, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('GC2 content')+
  xlim(c(0.25, 0.75))+
  myTheme+
  ggtitle(p_label)
pdf(file = "GC2_content_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

#GC3
t <- wilcox.test(data = GC3_data, GC3 ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

density_plot <- ggplot(data = GC3_data, aes(x = GC3, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('GC3 content')+
  xlim(c(0.25, 0.75))+
  myTheme+
  ggtitle(p_label)
pdf(file = "GC3_content_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()


#G4 predictions----
G4_data <- read_tsv(file = file.path(home, "raw_data/G4_RNA_screener/fpUTR_G4_screener.tsv"), col_names = T)

G4_data %>%
  select(-X1) %>%
  rename(transcript = description) %>%
  filter(transcript %in% merged_data$transcript) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_data, by = "gene") %>%
  group_by(transcript) %>%
  top_n(wt = G4NN, n = 1) -> merged_G4_data

#plot G4NN
t <- wilcox.test(data = merged_G4_data, G4NN ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

fourAdep_mean <- mean(merged_G4_data$G4NN[merged_G4_data$translation == "4A-dep"])
fourAindep_mean <- mean(merged_G4_data$G4NN[merged_G4_data$translation == "4A-indep"])

density_plot <- ggplot(data = merged_G4_data, aes(x = G4NN, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('G4NN score')+
  geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
  myTheme+
  ggtitle(p_label)
pdf(file = "fpUTR_G4NN_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()


