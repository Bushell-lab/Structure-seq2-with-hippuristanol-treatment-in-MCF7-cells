###This script was written by Joseph A. Waldron and produces panels 3B-E in Waldron et al. (2019) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#Imports----
library(grid)
library(gridExtra)
library(tidyverse)

#import variables----
source("Structure_seq_variables.R")

#write functions----
#makes a label from the output of either t.test or wilcox.test to include the p value and 95% confidence limits
myP <- function(x) {
  p <- as.numeric(x$p.value)
  if (p < 2.2e-16) {
    p_label <- "P < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("P =", rounded_p)
  }
  lower <- round(x$conf.int[1], digits = 3)
  upper <- round(x$conf.int[2], digits = 3)
  if (lower == 0 | upper == 0) {
    lower <- round(x$conf.int[1], digits = 4)
    upper <- round(x$conf.int[2], digits = 4)
    if (lower == 0 | upper == 0) {
      lower <- round(x$conf.int[1], digits = 5)
      upper <- round(x$conf.int[2], digits = 5)
      if (lower == 0 | upper == 0) {
        lower <- round(x$conf.int[1], digits = 6)
        upper <- round(x$conf.int[2], digits = 6)
      }
    }
  }
  return(paste0(p_label, "\n95% conf int:\n", lower, "   ", upper))
}

#exports just the legend of a plot
myLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#make theme----
myTheme <- theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#load common data
source("Structure_seq_common_data.R")

#G4 predictions
G4_data <- read_tsv(file = "fpUTR_G4_screener.tsv", col_names = T) #generate with Custom_scripts.sh

#sort data----
#the following pipe mutates G and C contenet into percentages by multiplying by 100, merges data, filters it to include only 4A-dep and 4A-indep transcripts
#and then selects the most abundant transcript per gene and then the maximum G4NN score per transcript
FASTA_compositions_list$fpUTR %>%
  mutate(G_content = G_content * 100,
         C_content = C_content * 100) %>%
  inner_join(G4_data, by = c("transcript" = "description")) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_data, by = "gene") %>%
  inner_join(abundance_data, by = "transcript") %>%
  filter(translation == "4A-dep" | translation == "4A-indep") %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  top_n(wt = G4NN, n = 1) %>%
  ungroup() -> merged_data

#calculate the number of 4A-dep transcripts and then select the same numebr of 4A-indep transcripts based on the lowest posterior probability
fourAdep_transcript_n <- n_distinct(merged_data$transcript[merged_data$translation == "4A-dep"])

merged_data %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = fourAdep_transcript_n) %>%
  ungroup() -> merged_data 

summary(merged_data$translation)

#plot data----
#length
t <- wilcox.test(data = merged_data, length ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

fourAdep_mean <- mean(merged_data$length[merged_data$translation == "4A-dep"])
fourAindep_mean <- mean(merged_data$length[merged_data$translation == "4A-indep"])

density_plot <- ggplot(data = merged_data, aes(x = length, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('length')+
  scale_x_log10(breaks=c(10,100,1000,10000),limits=c(10, 10000))+
  geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
  myTheme+
  ggtitle(p_label)
pdf(file = "fpUTR_length_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

#C_content
t <- wilcox.test(data = merged_data, C_content ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

fourAdep_mean <- mean(merged_data$C_content[merged_data$translation == "4A-dep"])
fourAindep_mean <- mean(merged_data$C_content[merged_data$translation == "4A-indep"])

density_plot <- ggplot(data = merged_data, aes(x = C_content, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('C content')+
  xlim(c(0,50))+
  geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
  myTheme+
  ggtitle(p_label)
pdf(file = "fpUTR_C_content_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

#G_content
t <- wilcox.test(data = merged_data, G_content ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

fourAdep_mean <- mean(merged_data$G_content[merged_data$translation == "4A-dep"])
fourAindep_mean <- mean(merged_data$G_content[merged_data$translation == "4A-indep"])

density_plot <- ggplot(data = merged_data, aes(x = G_content, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('G content')+
  xlim(c(0,50))+
  geom_vline(xintercept = fourAdep_mean, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = fourAindep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
  myTheme+
  ggtitle(p_label)
pdf(file = "fpUTR_G_content_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

#G4NN scores
t <- wilcox.test(data = merged_data, G4NN ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

fourAdep_mean <- mean(merged_data$G4NN[merged_data$translation == "4A-dep"])
fourAindep_mean <- mean(merged_data$G4NN[merged_data$translation == "4A-indep"])

density_plot <- ggplot(data = merged_data, aes(x = G4NN, colour = translation, fill = translation))+
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

