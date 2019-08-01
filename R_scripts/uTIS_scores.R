#Imports
library(tidyverse)
library(plyr)

#set working directory and home
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/uORF/panels')
home <- 'N:\\JWALDRON/Structure_seq/MCF7_2015'

#setwd('\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Paper/Figures/R/uORF/panels')
#home <- '\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/MCF7_2015'

#set posterior probability thresholds----
positive_change <- "0.2"
no_change <- "0.02"

#write a function to make a p value label from the output of either t.test or wilcox.test including 95% confidence intervals
myP <- function(x) {
  p <- as.numeric(x$p.value)
  if (p < 0.001) {
    rounded_p <- formatC(p, format = "e", digits = 2)
  }else{
    rounded_p <- round(p, digits = 3)
  }
  lower_int <- round(x$conf.int[1], digits = 5)
  lower_int <- formatC(lower_int, format = "f", digits = 5)
  upper_int <- round(x$conf.int[2], digits = 5)
  upper_int <- formatC(upper_int, format = "f", digits = 5)
  p_label = paste0('P = ', rounded_p, "\n95% conf int:\n", lower_int, "  ", upper_int)
  return(p_label)
}

#read in transcript ID lists used for positional changes plots
#fourAdep transcripts
IDs_list <- list()
for (pp in c(0.2, 0.5, 0.8)) {
  df <- read.table(file = file.path(home, paste0('raw_data/transcriptIDs/100nt_filtered_transcripts/fourAdep_pp', pp, '_transcripts.txt')))
  df$pp <- rep(pp, nrow(df))
  IDs_list[[as.character(pp)]] <- df
}

#fourAindep transcripts
for (pp in c(0.01, 0.02, 0.05)) {
  df <- read.table(file = file.path(home, paste0('raw_data/transcriptIDs/100nt_filtered_transcripts/fourAindep_pp', pp, '_transcripts.txt')))
  df$pp <- rep(pp, nrow(df))
  IDs_list[[as.character(pp)]] <- df
}

IDs <- do.call("rbind", IDs_list)

MCF7_IDs <- read_csv(file = "N:\\JWALDRON/Indexes/MCF7/2015/transcript_gene_maps/MCF7_2015_ensembl_IDs.csv", col_names = T)

uORF_scores <- read_csv(file = "N:\\JWALDRON/Structure_seq/refseq/playground/uORF_scores.csv")


#filter to include chosen posterior probability thresholds
IDs %>%
  as.tibble() %>%
  dplyr::rename(transcript = V1) %>%
  filter(pp == no_change | pp == positive_change) %>%
  mutate(translation = case_when(pp == positive_change ~ "4A-dep",
                                 pp == no_change ~ "4A-indep")) %>%
  inner_join(MCF7_IDs, by = "transcript") %>%
  inner_join(uORF_scores, by = "Ensembl_gene_symbol") -> plot_data


t <- wilcox.test(data = plot_data, uORF_score ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

violin_plot <- ggplot(data = plot_data, aes(x = translation, y = uORF_score, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab('uTIS score')+
  ylim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape = 16, size = 6)
pdf(file = "fourAdep_uORF_score_violin.pdf", width = 4, height = 4)
print(violin_plot)
dev.off()

density_plot <- ggplot(data = plot_data, aes(x = uORF_score, colour = factor(translation)))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('uTIS score')+
  xlim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)
pdf(file = "fourAdep_uORF_score_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

plot_data <- ddply(plot_data, .(translation), transform, ecd=ecdf(uORF_score)(uORF_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uORF_score))+
  stat_ecdf(aes(colour = translation), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#74add1", "#fdae61"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))
pdf(file = "fourAdep_uORF_score_cdf.pdf", width = 4, height = 4)
print(cdf_plot)
dev.off()

#ingolia data----
high_sensitivity <- read_csv(file = "N:\\JWALDRON/Ingolia/nature17978-s2A.csv", col_names = T)
low_sensitivity <- read_csv(file = "N:\\JWALDRON/Ingolia/nature17978-s2B.csv", col_names = T)

uORF_scores %>%
  dplyr::rename(Gene = Ensembl_gene_symbol) %>%
  inner_join(high_sensitivity, by = "Gene") %>%
  mutate(translation = rep("high\nsensitivity")) -> high_sensitivity_merged

uORF_scores %>%
  dplyr::rename(Gene = Ensembl_gene_symbol) %>%
  inner_join(low_sensitivity, by = "Gene") %>%
  mutate(translation = rep("low\nsensitivity")) -> low_sensitivity_merged

plot_data <- bind_rows(high_sensitivity_merged, low_sensitivity_merged)

t <- wilcox.test(data = plot_data, uORF_score ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

violin_plot <- ggplot(data = plot_data, aes(x = translation, y = uORF_score, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab('uTIS score')+
  ylim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape = 16, size = 6)
pdf(file = "ingolia_uORF_score_violin.pdf", width = 4, height = 4)
print(violin_plot)
dev.off()

density_plot <- ggplot(data = plot_data, aes(x = uORF_score, colour = factor(translation)))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('uTIS score')+
  xlim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)
pdf(file = "ingolia_uORF_score_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

plot_data <- ddply(plot_data, .(translation), transform, ecd=ecdf(uORF_score)(uORF_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uORF_score))+
  stat_ecdf(aes(colour = translation), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#74add1", "#fdae61"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))
pdf(file = "ingolia_uORF_score_cdf.pdf", width = 4, height = 4)
print(cdf_plot)
dev.off()

#low TE vs high TE----
low_TE <- read.table(file = file.path(home, 'raw_data/transcriptIDs/100nt_filtered_transcripts/low_TE_transcriptIDs.txt'))
low_TE %>%
  mutate(TE = rep("low TE", nrow(low_TE))) -> low_TE

high_TE <- read.table(file = file.path(home, 'raw_data/transcriptIDs/100nt_filtered_transcripts/high_TE_transcriptIDs.txt'))
high_TE %>%
  mutate(TE = rep("high TE", nrow(high_TE))) -> high_TE

low_TE %>%
  bind_rows(high_TE) %>%
  as.tibble() %>% 
  dplyr::rename(transcript = V1) %>%
  inner_join(MCF7_IDs, by = "transcript") %>%
  inner_join(uORF_scores, by = "Ensembl_gene_symbol") -> plot_data


t <- wilcox.test(data = plot_data, uORF_score ~ TE,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

violin_plot <- ggplot(data = plot_data, aes(x = factor(TE, levels = c("low TE", "high TE"), ordered = T), y = uORF_score, fill = factor(TE, levels = c("low TE", "high TE"), ordered = T)))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#7C71D8", "#FFE073"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab('uTIS score')+
  ylim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape = 16, size = 6)
pdf(file = "TE_uORF_score_violin.pdf", width = 4, height = 4)
print(violin_plot)
dev.off()

density_plot <- ggplot(data = plot_data, aes(x = uORF_score, colour = factor(TE)))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#7C71D8", "#FFE073"))+
  scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
  xlab('uTIS score')+
  xlim(c(0, 1))+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))+
  ggtitle(p_label)
pdf(file = "TE_uORF_score_density.pdf", width = 4, height = 4)
print(density_plot)
dev.off()

plot_data <- ddply(plot_data, .(TE), transform, ecd=ecdf(uORF_score)(uORF_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uORF_score))+
  stat_ecdf(aes(colour = TE), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#7C71D8", "#FFE073"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  theme_bw()+
  theme(legend.position='none', 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))
pdf(file = "TE_uORF_score_cdf.pdf", width = 4, height = 4)
print(cdf_plot)
dev.off()

