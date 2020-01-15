###This script was written by Joseph A. Waldron and produces panels S6B and S7E-F in Waldron et al. (2019) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#Imports----
library(tidyverse)
library(plyr)

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

#write theme----
cdf_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#load common data
source("Structure_seq_common_data.R")

#GTI scores
#generate with Custom_scripts.sh
GTI_uTIS_scores <- read_csv(file = "GTI_data_uTIS_scores.csv", col_names = T)

#ingolia data
#this is supplemental tables S2A/B from Iwasaki et al. (2016) Nature
high_sensitivity <- read_csv(file = "nature17978-s2A.csv", col_names = T)
low_sensitivity <- read_csv(file = "nature17978-s2B.csv", col_names = T)

#MCF7 gene names
#can be downloaded from the data folder of this repository
MCF7_IDs <- read_csv(file = "MCF7_2015_ensembl_IDs.csv", col_names = T)

#filtered transcript IDs that have 5'UTRs more than 100nt
#download from the data folder within this repository
filtered_transcripts <- read.table(file = "filtered_plus_100_transcripts.txt", header = F)
names(filtered_transcripts) <- "transcript"

#reformat, merged and plot data
#4A-dep transcripts
translation_data %>%
  inner_join(transcript_to_geneID, by = "gene") %>%
  filter(transcript %in% filtered_transcripts$transcript) %>%
  filter(translation == "4A-dep" | translation == "4A-indep") -> translation_df

n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"]) -> n_fourAdep_transcripts

translation_df %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
  ungroup() %>%
  select(transcript, translation) %>%
  inner_join(MCF7_IDs, by = "transcript") %>%
  inner_join(GTI_uTIS_scores, by = "Ensembl_gene_symbol") -> filtered_translation_data

t <- wilcox.test(data = filtered_translation_data, uTIS_score ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

plot_data <- ddply(filtered_translation_data, .(translation), transform, ecd=ecdf(uTIS_score)(uTIS_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uTIS_score))+
  stat_ecdf(aes(colour = translation), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#74add1", "#fdae61"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  cdf_theme

pdf(file = "fourAdep_uTIS_score_cdf.pdf", width = 5.5, height = 4)
print(cdf_plot)
dev.off()

#Ingolia data
high_sensitivity %>%
  bind_rows(low_sensitivity) %>%
  inner_join(GTI_uTIS_scores, by = c("Gene" = "Ensembl_gene_symbol")) %>%
  mutate(translation = factor(case_when(`Translation_fold_change_to_mean_[log2]` < 0 ~ "high\nsensitivity",
                                        `Translation_fold_change_to_mean_[log2]` > 0 ~ "low\nsensitivity"),
                              levels = c("high\nsensitivity", "low\nsensitivity"), ordered = T)) -> Ingolia_data

t <- wilcox.test(data = Ingolia_data, uTIS_score ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

plot_data <- ddply(Ingolia_data, .(translation), transform, ecd=ecdf(uTIS_score)(uTIS_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uTIS_score))+
  stat_ecdf(aes(colour = translation), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#74add1", "#fdae61"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  cdf_theme

pdf(file = "ingolia_uTIS_score_cdf.pdf", width = 5.5, height = 4)
print(cdf_plot)
dev.off()

#low TE vs high TE----
transcript_to_geneID %>%
  filter(transcript %in% filtered_transcripts$transcript) %>%
  inner_join(MCF7_IDs, by = "transcript")%>%
  inner_join(GTI_uTIS_scores, by = "Ensembl_gene_symbol")%>%
  inner_join(translation_data, by = "gene") %>%
  mutate(TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene,
         TE_2 = mu_PD2_MCF7.gene - mu_SD2_MCF7.gene,
         TE_3 = mu_PD3_MCF7.gene - mu_SD3_MCF7.gene,
         mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)),
         TE = factor(case_when(mean_TE > quantile(mean_TE, 2/3) ~ "high TE",
                        mean_TE < quantile(mean_TE, 1/3) ~ "low TE"),
                     levels = c("low TE", "high TE"), ordered = T)) %>%
  filter(TE == "high TE" | TE == "low TE") %>%
  select(transcript, TE, uTIS_score) -> TE_data

t <- wilcox.test(data = TE_data, uTIS_score ~ TE,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

plot_data <- ddply(TE_data, .(TE), transform, ecd=ecdf(uTIS_score)(uTIS_score))

cdf_plot<- ggplot(data = plot_data, aes(x = uTIS_score))+
  stat_ecdf(aes(colour = TE), size = 1)+ 
  xlab("uTIS score")+
  scale_color_manual(values=c("#7C71D8", "#FFE073"))+
  xlim(c(0, 1))+
  ylab("cumulative fraction")+
  ggtitle(p_label)+
  cdf_theme

pdf(file = "TE_uTIS_score_cdf.pdf", width = 5.5, height = 4)
print(cdf_plot)
dev.off()

