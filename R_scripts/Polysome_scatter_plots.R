###This script was written by Joseph A. Waldron and produces panels 2A and 3A in Waldron et al. (2020) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accession GSE134888 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888

#Imports----
library(tidyverse)

#import variables----
source("Structure_seq_variables.R")

#make theme----
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

#read in translation data----
translation_data <- read_tsv(file = 'GSE134888_penn-DOD-gene.mmdiffMCF7.tsv', col_names = T, skip = 1)

#4A-dep logFC plot----
#assign each gene as 4A-dep, 4A-indep, 4A-antidep and NA based on posterior probability and DOD
#eta1_1 = logFC in polysomes
#eta1_2 = logFC in sub-polysomes
#DOD = eta1_1 - eta1_2
translation_data %>%
  mutate(DOD = eta1_1 - eta1_2,
         translation = factor(case_when(posterior_probability > positive_change & DOD < 0 ~ "4A-dep",
                                        posterior_probability > positive_change & DOD > 0 ~ "4A-antidep",
                                        posterior_probability < no_change ~ "4A-indep",
                                        posterior_probability <= positive_change & posterior_probability >= no_change ~ "NA"),
                              levels = c("4A-dep", "4A-indep", "4A-antidep", "NA"),
                              labels = c("4A-dep", "4A-indep", "4A-antidep", "not assigned"), order = T),
         alpha_score = case_when(posterior_probability > positive_change & DOD < 0 ~ 1,
                                 posterior_probability > positive_change & DOD > 0 ~ 1,
                                 posterior_probability < no_change ~ 1,
                                 posterior_probability <= positive_change & posterior_probability >= no_change ~ 0.5)) %>%
  select(feature_id, eta1_1, eta1_2, translation, alpha_score) %>%
  rename(polys_logFC = eta1_1,
         subs_logFC = eta1_2) -> logFC_plot_data

#plot scatter
logFC_scatter <- ggplot(data = logFC_plot_data, aes(x = subs_logFC, y = polys_logFC, colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#74add1", "#fdae61", "#A468D5", "grey"))+
  scale_alpha(guide = F)+
  geom_abline(lty=2)+
  geom_hline(yintercept = 0, lty=2)+
  geom_vline(xintercept = 0, lty=2)+
  xlim(c(-1.5, 1.5))+
  ylim(c(-1.5, 1.5))+
  xlab("Sub-polysomes log-fold change")+
  ylab("Polysomes log-fold change")+
  myTheme

pdf(file = "subs_vs_polys_logFC.pdf", width = 6.5, height = 5)
print(logFC_scatter)
dev.off()

#TE plot----
#assign each gene as high_TE, low_TE or NA based on whether it is in the top, bottom or middle third respectively for the ratio of polysomal RNA over sub-polysomal RNA
translation_data %>%
  mutate(polysomes = rowMeans(cbind(mu_PD1_MCF7.gene, mu_PD2_MCF7.gene, mu_PD3_MCF7.gene)), #calculates the mean of the normalised log values for polysomal RNA
         sub_polysomes = rowMeans(cbind(mu_SD1_MCF7.gene, mu_SD2_MCF7.gene, mu_SD3_MCF7.gene)), #calculates the mean of the normalised log values for sub-polysomal RNA
         TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 1
         TE_2 = mu_PD2_MCF7.gene - mu_SD2_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 2
         TE_3 = mu_PD3_MCF7.gene - mu_SD3_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 3
         mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)), #calculates the mean of TE_1-3
         TE = factor(case_when(mean_TE > quantile(mean_TE, 2/3) ~ "high_TE",
                               mean_TE < quantile(mean_TE, 1/3) ~ "low_TE",
                               mean_TE >= quantile(mean_TE, 1/3) & mean_TE <= quantile(mean_TE, 2/3) ~ "NA"),
                     levels = c("low_TE", "high_TE", "NA"),
                     labels = c("low TE", "high TE", "not assigned"), ordered = T)) %>%
  select(feature_id, polysomes, sub_polysomes, TE) -> TE_plot_data

#plot scatter
TE_scatter <- ggplot(data = TE_plot_data, aes(x = sub_polysomes, y = polysomes, colour = TE))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#7C71D8", "#FFE073", "grey20"))+
  geom_abline(lty=2)+
  geom_hline(yintercept = 0, lty=2)+
  geom_vline(xintercept = 0, lty=2)+
  xlim(c(-1, 7.5))+
  ylim(c(-1, 7.5))+
  ylab("Polysomes log(e) FPKM")+
  xlab("Sub-polysomes log(e) FPKM")+
  myTheme
pdf(file = "subs_vs_polys_TE.pdf", width = 6.5, height = 5)
print(TE_scatter)
dev.off()
