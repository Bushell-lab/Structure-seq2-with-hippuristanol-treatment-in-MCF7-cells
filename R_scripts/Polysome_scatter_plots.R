###Imports
library(tidyverse)

#import variables----
source("N:\\JWALDRON/R_scripts/structure_seq_variables.R")

#set home and working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Polysomes/panels')

#make theme----
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

#read in translation data
translation_data <- read_tsv(file = file.path(home, 'raw_data/polysomes/penn-DOD-gene.mmdiffMCF7'), col_names = T, skip = 1)

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


translation_data %>%
  mutate(polysomes = rowMeans(cbind(mu_PD1_MCF7.gene, mu_PD2_MCF7.gene, mu_PD3_MCF7.gene)),
         sub_polysomes = rowMeans(cbind(mu_SD1_MCF7.gene, mu_SD2_MCF7.gene, mu_SD3_MCF7.gene)),
         TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene,
         TE_2 = mu_PD1_MCF7.gene - mu_SD2_MCF7.gene,
         TE_3 = mu_PD1_MCF7.gene - mu_SD3_MCF7.gene,
         mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)),
         TE = factor(case_when(mean_TE > quantile(mean_TE, 0.8) ~ "high_TE",
                        mean_TE < quantile(mean_TE, 0.2) ~ "low_TE",
                        mean_TE >= quantile(mean_TE, 0.2) & mean_TE <= quantile(mean_TE, 0.8) ~ "NA"),
                     levels = c("low_TE", "high_TE", "NA"),
                     labels = c("low TE", "high TE", "not assigned"), ordered = T),
         alpha_score = case_when(mean_TE > quantile(mean_TE, 0.8) ~ 1,
                                 mean_TE < quantile(mean_TE, 0.2) ~ 1,
                                 mean_TE >= quantile(mean_TE, 0.2) & mean_TE <= quantile(mean_TE, 0.8) ~ 0.5)) %>%
  select(feature_id, polysomes, sub_polysomes, TE, alpha_score) -> TE_plot_data

#plot scatters
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

TE_scatter <- ggplot(data = TE_plot_data, aes(x = sub_polysomes, y = polysomes, colour = TE, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#7C71D8", "#FFE073", "grey20"))+
  scale_alpha(guide = F)+
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
