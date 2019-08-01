#load packages
library(tidyverse)
library(scales)
library(viridis)

#import functions----
source("N:\\JWALDRON/R_scripts/functions.R")

#import variables----
source("N:\\JWALDRON/R_scripts/structure_seq_variables.R")

#set directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Stats/panels')

#write themes----
violin_theme <- theme_bw()+
  theme(legend.position='none',
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

scatter_theme <- theme_bw()+
  theme(legend.position = c(0.11,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))

#load data----
#common data
source("N:\\JWALDRON/R_scripts/Structure_seq_common_data.R")

translation_data %>%
  filter(translation == "4A-dep" | translation == "4A-indep" ) %>%
  select(gene, translation, posterior_probability) -> translation_list

#stats
stats_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (condition in c("control", "hippuristanol")) {
    if (region == "tpUTR") {
      df <- read_csv(file = file.path(home, paste0('raw_data/statistics/', condition, '_', region, '_', tp_trim, 'trim_20minlen_statistics.csv')), col_names = T)
    } else {
      df <- read_csv(file = file.path(home, paste0('raw_data/statistics/', condition, '_', region, '_0trim_20minlen_statistics.csv')), col_names = T)
    }
    names(df) <- c("transcript", "max", "average", "std", "gini")
    df$condition <- rep(condition, nrow(df))
    df$region <- rep(region, nrow(df))
    stats_list[[paste(condition, region, sep = "_")]] <- df
  }
}
stats_data <- do.call("rbind", stats_list)
stats_data$condition <- factor(stats_data$condition, levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)

#strandedness
strands_list <- list()
for (condition in c("control", "hippuristanol")) {
  df <- read_csv(file = file.path(home, paste0('raw_data/folding/fpUTRs/', condition, '_fpUTR_strands.csv')), col_names = T)
  df$condition <- rep(condition, nrow(df))
  strands_list[[condition]] <- df
}
strandedness <- do.call("rbind", strands_list)
strandedness$condition <- factor(strandedness$condition, levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)

#MFE
MFE_list <- list()
for (condition in c("control", "hippuristanol")) {
  df <- read_csv(file = file.path(home, paste0('raw_data/folding/fpUTRs/', condition, '_fpUTR_MFE.csv')), col_names = T)
  df$condition <- rep(condition, nrow(df))
  MFE_list[[condition]] <- df
}
MFE <- do.call("rbind", MFE_list)
MFE$condition <- factor(MFE$condition, levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)

#filter and plot average data----
#filter data
stats_data %>%
  filter(is.na(average)) %>%
  pull(transcript) -> remove_IDs # generates a list of transcript IDs which have an NA for average in either condition and any region

stats_data %>%
  filter(!(transcript %in% remove_IDs)) %>%
  select(transcript, average, condition, region) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_coverage > coverage,
         hippuristanol_plus_DMS_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_average_data
  
print(paste("average plots n = ", n_distinct(filtered_average_data$transcript)))

#plot
for (region in c("fpUTR", "CDS", "tpUTR")) {
  filtered_average_data[filtered_average_data$region == region,] %>%
    select(transcript, condition, average) -> df
  
  average_axis_lims <- myAxisLims(df$average, 0.1)
  
  t <- wilcox.test(data = df, average ~ condition,
                   paired = T,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  df %>%
    ggplot(aes(x = condition, y = average, fill = condition))+
    geom_violin(alpha = 0.5)+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab("average reactivity")+
    ylim(average_axis_lims$lower_lim, average_axis_lims$upper_lim)+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> average_violin
  
  pdf(file = paste0(region, '_average_reactivity_violin.pdf'), height = 4, width = 4)
  print(average_violin)
  dev.off()
  
  df %>%
    spread(key = condition, value = average) %>%
    ggplot(aes(Ctrl, Hipp))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    xlim(average_axis_lims$lower_lim, average_axis_lims$upper_lim)+
    ylim(average_axis_lims$lower_lim, average_axis_lims$upper_lim)+
    scatter_theme+
    xlab("average reactivity, Ctrl")+
    ylab("average reactivity, Hipp")+
    geom_abline(color="red", size=1) -> average_scatter

  pdf(file = paste0(region, '_average_reactivity_scatter.pdf'), height = 4, width = 4)
  print(average_scatter)
  dev.off()
  
  #fourAdep
  df %>%
    spread(key = condition, value = average) %>%
    mutate(delta = Hipp - Ctrl) %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(translation_list, by = "gene") -> translation_df
  
  n_fourAdep_transcripts <- n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"])
  print(paste("delta average n =", n_fourAdep_transcripts))
  
  translation_df %>%
    group_by(translation) %>%
    top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
    ungroup() -> delta_df
  
  delta_axis_lims <- myAxisLims(delta_df$delta, 0.1)
  
  t <- wilcox.test(data = delta_df, delta ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  delta_df %>%
    ggplot(aes(x = translation, y = delta, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab(expression(paste(Delta, " reactivity")))+
    ylim(c(delta_axis_lims$lower_lim, delta_axis_lims$upper_lim))+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> average_delta_violin
  
  pdf(file = paste0(region, '_delta_reactivity_violin.pdf'), height = 4, width = 4)
  print(average_delta_violin)
  dev.off()
  
  #TE
  df %>%
    spread(key = condition, value = average) %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(translation_data, by = "gene") %>% 
    mutate(TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene,
           TE_2 = mu_PD1_MCF7.gene - mu_SD2_MCF7.gene,
           TE_3 = mu_PD1_MCF7.gene - mu_SD3_MCF7.gene,
           mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)),
           translation = factor(case_when(mean_TE > quantile(mean_TE, 0.8) ~ "high TE",
                                 mean_TE < quantile(mean_TE, 0.2) ~ "low TE"), levels = c("low TE", "high TE"), ordered = T)) %>%
    filter(translation == "low TE" | translation == "high TE") %>%
    select(transcript, Ctrl, translation) -> TE_data
  
  print(paste(region, "high TE n =", n_distinct(TE_data$transcript[TE_data$translation == "high TE"])))
  print(paste(region, "low TE n =", n_distinct(TE_data$transcript[TE_data$translation == "low TE"])))
  
  t <- wilcox.test(data = TE_data, Ctrl ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  TE_data %>%
    ggplot(aes(x = translation, y = Ctrl, fill = translation))+
    geom_violin(alpha = 0.5)+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    scale_fill_manual(values=c("#7C71D8", "#FFE073"))+
    ylab("average reactivity")+
    ylim(average_axis_lims$lower_lim, average_axis_lims$upper_lim)+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> TE_violin
  
  pdf(file = paste0(region, '_TE_reactivity_violin.pdf'), height = 4, width = 4)
  print(TE_violin)
  dev.off()
}

#filter and plot gini data----
#filter data
stats_data %>%
  filter(is.na(gini)) %>%
  pull(transcript) -> remove_IDs # generates a list of transcript IDs which have an NA for gini in either condition and any region
  
stats_data %>%
  filter(!(transcript %in% remove_IDs)) %>%
  select(transcript, gini, condition, region) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_coverage > coverage,
         hippuristanol_plus_DMS_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_gini_data

print(paste("gini plots n = ", n_distinct(filtered_gini_data$transcript)))

#plot
for (region in c("fpUTR", "CDS", "tpUTR")) {
  filtered_gini_data[filtered_gini_data$region == region,] %>%
    select(transcript, condition, gini) -> df
  
  gini_axis_lims <- myAxisLims(df$gini, 0.1)
  
  t <- wilcox.test(data = df, gini ~ condition,
                   paired = T,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  df %>%
    ggplot(aes(x = condition, y = gini, fill = condition))+
    geom_violin(alpha = 0.5)+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab("Gini coefficient")+
    ylim(gini_axis_lims$lower_lim, gini_axis_lims$upper_lim)+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> gini_violin
  
  pdf(file = paste0(region, '_gini_violin.pdf'), height = 4, width = 4)
  print(gini_violin)
  dev.off()
  
  df %>%
    spread(key = condition, value = gini) %>%
    ggplot(aes(Ctrl, Hipp))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    xlim(gini_axis_lims$lower_lim, gini_axis_lims$upper_lim)+
    ylim(gini_axis_lims$lower_lim, gini_axis_lims$upper_lim)+
    scatter_theme+
    xlab("Gini coefficient, Ctrl")+
    ylab("Gini coefficient, Hipp")+
    geom_abline(color="red", size=1) -> gini_scatter
  
  pdf(file = paste0(region, '_gini_scatter.pdf'), height = 4, width = 4)
  print(gini_scatter)
  dev.off()
  
  #fourAdep
  df %>%
    spread(key = condition, value = gini) %>%
    mutate(delta = Hipp - Ctrl) %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(translation_list, by = "gene") -> translation_df
  
  n_fourAdep_transcripts <- n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"])
  print(paste("delta gini n =", n_fourAdep_transcripts))
  
  translation_df %>%
    group_by(translation) %>%
    top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
    ungroup() -> delta_df
  
  delta_axis_lims <- myAxisLims(delta_df$delta, 0.1)
  
  t <- wilcox.test(data = delta_df, delta ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  delta_df %>%
    ggplot(aes(x = translation, y = delta, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab(expression(paste(Delta, " Gini coefficient")))+
    ylim(c(delta_axis_lims$lower_lim, delta_axis_lims$upper_lim))+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> gini_delta_violin
  
  pdf(file = paste0(region, '_delta_gini_violin.pdf'), height = 4, width = 4)
  print(gini_delta_violin)
  dev.off()
  
  #TE
  df %>%
    spread(key = condition, value = gini) %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(translation_data, by = "gene") %>% 
    mutate(TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene,
           TE_2 = mu_PD1_MCF7.gene - mu_SD2_MCF7.gene,
           TE_3 = mu_PD1_MCF7.gene - mu_SD3_MCF7.gene) %>%
    mutate(mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)),
           translation = factor(case_when(mean_TE > quantile(mean_TE, 0.8) ~ "high TE",
                                   mean_TE < quantile(mean_TE, 0.2) ~ "low TE"), levels = c("low TE", "high TE"), ordered = T)) %>%
    filter(translation == "low TE" | translation == "high TE") %>%
    select(transcript, Ctrl, translation) -> TE_data
  
  print(paste(region, "high TE n =", n_distinct(TE_data$transcript[TE_data$translation == "high TE"])))
  print(paste(region, "low TE n =", n_distinct(TE_data$transcript[TE_data$translation == "low TE"])))
  
  t <- wilcox.test(data = TE_data, Ctrl ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  TE_data %>%
    ggplot(aes(x = translation, y = Ctrl, fill = translation))+
    geom_violin(alpha = 0.5)+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    scale_fill_manual(values=c("#7C71D8", "#FFE073"))+
    ylab("Gini coefficient")+
    ylim(gini_axis_lims$lower_lim, gini_axis_lims$upper_lim)+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> TE_violin
  
  pdf(file = paste0(region, '_TE_gini_violin.pdf'), height = 4, width = 4)
  print(TE_violin)
  dev.off()
}

#filter and plot MFE data----
MFE %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_coverage > coverage,
         hippuristanol_plus_DMS_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup()%>%
  select(transcript, condition, Free_Energy) %>%
  mutate(Free_Energy = Free_Energy - 0.1) -> filtered_MFE_data

print(paste("MFE plots n =", n_distinct(filtered_MFE_data$transcript)))

t <- wilcox.test(data = filtered_MFE_data, Free_Energy ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

filtered_MFE_data %>%
  ggplot(aes(x = condition, y = -Free_Energy, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste("MFE (", Delta, "G)")))+
  scale_y_continuous(trans=reverselog_trans(base=10), labels=trans_format("identity", function(x) -x))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> MFE_violin

pdf(file = 'MFE_violin.pdf', height = 4, width = 4)
MFE_violin
dev.off()

filtered_MFE_data %>%
  spread(key = condition, value = Free_Energy) %>%
  ggplot(aes(-Ctrl, -Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  scale_x_continuous(trans=reverselog_trans(base=10),
                     labels=trans_format("identity", function(x) -x))+
  scale_y_continuous(trans=reverselog_trans(base=10),
                     labels=trans_format("identity", function(x) -x))+
  theme_bw()+
  theme(legend.position = c(0.13,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))+
  xlab(expression(paste("MFE (", Delta, "G), Ctrl")))+
  ylab(expression(paste("MFE (", Delta, "G), Hipp")))+
  geom_abline(color="red", size=1) -> MFE_scatter

pdf(file = 'MFE_scatter.pdf', height = 4, width = 4)
MFE_scatter
dev.off()

filtered_MFE_data %>%
  inner_join(FASTA_compositions_list$fpUTR, by = "transcript") %>%
  mutate(normalised_free_energy = Free_Energy / length) %>%
  select(transcript, condition, normalised_free_energy) %>%
  spread(key = condition, value = normalised_free_energy) %>%
  mutate(delta = Hipp - Ctrl) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") -> MFE_df

n_fourAdep_transcripts <- n_distinct(MFE_df$transcript[MFE_df$translation == "4A-dep"])
print(paste("delta MFE n =", n_fourAdep_transcripts))

translation_df %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
  ungroup() -> MFE_delta

MFE_delta_axis_lims <- myAxisLims(MFE_delta$delta, 0.1)

t <- wilcox.test(data = MFE_delta, delta ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

MFE_delta %>%
  ggplot(aes(x = translation, y = delta, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste(Delta, " normalised MFE")))+
  ylim(c(MFE_delta_axis_lims$lower_lim, MFE_delta_axis_lims$upper_lim))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> MFE_delta_violin

pdf(file = file.path(paste0(getwd(), '/', 'MFE_delta_violin.pdf')), height = 4, width = 4)
MFE_delta_violin
dev.off()

#filter and plot strandedness data----
strandedness %>%
  mutate(strandedness = double_percent * 100) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_coverage > coverage,
         hippuristanol_plus_DMS_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup()%>%
  select(transcript, condition, strandedness) -> filtered_strandedness_data

print(paste("strandedness plots n =", n_distinct(filtered_strandedness_data$transcript)))

t <- wilcox.test(data = filtered_strandedness_data, strandedness ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

filtered_strandedness_data %>%
  ggplot(aes(x = condition, y = strandedness, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab("strandedness (%)")+
  ylim(c(10,90))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> strandedness_violin

pdf(file = 'strandedness_violin.pdf', height = 4, width = 4)
strandedness_violin
dev.off()

filtered_strandedness_data %>%
  spread(key = condition, value = strandedness) %>%
  ggplot(aes(Ctrl, Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  xlim(c(10,90))+
  ylim(c(10,90))+
  scatter_theme+
  xlab("strandedness (%), Ctrl")+
  ylab("strandedness (%), Hipp")+
  geom_abline(color="red", size=1) -> strandedness_scatter

pdf(file = 'strandedness_scatter.pdf', height = 4, width = 4)
strandedness_scatter
dev.off()

filtered_strandedness_data %>%
  spread(key = condition, value = strandedness) %>%
  mutate(delta = Hipp - Ctrl) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") -> strandedness_df

n_fourAdep_transcripts <- n_distinct(strandedness_df$transcript[strandedness_df$translation == "4A-dep"])
print(paste("delta strandedness n =", n_fourAdep_transcripts))

strandedness_df %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
  ungroup() -> strandedness_delta

strandedness_delta_axis_lims <- myAxisLims(strandedness_delta$delta, 0.1)

t <- wilcox.test(data = strandedness_delta, delta ~ translation,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

strandedness_delta_violin <- ggplot(data = strandedness_delta, aes(x = translation, y = delta, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste(Delta, " strandedness (%)")))+
  ylim(c(-7.5,7.5))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6)

pdf(file = 'strandedness_delta_violin.pdf', height = 4, width = 4)
strandedness_delta_violin
dev.off()

