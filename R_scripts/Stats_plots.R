###This script was written by Joseph A. Waldron and produces panels 2B-D, 4A-C, S4A-J and S7B-C in Waldron et al. (2020) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#load packages----
library(tidyverse)
library(viridis)

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

#calculates axis limits to plot data minus the top and bottom n% of values
myAxisLims <- function(x, n) {
  upper_quan <- as.numeric(quantile (x, prob = 1 - (n / 100), na.rm=TRUE))
  lower_quan <- as.numeric(quantile (x, prob = n / 100, na.rm=TRUE))
  return(list(upper_lim = upper_quan, lower_lim = lower_quan))
}

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
#load common data
source("Structure_seq_common_data.R")

#stats data
#uses a for loop to read in data for each region and each condition
#data can be generated with SF2_pipeline_3b_statistics.sh
stats_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (condition in c("control", "hippuristanol")) {
    if (region == "tpUTR") {
      df <- read_csv(file = paste0(condition, '_', region, '_', tp_trim, 'trim_20minlen_statistics.csv'), col_names = T)
    } else {
      df <- read_csv(file = paste0(condition, '_', region, '_0trim_20minlen_statistics.csv'), col_names = T)
    }
    names(df) <- c("transcript", "max", "average", "std", "gini")
    df$condition <- factor(rep(condition, nrow(df)), levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)
    df$region <- factor(rep(region, nrow(df)))
    stats_list[[paste(condition, region, sep = "_")]] <- df
  }
}
stats_data <- do.call("rbind", stats_list)

#MFE and strandedness
#uses a for loop to read in MFE data for each each condition for whole 5'UTRs (if 5'UTR <= 100) and windows (if 5'UTR > 100)
#data can be generated with SF2_pipeline_3d_folding.sh
structure_statistics_list <- list()
for (condition in c("control", "hippuristanol")) {
  whole_fpUTRs <- read_csv(file = paste(condition, 'fpUTR_structure_statistics.csv', sep = "_"), col_names = T)
  whole_fpUTRs %>%
    mutate(condition = factor(rep(condition), levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T),
           step = rep(0)) -> whole_fpUTRs
  
  windows <- read_csv(file = paste(condition, 'fpUTR_100win_10step_structure_statistics.csv', sep = "_"), col_names = T)
  windows %>%
    mutate(step = as.numeric(str_replace(transcript, ".+_", "")),
           transcript = str_replace(transcript, "_.+", ""),
           condition = factor(rep(condition), levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)) -> windows
  
  structure_statistics_list[[condition]] <- bind_rows(whole_fpUTRs, windows)
}
structure_statistics <- do.call("rbind", structure_statistics_list)


#reformat, filter and plot data----

#average reactivity----
#The following two pipes will remove transcripts which have an NA for average in either condition and any region
#filter by coverage and 5' coverage (fp_coverage) and then pick the most abundant transcript per gene

stats_data %>%
  filter(is.na(average)) %>%
  pull(transcript) -> remove_IDs

stats_data %>%
  filter(!(transcript %in% remove_IDs)) %>%
  select(transcript, average, condition, region) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_average_data
  
print(paste("average plots n = ", n_distinct(filtered_average_data$transcript)))

#plot
#the following for loop will create and export as pdfs, the following plots 
#violin and scatter plot for control vs hippuristanol average reactivity
#violin plot for 4A-dep vs 4A-indep delta reactivity
#violin plots of control reactivity in high TE vs low TE mRNAs (top and bottom third of mRNAs based on ratio of polysomal to sub-polysomal RNA)
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
    mutate(TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 1
           TE_2 = mu_PD2_MCF7.gene - mu_SD2_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 2
           TE_3 = mu_PD3_MCF7.gene - mu_SD3_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 3
           mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3)),
           translation = factor(case_when(mean_TE > quantile(mean_TE, 2/3) ~ "high TE",
                                 mean_TE < quantile(mean_TE, 1/3) ~ "low TE"), levels = c("low TE", "high TE"), ordered = T)) %>%
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

#gini----
#The following two pipes will remove transcripts which have an NA for gini in either condition and any region
#filter by coverage and 5' coverage (fp_coverage) and then pick the most abundant transcript per gene
stats_data %>%
  filter(is.na(gini)) %>%
  pull(transcript) -> remove_IDs # generates a list of transcript IDs which have an NA for gini in either condition and any region
  
stats_data %>%
  filter(!(transcript %in% remove_IDs)) %>%
  select(transcript, gini, condition, region) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_gini_data

print(paste("gini plots n = ", n_distinct(filtered_gini_data$transcript)))

#plot
#the following for loop will create and export as pdfs, violin and scatter plot for control vs hippuristanol gini

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
}

#MFE and strandedness----
#filter and summarise
#The following two pipes will remove transcripts which have an NA for DeltaG or strandedness in either condition
#filter by coverage and 5' coverage (fp_coverage), pick the most abundant transcript per gene
#and then calculate the minimum and average MFE and strandedness

structure_statistics %>%
  filter(is.na(CT_DeltaG) | is.na(CT_double_stranded)) %>%
  pull(transcript) -> remove_IDs

structure_statistics %>%
  filter(!(transcript %in% remove_IDs)) %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() %>%
  group_by(condition, transcript) %>%
  summarise(avg_MFE = mean(CT_DeltaG),
            min_MFE = min(CT_DeltaG),
            avg_strandedness = mean(CT_double_stranded) * 100,
            max_strandedness = max(CT_double_stranded) * 100) -> summarised_structure_statistics

print(paste("structure statistics plots n =", n_distinct(summarised_structure_statistics$transcript)))

#plot minimum MFE
t <- wilcox.test(data = summarised_structure_statistics, min_MFE ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

summarised_structure_statistics %>%
  ggplot(aes(x = condition, y = min_MFE, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol)")))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> min_MFE_violin

pdf(file = 'min_MFE_violin.pdf', height = 4, width = 4)
min_MFE_violin
dev.off()

summarised_structure_statistics %>%
  select(transcript, condition, min_MFE) %>%
  spread(key = condition, value = min_MFE) %>%
  ggplot(aes(Ctrl, Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  theme_bw()+
  theme(legend.position = c(0.13,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))+
  xlab(expression(paste("MFE (", Delta, "G kcal/mol), Ctrl")))+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol), Hipp")))+
  geom_abline(color="red", size=1) -> min_MFE_scatter

pdf(file = 'min_MFE_scatter.pdf', height = 4, width = 4)
min_MFE_scatter
dev.off()

#plot average free energy
t <- wilcox.test(data = summarised_structure_statistics, avg_MFE ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

summarised_structure_statistics %>%
  ggplot(aes(x = condition, y = avg_MFE, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol)")))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> avg_MFE_violin

pdf(file = 'avg_MFE_violin.pdf', height = 4, width = 4)
avg_MFE_violin
dev.off()

summarised_structure_statistics %>%
  select(transcript, condition, avg_MFE) %>%
  spread(key = condition, value = avg_MFE) %>%
  ggplot(aes(Ctrl, Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  theme_bw()+
  theme(legend.position = c(0.13,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))+
  xlab(expression(paste("MFE (", Delta, "G kcal/mol), Ctrl")))+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol), Hipp")))+
  geom_abline(color="red", size=1) -> avg_MFE_scatter

pdf(file = 'avg_MFE_scatter.pdf', height = 4, width = 4)
avg_MFE_scatter
dev.off()

#plot delta of the minimum MFE 
summarised_structure_statistics %>%
  select(transcript, condition, min_MFE) %>%
  spread(key = condition, value = min_MFE) %>%
  mutate(delta = Hipp - Ctrl) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") -> MFE_df

n_fourAdep_transcripts <- n_distinct(MFE_df$transcript[MFE_df$translation == "4A-dep"])
print(paste("delta MFE n =", n_fourAdep_transcripts))

MFE_df %>%
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

pdf(file = 'MFE_delta_violin.pdf', height = 4, width = 4)
MFE_delta_violin
dev.off()

#plot maximum strandedness
t <- wilcox.test(data = summarised_structure_statistics, max_strandedness ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

summarised_structure_statistics %>%
  ggplot(aes(x = condition, y = max_strandedness, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab("strandedness (%)")+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> max_strandedness_violin

pdf(file = 'max_strandedness_violin.pdf', height = 4, width = 4)
max_strandedness_violin
dev.off()

summarised_structure_statistics %>%
  select(transcript, condition, max_strandedness) %>%
  spread(key = condition, value = max_strandedness) %>%
  ggplot(aes(Ctrl, Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  theme_bw()+
  theme(legend.position = c(0.13,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))+
  xlab("strandedness (%), Ctrl")+
  ylab("strandedness (%), Hipp")+
  geom_abline(color="red", size=1) -> max_strandedness_scatter

pdf(file = 'max_strandedness_scatter.pdf', height = 4, width = 4)
max_strandedness_scatter
dev.off()

#plot average strandedness
t <- wilcox.test(data = summarised_structure_statistics, avg_strandedness ~ condition,
                 paired = T,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

summarised_structure_statistics %>%
  ggplot(aes(x = condition, y = avg_strandedness, fill = condition))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab("strandedness (%)")+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> avg_strandedness_violin

pdf(file = 'avg_strandedness_violin.pdf', height = 4, width = 4)
avg_strandedness_violin
dev.off()

summarised_structure_statistics %>%
  select(transcript, condition, avg_strandedness) %>%
  spread(key = condition, value = avg_strandedness) %>%
  ggplot(aes(Ctrl, Hipp))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  theme_bw()+
  theme(legend.position = c(0.13,0.78),
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))+
  xlab("strandedness (%), Ctrl")+
  ylab("strandedness (%), Hipp")+
  geom_abline(color="red", size=1) -> avg_strandedness_scatter

pdf(file = 'avg_strandedness_scatter.pdf', height = 4, width = 4)
avg_strandedness_scatter
dev.off()

#plot delta of the maximum strandedness
summarised_structure_statistics %>%
  select(transcript, condition, max_strandedness) %>%
  spread(key = condition, value = max_strandedness) %>%
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

