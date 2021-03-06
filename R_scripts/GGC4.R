###This script was written by Joseph A. Waldron and produces panels 3G-J in Waldron et al. (2019) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#load packages
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

#write themes----
violin_theme <- theme_bw()+
  theme(legend.position='none',
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#load common data
source("Structure_seq_common_data.R")

#motif data
#generate with SF2_pipeline_3e_GGC4.sh
motif_stats_list <- list()
for (motif in c("GGCGGCGGCGGC", "GCCGCCGCCGCC")) {
  for (condition in c("control", "hippuristanol")) {
    df = read_csv(file = file.path(paste0(condition, '_fpUTR_', motif, '_0fp_0tp_0trim_0minlen_statistics.csv')), col_names = T)
    names(df) <- c("transcript", "max", "average", "std", "gini")
    df$condition <- rep(condition, nrow(df))
    df$motif <- rep(motif, nrow(df))
    motif_stats_list[[paste(condition, motif, sep = "_")]] <- df[, c("transcript", "condition", "average", "motif")]
  }
}
motif_stats_data <- do.call("rbind", motif_stats_list)

#5'UTR stats data
#generate with SF2_pipeline_3b_statistics.sh
avg_ctrl <- read_csv(file = 'control_fpUTR_0trim_20minlen_statistics.csv', col_names = T)
avg_hipp <- read_csv(file = 'hippuristanol_fpUTR_0trim_20minlen_statistics.csv', col_names = T)

#merge and filter and normalise data----
#the following pipe splits the transcript and position into two variables, randomly selects one motif per transcript, normalises motif reactivity
#to the average reactivity of the whole 5'UTR, filters by coverage and selects the most abundant transcript per gene
motif_stats_data %>%
  spread(key = condition, value = average) %>%
  mutate(position = str_replace(transcript, ".+_", ""),
         transcript = str_replace(transcript, "_.+", "")) %>%
  group_by(transcript, motif) %>%
  sample_n(size = 1) %>% # randomly selects 1 of each motif per transcript
  ungroup() %>%
  inner_join(avg_ctrl[,c("transcript", "control_fpUTR_average")], by = "transcript") %>% #merges with the average control reactivity for whole 5'UTR
  inner_join(avg_hipp[,c("transcript", "hippuristanol_fpUTR_average")], by = "transcript") %>% #merges with the average hippuristanol reactivity for whole 5'UTR
  mutate(normalised_control = control - control_fpUTR_average,
         normalised_hippuristanol = hippuristanol - hippuristanol_fpUTR_average,
         delta = hippuristanol - control,
         normalised_delta = delta - (hippuristanol_fpUTR_average - control_fpUTR_average)) %>% #normalises ctrl, hipp and delta reactivity of the motif to the average ctrl, hipp and delta reactivity for the whole 5'UTR
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>% #filters by coverage
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  group_by(gene, motif) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_df

summary(as.factor(filtered_df$motif))

#plot normalised hipp----
t <- wilcox.test(data = filtered_df, normalised_hippuristanol ~ motif,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

filtered_df %>%
  ggplot(aes(x = factor(motif, levels = c("GGCGGCGGCGGC", "GCCGCCGCCGCC"), labels = c("(GGC)4", "(GCC)4")), y = normalised_hippuristanol, fill = factor(motif)))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab("normalised hipp reactivity")+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4) -> average_violin

pdf(file = 'GGC_average_reactivity_violin.pdf', height = 4, width = 4)
print(average_violin)
dev.off()

#plot normalised delta----
t <- wilcox.test(data = filtered_df, normalised_delta ~ motif,
                 paired = F,
                 alternative = "two.sided",
                 var.equal = F,
                 conf.int = T)
p_label <- myP(t)

filtered_df %>%
  ggplot(aes(x = factor(motif, levels = c("GGCGGCGGCGGC", "GCCGCCGCCGCC"), labels = c("(GGC)4", "(GCC)4")), y = normalised_delta, fill = factor(motif)))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste("normalised ", Delta, " reactivity")))+
  violin_theme+
  ylim(-0.25, 0.25)+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4) -> delta_violin

pdf(file = 'GGC_delta_reactivity_violin.pdf', height = 4, width = 4)
print(delta_violin)
dev.off()

#plot normalised delta between 4A-dep and 4A-indep----
#the following for loop plots the delta reactivity between 4A-dep and 4A-indep transcripts for each motif
#it selects the same number of 4A-indep motifs as 4A-dep, based on the lowest posterior probability
for (motif in c("GGCGGCGGCGGC", "GCCGCCGCCGCC")) {
  filtered_df[filtered_df$motif == motif,] %>%
    inner_join(translation_list, by = "gene") %>%
    filter(translation == "4A-dep" | translation == "4A-indep") %>%
    select(transcript, position, normalised_delta, translation, posterior_probability) -> translation_df
  
  fourAdep_transcript_n <- n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"])
  
  translation_df %>%
    group_by(translation) %>%
    top_n(wt = -posterior_probability, n = fourAdep_transcript_n) %>%
    ungroup() -> delta_df
  
  print(paste(motif, "4A-dep =", sum(delta_df$translation == "4A-dep")))
  print(paste(motif, "4A-indep =", sum(delta_df$translation == "4A-indep")))
  
  t <- wilcox.test(data = delta_df, normalised_delta ~ translation,
                   paired = F,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  fourAdep_delta_violin <- ggplot(data = delta_df, aes(x = translation, y = normalised_delta, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab(expression(paste("normalised ", Delta, " reactivity")))+
    ylim(c(-0.2, 0.2))+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=4)
  
  pdf(file = paste0(motif, '_delta_reactivity_violin.pdf'), height = 4, width = 4)
  print(fourAdep_delta_violin)
  dev.off()
}

