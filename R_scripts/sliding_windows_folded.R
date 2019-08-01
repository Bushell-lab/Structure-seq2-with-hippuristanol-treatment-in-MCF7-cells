###Imports
library(gridExtra)
library(grid)
library(tidyverse)

#import functions----
source("N:\\JWALDRON/R_scripts/functions.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/functions.R")

#import variables----
source("N:\\JWALDRON/R_scripts/structure_seq_variables.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/structure_seq_variables.R")

min_length <- 100 #sets minimum 5'UTR length
wStep <- 10 #sets window step size
wLen <- 25 #sets window size
bin_n <- 25 #sets number of bins

#create theme----
myTheme <- theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

#set working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Sliding_windows/panels')
#setwd('\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/Paper/Figures/R/Sliding_windows/panels')

#load common data----
source("N:\\JWALDRON/R_scripts/Structure_seq_common_data.R")
#home <- '\\\\data.beatson.gla.ac.uk/data/JWALDRON/Structure_seq/MCF7_2015'
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/R_scripts/Structure_seq_common_data.R")

#load windows data----
windows_FASTA <- read_csv(file = file.path(home, paste0('raw_data/sliding_windows/composition_files/control_fpUTR_hippuristanol_fpUTR_', wLen, 'win_', wStep, 'step_composition.csv')), col_names = T)
windows_FASTA %>%
  mutate(step = as.integer(str_replace(transcript, ".+\\_", "")),
         transcript = str_replace(transcript, "\\_.+", "")) -> windows_FASTA

MFE_data_list <- list()
for (condition in c("control", "hipp")) {
  for (i in letters[1:21]) {
    df <- read_csv(file = file.path(home, 'raw_data/folding/all_windows/output_files_fpUTR_25win_5step_310.15_RNAstructure-mfe_md_99999/', paste0(condition, '_MFE_', i, '.csv')), col_names = T)
    
    df %>%
      mutate(step = as.integer(str_replace(transcript, ".+\\_", "")),
             transcript = str_replace(transcript, "\\_.+", ""),
             condition = factor(rep(condition))) %>%
      inner_join(coverage_data, by = "transcript") %>%
      filter(control_plus_DMS_coverage > coverage,
             hippuristanol_plus_DMS_coverage > coverage,
             control_minus_DMS_fp_10_coverage > fp_coverage,
             hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) -> MFE_data_list[[paste(condition, i, sep = "_")]]
  }
}

MFE_data <- do.call("rbind", MFE_data_list)

#merged and filter data----
MFE_data %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(FASTA_compositions_list$fpUTR[, c("transcript", "length")], by = "transcript") %>%
  rename(fpUTR_length = length) %>%
  inner_join(windows_FASTA, by = c("transcript", "step")) %>%
  mutate(GC_content = GC_content * 100,
         GA_content = GA_content * 100) %>%
  rename(windows_GC_content = GC_content,
         windows_GA_content = GA_content,
         wLen = length) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  select(transcript, Free_Energy, step, condition, fpUTR_length, wLen, windows_GC_content, windows_GA_content) -> merged_data

#bin data
merged_data %>%
  mutate(position = 1 + (step * wStep) + round((wLen / 2)),
         normalised_position = position / fpUTR_length,
         bin = calculate_bins(normalised_position, bin_n)) -> binned_data

#plot GC content----
summarised_GC_data_list <- list()
for (bin in 1:bin_n) {
  df <- binned_data[binned_data$condition == "control" & binned_data$bin == bin,]
  
  df %>%
    summarise(GC_mean = mean(windows_GC_content),
              GC_sd = sd(windows_GC_content),
              observations = nrow(df),
              GC_error = qt(0.975, df = observations - 1) * GC_sd / sqrt(observations),
              GC_lower = GC_mean - GC_error,
              GC_upper = GC_mean + GC_error) %>%
    mutate(bin = rep(bin)) -> summarised_GC_data_list[[bin]]
}
summarised_GC_data <- do.call("rbind", summarised_GC_data_list)

GC_content_plot <- ggplot(data = summarised_GC_data, aes(x = bin))+
  geom_line(aes(y = GC_mean), size = 1, colour = "grey")+
  scale_colour_manual(values=c("grey"))+
  geom_ribbon(aes(ymin = GC_lower, ymax = GC_upper), alpha = 0.3, fill = "grey")+
  ylab("GC content (%)")+
  myTheme

pdf(file = "sliding_windows_fpUTR_GC_content_all_transcripts.pdf", width = 4, height = 2)
print(GC_content_plot)
dev.off()

#plot MFE----
summarised_MFE_data_list <- list()
for (condition in c("control", "hipp")) {
  for (bin in 1:bin_n) {
    df <- binned_data[binned_data$condition == condition & binned_data$bin == bin,]
    df %>%
      summarise(Free_Energy_mean = mean(Free_Energy),
                Free_Energy_sd = sd(Free_Energy),
                observations = nrow(df),
                Free_Energy_error = qt(0.975, df = observations - 1) * Free_Energy_sd/ sqrt(observations),
                Free_Energy_lower = Free_Energy_mean - Free_Energy_error,
                Free_Energy_upper = Free_Energy_mean + Free_Energy_error) %>%
      mutate(condition = rep(condition),
             bin = rep(bin)) -> summarised_MFE_data_list[[paste(condition, bin, sep = "_")]]
  }
}
summarised_MFE_data <- do.call("rbind", summarised_MFE_data_list)

MFE_plot <- ggplot(data = summarised_MFE_data, aes(x = bin))+
  geom_line(aes(y = Free_Energy_mean, colour = factor(condition, levels = c("control", "hipp"), ordered = T)), size = 1)+
  geom_ribbon(aes(ymin = Free_Energy_lower, ymax = Free_Energy_upper, fill = condition), alpha = 0.3, colour = NA)+
  ylab(expression(paste("MFE (", Delta, "G)")))+
  myTheme

pdf(file = "sliding_windows_fpUTR_MFE_all_transcripts.pdf", width = 4, height = 2)
print(MFE_plot)
dev.off()

#export legend
legend_plot <- ggplot(data = summarised_MFE_data, aes(x = bin, y = Free_Energy_mean, colour = factor(condition, levels = c("control", "hipp"),
                                                                                                             labels = c("Ctrl", "Hipp"), ordered = T)))+
  geom_line()+
  theme_bw()+
  theme(legend.title = element_blank())

ctrl_hipp_legend <- myLegend(legend_plot)
pdf(file = 'ctrl_hipp_legend.pdf', height = 1, width = 1)
grid.arrange(ctrl_hipp_legend)
dev.off()

#compare between 4A-dep and indep transcripts
translation_data %>%
  inner_join(transcript_to_geneID, by = "gene") %>%
  filter(transcript %in% binned_data$transcript) %>%
  filter(translation == "4A-dep" | translation == "4A-indep") -> translation_df

n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"]) -> n_fourAdep_transcripts

translation_df %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
  ungroup() %>%
  select(transcript, translation) %>%
  inner_join(binned_data, by = "transcript") -> binned_translation_data

summarised_data_list <- list()
for (condition in c("silico")) {
  for (bin in 1:bin_n) {
    df <- binned_translation_data[binned_translation_data$condition == condition & binned_translation_data$bin == bin,]
    
    df %>%
      group_by(translation) %>%
      summarise(Free_Energy_mean = mean(Free_Energy),
                GC_mean = mean(windows_GC_content),
                GA_mean = mean(windows_GA_content),
                Free_Energy_sd = sd(Free_Energy),
                GC_sd = sd(windows_GC_content),
                GA_sd = sd(windows_GA_content),
                observations = nrow(df) / n_distinct(binned_translation_data$translation),
                Free_Energy_error = qt(0.975, df = observations - 1) * Free_Energy_sd/ sqrt(observations),
                Free_Energy_lower = Free_Energy_mean - Free_Energy_error,
                Free_Energy_upper = Free_Energy_mean + Free_Energy_error,
                GC_error = qt(0.975, df = observations - 1) * GC_sd/ sqrt(observations),
                GC_lower = GC_mean - GC_error,
                GC_upper = GC_mean + GC_error,
                GA_error = qt(0.975, df = observations - 1) * GA_sd/ sqrt(observations),
                GA_lower = GA_mean - GA_error,
                GA_upper = GA_mean + GA_error) %>%
      mutate(condition = rep(condition),
             bin = rep(bin)) -> summarised_data_list[[paste(condition, bin, sep = "_")]]
  }
}
summarised_data <- do.call("rbind", summarised_data_list)

GC_content_plot <- ggplot(data = summarised_data, aes(x = bin))+
  #geom_line(size = 1)+
  #geom_ribbon(aes(ymin = GC_lower, ymax = GC_upper), alpha = 0.2, colour = NA)+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_line(aes(y = GC_mean, colour = translation), size = 1)+
  geom_ribbon(aes(ymin = GC_lower, ymax = GC_upper, fill = translation), alpha = 0.2, colour = NA)+
  ylab("GC content (%)")+
  myTheme

pdf(file = "sliding_windows_fpUTR_GC_content_fourAdep.pdf", width = 4, height = 2)
print(GC_content_plot)
dev.off()

GA_content_plot <- ggplot(data = summarised_data, aes(x = bin))+
  #geom_line(size = 1)+
  #geom_ribbon(aes(ymin = GC_lower, ymax = GC_upper), alpha = 0.2, colour = NA)+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_line(aes(y = GA_mean, colour = translation), size = 1)+
  geom_ribbon(aes(ymin = GA_lower, ymax = GA_upper, fill = translation), alpha = 0.2, colour = NA)+
  ylab("GA content (%)")+
  myTheme

pdf(file = "sliding_windows_fpUTR_GA_content_fourAdep.pdf", width = 4, height = 2)
print(GA_content_plot)
dev.off()
