###This script was written by Joseph A. Waldron and produces all panels from Figure 6 in Waldron et al. (2019) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)
###dStruct package can be downloaded from https://github.com/AviranLab/dStruct

#load packages----
library(tidyverse)
library(dStruct)
library(grid)
library(gridExtra)
library(viridis)

#import variables----
source("Structure_seq_variables.R")

#set the significance threshold for the windows identified by dStruct
sig_threshold <- 0.25

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

#bins data
calculate_bins <- function(x, n) {
  bins <- n
  cut_size <- 1 / bins
  breaks <- seq(0, 1, cut_size)
  start <- cut_size / 2
  stop <- 1 - start
  bin_centres <- seq(start, stop, cut_size)
  bin <- as.numeric(cut(x, breaks, include.lowest = F, labels = bin_centres))
  return(bin)
}

#write themes----
violin_theme <- theme_bw()+
  theme(legend.position='none',
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

my_theme <- theme_bw()+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())

x_and_y_theme <- my_theme+
  theme(axis.text = element_text(size = 18))

y_only_theme <- my_theme+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank())

no_labels_theme <- my_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

#load data----
#load common data
source("Structure_seq_common_data.R")

#make a list of filtered transcripts----
#the following pipe makes a vector of all transcript IDs that have a 5'UTR, CDS and 3'UTR,
#have coverage and 5' coverage above thresholds, and are the most abundant transcript
FASTA_compositions %>%
  select(transcript, length, region) %>%
  spread(key = region, value = length) %>%
  rename(fpUTR_length = fpUTR,
         CDS_length = CDS,
         tpUTR_length = tpUTR) %>%
  filter(!(is.na(fpUTR_length)) & !(is.na(CDS_length)) & !(is.na(tpUTR_length))) %>% #ensures all transcripts have a 5'UTR, CDS and 'UTR
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(rowSums(cbind(control_plus_DMS_1_coverage, control_plus_DMS_2_coverage, control_plus_DMS_3_coverage)) > coverage,
         rowSums(cbind(hippuristanol_plus_DMS_1_coverage, hippuristanol_plus_DMS_2_coverage, hippuristanol_plus_DMS_3_coverage)) > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  pull(transcript) -> filtered_transcripts

#read in csvs----
#the following for loops reads in a reactivity csv for each replicate in each condition for each transcript and makes a data frame for each transcript with 6 columns;
#A1-3 being the reactivity in control 1-3 and B1-3 being the reactivity in hippuristanol 1-3. It then saves each data frame in a list.
#Directories containing all transcript csv files for every replicate within each condition can be created by running SF2_pipeline_3c_react_CSVs.sh

reactivity_list <- list()
for (transcript in filtered_transcripts) {
  alist <- list()
  for (condition in c("control", "hippuristanol")) {
    for (i in 1:3) {
      fyle <- paste0(paste(transcript, condition, i, sep = "_"), ".csv")
      path <- paste(condition, i, "all_csvs", sep = "_")
      df <- read.csv(file = file.path(path, fyle), header = T)
      df$condition <- rep(paste0(condition, i), nrow(df))
      alist[[paste(condition, i, sep = "_")]] <- df
    }
  }
  df <- do.call("rbind", alist)
  df %>%
    spread(key = condition, value = Reactivity) %>%
    select(-Position, -Nucleotide) %>%
    as.data.frame() -> df
  
  colnames(df) <- gsub("control", "A", colnames(df))
  colnames(df) <- gsub("hippuristanol", "B", colnames(df))
  
  reactivity_list[[transcript]] <- df
}

#perform dStruc----
dStruct_output <- dStructome(reactivity_list, reps_A = 3, reps_B = 3,
                     min_length = 10, 
                     batches = T,
                     check_signal_strength = T,
                     check_nucs = T,
                     check_quality = T,
                     processes = 1) #'mc.cores' > 1 is not supported on Windows

#assess dStruct windows----
#the following pipe makes a tibble with the length of each UTR and the CDS for each transcript
FASTA_compositions_list$fpUTR[,c("transcript", "length")] %>%
  rename(fpUTR_length = length) %>%
  inner_join(FASTA_compositions_list$CDS[,c("transcript", "length")], by = "transcript") %>%
  rename(CDS_length = length) %>%
  inner_join(FASTA_compositions_list$tpUTR[,c("transcript", "length")], by = "transcript") %>%
  rename(tpUTR_length = length) -> region_lengths

#the following pipe assigns each window to a location based on the positions of the window and the UTR/CDS lengths of that transcript
#if the whole window is in one region, it is assigned to that region. If it starts in one region and finishes in another it is assigned to that UTR/CDS junction
#the number of nts within each region is also calcualted
dStruct_output %>%
  as.tibble() %>%
  rename(transcript = t) %>%
  inner_join(region_lengths, by = "transcript") %>%
  mutate(Start_region = case_when(Start <= fpUTR_length ~ "fpUTR",
                                  Start > fpUTR_length & Start <= (fpUTR_length + CDS_length) ~ "CDS",
                                  Start > (fpUTR_length + CDS_length) ~ "tpUTR"),
         Stop_region = case_when(Stop <= fpUTR_length ~ "fpUTR",
                                 Stop > fpUTR_length & Stop <= (fpUTR_length + CDS_length) ~ "CDS",
                                 Stop > (fpUTR_length + CDS_length) ~ "tpUTR"),
         location = case_when(Start_region == "fpUTR" & Stop_region == "fpUTR" ~ "fpUTR",
                              Start_region == "fpUTR" & Stop_region == "CDS" ~ "fpUTR_CDS",
                              Start_region == "CDS" & Stop_region == "CDS" ~ "CDS",
                              Start_region == "CDS" & Stop_region == "tpUTR" ~ "CDS_tpUTR",
                              Start_region == "tpUTR" & Stop_region == "tpUTR" ~ "tpUTR"),
         window_length = (Stop - Start) + 1,
         fpUTR_count = case_when(location == "fpUTR" ~ window_length,
                                 location == "fpUTR_CDS" ~ (fpUTR_length - Start) + 1,
                                 location != "fpUTR" & location != "fpUTR_CDS" ~ 0),
         CDS_count = case_when(location == "CDS" ~ window_length,
                               location == "fpUTR_CDS" ~ as.numeric(Stop) - fpUTR_length,
                               location == "CDS_tpUTR" ~ ((fpUTR_length + CDS_length) - Start) + 1,
                               location != "CDS" & location != "CDS_tpUTR" ~ 0),
         tpUTR_count = case_when(location == "tpUTR" ~ window_length,
                                 location == "CDS_tpUTR" ~ as.numeric(Stop) - (fpUTR_length + CDS_length),
                                 location != "tpUTR" & location != "CDS_tpUTR" ~ 0)) %>%
  as.data.frame() -> merged_data

#plot FDRs----
fdr_plot <- ggplot(data = merged_data, aes(x = FDR, fill = factor(location, levels = c("fpUTR", "fpUTR_CDS", "CDS", "CDS_tpUTR", "tpUTR"), labels = c("5\'UTR", "5\'UTR / CDS", "CDS", "CDS / 3\'UTR", "3\'UTR"), ordered = T),
                                            colour = factor(location, levels = c("fpUTR", "fpUTR_CDS", "CDS", "CDS_tpUTR", "tpUTR"), labels = c("5\'UTR", "5\'UTR / CDS", "CDS", "CDS / 3\'UTR", "3\'UTR"), ordered = T)))+
  geom_density(alpha = 0.5)+
  geom_vline(xintercept = sig_threshold, linetype="dashed", size = 1)+
  xlab("FDR adjusted p-values")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))
pdf(file = "dStruc_fdr.pdf", width = 6, height = 4)
fdr_plot
dev.off()

#filter to include only windows below a significance threshold----
merged_data %>%
  filter(FDR < sig_threshold) -> sig_data

print(paste("n sig windows =", nrow(sig_data)))
print(paste("n transcripts =", n_distinct(sig_data$transcript)))

#plot window lengths----
win_len_plot <- ggplot(data = sig_data, aes(x = window_length, fill = factor(location, levels = c("fpUTR", "fpUTR_CDS", "CDS", "CDS_tpUTR", "tpUTR"), labels = c("5\'UTR", "5\'UTR / CDS", "CDS", "CDS / 3\'UTR", "3\'UTR"), ordered = T),
                                           colour = factor(location, levels = c("fpUTR", "fpUTR_CDS", "CDS", "CDS_tpUTR", "tpUTR"), labels = c("5\'UTR", "5\'UTR / CDS", "CDS", "CDS / 3\'UTR", "3\'UTR"), ordered = T)))+
  geom_density(alpha = 0.5)+
  xlab("window length")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))
pdf(file = "dStruc_window_length.pdf", width = 6, height = 4)
win_len_plot
dev.off()

#calculate reactivity of windows----
#the following pipe calculates the average control, hippuristanol and delta reactivity for each window
windows_reactivity_list <- list()
for (win in 1:nrow(sig_data)) {
  transcript <- sig_data[win, "transcript"]
  Start <- sig_data[win, "Start"]
  Stop <- sig_data[win, "Stop"]
  
  reactivity_list[[transcript]][Start:Stop,] %>%
    mutate(D1 = B1 - A1,
           D2 = B2 - A2,
           D3 = B3 - A3,
           mean_delta = rowMeans(cbind(D1, D2, D3)),
           mean_ctrl = rowMeans(cbind(A1, A2, A3)),
           mean_hipp = rowMeans(cbind(B1, B2, B3))) %>%
    summarise(delta = mean(mean_delta, na.rm = T),
              ctrl = mean(mean_ctrl, na.rm = T),
              hipp = mean(mean_hipp, na.rm = T)) %>%
    mutate(transcript = transcript,
           Start = Start,
           Stop = Stop) -> windows_reactivity_list[[paste(transcript, Start, Stop, sep = "_")]]
}
windows_reactivity <- do.call("rbind", windows_reactivity_list)

windows_reactivity %>%
  inner_join(sig_data, by = c("transcript", "Start", "Stop")) -> windows_reactivity

#plot average reactivity of windows----
plot_list <- list()
for (location in c("fpUTR", "fpUTR_CDS", "CDS", "CDS_tpUTR", "tpUTR")) {
  df <- windows_reactivity[windows_reactivity$location == location,]
  
  print(paste(location, "n =", nrow(df)))
  
  axis_lims <- myAxisLims(c(df$ctrl, df$hipp), 0.1)
  
  t <- wilcox.test(df$ctrl, df$hipp,
                   paired = T,
                   alternative = "two.sided",
                   var.equal = F,
                   conf.int = T)
  p_label <- myP(t)
  
  df %>%
    select(transcript, Start, Stop, location, ctrl, hipp) %>%
    gather(key = condition, value = reactivity, ctrl, hipp) %>%
    ggplot(aes(x = factor(condition), y = reactivity, fill = factor(condition)))+
    geom_violin(alpha = 0.5)+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab("window reactivity")+
    ylim(0, 1.75)+
    violin_theme+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> plot_list[[paste(location, "reactivity_violin", sep = "_")]]
  
  df %>%
    ggplot(aes(ctrl, hipp))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    xlim(0, 1.75)+
    ylim(0, 1.75)+
    theme_bw()+
    theme(legend.position = c(0.87, 0.23),
          legend.title = element_blank(),
          axis.text = element_text(size=18), 
          axis.title = element_text(size=20))+
    xlab("window reactivity, Ctrl")+
    ylab("window reactivity, Hipp")+
    geom_abline(color="red", size=1) -> plot_list[[paste(location, "reactivity_scatter", sep = "_")]]
  
  #fourAdep
  df %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(translation_list, by = "gene") -> delta_df
  
  n_fourAdep_transcripts <- n_distinct(delta_df$transcript[delta_df$translation == "4A-dep"])
  n_fourAindep_transcripts <- n_distinct(delta_df$transcript[delta_df$translation == "4A-indep"])
  print(paste(location, "4A-dep n =", n_fourAdep_transcripts))
  print(paste(location, "4A-indep n =", n_fourAindep_transcripts))
  
  delta_axis_lims <- myAxisLims(delta_df$delta, 0.5)
  
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
    stat_summary(fun.y=mean, geom='point', shape=16, size=6) -> plot_list[[paste(location, "delta_reactivity_violin", sep = "_")]]
}

pdf(file = 'dStruct_window_reactivities.pdf', height = 8, width = 20)
grid.arrange(plot_list$fpUTR_reactivity_violin, plot_list$fpUTR_CDS_reactivity_violin, plot_list$CDS_reactivity_violin, plot_list$CDS_tpUTR_reactivity_violin, plot_list$tpUTR_reactivity_violin,
             plot_list$fpUTR_reactivity_scatter, plot_list$fpUTR_CDS_reactivity_scatter, plot_list$CDS_reactivity_scatter, plot_list$CDS_tpUTR_reactivity_scatter, plot_list$tpUTR_reactivity_scatter,
             nrow = 2)
dev.off()

pdf(file = 'dStruct_window_delta_reactivities.pdf', height = 4, width = 20)
grid.arrange(plot_list$fpUTR_delta_reactivity_violin, plot_list$fpUTR_CDS_delta_reactivity_violin, plot_list$CDS_delta_reactivity_violin, plot_list$CDS_tpUTR_delta_reactivity_violin, plot_list$tpUTR_delta_reactivity_violin,
             nrow = 1)
dev.off()

#plot binned delta reactivities-----
#first split any windows that overlap junctions into each region, but only if a window has at least 10nt in a region is it assigned to that region
#fpUTR/CDS
windows_reactivity %>%
  filter(location == "fpUTR_CDS") %>%
  select(transcript, Start, Stop, fpUTR_count, CDS_count) -> fpUTR_CDS_junctions

fpUTR_CDS_junctions_reactivity_list <- list()
for (n in 1:nrow(fpUTR_CDS_junctions)) {
  df <- fpUTR_CDS_junctions[n,]
  transcript <- df$transcript
  
  fpUTR_start <- df$Start
  fpUTR_end <- df$Start + df$fpUTR_count - 1
  if (fpUTR_end > (fpUTR_start + 9)) {
    reactivity_list[[df$transcript]][fpUTR_start:fpUTR_end,] %>%
      mutate(D1 = B1 - A1,
             D2 = B2 - A2,
             D3 = B3 - A3,
             mean_delta = rowMeans(cbind(D1, D2, D3)),
             mean_ctrl = rowMeans(cbind(A1, A2, A3)),
             mean_hipp = rowMeans(cbind(B1, B2, B3))) %>%
      summarise(delta = mean(mean_delta, na.rm = T),
                ctrl = mean(mean_ctrl, na.rm = T),
                hipp = mean(mean_hipp, na.rm = T)) %>%
      mutate(transcript = transcript,
             Start = fpUTR_start,
             Stop = fpUTR_end,
             location = "fpUTR") -> fpUTR_CDS_junctions_reactivity_list[[paste(transcript, fpUTR_start, fpUTR_end, sep = "_")]]
  }
  
  CDS_start <- df$Start + df$fpUTR_count
  CDS_end <- df$Stop
  if (CDS_end > (CDS_start + 9)) {
    reactivity_list[[df$transcript]][CDS_start:CDS_end,] %>%
      mutate(D1 = B1 - A1,
             D2 = B2 - A2,
             D3 = B3 - A3,
             mean_delta = rowMeans(cbind(D1, D2, D3)),
             mean_ctrl = rowMeans(cbind(A1, A2, A3)),
             mean_hipp = rowMeans(cbind(B1, B2, B3))) %>%
      summarise(delta = mean(mean_delta, na.rm = T),
                ctrl = mean(mean_ctrl, na.rm = T),
                hipp = mean(mean_hipp, na.rm = T)) %>%
      mutate(transcript = transcript,
             Start = CDS_start,
             Stop = CDS_end,
             location = "CDS") -> fpUTR_CDS_junctions_reactivity_list[[paste(transcript, CDS_start, CDS_end, sep = "_")]]
  }
}
fpUTR_CDS_junctions_reactivity <- do.call("rbind", fpUTR_CDS_junctions_reactivity_list)
fpUTR_CDS_junctions_reactivity <- inner_join(fpUTR_CDS_junctions_reactivity, region_lengths, by = "transcript")

#CDS/tpUTR
windows_reactivity %>%
  filter(location == "CDS_tpUTR") %>%
  select(transcript, Start, Stop, CDS_count, tpUTR_count) -> CDS_tpUTR_junctions

CDS_tpUTR_junctions_reactivity_list <- list()
for (n in 1:nrow(CDS_tpUTR_junctions)) {
  df <- CDS_tpUTR_junctions[n,]
  transcript <- df$transcript
  
  CDS_start <- df$Start
  CDS_end <- df$Start + df$CDS_count - 1
  if (CDS_end > (CDS_start + 9)){
    reactivity_list[[df$transcript]][CDS_start:CDS_end,] %>%
      mutate(D1 = B1 - A1,
             D2 = B2 - A2,
             D3 = B3 - A3,
             mean_delta = rowMeans(cbind(D1, D2, D3)),
             mean_ctrl = rowMeans(cbind(A1, A2, A3)),
             mean_hipp = rowMeans(cbind(B1, B2, B3))) %>%
      summarise(delta = mean(mean_delta, na.rm = T),
                ctrl = mean(mean_ctrl, na.rm = T),
                hipp = mean(mean_hipp, na.rm = T)) %>%
      mutate(transcript = transcript,
             Start = CDS_start,
             Stop = CDS_end,
             location = "CDS") -> CDS_tpUTR_junctions_reactivity_list[[paste(transcript, CDS_start, CDS_end, sep = "_")]]
  }
  
  tpUTR_start <- df$Start + df$CDS_count
  tpUTR_end <- df$Stop
  if (tpUTR_end > (tpUTR_start + 9)) {
    reactivity_list[[df$transcript]][tpUTR_start:tpUTR_end,] %>%
      mutate(D1 = B1 - A1,
             D2 = B2 - A2,
             D3 = B3 - A3,
             mean_delta = rowMeans(cbind(D1, D2, D3)),
             mean_ctrl = rowMeans(cbind(A1, A2, A3)),
             mean_hipp = rowMeans(cbind(B1, B2, B3))) %>%
      summarise(delta = mean(mean_delta, na.rm = T),
                ctrl = mean(mean_ctrl, na.rm = T),
                hipp = mean(mean_hipp, na.rm = T)) %>%
      mutate(transcript = transcript,
             Start = tpUTR_start,
             Stop = tpUTR_end,
             location = "tpUTR") -> CDS_tpUTR_junctions_reactivity_list[[paste(transcript, tpUTR_start, tpUTR_end, sep = "_")]]
  }
}
CDS_tpUTR_junctions_reactivity <- do.call("rbind", CDS_tpUTR_junctions_reactivity_list)
CDS_tpUTR_junctions_reactivity <- inner_join(CDS_tpUTR_junctions_reactivity, region_lengths, by = "transcript")

#bin data
#the following three pipes take the reactivity data from each region, bind rows with the reactivity data for that region from the junctions
#and then bins the data

#fpUTR
windows_reactivity %>%
  filter(location == "fpUTR") %>%
  select(delta, ctrl, hipp, transcript, Start, Stop, location, fpUTR_length, CDS_length, tpUTR_length) %>%
  bind_rows(fpUTR_CDS_junctions_reactivity[fpUTR_CDS_junctions_reactivity$location == "fpUTR",]) %>%
  mutate(position = Start + ((Stop - Start) / 2),
         normalised_position = position / fpUTR_length,
         bin = calculate_bins(normalised_position, 25)) %>%
  select(transcript, delta, bin, location, ctrl, hipp) -> fpUTR_binned

#CDS
windows_reactivity %>%
  filter(location == "CDS") %>%
  select(delta, ctrl, hipp, transcript, Start, Stop, location, fpUTR_length, CDS_length, tpUTR_length) %>%
  bind_rows(fpUTR_CDS_junctions_reactivity[fpUTR_CDS_junctions_reactivity$location == "CDS",]) %>%
  bind_rows(CDS_tpUTR_junctions_reactivity[CDS_tpUTR_junctions_reactivity$location == "CDS",]) %>%
  mutate(position = Start + ((Stop - Start) / 2),
         normalised_position = (position - fpUTR_length) / CDS_length,
         bin = calculate_bins(normalised_position, 50)) %>%
  select(transcript, delta, bin, location, ctrl, hipp) -> CDS_binned

#tpUTR
windows_reactivity %>%
  filter(location == "tpUTR") %>%
  select(delta, ctrl, hipp, transcript, Start, Stop, location, fpUTR_length, CDS_length, tpUTR_length) %>%
  bind_rows(CDS_tpUTR_junctions_reactivity[CDS_tpUTR_junctions_reactivity$location == "tpUTR",]) %>%
  mutate(position = Start + ((Stop - Start) / 2),
         normalised_position = (position - (fpUTR_length + CDS_length)) / tpUTR_length,
         bin = calculate_bins(normalised_position, 25)) %>%
  select(transcript, delta, bin, location, ctrl, hipp) -> tpUTR_binned

binned_data <- bind_rows(fpUTR_binned, CDS_binned, tpUTR_binned)

#calculate 95% confidence intervals
binned_intervals_list <- list()
for (location in c("fpUTR", "CDS", "tpUTR")) {
  df <- binned_data[binned_data$location == location,]
  for (n in 1:max(df$bin)) {
    bin_n <- df[df$bin == n,]
    
    t <- t.test(bin_n$hipp, bin_n$ctrl, paired = T, conf.int = T)
    upper <- t$conf.int[[1]]
    lower <- t$conf.int[[2]]
    
    binned_intervals_list[[paste(location, n, sep = "_")]] <- data.frame(location = location, bin = n, upper = upper, lower =lower)
  }
}
binned_intervals_data <- do.call("rbind", binned_intervals_list)

#summarise data
binned_data %>%
  group_by(location, bin) %>%
  summarise(average_delta = mean(delta)) %>%
  inner_join(binned_intervals_data, by = c("location", "bin")) %>%
  ungroup() -> summarised_binned_data

#calculate axis limits
lower_delta_ylim <- min(c(summarised_binned_data$average_delta, summarised_binned_data$upper))
upper_delta_ylim <- max(c(summarised_binned_data$average_delta,summarised_binned_data$lower))

#plot
binned_delta_plot_list <- list()
for (location in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_binned_data[summarised_binned_data$location == location,]
  xlim <- max(df$bin)
  
  if(location == "fpUTR") {
    binned_delta_plot_list[[location]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    binned_delta_plot_list[[location]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_binned_data[summarised_binned_data$location == "fpUTR",]$bin)
bwidth <- max(summarised_binned_data[summarised_binned_data$location == "CDS",]$bin)
cwidth <- max(summarised_binned_data[summarised_binned_data$location == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'dStruct_windows_binned.pdf', width = 20, height = 3)
grid.arrange(binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()

#plot binned fourAdep vs 4A-indep
binned_data %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") -> fourAdep_binned

#summarise data
fourAdep_binned %>%
  group_by(location, bin, translation) %>%
  summarise(average_delta = mean(delta)) %>%
  ungroup() -> summarised_fourAdep_binned

#calculate axis limits
lower_delta_ylim <- min(c(summarised_fourAdep_binned$average_delta))
upper_delta_ylim <- max(c(summarised_fourAdep_binned$average_delta))

#plot
binned_delta_plot_list <- list()
for (location in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_fourAdep_binned[summarised_fourAdep_binned$location == location,]
  xlim <- max(df$bin)
  
  if(location == "fpUTR") {
    binned_delta_plot_list[[location]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#74add1", "#fdae61"))+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    binned_delta_plot_list[[location]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#74add1", "#fdae61"))+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$location == "fpUTR",]$bin)
bwidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$location == "CDS",]$bin)
cwidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$location == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'dStruct_fourAdep_binned.pdf', width = 20, height = 3)
grid.arrange(binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()

#fpUTR/CDS junction----
#plot the average delta reactivity for 4A-dep and 4A-indep transcripts at each nt, 20nt either side of the translation start site,
#using only those windows overlapping this junction
windows_reactivity %>%
  filter(location == "fpUTR_CDS") %>%
  select(transcript, Start, Stop, fpUTR_count, CDS_count) -> fpUTR_CDS_junctions

#the following for loop calculates the positions relative to the translation start site
fpUTR_CDS_junctions_reactivity_list <- list()
for (n in 1:nrow(fpUTR_CDS_junctions)) {
  df <- fpUTR_CDS_junctions[n,]
  window <- reactivity_list[[df$transcript]][df$Start:df$Stop,]
  
  window$position <- c(-df$fpUTR_count:-1, 1:df$CDS_count) #this is the line that calculates the new positions relative to the start site
  window$transcript <- rep(df$transcript, nrow(window))
  fpUTR_CDS_junctions_reactivity_list[[df$transcript]] <- window
}
fpUTR_CDS_junctions_reactivity <- do.call("rbind", fpUTR_CDS_junctions_reactivity_list)

#filter to include only 20nt either side of the start site, then average the reactivity across replicates at each position in each window
#and then merge with the translation data and then average the delta for each position within each group of mRNAs
fpUTR_CDS_junctions_reactivity %>%
  filter(position >= -20 & position <= 20) %>%
  mutate(D1 = B1 - A1,
         D2 = B2 - A2,
         D3 = B3 - A3,
         average_delta = rowMeans(cbind(D1, D2, D3))) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") %>%
  group_by(position, translation) %>%
  summarise(delta = mean(average_delta, na.rm = T)) %>%
  ungroup() %>%
  filter(position != 2 & position != 3) -> fpUTR_CDS_junctions_fourAdep #this final line removes position 2 and 3 as these are UG from the AUG codon (DMS only methylates As and Cs)

#calculate axis limits
lower_delta_ylim <- min(fpUTR_CDS_junctions_fourAdep$delta)
upper_delta_ylim <- max(fpUTR_CDS_junctions_fourAdep$delta)

#plot
fpUTR_CDS_junction_fourAdep_plot <- ggplot(data = fpUTR_CDS_junctions_fourAdep, aes(x = position, y = delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
  geom_bar(position = position_dodge(), stat='identity')+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_vline(xintercept = 0, size = 3)+
  x_and_y_theme+
  scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
  coord_trans(limx = c(-20.5, 20.5))

pdf(file = "dStruct_fpUTR_CDS_junction_fourAdep.pdf", height = 3, width = 4)
print(fpUTR_CDS_junction_fourAdep_plot)
dev.off()

#CDS/tpUTR junction----
#plot the average delta reactivity for 4A-dep and 4A-indep transcripts at each nt, 20nt either side of the translation stop site,
#using only those windows overlapping this junction
windows_reactivity %>%
  filter(location == "CDS_tpUTR") %>%
  select(transcript, Start, Stop, CDS_count, tpUTR_count) -> CDS_tpUTR_junctions

#the following for loop calculates the positions relative to the translation stop site
CDS_tpUTR_junctions_reactivity_list <- list()
for (n in 1:nrow(CDS_tpUTR_junctions)) {
  df <- CDS_tpUTR_junctions[n,]
  window <- reactivity_list[[df$transcript]][df$Start:df$Stop,]
  
  window$position <- c(-df$CDS_count:-1, 1:df$tpUTR_count) #this is the line that calculates the new positions relative to the stop site
  window$transcript <- rep(df$transcript, nrow(window))
  CDS_tpUTR_junctions_reactivity_list[[df$transcript]] <- window
}
CDS_tpUTR_junctions_reactivity <- do.call("rbind", CDS_tpUTR_junctions_reactivity_list)

#filter to include only 20nt either side of the stop site, then average the reactivity across replicates at each position in each window
#and then merge with the translation data and then average the delta for each position within each group of mRNAs
CDS_tpUTR_junctions_reactivity %>%
  filter(position >= -20 & position <= 20) %>%
  mutate(D1 = B1 - A1,
         D2 = B2 - A2,
         D3 = B3 - A3,
         average_delta = rowMeans(cbind(D1, D2, D3))) %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(translation_list, by = "gene") %>%
  group_by(position, translation) %>%
  summarise(delta = mean(average_delta, na.rm = T)) %>%
  ungroup() %>%
  filter(position != -3) -> CDS_tpUTR_junctions_fourAdep #this final line removes position -3 as this is always a U in the stop codons (UAA,UAG,UGA) (DMS only methylates As and Cs)

#calculate axis limits
lower_delta_ylim <- min(CDS_tpUTR_junctions_fourAdep$delta)
upper_delta_ylim <- max(CDS_tpUTR_junctions_fourAdep$delta)

#plot
CDS_tpUTR_junction_fourAdep_plot <- ggplot(data = CDS_tpUTR_junctions_fourAdep, aes(x = position, y = delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
  geom_bar(position = position_dodge(), stat='identity')+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_vline(xintercept = 0, size = 3)+
  x_and_y_theme+
  scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
  coord_trans(limx = c(-20.5, 20.5))

pdf(file = "dStruct_CDS_tpUTR_junction_fourAdep.pdf", height = 3, width = 4)
print(CDS_tpUTR_junction_fourAdep_plot)
dev.off()

