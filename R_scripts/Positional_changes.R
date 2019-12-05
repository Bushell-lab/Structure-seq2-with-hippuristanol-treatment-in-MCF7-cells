###This script was written by Joseph A. Waldron and produces panels 1C, 2E, 4D, S5A-C, G-H, S6A and S7D in Waldron et al. (2020) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#set your working directory here


#Save the working directory as a string so that we can return here after reading in the csv files
working_dir <- getwd()

#load packages----
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#import variables----
source("Structure_seq_variables.R")

#min length of each UTR/CDS
min_length <- 100

#set the number of bins for each region
fpUTR_bins <- 25
CDS_bins <- 50
tpUTR_bins <- 25

#write functions----
#exports just the legend of a plot
myLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#reads in a csv file (for use with parLapply)
read_react_csv <- function(k){
  df <- read.csv(file = k, header = T)
  df$transcript <- rep(k, nrow(df))
  return(df)
}

#Returns transcript ID (for use with parLapply)
get_transcript_ID <- function(x) {
  transcript <- str_replace(x, "_.+", "")
  return(transcript)
}

#calculate bins
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

#create themes----
my_theme <- theme_bw()+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())
  
x_and_y_theme <- my_theme+
  theme(axis.text = element_text(size = 18))

x_only_theme <- my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank())
  
y_only_theme <- my_theme+
  theme(axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank())

no_labels_theme <- my_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

#load data----
#load common data----
source("Structure_seq_common_data.R")

#make a list of filtered transcripts----
#the following pipe makes a vector of all transcript IDs that have a 5'UTR, CDS and 3'UTR more than 100nt in length,
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
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage,
         fpUTR_length > min_length,
         CDS_length > min_length,
         tpUTR_length > (min_length + tp_trim)) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  pull(transcript) -> filtered_transcripts

#read in csvs with parLapply----
#uses a for loop to read in the csv for every transcript in filtered_transcripts for each condition and region
#a directory containing all transcript csv files for every condition and region can be created by running SF2_pipeline_3c_react_CSVs.sh

#Calculate the number of cores
no_cores <- detectCores() - 1

reactivity_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (condition in c("control", "hippuristanol")) {
    setwd(file.path(working_dir, paste(condition, region, 'all_csvs', sep = "_")))
    csv_list <- dir(pattern = "*.csv") #creates a list of all the csv files in the directory
    filtered_list <- csv_list[lapply(csv_list, get_transcript_ID) %in% filtered_transcripts] #filters the csv list to include only the filtered transcript IDs
    cl <- makeCluster(no_cores) #Initiates cluster
    alist <- parLapply(cl, filtered_list, read_react_csv) #reads in all the csvs from the filtered list and stores them in alist
    stopCluster(cl) #Stops cluster
    df <- do.call("rbind", alist) #combines the list of csv files into a data frame
    df$condition <- rep(condition, nrow(df))
    df$region <- rep(region, nrow(df))
    df$transcript <- str_replace(df$transcript, "_.+", "")
    reactivity_list[[paste(condition, region, sep = "_")]] <- df
  }
}
reactivity_data <- do.call("rbind", reactivity_list)
rm(reactivity_list)

#return to your original working directory
setwd(working_dir)

#bin data----
#the following three pipes bin the data for each region and calculate the
#mean control, hippuristanol and delta reactivity for each bin within each transcript
#fpUTR
reactivity_data[reactivity_data$region == "fpUTR",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>% #calculates 5'UTR length
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, fpUTR_bins)) %>% #calculates bins
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("fpUTR")) -> fpUTR_binned

#CDS
reactivity_data[reactivity_data$region == "CDS",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>% #calculates CDS length
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, CDS_bins)) %>% #calculates bins
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("CDS")) -> CDS_binned

#tpUTR
reactivity_data[reactivity_data$region == "tpUTR",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>% #calculates 3'UTR length
  filter(Position <= (length - tp_trim)) %>% #trims last n nt
  mutate(trimmed_length = rep(max(Position))) %>% # recalculates length
  ungroup() %>%
  mutate(normalised_position = Position / trimmed_length,
         bin = calculate_bins(normalised_position, tpUTR_bins)) %>% #calculates bins
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("tpUTR")) -> tpUTR_binned

binned_data <- bind_rows(fpUTR_binned, CDS_binned, tpUTR_binned)
print(paste("all transcripts binned n =", n_distinct(binned_data$transcript)))

#calculate single nt reactivity----
#the following for loop extracts the first (FP) and last (TP) 60nt from each region,
#after removing the final 125nt from the 3'UTR
single_nt_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  reactivity_data[reactivity_data$region == region,] %>%
    spread(key = condition, value = Reactivity) %>%
    mutate(delta = hippuristanol - control) %>%
    filter(Position <= 60) %>%
    arrange(transcript, Position) %>% 
    mutate(window = rep(seq(2,59, by = 3), each = 3, times = n_distinct(transcript))) -> single_nt_list[[paste0(region, "_FP")]]
  
  if (region == "tpUTR") {#the last n nt of each transcript is trimmed prior to any analysis
    reactivity_data[reactivity_data$region == region,] %>%
      spread(key = condition, value = Reactivity) %>%
      mutate(delta = hippuristanol - control) %>%
      group_by(transcript) %>%
      top_n(n = tp_trim + 60, wt = Position) %>% #extracts the last 60nt + 125nt to trim
      top_n(n = -60, wt = Position) %>% #removes the last 125nt
      arrange(transcript, Position) %>%
      mutate(window = rep(seq(-59,-2, by = 3), each = 3, times = n_distinct(transcript))) %>%
      ungroup() -> single_nt_list[[paste0(region, "_TP")]]
  } else {
    reactivity_data[reactivity_data$region == region,] %>%
      spread(key = condition, value = Reactivity) %>%
      mutate(delta = hippuristanol - control) %>%
      group_by(transcript) %>%
      top_n(n = 60, wt = Position) %>% #extracts the last 60nt
      arrange(transcript, Position) %>%
      mutate(window = rep(seq(-59,-2, by = 3), each = 3, times = n_distinct(transcript))) %>%
      ungroup() -> single_nt_list[[paste0(region, "_TP")]]
  }
}

single_nt_data <- do.call("rbind", single_nt_list)
print(paste("all transcripts single_nt n =", n_distinct(single_nt_data$transcript)))

#plot binned reactivity for all transcripts----

#calculate 95% confidence intervals
binned_intervals_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- binned_data[binned_data$region == region,]
  for (n in 1:max(df$bin)) {
    bin_n <- df[df$bin == n,]
    
    t <- t.test(bin_n$hippuristanol, bin_n$control, paired = T, conf.int = T)
    upper <- t$conf.int[[1]]
    lower <- t$conf.int[[2]]
    
    binned_intervals_list[[paste(region, n, sep = "_")]] <- data.frame(region = region, bin = n, upper = upper, lower =lower)
  }
}
binned_intervals_data <- do.call("rbind", binned_intervals_list)

#summarise data
binned_data %>%
  group_by(region, bin) %>%
  summarise(average_delta = mean(net_change, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  inner_join(binned_intervals_data, by = c("region", "bin")) %>%
  ungroup() -> summarised_binned_data

#calculate axis limits
raw_ylim <- max(c(summarised_binned_data$average_control,summarised_binned_data$average_hippuristanol), na.rm = T)

lower_delta_ylim <- min(c(summarised_binned_data$average_delta, summarised_binned_data$upper))
upper_delta_ylim <- max(c(summarised_binned_data$average_delta,summarised_binned_data$lower))

#plot
binned_raw_plot_list <- list()
binned_delta_plot_list <- list()

for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_binned_data[summarised_binned_data$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      y_only_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      no_labels_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_binned_data[summarised_binned_data$region == "fpUTR",]$bin)
bwidth <- max(summarised_binned_data[summarised_binned_data$region == "CDS",]$bin)
cwidth <- max(summarised_binned_data[summarised_binned_data$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'all_transcripts_binned.pdf', width = 20, height = 6)
grid.arrange(binned_raw_plot_list$fpUTR, binned_raw_plot_list$CDS, binned_raw_plot_list$tpUTR,
             binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 2, widths = plot_widths, heights = c(3,2))
dev.off()

#plot single nt reactivity for all transcripts----

#calculate and 95% confidence intervals
single_nt_intervals_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (n in c(seq(-59,-2, by = 3), seq(2,59, by = 3))) {
    window_n <- single_nt_data[single_nt_data$region == region & single_nt_data$window == n,]
    
    t <- t.test(window_n$hippuristanol, window_n$control, paired = T, conf.int = T)
    upper <- t$conf.int[[1]]
    lower <- t$conf.int[[2]]
    
    single_nt_intervals_list[[paste(region, n, sep = "_")]] <- data.frame(region = region, window = n, upper = upper, lower =lower)
  }
}
single_nt_intervals_data <- do.call("rbind", single_nt_intervals_list)

#summarise data
single_nt_data %>%
  group_by(region, window) %>%
  summarise(average_delta = mean(delta, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  inner_join(single_nt_intervals_data, by = c("region", "window")) %>%
  ungroup() %>%
  mutate(location = case_when(window > 0 ~ "Start",
                              window < 0 ~ "End"))-> summarised_single_nt_data

#calculate axis limits
raw_ylim <- max(c(summarised_single_nt_data$average_control, summarised_single_nt_data$average_hippuristanol), na.rm = T)

lower_delta_ylim <- min(c(summarised_single_nt_data$average_delta, summarised_single_nt_data$upper))
upper_delta_ylim <- max(c(summarised_single_nt_data$average_delta, summarised_single_nt_data$lower))

single_nt_raw_plot_list <- list()
single_nt_delta_plot_list <- list()

for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (location in c("Start", "End")) {
    df <- summarised_single_nt_data[summarised_single_nt_data$region == region & summarised_single_nt_data$location == location,]
    
    lower_xlim <- min(df$window) - 1
    upper_xlim <- max(df$window) + 1
    
    if (region == "fpUTR" & location == "Start") {
      df %>%
        select(window, average_control, average_hippuristanol) %>%
        gather(key = Condition, value = Reactivity, average_control, average_hippuristanol) %>%
        ggplot(aes(x = window, y = Reactivity, colour = factor(Condition)))+
        geom_line(size = 1)+
        y_only_theme+
        scale_y_continuous(limits = c(0,raw_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5)) -> single_nt_raw_plot_list[[paste(region, location, sep = "_")]]
      
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window))+
        geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
        geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
        x_and_y_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
      
    }else{
      df %>%
        select(window, average_control, average_hippuristanol) %>%
        gather(key = Condition, value = Reactivity, average_control, average_hippuristanol) %>%
        ggplot(aes(x = window, y = Reactivity, colour = factor(Condition)))+
        geom_line(size = 1)+
        no_labels_theme+
        scale_y_continuous(limits = c(0,raw_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5)) -> single_nt_raw_plot_list[[paste(region, location, sep = "_")]]
      
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window))+
        geom_bar(aes(y = average_delta), fill = "grey", position = position_dodge(), stat='identity')+
        geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
        x_only_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
    }
  }
}

#export figure
start_margin = theme(plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"))
end_margin = theme(plot.margin = unit(c(0.5 ,0, 0.5, 1), "cm"))

fpUTR_Start <- list(single_nt_raw_plot_list[[1]], single_nt_delta_plot_list[[1]])
fpUTR_End <- list(single_nt_raw_plot_list[[2]], single_nt_delta_plot_list[[2]])

CDS_Start <- list(single_nt_raw_plot_list[[3]], single_nt_delta_plot_list[[3]])
CDS_End <- list(single_nt_raw_plot_list[[4]], single_nt_delta_plot_list[[4]])

tpUTR_Start <- list(single_nt_raw_plot_list[[5]], single_nt_delta_plot_list[[5]])
tpUTR_End <- list(single_nt_raw_plot_list[[6]], single_nt_delta_plot_list[[6]])

pdf(file = 'all_transcripts_single_nt.pdf', width = 20, height = 6)
grid.arrange(grobs = c(lapply(fpUTR_Start, "+", start_margin), lapply(fpUTR_End, "+", end_margin),
                       lapply(CDS_Start, "+", start_margin), lapply(CDS_End, "+", end_margin),
                       lapply(tpUTR_Start, "+", start_margin), lapply(tpUTR_End, "+", end_margin)),
             nrow = ,
             layout_matrix = matrix(data = 1:12, nrow = 2),
             heights = c(3,2))
dev.off()

#export legend
summarised_binned_data[summarised_binned_data$region == "fpUTR",] %>%
  select(bin, average_control, average_hippuristanol) %>%
  gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
  ggplot(aes(x = bin, y = Reactivity, colour = factor(condition, levels = c("average_control", "average_hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T)))+
  geom_line(size = 1)+
  theme(legend.title = element_blank()) -> legend_plot

ctrl_hipp_legend <- myLegend(legend_plot)
pdf(file = 'ctrl_hipp_legend.pdf', height = 1, width = 1)
grid.arrange(ctrl_hipp_legend)
dev.off()

#TE data----
#calculate TE for all filtered transcripts
transcript_to_geneID %>%
  filter(transcript %in% filtered_transcripts) %>%
  inner_join(translation_data, by = "gene") %>%
  mutate(TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 1
         TE_2 = mu_PD2_MCF7.gene - mu_SD2_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 2
         TE_3 = mu_PD3_MCF7.gene - mu_SD3_MCF7.gene, #subtracts normalised log values for sub-polysomal RNA from polysomal RNA in sample 3
         mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3))) -> TE_data #calculates the mean of TE_1-3

#get transcript IDs of the top and bottom third of transcripts based on TE
TE_data %>%
  filter(mean_TE > quantile(TE_data$mean_TE, 2/3)) %>%
  pull(transcript) -> high_TE_transcripts

TE_data %>%
  filter(mean_TE < quantile(TE_data$mean_TE, 1/3)) %>%
  pull(transcript) -> low_TE_transcripts

#merge with binned reactivity data
binned_data %>%
  filter(transcript %in% high_TE_transcripts | transcript %in% low_TE_transcripts) %>%
  select(-hippuristanol, -net_change) %>%
  mutate(translation = factor(case_when(transcript %in% low_TE_transcripts ~ "low_TE",
                                 transcript %in% high_TE_transcripts ~ "high_TE"), levels = c("low_TE", "high_TE"), ordered = T)) -> TE_binned

print(paste("low TE binned n =", n_distinct(TE_binned$transcript[TE_binned$translation == "low_TE"])))
print(paste("high TE binned n =", n_distinct(TE_binned$transcript[TE_binned$translation == "high_TE"])))

#merge with single_nt reactivity data
single_nt_data %>%
  filter(transcript %in% high_TE_transcripts | transcript %in% low_TE_transcripts) %>%
  select(-hippuristanol, -delta) %>%
  mutate(translation = factor(case_when(transcript %in% low_TE_transcripts ~ "low_TE",
                                        transcript %in% high_TE_transcripts ~ "high_TE"), levels = c("low_TE", "high_TE"), ordered = T)) -> TE_single_nt

print(paste("low TE single_nt n =", n_distinct(TE_single_nt$transcript[TE_single_nt$translation == "low_TE"])))
print(paste("high TE single_nt n =", n_distinct(TE_single_nt$transcript[TE_single_nt$translation == "high_TE"])))

#plot binned reactivity for TE transcripts----

#calculate 95% confidence intervals
binned_intervals_list <- list()

for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- TE_binned[TE_binned$region == region,]
  for (n in 1:max(df$bin)) {
    bin_n <-  df[df$bin == n,]
    
    t <- t.test(data = bin_n, control ~ translation, paired = F, conf.int = T)
    upper <- t$conf.int[[1]]
    lower <- t$conf.int[[2]]
    
    binned_intervals_list[[paste(region, n, sep = "_")]] <- data.frame(region = region, bin = n, upper = upper, lower =lower)
  }
}
binned_intervals_data <- do.call("rbind", binned_intervals_list)

#summarise data
TE_binned %>%
  group_by(region, bin, translation) %>%
  summarise(average_control = mean(control, na.rm = T))  %>%
  spread(key = translation, value = average_control) %>%
  mutate(delta = low_TE - high_TE) %>%
  inner_join(binned_intervals_data, by = c("region", "bin")) %>%
  ungroup() -> summarised_TE_binned

#calculate axis limits
raw_ylim <- max(c(summarised_TE_binned$high_TE, summarised_TE_binned$low_TE), na.rm = T)
lower_delta_ylim <- min(c(summarised_TE_binned$delta, summarised_TE_binned$upper))
upper_delta_ylim <- max(c(summarised_TE_binned$delta,summarised_TE_binned$lower))

#plot
binned_raw_plot_list <- list()
binned_delta_plot_list <- list()

for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_TE_binned[summarised_TE_binned$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    df %>%
      select(bin, low_TE, high_TE) %>%
      gather(key = translation, value = Reactivity, low_TE, high_TE) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(translation, levels = c("low_TE", "high_TE"), ordered = T)))+
      scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
      geom_line(size = 1)+
      y_only_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    df %>%
      select(bin, low_TE, high_TE) %>%
      gather(key = translation, value = Reactivity, low_TE, high_TE) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(translation, levels = c("low_TE", "high_TE"), ordered = T)))+
      scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
      geom_line(size = 1)+
      no_labels_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin))+
      geom_bar(aes(y = delta), fill = "grey", position = position_dodge(), stat='identity')+
      geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_TE_binned[summarised_TE_binned$region == "fpUTR",]$bin)
bwidth <- max(summarised_TE_binned[summarised_TE_binned$region == "CDS",]$bin)
cwidth <- max(summarised_TE_binned[summarised_TE_binned$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'TE_binned.pdf', width = 20, height = 6)
grid.arrange(binned_raw_plot_list$fpUTR, binned_raw_plot_list$CDS, binned_raw_plot_list$tpUTR,
             binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 2, widths = plot_widths, heights = c(3,2))
dev.off()

#plot single nt reactivity for TE transcripts----

#calculate 95% confidence intervals
single_nt_intervals_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (n in c(seq(-59,-2, by = 3), seq(2,59, by = 3))) {
    window_n <- TE_single_nt[TE_single_nt$region == region & TE_single_nt$window == n,]
    
    t <- t.test(data = window_n, control ~ translation, paired = F, conf.int = T)
    upper <- t$conf.int[[1]]
    lower <- t$conf.int[[2]]
    
    single_nt_intervals_list[[paste(region, n, sep = "_")]] <- data.frame(region = region, window = n, upper = upper, lower =lower)
  }
}
single_nt_intervals_data <- do.call("rbind", single_nt_intervals_list)

#summarise data
TE_single_nt %>%
  group_by(region, window, translation) %>%
  summarise(average_control = mean(control, na.rm = T)) %>%
  spread(key = translation, value = average_control) %>%
  mutate(delta = low_TE - high_TE) %>%
  inner_join(single_nt_intervals_data, by = c("region", "window")) %>%
  ungroup() %>%
  mutate(location = case_when(window > 0 ~ "Start",
                              window < 0 ~ "End")) -> summarised_TE_single_nt

#calculate axis limits
raw_ylim <- max(c(summarised_TE_single_nt$high_TE, summarised_TE_single_nt$low_TE), na.rm = T)

lower_delta_ylim <- min(c(summarised_TE_single_nt$delta, summarised_TE_single_nt$upper))
upper_delta_ylim <- max(c(summarised_TE_single_nt$delta, summarised_TE_single_nt$lower))

#plot
single_nt_raw_plot_list <- list()
single_nt_delta_plot_list <- list()

for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (location in c("Start", "End")) {
    df <- summarised_TE_single_nt[summarised_TE_single_nt$region == region & summarised_TE_single_nt$location == location,]
    
    lower_xlim <- min(df$window) - 1
    upper_xlim <- max(df$window) + 1
    
    if (region == "fpUTR" & location == "Start") {
      df %>%
        select(window, low_TE, high_TE) %>%
        gather(key = translation, value = Reactivity, low_TE, high_TE) %>%
        ggplot(aes(x = window, y = Reactivity, colour = factor(translation, levels = c("low_TE", "high_TE"), ordered = T)))+
        geom_line(size = 1)+
        scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
        y_only_theme+
        scale_y_continuous(limits = c(0,raw_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5)) -> single_nt_raw_plot_list[[paste(region, location, sep = "_")]]
      
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window))+
        geom_bar(aes(y = delta), fill = "grey", position = position_dodge(), stat='identity')+
        geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
        x_and_y_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
      
    }else{
      df %>%
        select(window, low_TE, high_TE) %>%
        gather(key = translation, value = Reactivity, low_TE, high_TE) %>%
        ggplot(aes(x = window, y = Reactivity, colour = factor(translation, levels = c("low_TE", "high_TE"), ordered = T)))+
        geom_line(size = 1)+
        scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
        no_labels_theme+
        scale_y_continuous(limits = c(0,raw_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5)) -> single_nt_raw_plot_list[[paste(region, location, sep = "_")]]
      
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window))+
        geom_bar(aes(y = delta), fill = "grey", position = position_dodge(), stat='identity')+
        geom_ribbon(aes(ymin = upper, ymax = lower), alpha = 0.3, colour = NA)+
        x_only_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
    }
  }
}

#export figure
start_margin = theme(plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"))
end_margin = theme(plot.margin = unit(c(0.5 ,0, 0.5, 1), "cm"))

fpUTR_Start <- list(single_nt_raw_plot_list[[1]], single_nt_delta_plot_list[[1]])
fpUTR_End <- list(single_nt_raw_plot_list[[2]], single_nt_delta_plot_list[[2]])

CDS_Start <- list(single_nt_raw_plot_list[[3]], single_nt_delta_plot_list[[3]])
CDS_End <- list(single_nt_raw_plot_list[[4]], single_nt_delta_plot_list[[4]])

tpUTR_Start <- list(single_nt_raw_plot_list[[5]], single_nt_delta_plot_list[[5]])
tpUTR_End <- list(single_nt_raw_plot_list[[6]], single_nt_delta_plot_list[[6]])

pdf(file = 'TE_single_nt.pdf', width = 20, height = 6)
grid.arrange(grobs = c(lapply(fpUTR_Start, "+", start_margin), lapply(fpUTR_End, "+", end_margin),
                       lapply(CDS_Start, "+", start_margin), lapply(CDS_End, "+", end_margin),
                       lapply(tpUTR_Start, "+", start_margin), lapply(tpUTR_End, "+", end_margin)),
             nrow = ,
             layout_matrix = matrix(data = 1:12, nrow = 2),
             heights = c(3,2))
dev.off()

#export legend
summarised_TE_binned[summarised_TE_binned$region == "fpUTR",] %>%
  select(bin, low_TE, high_TE) %>%
  gather(key = translation, value = Reactivity, low_TE, high_TE) %>%
  ggplot(aes(x = bin, y = Reactivity, colour = factor(translation, levels = c("low_TE", "high_TE"), labels = c("low TE", "high TE"), ordered = T)))+
  geom_line(size = 1)+
  scale_colour_manual(values=c("#7C71D8", "#FFE073"))+
  theme(legend.title = element_blank()) -> legend_plot

TE_legend <- myLegend(legend_plot)
pdf(file = 'TE_legend_legend.pdf', height = 1, width = 1)
grid.arrange(TE_legend)
dev.off()

#aTIS transcripts----
#read in gene names for MCF7 transcripts
MCF7_IDs <- read_csv(file = "MCF7_2015_ensembl_IDs.csv", col_names = T) #download from the data folder of this repository

#read in uTIS scores
uTIS_scores <- read_csv(file = "GTI_data_uTIS_scores.csv", col_names = T) #generate with Custom_scripts.sh

#extract filtered transcript IDs with a uTIS score of 0
uTIS_scores %>%
  inner_join(MCF7_IDs, by = "Ensembl_gene_symbol") %>%
  filter(transcript %in% filtered_transcripts & uTIS_score == 0) %>%
  pull(transcript) -> aTIS_transcripts

#merge with reactivity data
binned_data %>%
  mutate(uTIS = case_when(transcript %in% aTIS_transcripts ~ "aTIS",
                          !(transcript %in% aTIS_transcripts) ~ "all_transcripts")) -> aTIS_binned

print(paste("aTIS binned n = ", n_distinct(aTIS_binned$transcript[aTIS_binned$uTIS == "aTIS"])))
print(paste("aTIS all transcripts binned n = ", n_distinct(aTIS_binned$transcript[aTIS_binned$uTIS == "all_transcripts"])))

single_nt_data %>%
  mutate(uTIS = case_when(transcript %in% aTIS_transcripts ~ "aTIS",
                          !(transcript %in% aTIS_transcripts) ~ "all_transcripts"))  -> aTIS_single_nt

print(paste("aTIS single_nt n = ", n_distinct(aTIS_single_nt$transcript[aTIS_single_nt$uTIS == "aTIS"])))
print(paste("aTIS all transcripts single_nt n = ", n_distinct(aTIS_single_nt$transcript[aTIS_single_nt$uTIS == "all_transcripts"])))

#plot binned reactivity for aTIS transcripts----

#summarise data
aTIS_binned %>%
  group_by(region, bin, uTIS) %>%
  summarise(average_delta = mean(net_change, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() -> summarised_aTIS_binned

#calculate axis limits
lower_delta_ylim <- min(c(summarised_aTIS_binned$average_delta))
upper_delta_ylim <- max(c(summarised_aTIS_binned$average_delta))

#plot
binned_delta_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_aTIS_binned[summarised_aTIS_binned$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(uTIS, levels = c("all_transcripts", "aTIS"),
                                                                                                        ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#5D1882", "#A6A61C"))+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(uTIS, levels = c("all_transcripts", "aTIS"),
                                                                                                        ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#5D1882", "#A6A61C"))+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_aTIS_binned[summarised_aTIS_binned$region == "fpUTR",]$bin)
bwidth <- max(summarised_aTIS_binned[summarised_aTIS_binned$region == "CDS",]$bin)
cwidth <- max(summarised_aTIS_binned[summarised_aTIS_binned$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'aTIS_binned.pdf', width = 20, height = 3)
grid.arrange(binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()

#plot single nt reactivity for aTIS transcripts----

#summarise data
aTIS_single_nt %>%
  group_by(region, window, uTIS) %>%
  summarise(average_delta = mean(delta, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() %>%
  mutate(location = case_when(window > 0 ~ "Start",
                              window < 0 ~ "End"))-> summarised_aTIS_single_nt

#calculate axis limits
lower_delta_ylim <- min(c(summarised_aTIS_single_nt$average_delta))
upper_delta_ylim <- max(c(summarised_aTIS_single_nt$average_delta))

#plot
single_nt_delta_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (location in c("Start", "End")) {
    df <- summarised_aTIS_single_nt[summarised_aTIS_single_nt$region == region & summarised_aTIS_single_nt$location == location,]
    
    lower_xlim <- min(df$window) - 1
    upper_xlim <- max(df$window) + 1
    
    if (region == "fpUTR" & location == "Start") {
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window, y = average_delta, fill = factor(uTIS, levels = c("all_transcripts", "aTIS"),
                                                                                                                                            ordered = T)))+
        geom_bar(position = position_dodge(), stat='identity')+
        scale_fill_manual(values=c("#5D1882", "#A6A61C"))+
        x_and_y_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
      
    }else{
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window, y = average_delta, fill = factor(uTIS, levels = c("all_transcripts", "aTIS"),
                                                                                                                                            ordered = T)))+
        geom_bar(position = position_dodge(), stat='identity')+
        scale_fill_manual(values=c("#5D1882", "#A6A61C"))+
        x_only_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
    }
  }
}

#export figure
start_margin = theme(plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"))
end_margin = theme(plot.margin = unit(c(0.5 ,0, 0.5, 1), "cm"))

fpUTR_Start <- list(single_nt_delta_plot_list[[1]])
fpUTR_End <- list(single_nt_delta_plot_list[[2]])

CDS_Start <- list(single_nt_delta_plot_list[[3]])
CDS_End <- list(single_nt_delta_plot_list[[4]])

tpUTR_Start <- list(single_nt_delta_plot_list[[5]])
tpUTR_End <- list(single_nt_delta_plot_list[[6]])

pdf(file = 'aTIS_single_nt.pdf', width = 20, height = 3)
grid.arrange(grobs = c(lapply(fpUTR_Start, "+", start_margin), lapply(fpUTR_End, "+", end_margin),
                       lapply(CDS_Start, "+", start_margin), lapply(CDS_End, "+", end_margin),
                       lapply(tpUTR_Start, "+", start_margin), lapply(tpUTR_End, "+", end_margin)),
             nrow = 1)
dev.off()

#export legend
legend_plot <- ggplot(data = summarised_aTIS_single_nt, aes(x = window, y = average_delta, fill = factor(uTIS, levels = c("all_transcripts", "aTIS"),
                                                                                      ordered = T, labels = c("all\ntranscripts", "aTIS"))))+
  geom_bar(position = position_dodge(), stat='identity')+
  scale_fill_manual(values=c("#5D1882", "#A6A61C"))+
  theme(legend.title = element_blank())

aTIS_legend <- myLegend(legend_plot)
pdf(file = 'aTIS_legend.pdf', height = 1, width = 1)
grid.arrange(aTIS_legend)
dev.off()

#fourAdep----
#filter translation data to include only 4A-dep and 4A-indep transcripts in filtered_transcripts
translation_data %>%
  inner_join(transcript_to_geneID, by = "gene") %>%
  filter(transcript %in% filtered_transcripts,
         translation == "4A-dep" | translation == "4A-indep") -> translation_df

#calculate the number of 4A-dep transcripts and then select the same number of 4A-indep transcripts based on the lowest posterior probability
n_fourAdep_transcripts <- n_distinct(translation_df$transcript[translation_df$translation == "4A-dep"])

translation_df %>%
  group_by(translation) %>%
  top_n(wt = -posterior_probability, n = n_fourAdep_transcripts) %>%
  ungroup() -> filtered_translation_data

#merge with reactivity data
filtered_translation_data %>%
  inner_join(binned_data, by = "transcript") -> fourAdep_binned

print(paste("4A-dep binned =", n_distinct(fourAdep_binned$transcript[fourAdep_binned$translation == "4A-dep"])))
print(paste("4A-indep binned =", n_distinct(fourAdep_binned$transcript[fourAdep_binned$translation == "4A-indep"])))

filtered_translation_data %>%
  inner_join(single_nt_data, by = "transcript") -> fourAdep_single_nt

print(paste("4A-dep single_nt =", n_distinct(fourAdep_single_nt$transcript[fourAdep_single_nt$translation == "4A-dep"])))
print(paste("4A-indep single_nt =", n_distinct(fourAdep_single_nt$transcript[fourAdep_single_nt$translation == "4A-indep"])))

#plot binned reactivity for fourAdep----

#summarise data
fourAdep_binned %>%
  group_by(region, bin, translation) %>%
  summarise(average_delta = mean(net_change, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() -> summarised_fourAdep_binned

#calculate axis limits
lower_delta_ylim <- min(c(summarised_fourAdep_binned$average_delta))
upper_delta_ylim <- max(c(summarised_fourAdep_binned$average_delta))

#plot
binned_delta_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_fourAdep_binned[summarised_fourAdep_binned$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#74add1", "#fdae61"))+
      y_only_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
    
  }else{
    binned_delta_plot_list[[region]] <- ggplot(data = df, aes(x = bin, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
      geom_bar(position = position_dodge(), stat='identity')+
      scale_fill_manual(values=c("#74add1", "#fdae61"))+
      no_labels_theme+
      scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5))
  }
}

#export figures
awidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$region == "fpUTR",]$bin)
bwidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$region == "CDS",]$bin)
cwidth <- max(summarised_fourAdep_binned[summarised_fourAdep_binned$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'fourAdep_binned.pdf', width = 20, height = 3)
grid.arrange(binned_delta_plot_list$fpUTR, binned_delta_plot_list$CDS, binned_delta_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()

#plot single nt reactivity for fourAdep transcripts----

#summarise data
fourAdep_single_nt %>%
  group_by(region, window, translation) %>%
  summarise(average_delta = mean(delta, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() %>%
  mutate(location = case_when(window > 0 ~ "Start",
                              window < 0 ~ "End"))-> summarised_fourAdep_single_nt

#calculate axis limits
lower_delta_ylim <- min(c(summarised_fourAdep_single_nt$average_delta))
upper_delta_ylim <- max(c(summarised_fourAdep_single_nt$average_delta))

#plot
single_nt_delta_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (location in c("Start", "End")) {
    df <- summarised_fourAdep_single_nt[summarised_fourAdep_single_nt$region == region & summarised_fourAdep_single_nt$location == location,]
    
    lower_xlim <- min(df$window) - 1
    upper_xlim <- max(df$window) + 1
    
    if (region == "fpUTR" & location == "Start") {
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
        geom_bar(position = position_dodge(), stat='identity')+
        scale_fill_manual(values=c("#74add1", "#fdae61"))+
        x_and_y_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
      
    }else{
      single_nt_delta_plot_list[[paste(region, location, sep = "_")]] <- ggplot(data = df, aes(x = window, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
        geom_bar(position = position_dodge(), stat='identity')+
        scale_fill_manual(values=c("#74add1", "#fdae61"))+
        x_only_theme+
        scale_y_continuous(limits = c(lower_delta_ylim, upper_delta_ylim))+
        coord_trans(limx = c(lower_xlim - 0.5, upper_xlim + 0.5))
    }
  }
}

#export figure
start_margin = theme(plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"))
end_margin = theme(plot.margin = unit(c(0.5 ,0, 0.5, 1), "cm"))

fpUTR_Start <- list(single_nt_delta_plot_list[[1]])
fpUTR_End <- list(single_nt_delta_plot_list[[2]])

CDS_Start <- list(single_nt_delta_plot_list[[3]])
CDS_End <- list(single_nt_delta_plot_list[[4]])

tpUTR_Start <- list(single_nt_delta_plot_list[[5]])
tpUTR_End <- list(single_nt_delta_plot_list[[6]])

pdf(file = 'fourAdep_single_nt.pdf', width = 20, height = 3)
grid.arrange(grobs = c(lapply(fpUTR_Start, "+", start_margin), lapply(fpUTR_End, "+", end_margin),
                       lapply(CDS_Start, "+", start_margin), lapply(CDS_End, "+", end_margin),
                       lapply(tpUTR_Start, "+", start_margin), lapply(tpUTR_End, "+", end_margin)),
             nrow = 1)
dev.off()

#export legend
legend_plot <- ggplot(data = summarised_fourAdep_single_nt, aes(x = window, y = average_delta, fill = factor(translation, levels = c("4A-dep", "4A-indep"), ordered = T)))+
  geom_bar(position = position_dodge(), stat='identity')+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  theme(legend.title = element_blank())

delta_legend <- myLegend(legend_plot)
pdf(file = 'forAdep_delta_legend.pdf', height = 1, width = 1)
grid.arrange(delta_legend)
dev.off()

#cytosines only----
#bin data as above but for cytosines only

#fpUTR
reactivity_data[reactivity_data$region == "fpUTR" & reactivity_data$Nucleotide == "C",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>% 
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, fpUTR_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("fpUTR")) -> fpUTR_binned_cytosines

#CDS
reactivity_data[reactivity_data$region == "CDS" & reactivity_data$Nucleotide == "C",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>%
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, CDS_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("CDS")) -> CDS_binned_cytosines

#tpUTR
reactivity_data[reactivity_data$region == "tpUTR" & reactivity_data$Nucleotide == "C",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>%
  filter(Position <= (length - tp_trim)) %>%
  mutate(trimmed_length = rep(max(Position))) %>%
  ungroup() %>%
  mutate(normalised_position = Position / trimmed_length,
         bin = calculate_bins(normalised_position, tpUTR_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("tpUTR")) -> tpUTR_binned_cytosines

binned_cytosines <- bind_rows(fpUTR_binned_cytosines, CDS_binned_cytosines, tpUTR_binned_cytosines)

#summarise data
binned_cytosines %>%
  group_by(region, bin) %>%
  summarise(average_delta = mean(net_change, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() -> summarised_binned_cytosines

#calculate axis limits
raw_ylim <- max(c(summarised_binned_cytosines$average_control,summarised_binned_cytosines$average_hippuristanol), na.rm = T)

#plot
binned_raw_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_binned_cytosines[summarised_binned_cytosines$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      y_only_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
  }else{
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      no_labels_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
  }
}

#export figures
awidth <- max(summarised_binned_cytosines[summarised_binned_cytosines$region == "fpUTR",]$bin)
bwidth <- max(summarised_binned_cytosines[summarised_binned_cytosines$region == "CDS",]$bin)
cwidth <- max(summarised_binned_cytosines[summarised_binned_cytosines$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'all_transcripts_binned_C_only.pdf', width = 20, height = 4)
grid.arrange(binned_raw_plot_list$fpUTR, binned_raw_plot_list$CDS, binned_raw_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()

#adenines only----
#bin data as above but for adenines only

#fpUTR
reactivity_data[reactivity_data$region == "fpUTR" & reactivity_data$Nucleotide == "A",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>%
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, fpUTR_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("fpUTR")) -> fpUTR_binned_adenines

#CDS
reactivity_data[reactivity_data$region == "CDS" & reactivity_data$Nucleotide == "A",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>%
  ungroup() %>%
  mutate(normalised_position = Position / length,
         bin = calculate_bins(normalised_position, CDS_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("CDS")) -> CDS_binned_adenines

#tpUTR
reactivity_data[reactivity_data$region == "tpUTR" & reactivity_data$Nucleotide == "A",] %>%
  spread(key = condition, value = Reactivity) %>%
  mutate(delta = hippuristanol - control) %>%
  group_by(transcript) %>%
  mutate(length = rep(max(Position))) %>%
  filter(Position <= (length - tp_trim)) %>%
  mutate(trimmed_length = rep(max(Position))) %>%
  ungroup() %>%
  mutate(normalised_position = Position / trimmed_length,
         bin = calculate_bins(normalised_position, tpUTR_bins)) %>%
  group_by(transcript, bin) %>%
  summarise(control = mean(control, na.rm = T),
            hippuristanol = mean(hippuristanol, na.rm = T),
            net_change = mean(delta, na.rm = T)) %>%
  ungroup() %>%
  mutate(region = rep("tpUTR")) -> tpUTR_binned_adenines

binned_adenines <- bind_rows(fpUTR_binned_adenines, CDS_binned_adenines, tpUTR_binned_adenines)

#summarise data
binned_adenines %>%
  group_by(region, bin) %>%
  summarise(average_delta = mean(net_change, na.rm = T),
            average_control = mean(control, na.rm = T),
            average_hippuristanol = mean(hippuristanol, na.rm = T)) %>%
  ungroup() -> summarised_binned_adenines

#calculate axis limits
raw_ylim <- max(c(summarised_binned_adenines$average_control,summarised_binned_adenines$average_hippuristanol), na.rm = T)

#plot
binned_raw_plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_binned_adenines[summarised_binned_adenines$region == region,]
  xlim <- max(df$bin)
  
  if(region == "fpUTR") {
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      y_only_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
    
  }else{
    df %>%
      select(bin, average_control, average_hippuristanol) %>%
      gather(key = condition, value = Reactivity, average_control, average_hippuristanol) %>%
      ggplot(aes(x = bin, y = Reactivity, colour = factor(condition)))+
      geom_line(size = 1)+
      no_labels_theme+
      scale_y_continuous(limits = c(0,raw_ylim))+
      coord_trans(limx = c(0.5,xlim + 0.5)) -> binned_raw_plot_list[[region]]
  }
}

#export figures
awidth <- max(summarised_binned_adenines[summarised_binned_adenines$region == "fpUTR",]$bin)
bwidth <- max(summarised_binned_adenines[summarised_binned_adenines$region == "CDS",]$bin)
cwidth <- max(summarised_binned_adenines[summarised_binned_adenines$region == "tpUTR",]$bin)
plot_widths <- c(awidth, bwidth, cwidth)

pdf(file = 'all_transcripts_binned_A_only.pdf', width = 20, height = 4)
grid.arrange(binned_raw_plot_list$fpUTR, binned_raw_plot_list$CDS, binned_raw_plot_list$tpUTR,
             nrow = 1, widths = plot_widths)
dev.off()
