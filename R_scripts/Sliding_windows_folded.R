###This script was written by Joseph A.Waldron and produces panels 1F-G in Waldron et al. (2019) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accessions GSE134865 and GSE134888 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888

###Imports
library(gridExtra)
library(grid)
library(tidyverse)

#set home directory----
home <- '' #this needs to be set to the directory containing the data

#set variables----
#filter thresholds
coverage <- 1
fp_coverage <- 1.5
positive_change <- 0.25
min_length <- 100 #sets minimum 5'UTR length
wStep <- 10 #sets window step size
wLen <- 50 #sets window size
bin_n <- 10 #sets number of bins

#write functions----
#calculate_bins
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

#exports just the legend of a plot
myLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#create theme----
myTheme <- theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

#load data----
#coverage data
coverage_data <- read_csv(file = file.path(home, 'plus_DMS_coverage.csv'), col_names = T) #download from GSE134865
ctrl_fp_coverage_data <- read_csv(file = file.path(home, 'control_minus_DMS_fp_10_coverage.csv'), col_names = T) #download from GSE134865
hipp_fp_coverage_data <- read_csv(file = file.path(home, 'hippuristanol_minus_DMS_fp_10_coverage.csv'), col_names = T) #download from GSE134865
fp_coverage_data <- inner_join(ctrl_fp_coverage_data, hipp_fp_coverage_data, by = "transcript")
rm(ctrl_fp_coverage_data, hipp_fp_coverage_data)

#totals data
totals_data <- read_tsv(file = file.path(home, 'penn-DE.mmdiffMCF7'), col_names = T, skip = 1) #download from GSE134888
totals_data %>%
  mutate(abundance = case_when(posterior_probability > positive_change ~ alpha1,
                               posterior_probability < positive_change ~ alpha0)) %>%
  rename(transcript = feature_id) %>%
  select(transcript, abundance) -> abundance_data
rm(totals_data)

#transcript to gene ID
transcript_to_geneID <- read_tsv(file = file.path(home, 'MCF7_2015_transcript_to_gene_map.txt'), col_names = T) #download from GSE134865

#5'UTR FASTA composition data
fpUTR_fasta_composition <- read_csv(file = file.path(home, 'MCF7_2015_fpUTRs_composition.csv'), col_names = T) #download data from GSE134865

#windows data
#read in fasta composition of all window sequences and seperate transcript and step into two different variables
windows_FASTA <- read_csv(file = file.path(home, paste('fpUTR', wLen, 'win', wStep, 'step_composition.csv', sep = "_")), col_names = T) #download data from GSE134865
windows_FASTA %>%
  mutate(step = as.integer(str_replace(transcript, ".+\\_", "")),
         transcript = str_replace(transcript, "\\_.+", "")) -> windows_FASTA

#use a for loop to read in MFEs of all window sequences under each condition, seperate transcript and step into two different variables
#and filter by coverage and 5' coverage and select the most abundant transcript per gene
MFE_data_list <- list()
for (condition in c("control", "hippuristanol")) {
  df <- read_csv(file = file.path(home, paste(condition, 'fpUTR', wLen, 'win', wStep, 'step_MFE.csv', sep = "_")), col_names = T)#download data from GSE134865
  
  df %>%
    mutate(step = as.integer(str_replace(transcript, ".+\\_", "")),
           transcript = str_replace(transcript, "\\_.+", ""),
           condition = factor(rep(condition))) %>%
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
    select(transcript, Free_Energy, step, condition) -> MFE_data_list[[condition]]
}
MFE_data <- do.call("rbind", MFE_data_list)

#merge data----
MFE_data %>%
  inner_join(fpUTR_fasta_composition[, c("transcript", "length")], by = "transcript") %>%
  rename(fpUTR_length = length) %>%
  inner_join(windows_FASTA, by = c("transcript", "step")) %>%
  mutate(GC_content = GC_content * 100) %>%
  rename(windows_GC_content = GC_content,
         wLen = length) %>%
  select(transcript, Free_Energy, step, condition, fpUTR_length, wLen, windows_GC_content) -> merged_data

#bin data----
merged_data %>%
  mutate(position = 1 + (step * wStep) + round((wLen / 2)),
         normalised_position = position / fpUTR_length,
         bin = calculate_bins(normalised_position, bin_n)) %>%
  group_by(condition, transcript, bin) %>%
  summarise(Free_Energy = mean(Free_Energy),
            windows_GC_content = mean(windows_GC_content)) %>%
  ungroup() -> binned_data

#plot GC content----
#summarise data and calculate 95% confidence intervals
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

#plot
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
#summarise data and calculate 95% confidence intervals
summarised_MFE_data_list <- list()
for (condition in c("control", "hippuristanol")) {
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

#plot
MFE_plot <- ggplot(data = summarised_MFE_data, aes(x = bin))+
  geom_line(aes(y = Free_Energy_mean, colour = factor(condition, levels = c("control", "hippuristanol"), ordered = T)), size = 1)+
  geom_ribbon(aes(ymin = Free_Energy_lower, ymax = Free_Energy_upper, fill = condition), alpha = 0.3, colour = NA)+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol)")))+
  myTheme

pdf(file = "sliding_windows_fpUTR_MFE_all_transcripts.pdf", width = 4, height = 2)
print(MFE_plot)
dev.off()

#export legend
legend_plot <- ggplot(data = summarised_MFE_data, aes(x = bin, y = Free_Energy_mean, colour = factor(condition, levels = c("control", "hippuristanol"),
                                                                                                             labels = c("Ctrl", "Hipp"), ordered = T)))+
  geom_line()+
  theme_bw()+
  theme(legend.title = element_blank())

ctrl_hipp_legend <- myLegend(legend_plot)
pdf(file = 'ctrl_hipp_legend.pdf', height = 1, width = 1)
grid.arrange(ctrl_hipp_legend)
dev.off()
