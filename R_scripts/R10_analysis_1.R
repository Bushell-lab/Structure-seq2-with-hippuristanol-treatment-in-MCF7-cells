###This script was written by Joseph A. Waldron and in combinatin with R10_analysis_2.R produces Figures 4e-g and Extended data Figures 6e-j in Schmidt et al. (in prep)
###Input data first needs to be generated using the R10_analysis_1.sh shell script from this repository (see README file)

#load packages
library(tidyverse)
library(gridExtra)
library(grid)

#import common variables----
source("Structure_seq_variables.R")

#write functions----
#makes a label with a numeric p value from the output of either t.test or wilcox.test
myP_numeric <- function(p) {
  if (p >= 2.2e-16) {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    }else{
      rounded_p <- round(p, digits = 3)
    }
    p_label = paste0("P = ", rounded_p)
  }else{
    p_label = "P < 2.2e-16"
  }
  return(p_label)
}

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
bars_theme <- theme(axis.text = element_text(size = 18),
                    legend.position="none",
                    axis.title.y = element_text(size = 18),
                    axis.title.x = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(hjust = 0.5, colour="grey20", size=24, face="bold"))

boxplot_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.text = element_text(size=16),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank())

violin_theme <- theme_bw()+
  theme(legend.position='none',
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

density_theme <- theme_bw()+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#common data
source("Structure_seq_common_data.R")

#motif data
#the following lines create a vector to use for the col_types for the read_csv function based on the properties used in the react_composition.py script
size <- 10 #size of motif
fp <- 50 #number of nts included 5' of motif
tp <- 50 #number of nts included 3' of motif
n <- fp + size + tp
c_length <- paste(rep("c", n), collapse = "")
n <- 3 * (fp + size + tp)
d_length <- paste(rep("d", n), collapse = "")
col_types <- paste0("ccdd", c_length, d_length, "dddddddddddd", collapse = "")

#the following for loop reads in the react_composition for non-overlapping 10nt AG motifs within each region and saves the data in a list
data_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- read_csv(file = paste("control", region, "hippuristanol", region, "10_AG_1.0_unique_50fp_50tp.csv", sep = "_"), col_names = T, col_types = col_types)
  df$region <- rep(region, nrow(df))
  data_list[[region]] <- df
}

#sliding windows
#the following for loop reads in the sliding windows reactivity data and fasta composition for each region
#The transcript and step are combined for the reactivity ID
#The fasta composition already has transcript_step as the transcript ID. As only GC content is needed this is all that is kept (as a %)
windows_list <- list()
windows_composition_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  windows_list[[region]] <- read_csv(file = paste("control", region, "hippuristanol", region, "20win_10step.csv", sep = "_"), col_names = T)
  windows_list[[region]]$ID <- str_c(windows_list[[region]]$transcript, windows_list[[region]]$step, sep = "_")
  
  windows_composition_list[[region]] <- read_csv(file = paste(region, "20win_10step_composition.csv", sep = "_"), col_names = T)
  windows_composition_list[[region]] %>%
    select(transcript, GC_content) %>%
    mutate(GC_content = GC_content * 100) -> windows_composition_list[[region]]
}

#filter data----
#the following for loop filters by coverage and 5' coverage (fp_coverage)
#it removes any motif that is positioned less than 50 nt from a UTR/CDS boundary or the 5' or 3' end of the transcript
#it removes any motifs that have less than 3 adenines or guanines (ensuring the R10 motifs are mixtures of both purines)
#it picks the most abundant transcript per gene
#it then calclates the number of motifs and transcripts and saves the data in lists
filtered_list <- list()
n_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  data_list[[region]] %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(coverage_data, by = "transcript") %>%
    inner_join(fp_coverage_data, by = "transcript") %>%
    inner_join(abundance_data, by = "transcript") %>%
    mutate(motif = str_c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)) %>%
    filter(control_plus_DMS_1_coverage > coverage,
           control_plus_DMS_2_coverage > coverage,
           control_plus_DMS_3_coverage > coverage,
           hippuristanol_plus_DMS_1_coverage > coverage,
           hippuristanol_plus_DMS_2_coverage > coverage,
           hippuristanol_plus_DMS_3_coverage > coverage,
           control_minus_DMS_fp_10_coverage > fp_coverage,
           hippuristanol_minus_DMS_fp_10_coverage > fp_coverage,
           fp1 != "-" & tp50 != "-", # this line removes any motif that is positioned less than 50 nt from a UTR/CDS boundary or the 5' or 3' end of the transcript
           str_count(motif, "A") > 2 & str_count(motif, "G") > 2) %>% #this line ensures every motif has at least 3 adenines and 3 guanines
    group_by(gene) %>%
    top_n(n = 1, wt = abundance) %>% #this line picks the most abundant transcript per gene
    ungroup() -> df
  
  n_list[[paste0(region, "_motif_n")]] <- paste0("n=", as.character(nrow(df)), " motifs")
  n_list[[paste0(region, "_transcript_n")]] <- paste0("n=", as.character(n_distinct(df$transcript)), " transcripts")
  
  filtered_list[[region]] <- df
}
print(n_list)

#quantify----
#R10 windows
#the following for loop quantifies the average delta reactivity and GC content for every 20nt window (10nt steps) either side of each motif within each region
#it then saves the data as a data frame in a list. The do.call function then converts the list of data frames into one data frame
fp_vs_tp_R10_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- filtered_list[[region]]
  D1 <- which(colnames(df) == "D1") #identifies the first delta reactivity column in the data
  nt1 <- which(colnames(df) == "fp1") #identifies the first nucleotide column in the data
  
  for (i in c(seq(0, 30, by = 10))) {
    
    #extract delta reactivity
    fp_start <- D1 + i #identifies the first column of the 5' window
    fp_end <- fp_start + 19 #identifies the last column of the 5' window
    
    tp_start <- D1 + fp + size + (30 - i) #identifies the first column of the 3' window
    tp_end <- tp_start + 19 #identifies the last column of the 3' window
    
    fp_delta <- rowMeans(df[fp_start:fp_end], na.rm = T) #calculates the mean reactivity of the 5' window
    tp_delta <- rowMeans(df[tp_start:tp_end], na.rm = T) #calculates the mean reactivity of the 3' window
    
    #extract GC content
    fp_start <- nt1 + i #identifies the first column of the 5' window
    fp_end <- fp_start + 19 #identifies the last column of the 5' window
    
    tp_start <- nt1 + fp + size + (30 - i) #identifies the first column of the 3' window
    tp_end <- tp_start + 19 #identifies the last column of the 3' window
    
    fp_nts <- df[fp_start:fp_end]
    tp_nts <- df[tp_start:tp_end]
    
    fp_nts$motif <-  do.call(paste, c(fp_nts, sep="")) #collapses the columns into one
    tp_nts$motif <-  do.call(paste, c(tp_nts, sep="")) #collapses the columns into one
    
    fp_nts %>%
      mutate(G_count = str_count(motif, "G"),
             C_count = str_count(motif, "C"),
             GC_content = (rowSums(cbind(G_count, C_count)))/20 * 100) %>% #calculates the GC content (as a %)
      pull(GC_content) -> fp_GC_content
    
    tp_nts %>%
      mutate(G_count = str_count(motif, "G"),
             C_count = str_count(motif, "C"),
             GC_content = (rowSums(cbind(G_count, C_count)))/20 * 100) %>% #calculates the GC content (as a %)
      pull(GC_content) -> tp_GC_content
    
    fp_vs_tp_R10_list[[paste(region,as.character(i),sep = "_")]] <- data.frame(transcript = df$transcript,
                                                                            start = df$start,
                                                                            fp_delta = fp_delta,
                                                                            tp_delta = tp_delta,
                                                                            fp_GC_content = fp_GC_content,
                                                                            tp_GC_content = tp_GC_content,
                                                                            region = rep(region),
                                                                            motif = rep("AG"),
                                                                            distance = rep(paste0(as.character(31-i), ":", as.character(50-i))))
  }
}
fp_vs_tp_R10_data <- do.call("rbind", fp_vs_tp_R10_list)

#random windows
#the following for loop uses the sliding windows data to randomly select the same number of windows per transcript as R10 motifs
#and then extract the delta reactivity and GC content of the surrounding windows and save as a data frame in a list
fp_vs_tp_random_list <- list()
random_motifs_list <- list()
set.seed(020588) #set the seed for random window selection
for (region in c("fpUTR", "CDS", "tpUTR")) {
  
  #calculate the number of motifs per transcript
  fp_vs_tp_R10_list[[paste0(region, "_0")]] %>%
    group_by(transcript) %>%
    count() -> motif_counts
  
  windows_list[[region]] %>%
    group_by(transcript) %>%
    mutate(max_step = max(step)) %>%
    filter(transcript %in% fp_vs_tp_R10_list[[paste0(region, "_0")]]$transcript, #ensures only those transcripts with R10 motifs are included
           step > 4 & step < (max_step - 4)) %>% #ensures the window is at least 50 nt away from the 5' and 3' of the UTR/CDS boundary or the 5' or 3 end of the transcript
    inner_join(motif_counts, by = "transcript") %>%
    group_by(transcript) %>%
    sample_n(size = n) %>% #selects the same number of random windows as the number or R10 motifs per transcript
    ungroup() -> random_motifs
  
  random_motifs_list[[region]] <- random_motifs[,c("transcript", "step")] #saves the random motif coordinates in a list
  
  for (i in 1:4) {
    #calculate the adjacent windows
    random_motifs %>%
      mutate(fp_step_ID = str_c(transcript, step - i, sep = "_"),
             tp_step_ID = str_c(transcript, step + i, sep = "_")) -> adjacent_motifs
    
    #extract delta reactivity
    windows_list[[region]] %>%
      filter(ID %in% adjacent_motifs$fp_step_ID) %>%
      mutate(delta = net_change / 20) %>% #the net change is the sum of the delta reactivity at each position, so needs to be divided by window length (20)
      pull(delta) -> fp_delta
    
    windows_list[[region]] %>%
      filter(ID %in% adjacent_motifs$tp_step_ID) %>%
      mutate(delta = net_change / 20) %>% #the net change is the sum of the delta reactivity at each position, so needs to be divided by window length (20)
      pull(delta) -> tp_delta
    
    #extract GC content
    windows_composition_list[[region]] %>%
      filter(transcript %in% adjacent_motifs$fp_step_ID) %>%
      pull(GC_content) -> fp_GC_content
    
    windows_composition_list[[region]] %>%
      filter(transcript %in% adjacent_motifs$tp_step_ID) %>%
      pull(GC_content) -> tp_GC_content
    
    fp_vs_tp_random_list[[paste(region,as.character(i),sep = "_")]] <- data.frame(transcript = random_motifs$transcript,
                                                                                  start = random_motifs$step * 10,
                                                                                  fp_delta = fp_delta,
                                                                                  tp_delta = tp_delta,
                                                                                  fp_GC_content = fp_GC_content,
                                                                                  tp_GC_content = tp_GC_content,
                                                                                  region = rep(region),
                                                                                  motif = rep("random"),
                                                                                  distance = rep(paste0(as.character(((i-1) * 10) + 1), ":", as.character(((i-1) * 10) + 20))))
  }
}
fp_vs_tp_random_data <- do.call("rbind", fp_vs_tp_random_list)
fp_vs_tp_data <- bind_rows(fp_vs_tp_R10_data, fp_vs_tp_random_data)

#plot boxplots----
#sliding windows for AG10 motifs within 5'UTRs delta reactivity
#calculate p-values
p_label_list <- list()
for (distance in c("1:20", "11:30", "21:40", "31:50")) {
  df <- fp_vs_tp_data[fp_vs_tp_data$region == "fpUTR" & fp_vs_tp_data$motif == "AG" & fp_vs_tp_data$distance == distance,]
  
  t <- wilcox.test(df$tp_delta, df$fp_delta, paired = T)
  
  p_label_list[[distance]] <- myP_numeric(t$p.value)
}

#plot
fp_vs_tp_data %>%
  filter(region == "fpUTR" & motif == "AG") %>%
  gather(key = position, value = delta, fp_delta, tp_delta ) %>%
  ggplot(aes(x = factor(distance, levels = c("1:20", "11:30", "21:40", "31:50"), ordered = T),
             y = delta,
             fill = factor(position, levels = c("fp_delta", "tp_delta"), labels = c("5\'", "3\'"), ordered = T)))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  geom_boxplot(outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(-0.2, 0.20))+
  xlab('distance from motif (nt)')+
  boxplot_theme+
  theme(axis.title.x = element_text(size = 20))+
  annotate("text", x = 1, y = 0.2, size=4, label = p_label_list[[1]])+
  annotate("text", x = 2, y = 0.2, size=4, label = p_label_list[[2]])+
  annotate("text", x = 3, y = 0.2, size=4, label = p_label_list[[3]])+
  annotate("text", x = 4, y = 0.2, size=4, label = p_label_list[[4]]) -> fp_vs_tp_delta_box_plot

pdf(file = "Extended_data_Fig_6e_fpUTR_distances_boxplot.pdf", width = 6, height = 4)
print(fp_vs_tp_delta_box_plot)
dev.off()

#5' vs 3' 31:50nts delta reactivity for AG vs random for each region
#calculate p-values
p_label_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (motif  in c("AG", "random")) {
    df <- fp_vs_tp_data[fp_vs_tp_data$distance == "31:50" & fp_vs_tp_data$region == region & fp_vs_tp_data$motif == motif,]
    
    t <- wilcox.test(df$tp_delta, df$fp_delta, paired = T)
    
    p_label_list[[paste(region, motif, sep = "_")]] <- myP_numeric(t$p.value)
  }
}

#plot
fp_vs_tp_data %>%
  filter(distance == "31:50") %>%
  gather(key = position, value = delta, fp_delta, tp_delta ) %>%
  mutate(position = factor(position, levels = c("fp_delta", "tp_delta"), labels = c("5\'", "3\'"), ordered = T),
         x_label = str_c(region, motif, sep = "_"),
         x_label = factor(x_label, levels = c("fpUTR_AG", "fpUTR_random", "CDS_AG", "CDS_random", "tpUTR_AG", "tpUTR_random"),
                          labels = c("5\'UTR\nAG", "5\'UTR\nrandom", "CDS\nAG", "CDS\nrandom", "3\'UTR\nAG", "3\'UTR\nrandom"), ordered = T)) %>%
  ggplot(aes(x = x_label, y = delta, fill = position))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  geom_boxplot(outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(-0.2, 0.2))+
  boxplot_theme+
  annotate("text", x = 1, y = 0.2, size=4, label = p_label_list[[1]])+
  annotate("text", x = 2, y = 0.2, size=4, label = p_label_list[[2]])+
  annotate("text", x = 3, y = 0.2, size=4, label = p_label_list[[3]])+
  annotate("text", x = 4, y = 0.2, size=4, label = p_label_list[[4]])+
  annotate("text", x = 5, y = 0.2, size=4, label = p_label_list[[5]])+
  annotate("text", x = 6, y = 0.2, size=4, label = p_label_list[[6]]) -> fp_vs_tp_delta_box_plot

pdf(file = "Fig_4e_31_50_boxplot.pdf", width = 8, height = 4)
print(fp_vs_tp_delta_box_plot)
dev.off()

#1:20 vs 31:50nts for AG and random for each region
#calculate p-values
p_label_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  for (motif  in c("AG", "random")) {
    df <- fp_vs_tp_data[fp_vs_tp_data$region == region & fp_vs_tp_data$motif == motif,]
    
    df %>%
      select(transcript, start, tp_delta, distance) %>%
      filter(distance == "1:20" | distance == "31:50") %>%
      spread(key = distance, value = tp_delta, ) -> df
    
    t <- wilcox.test(df$`1:20`, df$`31:50`, paired = T)
    
    p_label_list[[paste(region, motif, sep = "_")]] <- myP_numeric(t$p.value)
  }
}

#plot
fp_vs_tp_data %>%
  filter(distance == "1:20" | distance == "31:50") %>%
  mutate(distance = factor(distance, levels = c("1:20", "31:50"), ordered = T),
         x_label = str_c(region, motif, sep = "_"),
         x_label = factor(x_label, levels = c("fpUTR_AG", "fpUTR_random", "CDS_AG", "CDS_random", "tpUTR_AG", "tpUTR_random"),
                          labels = c("5\'UTR\nAG", "5\'UTR\nrandom", "CDS\nAG", "CDS\nrandom", "3\'UTR\nAG", "3\'UTR\nrandom"), ordered = T)) %>%
  ggplot(aes(x = x_label, y = tp_delta, fill = distance))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  geom_boxplot(outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(-0.2, 0.2))+
  boxplot_theme+
  annotate("text", x = 1, y = 0.2, size=4, label = p_label_list[[1]])+
  annotate("text", x = 2, y = 0.2, size=4, label = p_label_list[[2]])+
  annotate("text", x = 3, y = 0.2, size=4, label = p_label_list[[3]])+
  annotate("text", x = 4, y = 0.2, size=4, label = p_label_list[[4]])+
  annotate("text", x = 5, y = 0.2, size=4, label = p_label_list[[5]])+
  annotate("text", x = 6, y = 0.2, size=4, label = p_label_list[[6]]) -> fp_vs_tp_delta_box_plot

pdf(file = "Extended_data_Fig_6f_prox_vs_distal_tp_boxplot.pdf", width = 8, height = 4)
print(fp_vs_tp_delta_box_plot)
dev.off()

#bars----
#summarise
#the following for loop calculates the mean delta reactivity at each position from -50:+50 of every AG motif in each region
#the mean delta reactivity for the whole motif is collapsed into a single value (position 0)
summarised_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- filtered_list[[region]]
  
  col_d_start <- which(colnames(df) == "D1") #identifies the first delta reactivity column (-50)
  col_d_end <- col_d_start + fp + size + tp - 1 #identifies the last delta reactivity column (+50)
  
  df %>%
    select(col_d_start:col_d_end) %>%
    summarise_all(funs(mean(., na.rm = TRUE))) %>%
    gather() %>%
    mutate(position = c(-fp:-1, rep(0, size), 1:tp)) %>%
    group_by(position) %>%
    summarise(value = mean(value, na.rm = T)) %>%
    mutate(region = rep(region)) -> df
  
  summarised_list[[region]] <- df
}
summarised_data <- do.call("rbind", summarised_list)

#plot
upper_axis_lim <- (max(summarised_data$value, na.rm = T))
lower_axis_lim <- (min(summarised_data$value, na.rm = T))

plot_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- summarised_list[[region]]
  
  plot_list[[region]] <- ggplot(data = df, aes(x = position, y = value))+
    geom_bar(stat = "identity", aes(fill = value))+
    scale_fill_gradient2(low = 'blue', mid = 'grey88', high = 'red', space = 'Lab', limits=c(lower_axis_lim, upper_axis_lim))+
    bars_theme+
    scale_y_continuous(limits = c(lower_axis_lim, upper_axis_lim))+
    coord_trans(limy = c(lower_axis_lim, upper_axis_lim),
                limx = c(-50.5, 50.5))+
    geom_vline(xintercept = -0.5, colour="black", size=0.05)+
    geom_vline(xintercept = 0.5, colour="black", size=0.05)+
    ylab(expression(paste(Delta, " reactivity")))
}

pdf(file = "Extended_data_Fig_6g_bars.pdf", width = 20, height = 4)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
             nrow = 1)
dev.off()

#windows GC content----
#31:50 downstream of AG10 vs random motifs in 5'UTRs
fp_vs_tp_data %>%
  filter(distance == "31:50",
         region == "fpUTR") -> df

t <- wilcox.test(data = df, tp_GC_content ~ motif,
                 paired = F, alternative = "two.sided", var.equal = F, conf.int = T)
p_label <- myP(t)

df %>%
  ggplot(aes(x = motif, y = tp_GC_content, fill = motif))+
  geom_violin(alpha = 0.5, fill = "#00BFC4")+
  geom_boxplot(width = 0.2, outlier.shape=NA, fill = "#00BFC4")+
  ylab("GC content (%)")+
  ylim(c(0, 100))+
  violin_theme+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4)+
  ggtitle(p_label) -> windows_GC_content_violin

#whole 5'UTR GC content----
#the following pipe identifies all transcripts with R10 motifs (with at least 3 adenines and guanines) in their 5'UTRs
data_list$fpUTR %>%
  mutate(motif = str_c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)) %>%
  filter(str_count(motif, "A") > 2 & str_count(motif, "G") > 2) %>%
  group_by(transcript) %>%
  sample_n(size = 1) %>% #ensures each transcript is only in the list once (some transcripts have more than 1 R10 motif)
  pull(transcript) -> R10_IDs

#the following pipe takes all transcripts with 5'UTRs at least 100nt in length, selects the most abundant transcript and
#creates a factor depending on whether they have an R10 motif in their 5'UTR or not
FASTA_compositions_list$fpUTR %>%
  filter(length > 100) %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  mutate(GC_content = GC_content * 100,
         AG10 = factor(case_when(transcript %in% R10_IDs ~ "R10\n5'UTRs",
                                 !(transcript %in% R10_IDs) ~ "non-R10\n5'UTRs"),
                       levels = c("R10\n5'UTRs", "non-R10\n5'UTRs"), ordered = T)) %>%
  select(transcript, AG10, GC_content) -> fpUTR_GC_content_data

#make a size matched group of non-R10 mRNAs
R10_n <- nrow(fpUTR_GC_content_data[fpUTR_GC_content_data$AG10 == "R10\n5'UTRs",])

set.seed(020588) #sets seed for random selection of size matched group

fpUTR_GC_content_data %>%
  group_by(AG10) %>%
  sample_n(size = R10_n) %>%
  ungroup() -> fpUTR_GC_content_size_matched

#print group sizes and output from wilcoxon test
summary(fpUTR_GC_content_size_matched)

t <- wilcox.test(data = fpUTR_GC_content_size_matched, GC_content ~ AG10,
                 paired = F, alternative = "two.sided", var.equal = F, conf.int = T)
p_label <- myP(t)

#plot
fpUTR_GC_content_size_matched %>%
  ggplot(aes(x = AG10, y = GC_content, fill = AG10))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  ylab("GC content (%)")+
  ylim(c(0, 100))+
  violin_theme+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4)+
  ggtitle(p_label) -> whole_fpUTR_GC_content_violin

#windows MFE----
#generate transcript IDs
#the following pipe generates a list of transcript IDs that can be used to fold windows +16:65nt downstream of the motifs.
#the ID will be the same whatever the -fp and -tp settings used in the react_composition.py script
filtered_list[["fpUTR"]] %>%
  mutate(ID = str_c(transcript, query, start, end, sep = "_")) %>%
  pull(ID) -> filtered_R10_IDs
write.table(filtered_R10_IDs, file = "filtered_10_AG_1.0_IDs.txt", col.names = F, row.names = F, quote = F)

#the following pipe converts the random IDs into the ID that would be associated with the +16:65nt window from the output of the react_windows.py script with wlen=50 and 2step=5
random_motifs_list$fpUTR %>%
  mutate(new_step = (step * 2) + 5, #this will generate the step for the +16:65 downstream window of the random motif with wlen=50 and wstep=5
         min_length = 50 + (new_step * 5), #this calculates the min length of a 5'UTR for the new downstream +16:65nt step
         ID = str_c(transcript, new_step, sep = "_")) %>%
  inner_join(FASTA_compositions_list$fpUTR, by = "transcript") %>%
  filter(length >= min_length) %>% #this removes any random motifs that do not have 65nt downstream
  pull(ID) -> random_IDs
write.table(random_IDs, file = "filtered_random_IDs.txt", col.names = F, row.names = F, quote = F)

#run R10_analysis_2.sh then proceed to R10_analysis_2.R but keep all variables in the environment

