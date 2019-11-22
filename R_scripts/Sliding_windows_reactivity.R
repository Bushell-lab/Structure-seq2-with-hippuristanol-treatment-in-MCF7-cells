###This script was written by Joseph A.Waldron and produces panels 5B-H in Waldron et al. (2019) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accessions GSE134865 and GSE134888 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888

#Imports----
library(gridExtra)
library(grid)
library(tidyverse)
library(BBmisc)
library(viridis)

#set home directory----
home <- '' #this needs to be set to the directory containing the data

#set variables----
#posterior probability thresholds
positive_change <- 0.25
no_change <- 0.02

#filter thresholds
coverage <- 1
fp_coverage <- 1.5

#min 5'UTR length
min_length <- 50

#step size
wStep <- 3

#write functions----
#makes a label with a numeric p value  from the output of either t.test or wilcox.test
myP <- function(p) {
  if (p >= 2.2e-16) {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    }else{
      rounded_p <- round(p, digits = 3)
    }
    p_label = paste0("P = ", rounded_p)
  }else{
    p_label = paste0("P < 2.2e-16")
  }
  return(p_label)
}

#makes a label from the output of a correlation test to include r and p values
myR <- function(x) {
  r <- round(as.numeric(x$estimate), digits = 2)
  p <- as.numeric(x$p.value)
  if (p > 0.001) {
    rounded_p <- round(p, digits = 3)
    p_label <- paste('P =', rounded_p)
  }else{
    if (p<2.2E-16) {
      p_label <- 'P < 2.2e-16'
    } else {
      rounded_p <- formatC(p, format = "e", digits = 2)
      p_label <- paste('P =', rounded_p)
    }
  }
  return(paste0('r = ', r, '\n', p_label))
}

#calculates axis limits to plot data minus the top and bottom n% of values
myAxisLims <- function(x, n) {
  upper_quan <- as.numeric(quantile (x, prob = 1 - (n / 100), na.rm=TRUE))
  lower_quan <- as.numeric(quantile (x, prob = n / 100, na.rm=TRUE))
  return(list(upper_lim = upper_quan, lower_lim = lower_quan))
}

#write themes----
boxplot_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))

flipped_violin_theme <- theme_bw()+
  theme(legend.position='none', 
        axis.text.y = element_text(size=18), 
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 1, vjust = 1, size = 14, face = "bold"))

violin_theme <- theme_bw()+
  theme(legend.position='none', 
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 1, vjust = 1, size = 14, face = "bold"))

scatter_theme <- theme_bw()+
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

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

#translation data
translation_data <- read_tsv(file = file.path(home, 'penn-DOD-gene.mmdiffMCF7'), col_names = T, skip = 1) #download from GSE134888
translation_data %>%
  rename(gene = feature_id) %>%
  mutate(DOD = eta1_1 - eta1_2,
         translation = factor(case_when(posterior_probability > positive_change & DOD < 0 ~ "4A-dep",
                                        posterior_probability < no_change ~ "4A-indep"), levels = c("4A-dep", "4A-indep"), ordered = T)) %>%
  filter(translation == "4A-dep" | translation == "4A-indep" ) %>%
  select(gene, translation, posterior_probability) -> translation_list

#transcript to gene ID
transcript_to_geneID <- read_tsv(file = file.path(home, 'MCF7_2015_transcript_to_gene_map.txt'), col_names = T) #download from GSE134865

#5'UTR FASTA composition data
fpUTR_fasta_composition <- read_csv(file = file.path(home, 'MCF7_2015_fpUTRs_composition.csv'), col_names = T) #download data from GSE134865

#windows data
#uses a for loop to read in sliding window data with window sizes of 15, 20, 25, 30 and 40nt and filter by coverage and 5' coverage
#download data from GSE134865
windows_data <- list()
for (wLen in c(15, 20, 25, 30, 40)) {
  windows <- read_csv(file = file.path(home, paste0('control_fpUTR_hippuristanol_fpUTR_', wLen, 'win_', wStep, 'step.csv')), col_names = T)
  windows %>%
    mutate(wLen = rep(wLen, nrow(windows))) %>%
    inner_join(coverage_data, by = "transcript") %>%
    inner_join(fp_coverage_data, by = "transcript") %>%
    filter(control_plus_DMS_1_coverage > coverage,
           control_plus_DMS_2_coverage > coverage,
           control_plus_DMS_3_coverage > coverage,
           hippuristanol_plus_DMS_1_coverage > coverage,
           hippuristanol_plus_DMS_2_coverage > coverage,
           hippuristanol_plus_DMS_3_coverage > coverage,
           control_minus_DMS_fp_10_coverage > fp_coverage,
           hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) -> windows_data[[wLen]]
}

#wLen comparisons----
#the following for loop calculates the biggest decrease (min) and biggest increase (max) delta reactivity for each gene,
#and then selects the 4A-dep genes, calcualtes the number of 4A-dep genes and then selects the same number of 4A-indep genes based on the lowest posterior propability,
#for each window size 
wLen_comparisions_list <- list()
for (wLen in c(15, 20, 25, 30, 40)) {
  windows_data[[wLen]] %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    group_by(gene) %>%
    summarise(min_net_change = min(net_change, na.rm = T),
              max_net_change = max(net_change, na.rm = T),
              wLen = mean(wLen)) %>%
    inner_join(translation_list, by = "gene") -> df
  
  fourAdep_transcripts <- n_distinct(df$gene[df$translation=="4A-dep"])
  
  df %>%
    group_by(translation) %>%
    top_n(wt = -posterior_probability, n = fourAdep_transcripts) %>%
    ungroup() -> wLen_comparisions_list[[wLen]]
}
wLen_comparisions <- do.call("rbind", wLen_comparisions_list)

#calculate P values and store in list
biggest_decrease_p_values_list <- list()
biggest_increase_p_values_list <- list()
for (wLen in c(15, 20, 25, 30, 40)) {
  t <- wilcox.test(data = wLen_comparisions_list[[wLen]], min_net_change ~ translation, paired = F)
  biggest_decrease_p_values_list[[wLen]] <- data.frame(wLen = wLen, p = t$p.value)
  t <- wilcox.test(data = wLen_comparisions_list[[wLen]], max_net_change ~ translation, paired = F)
  biggest_increase_p_values_list[[wLen]] <- data.frame(wLen = wLen, p = t$p.value)
}
biggest_decrease_p_values <- do.call("rbind", biggest_decrease_p_values_list)
biggest_increase_p_values <- do.call("rbind", biggest_increase_p_values_list)

biggest_decrease_p_labels <- paste("P =", round(biggest_decrease_p_values$p, digits = 3))
biggest_increase_p_labels <- paste("P =", round(biggest_increase_p_values$p, digits = 3))

#plot boxplots
axis_lims <- myAxisLims(wLen_comparisions$min_net_change, 1)
biggest_decrease_wLen_boxplots <- ggplot(data = wLen_comparisions, aes(x = factor(wLen), y = min_net_change, fill = translation))+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
  xlab('window length')+
  boxplot_theme+
  annotate("text", x = 1, y = axis_lims$upper_lim, size=4, label = biggest_decrease_p_labels[1])+
  annotate("text", x = 2, y = axis_lims$upper_lim, size=4, label = biggest_decrease_p_labels[2])+
  annotate("text", x = 3, y = axis_lims$upper_lim, size=4, label = biggest_decrease_p_labels[3])+
  annotate("text", x = 4, y = axis_lims$upper_lim, size=4, label = biggest_decrease_p_labels[4])+
  annotate("text", x = 5, y = axis_lims$upper_lim, size=4, label = biggest_decrease_p_labels[5])

#biggest increase
axis_lims <- myAxisLims(wLen_comparisions$max_net_change, 1)
biggest_increase_wLen_boxplots <- ggplot(data = wLen_comparisions, aes(x = factor(wLen), y = max_net_change, fill = translation))+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
  xlab('window length')+
  boxplot_theme+
  annotate("text", x = 1, y = axis_lims$upper_lim, size=4, label = biggest_increase_p_labels[1])+
  annotate("text", x = 2, y = axis_lims$upper_lim, size=4, label = biggest_increase_p_labels[2])+
  annotate("text", x = 3, y = axis_lims$upper_lim, size=4, label = biggest_increase_p_labels[3])+
  annotate("text", x = 4, y = axis_lims$upper_lim, size=4, label = biggest_increase_p_labels[4])+
  annotate("text", x = 5, y = axis_lims$upper_lim, size=4, label = biggest_increase_p_labels[5])

pdf(file = 'wLen_comparisons.pdf', height = 10, width = 7)
grid.arrange(biggest_decrease_wLen_boxplots, biggest_increase_wLen_boxplots,
             nrow = 2)
dev.off()

#window positions and delta reactivity----
#of 5'UTR length matched groups

#write a function that will create a matched controls group for one feature
myMatchedControl <- function(possible_matches_df, to_match, feature){
  #feature should be a character string with the name of the column from possible_matches_df to match
  #possible_matches_df should be a data.frame of the group with all possible matches
  #to_match should be a numerical value needing to be matched
  #it should be used within a for loop to obtain a match for each row of an existing data frame
  
  #first find the absolute differences between the to_match value and every row of the possible_matches_df
  possible_matches <- possible_matches_df
  possible_matches$feature_diff <- abs(possible_matches[,feature] - to_match)
  
  #rank to find those with the least difference in the variable to your gene
  possible_matches$feature_rank <- rank(possible_matches$feature_diff, ties.method = 'average')
  
  #Sort to get best ranked
  possible_matches <- sortByCol(possible_matches, 'feature_rank', asc = T)
  
  #keep the best match
  closest_match <- possible_matches[1,]
  closest_match_ID <- possible_matches[1, "gene"]
  
  #remove chosen match from possible list
  remaining_df <- possible_matches_df[possible_matches_df$"gene" != closest_match_ID,]
  
  return(list(closest_match, remaining_df))
}

#write a function to run the myMatchedControl function on a df of 4A-dep vs 4A-indep matching by length
makeMatchedGroup <- function(df) {
  df %>%
      filter(translation == "4A-dep") -> fourAdep
  
  df %>%
      filter(translation == "4A-indep") -> fourAindep
    
  #make 4A-dep and indep groups matched for 5'UTR length
  for(i in 1:nrow(fourAdep)){
    if(i==1){
      to_match <- as.data.frame(fourAdep)[i,'length']
      matched_mRNA <- myMatchedControl(to_match = to_match, possible_matches_df = as.data.frame(fourAindep), feature = "length")
      matched_mRNAs_df <- matched_mRNA[[1]]
    }else{
      to_match <- as.data.frame(fourAdep)[i,'length']
      matched_mRNA <- myMatchedControl(to_match = to_match, possible_matches_df = matched_mRNA[[2]], feature = "length")
      matched_mRNAs_df <- rbind(matched_mRNAs_df, matched_mRNA[[1]])
    }
  }
  matched_mRNAs_df %>%
    select(names(fourAdep)) -> matched_fourAindep
  
  return(rbind(fourAdep, matched_fourAindep))
}

#for a window length of 20nt and only 5'UTRs more than the minimum length threshold set at the start,
#filter data to include the window per gene which has the biggest increase or decrease in delta reactivity
#and calculate a relative position of that step within the 5'UTR normalised by 5'UTR length
#as some genes will have more than one window with the same delta reactivity, the average normalised position is calculated

#set constant window size
wLen <- 20

windows_data[[wLen]] %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(fpUTR_fasta_composition, by = "transcript") %>%
  filter(length > min_length) %>%
  group_by(gene) %>%
  top_n(n = -1, wt = net_change) %>%
  mutate(max_steps = ceiling(((length - wLen) + 1) / wStep),
         position = (step + 1) / max_steps) %>%
  summarise(avg_position = mean(position),
            length = mean(length),
            net_change = mean(net_change)) %>%
  inner_join(translation_list, by = "gene") -> biggest_decrease

windows_data[[wLen]] %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(fpUTR_fasta_composition, by = "transcript") %>%
  filter(length > min_length) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = net_change) %>%
  mutate(max_steps = ceiling(((length - wLen) + 1) / wStep),
         position = (step + 1) / max_steps) %>%
  summarise(avg_position = mean(position),
            length = mean(length),
            net_change = mean(net_change)) %>%
  inner_join(translation_list, by = "gene") -> biggest_increase

#match data
matched_data <- lapply(list(biggest_decrease, biggest_increase), makeMatchedGroup)

#plot matched biggest decrease
#positions
t <- wilcox.test(data = matched_data[[1]], avg_position ~ translation, paired = F)
p_label <- myP(t$p.value)

biggest_decrease_positions <- ggplot(data = matched_data[[1]], aes(x = translation, y = avg_position, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab('position')+
  ylim(c(0, 1.05))+
  flipped_violin_theme+
  ggtitle(p_label)+
  coord_flip()

#delta reactivity
t <- wilcox.test(data = matched_data[[1]], net_change ~ translation, paired = F)
p_label <- myP(t$p.value)
axis_lims <- myAxisLims(matched_data[[1]]$net_change, 1)

biggest_decrease_delta <- ggplot(data = matched_data[[1]], aes(x = factor(translation, levels = c("4A-dep", "4A-indep"), 
                                         labels = c("4A-dep", "4A-indep\n(matched)"), ordered = T), y = net_change, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab(expression(paste(Delta, " reactivity")))+
  ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
  violin_theme+
  ggtitle(p_label)+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4)

#plot matched biggest increase
#positions
t <- wilcox.test(data = matched_data[[2]], avg_position ~ translation, paired = F)
p_label <- myP(t$p.value)

biggest_increase_positions <- ggplot(data = matched_data[[2]], aes(x = translation, y = avg_position, fill = translation))+
  geom_violin(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  ylab('position')+
  ylim(c(0, 1.05))+
  flipped_violin_theme+
  ggtitle(p_label)+
  coord_flip()

pdf(file = 'position_violins.pdf', height = 10, width = 7)
grid.arrange(biggest_decrease_positions,
             biggest_increase_positions,
             nrow = 2)
dev.off()

pdf(file = 'matched_delta_violin.pdf', height = 4, width = 4)
grid.arrange(biggest_decrease_delta)
dev.off()

#window delta reactivity vs length and GC content----
#take all windows of length 20 and take the window with the minimum delta reactivity per gene
#in case there are two transcripts from the same gene with the same window, the most abundant is selected
#as some transcripts will have more than one window with the same delta reactivity a random window is then selected
windows_data[[wLen]] %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(fpUTR_fasta_composition, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  group_by(gene) %>%
  top_n(n = -1, wt = net_change) %>%
  top_n(n = 1, wt = abundance) %>%
  sample_n(size = 1) -> biggest_decrease

#length
r <- cor.test(x = biggest_decrease$net_change, y = biggest_decrease$length)
r_label <- myR(x = r)
y_axis_lims <- myAxisLims(biggest_decrease$net_change, 1)

length_vs_delta_scatter <- ggplot(data = biggest_decrease, aes(x = length, y = net_change))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  ylim(c(y_axis_lims$lower_lim, y_axis_lims$upper_lim))+
  scale_x_log10(limits = c(50, 1000), breaks = c(100, 1000))+
  scatter_theme+
  xlab("5\'UTR length")+
  ylab(expression(paste(Delta, " reactivity")))+
  ggtitle(r_label)

pdf(file = 'biggest_decrease_vs_length_scatter.pdf', height = 4, width = 4.5)
print(length_vs_delta_scatter)
dev.off()

#GC_content
r <- cor.test(x = biggest_decrease$net_change, y = biggest_decrease$GC_content)
r_label <- myR(x = r)
y_axis_lims <- myAxisLims(biggest_decrease$net_change, 1)

GC_content_vs_delta_scatter<- ggplot(data = biggest_decrease, aes(x = GC_content, y = net_change))+
  stat_binhex(bins=40)+
  scale_fill_viridis('Transcripts')+
  ylim(y_axis_lims$lower_lim, y_axis_lims$upper_lim)+
  xlim(c(0.35, 1))+
  scatter_theme+
  xlab("5\'UTR GC content")+
  ylab(expression(paste(Delta, " reactivity")))+
  ggtitle(r_label)

pdf(file = 'biggest_decrease_vs_GC_content_scatter.pdf', height = 4, width = 4.5)
print(GC_content_vs_delta_scatter)
dev.off()
