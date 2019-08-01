#Imports----
library(gridExtra)
library(grid)
library(tidyverse)
library(BBmisc)
library(viridis)

#import functions----
source("N:\\JWALDRON/R_scripts/functions.R")

#import variables----
source("N:\\JWALDRON/R_scripts/structure_seq_variables.R")
wStep <- 3

#set home and working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Sliding_windows/panels')

#write themes----
boxplot_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20))

scatter_theme <- theme_bw()+
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text = element_text(size=18), 
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#load data----
#common data
source("N:\\JWALDRON/R_scripts/Structure_seq_common_data.R")

translation_data %>%
  filter(translation == "4A-dep" | translation == "4A-indep" ) %>%
  select(gene, translation, posterior_probability) -> translation_list

#windows data
windows_data <- list()
for (region in c("fpUTR")) {
  for (wLen in c(15, 20, 25, 30, 40)) {
    windows <- read_csv(file = file.path(home, paste0('raw_data/sliding_windows/csv_files/control_', region, '_hippuristanol_', region, '_', wLen, 'win_', wStep, 'step.csv')), col_names = T)
    windows %>%
      mutate(region = rep(region, nrow(windows)),
             wLen = rep(wLen, nrow(windows))) %>%
      inner_join(coverage_data, by = "transcript") %>%
      filter(control_plus_DMS_coverage > coverage,
             hippuristanol_plus_DMS_coverage > coverage,
             control_minus_DMS_fp_10_coverage > fp_coverage,
             hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) -> windows_data[[paste(region, wLen, sep = "_")]]
  }
}

#wLen comparisons----
#calculate the biggest decrease (min) and biggest increase (max) delta reactivity for each gene
wLen_comparisions_list <- list()
for (region in c("fpUTR")) {
  for (wLen in c(15, 20, 25, 30, 40)) {
    windows <- windows_data[[paste(region, wLen, sep = "_")]]
    
    windows %>%
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
      ungroup() %>%
      mutate(region = rep(region)) -> wLen_comparisions_list[[paste(region, wLen, sep = "_")]]
  }
}
wLen_comparisions <- do.call("rbind", wLen_comparisions_list)

#calculate P values and store in list
biggest_decrease_p_values_list <- list()
biggest_increase_p_values_list <- list()
for (region in c("fpUTR")) {
  for (wLen in c(15, 20, 25, 30, 40)) {
    t <- wilcox.test(data = wLen_comparisions_list[[paste(region, wLen, sep = "_")]], min_net_change ~ translation,
                     conf.int = T, paired = F)
    biggest_decrease_p_values_list[[paste(region, wLen, sep = "_")]] <- data.frame(region = region, wLen = wLen, p = t$p.value, lower =  t$conf.int[1], upper = t$conf.int[2])
    t <- wilcox.test(data = wLen_comparisions_list[[paste(region, wLen, sep = "_")]], max_net_change ~ translation,
                     conf.int = T, paired = F)
    biggest_increase_p_values_list[[paste(region, wLen, sep = "_")]] <- data.frame(region = region, wLen = wLen, p = t$p.value, lower =  t$conf.int[1], upper = t$conf.int[2])
  }
}
biggest_decrease_p_values <- do.call("rbind", biggest_decrease_p_values_list)
biggest_increase_p_values <- do.call("rbind", biggest_increase_p_values_list)

biggest_decrease_adjusted_p_labels <- paste("P =", round(biggest_decrease_p_values$p, digits = 3))
biggest_increase_adjusted_p_labels <- paste("P =", round(biggest_increase_p_values$p, digits = 3))

#plot boxplots
plot_list <- list()
for (region in c("fpUTR")) {
  df <- wLen_comparisions[wLen_comparisions$region == region,]
  
  axis_lims <- myAxisLims(df$min_net_change, 1)
  plot_list[[paste0(region, "_biggest_decrease_wLen_comparisons")]] <- ggplot(data = df, aes(x = factor(wLen), y = min_net_change, fill = translation))+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(outlier.shape=NA)+
    ylab(expression(paste(Delta, " reactivity")))+
    ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
    xlab('window length')+
    boxplot_theme+
    annotate("text", x = 1, y = axis_lims$upper_lim, size=4, label = biggest_decrease_adjusted_p_labels[1])+
    annotate("text", x = 2, y = axis_lims$upper_lim, size=4, label = biggest_decrease_adjusted_p_labels[2])+
    annotate("text", x = 3, y = axis_lims$upper_lim, size=4, label = biggest_decrease_adjusted_p_labels[3])+
    annotate("text", x = 4, y = axis_lims$upper_lim, size=4, label = biggest_decrease_adjusted_p_labels[4])+
    annotate("text", x = 5, y = axis_lims$upper_lim, size=4, label = biggest_decrease_adjusted_p_labels[5])

  #biggest increase
  axis_lims <- myAxisLims(df$max_net_change, 1)
  plot_list[[paste0(region, "_biggest_increase_wLen_comparisons")]] <- ggplot(data = df, aes(x = factor(wLen), y = max_net_change, fill = translation))+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(outlier.shape=NA)+
    ylab(expression(paste(Delta, " reactivity")))+
    ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
    xlab('window length')+
    boxplot_theme+
    annotate("text", x = 1, y = axis_lims$upper_lim, size=4, label = biggest_increase_adjusted_p_labels[1])+
    annotate("text", x = 2, y = axis_lims$upper_lim, size=4, label = biggest_increase_adjusted_p_labels[2])+
    annotate("text", x = 3, y = axis_lims$upper_lim, size=4, label = biggest_increase_adjusted_p_labels[3])+
    annotate("text", x = 4, y = axis_lims$upper_lim, size=4, label = biggest_increase_adjusted_p_labels[4])+
    annotate("text", x = 5, y = axis_lims$upper_lim, size=4, label = biggest_increase_adjusted_p_labels[5])
  
  pdf(file = paste0(region, '/wLen_comparisons.pdf'), height = 10, width = 7)
  grid.arrange(plot_list[[paste0(region, "_biggest_decrease_wLen_comparisons")]],
               plot_list[[paste0(region, "_biggest_increase_wLen_comparisons")]],
               nrow = 2)
  dev.off()
}

#windows positions----
#set constant window size and minimum UTR/CDS length
wLen <- 20
min_length <- 50

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

#write a function to plot delta reactivity violins
writeDeltaViolin <- function(df) {
  t <- wilcox.test(data = df, net_change ~ translation, paired = F)
  p_label <- myP_numeric(t)
  axis_lims <- myAxisLims(df$net_change, 1)
  plot <- ggplot(data = df, aes(x = factor(translation, levels = c("4A-dep", "4A-indep"), 
                                           labels = c("4A-dep", "4A-indep\n(matched)"), ordered = T), y = net_change, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab(expression(paste(Delta, " reactivity")))+
    ylim(c(axis_lims$lower_lim, axis_lims$upper_lim))+
    theme_bw()+
    theme(legend.position='none', 
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(hjust = 1, vjust = 1, size = 14, face = "bold"))+
    ggtitle(p_label)+
    stat_summary(fun.y=mean, geom='point', shape=16, size=4)
  return(plot)
}

#write a function to plot position violins
writePositionViolin <- function(df) {
  t <- wilcox.test(data = df, avg_position ~ translation, paired = F)
  p_label <- myP_numeric(t)
  
  plot <- ggplot(data = df, aes(x = translation, y = avg_position, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab('position')+
    ylim(c(0, 1.05))+
    theme_bw()+
    theme(legend.position='none', 
          axis.text.y = element_text(size=18), 
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 1, vjust = 1, size = 14, face = "bold"))+
    ggtitle(p_label)+
    coord_flip()
  return(plot)
}

for (region in c("fpUTR")) {
  windows <- windows_data[[paste(region, wLen, sep = "_")]]
  
  windows %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    filter(length > min_length) %>%
    group_by(gene) %>%
    top_n(n = -1, wt = net_change) %>%
    mutate(max_steps = ceiling(((length - wLen) + 1) / wStep),
           position = (step + 1) / max_steps) %>%
    summarise(avg_position = mean(position, na.rm = T),
              length = mean(length, na.rm = T),
              net_change = mean(net_change)) %>%
    inner_join(translation_list, by = "gene") %>%
    mutate(region = rep(region)) -> biggest_decrease
  
  windows %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    filter(length > min_length) %>%
    group_by(gene) %>%
    top_n(n = 1, wt = net_change) %>%
    mutate(max_steps = ceiling(((length - wLen) + 1) / wStep),
           position = (step + 1) / max_steps) %>%
    summarise(avg_position = mean(position, na.rm = T),
              length = mean(length, na.rm = T),
              net_change = mean(net_change)) %>%
    inner_join(translation_list, by = "gene") %>%
    mutate(region = rep(region)) -> biggest_increase
  
  matched_data <- lapply(list(biggest_decrease, biggest_increase), makeMatchedGroup)
  position_plots <- lapply(matched_data, writePositionViolin)
  delta_plots <- lapply(matched_data, writeDeltaViolin)
  
  pdf(file = paste0(region, '/position_violins.pdf'), height = 10, width = 7)
  grid.arrange(position_plots[[1]],
               position_plots[[2]],
               nrow = 2)
  dev.off()
  
  pdf(file = paste0(region, '/matched_delta_violins.pdf'), height = 4, width = 4)
  grid.arrange(delta_plots[[1]])
  dev.off()
  
}

#window delta reactivity vs length----
#write function to make scatter plots for length vs delta reactivity
writeLengthScatter <- function(df) {
  r <- cor.test(x = df$avg_net_change, y = df$length)
  r_label <- myR(x = r)
  y_axis_lims <- myAxisLims(df$avg_net_change, 1)
  
  if (region == "fpUTR") {
    plot <- ggplot(data = df, aes(x = length, y = avg_net_change))+
      stat_binhex(bins=40)+
      scale_fill_viridis('Transcripts')+
      ylim(c(y_axis_lims$lower_lim, y_axis_lims$upper_lim))+
      scale_x_log10(limits = c(50, 1000), breaks = c(100, 1000))+
      scatter_theme+
      xlab("5\'UTR length")+
      ylab(expression(paste(Delta, " reactivity")))+
      ggtitle(r_label)
  } else {
    plot <- ggplot(data = df, aes(x = length, y = avg_net_change))+
      stat_binhex(bins=40)+
      scale_fill_viridis('Transcripts')+
      ylim(c(y_axis_lims$lower_lim, y_axis_lims$upper_lim))+
      scale_x_log10(limits = c(50, 10000), breaks = c(100, 1000, 10000))+
      scatter_theme+
      xlab(paste0(region, " length"))+
      ylab(expression(paste(Delta, " reactivity")))+
      ggtitle(r_label)
  }
  return(plot)
}

#write function to make scatter plots for GC_content vs delta reactivity
writeGC_contentScatter <- function(df) {
  r <- cor.test(x = df$avg_net_change, y = df$GC_content)
  r_label <- myR(x = r)
  y_axis_lims <- myAxisLims(df$avg_net_change, 1)
  
  plot<- ggplot(data = df, aes(x = GC_content, y = avg_net_change))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    ylim(y_axis_lims$lower_lim, y_axis_lims$upper_lim)+
    xlim(c(0.35, 1))+
    scatter_theme+
    xlab("5\'UTR GC content")+
    ylab(expression(paste(Delta, " reactivity")))+
    ggtitle(r_label)
}


for (region in c("fpUTR")) {
  windows <- windows_data[[paste(region, wLen, sep = "_")]]
  
  windows %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    group_by(gene) %>%
    top_n(n = -1, wt = net_change) %>%
    summarise(avg_net_change = mean(net_change, na.rm = T),
              length = mean(length, na.rm = T),
              GC_content = mean(GC_content, na.rm = T)) %>%
    mutate(direction = rep("biggest_decrease")) -> biggest_decrease_data
  
  windows %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    group_by(gene) %>%
    top_n(n = 1, wt = net_change) %>%
    summarise(avg_net_change = mean(net_change, na.rm = T),
              length = mean(length, na.rm = T),
              GC_content = mean(GC_content, na.rm = T)) %>%
    mutate(direction = rep("biggest_increase")) -> biggest_increase_data
  
  length_plots <- lapply(list(biggest_decrease_data, biggest_increase_data), writeLengthScatter)
  GC_content_plots <- lapply(list(biggest_decrease_data, biggest_increase_data), writeGC_contentScatter)
  
  #write figure
  pdf(file = paste0(region, '/biggest_decrease_vs_length_scatter.pdf'), height = 4, width = 4)
  print(length_plots[[1]])
  dev.off()
  
  pdf(file = paste0(region, '/biggest_increase_vs_length_scatter.pdf'), height = 4, width = 4)
  print(length_plots[[2]])
  dev.off()
  
  pdf(file = paste0(region, '/biggest_decrease_vs_GC_content_scatter.pdf'), height = 4, width = 4)
  print(GC_content_plots[[1]])
  dev.off()
  
  pdf(file = paste0(region, '/biggest_increase_vs_GC_content_scatter.pdf'), height = 4, width = 4)
  print(GC_content_plots[[2]])
  dev.off()
  
  margin = theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  pdf(file = paste0(region, '/length_and_GC_content_scatters.pdf'), height = 8, width = 9)
  grid.arrange(grobs = c(lapply(list(length_plots[[1]],
                                     length_plots[[2]],
                                     GC_content_plots[[1]],
                                     GC_content_plots[[2]]), "+", margin)),
               
               nrow = 2)
  dev.off()
}

#biggest decrease vs biggest increase----
for (region in c("fpUTR", "CDS", "tpUTR")) {
  windows <- windows_data[[paste(region, wLen, sep = "_")]]
  
  windows %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    group_by(gene) %>%
    top_n(n = -1, wt = net_change) %>%
    summarise(avg_net_change = mean(net_change, na.rm = T)) -> biggest_decrease_data
  
  windows %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    group_by(gene) %>%
    top_n(n = 1, wt = net_change) %>%
    summarise(avg_net_change = mean(net_change, na.rm = T)) -> biggest_increase_data
  
  biggest_decrease_data %>%
    inner_join(biggest_increase_data, by = "gene") -> merged_data
  
  r <- cor.test(x = merged_data$avg_net_change.x, merged_data$avg_net_change.y)
  r_label <- myR(x = r)
  x_axis_lims <- myAxisLims(merged_data$avg_net_change.x, 1)
  y_axis_lims <- myAxisLims(merged_data$avg_net_change.y, 1)
  
  biggest_increase_vs_decrease_scatter <- ggplot(data = merged_data, aes(x = avg_net_change.x, y = avg_net_change.y))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    xlim(x_axis_lims$lower_lim, x_axis_lims$upper_lim)+
    ylim(y_axis_lims$lower_lim, y_axis_lims$upper_lim)+
    theme_bw()+
    theme(legend.position='none', 
          axis.text = element_text(size=14), 
          axis.title = element_text(size=18), 
          plot.title = element_text(hjust = 1, vjust = 0, size=12, face="bold"))+
    xlab(expression(paste("biggest decrease ", Delta, " reactivity")))+
    ylab(expression(paste("biggest increase ", Delta, " reactivity")))+
    ggtitle(r_label)
  pdf(file = paste0(region, '/biggest_increase_vs_decrease_scatter.pdf'), height = 4, width = 4)
  print(biggest_increase_vs_decrease_scatter)
  dev.off()

}


#GC content ratio----
#write a function to plot violins
writeGC_contentRatioViolins <- function(df) {
  plot <- ggplot(data = df, aes(x = translation, y = GC_content_ratio, fill = translation))+
    geom_violin(alpha = 0.5)+
    scale_fill_manual(values=c("#74add1", "#fdae61"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab("GC content ratio")+
    ylim(c(0.5, 1.5))+
    theme_bw()+
    theme(legend.position='none', 
          axis.text.x = element_text(size=18), 
          axis.text.y = element_text(size=14), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=18))+
    stat_summary(fun.y=mean, geom='point', shape=16, size=4)
  return(plot)
}

for (region in c("fpUTR", "CDS", "tpUTR")) {
  windows <- windows_data[[paste(region, wLen, sep = "_")]]
  
  windows_composition <- read_csv(file = file.path(home, paste0('raw_data/sliding_windows/composition_files/control_', region, '_hippuristanol_', region, '_', wLen, 'win_', wStep, 'step_composition.csv')), col_names = T)
  windows_composition %>%
    mutate(step = as.numeric(str_replace(transcript, ".+\\_", "")),
           transcript = str_replace(transcript, "\\_.+", "")) %>%
    rename(window_GC_content = GC_content) %>%
    select(transcript, window_GC_content, step) -> windows_composition
  
  windows %>%
    inner_join(windows_composition, by = c("transcript", "step")) %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    select(gene, transcript, net_change, step, wLen, window_GC_content, length, GC_content) -> merged_data
  
  merged_data %>%
    group_by(gene) %>%
    top_n(n = -1, wt = net_change) %>%
    summarise(window_GC_content = mean(window_GC_content, na.rm = T),
              GC_content = mean(GC_content, na.rm = T),
              GC_content_ratio = window_GC_content/ GC_content) %>%
    inner_join(translation_list, by = "gene") -> biggest_decrease_data
  
  merged_data  %>%
    group_by(gene) %>%
    top_n(n = 1, wt = net_change) %>%
    summarise(window_GC_content = mean(window_GC_content, na.rm = T),
              GC_content = mean(GC_content, na.rm = T),
              GC_content_ratio = window_GC_content/ GC_content) %>%
    inner_join(translation_list, by = "gene") -> biggest_increase_data
  
  plots <- lapply(list(biggest_decrease_data, biggest_increase_data), writeGC_contentRatioViolins)
  
  
  #write figure
  pdf(file = paste0(region, '/GC_ratio_violins_pp', positive_change, '_', no_change, '.pdf'), height = 4, width = 8)
  grid.arrange(plots[[1]], plots[[2]], nrow = 1)
  dev.off()
}

#compare GC and length to average 5'UTR delta reactivity
for (region in c("fpUTR", "CDS", "tpUTR")) {
  if (region != "tpUTR") {
    control_stats = read_csv(file = file.path(home, paste0('raw_data/statistics/control_', region, '_0trim_20minlen_statistics.csv')), col_names = T)
    hippuristanol_stats = read_csv(file = file.path(home, paste0('raw_data/statistics/hippuristanol_', region, '_0trim_20minlen_statistics.csv')), col_names = T)
  } else {
    control_stats = read_csv(file = file.path(home, paste0('raw_data/statistics/control_', region, '_50trim_20minlen_statistics.csv')), col_names = T)
    hippuristanol_stats = read_csv(file = file.path(home, paste0('raw_data/statistics/hippuristanol_', region, '_50trim_20minlen_statistics.csv')), col_names = T)
  }
  
  control_stats %>%
    select(transcript, paste("control", region, "average", sep = "_")) %>%
    rename(control_average = paste("control", region, "average", sep = "_")) -> control_stats
  
  hippuristanol_stats %>%
    select(transcript, paste("hippuristanol", region, "average", sep = "_")) %>%
    rename(hippuristanol_average = paste("hippuristanol", region, "average", sep = "_")) -> hippuristanol_stats
  
  control_stats %>%
    inner_join(hippuristanol_stats, by = "transcript") %>%
    mutate(delta = hippuristanol_average - control_average) %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(FASTA_compositions_list[[region]], by = "transcript") %>%
    inner_join(coverage_data, by = "transcript") %>%
    inner_join(abundance_data, by = "transcript") %>%
    filter(control_plus_DMS_coverage > coverage,
           hippuristanol_plus_DMS_coverage > coverage,
           control_minus_DMS_fp_10_coverage > fp_coverage,
           hippuristanol_minus_DMS_fp_10_coverage >fp_coverage) %>%
    group_by(gene) %>%
    top_n(n = 1, wt = abundance) %>%
    ungroup() -> delta_average
  
  #plot mean delta vs length
  r <- cor.test(x = delta_average$delta, delta_average$length)
  r_label <- myR(x = r)
  axis_lims <- myAxisLims(delta_average$delta, 1)
  
  if(region == "fpUTR") {
    delta_length_scatter <- ggplot(data = delta_average, aes(x = length, y = delta))+
      stat_binhex(bins=40)+
      scale_fill_viridis('Transcripts')+
      ylim(axis_lims$lower_lim, axis_lims$upper_lim)+
      scale_x_log10(limits = c(50, 1000), breaks = c(100, 1000))+
      theme_bw()+
      theme(legend.position='none', 
            axis.text = element_text(size=14), 
            axis.title = element_text(size=18), 
            plot.title = element_text(hjust = 1, vjust = 0, size=12, face="bold"))+
      xlab(paste0(region, " length"))+
      ylab(expression(paste("mean ", Delta, " reactivity")))+
      ggtitle(r_label)
  } else {
    delta_length_scatter <- ggplot(data = delta_average, aes(x = length, y = delta))+
      stat_binhex(bins=40)+
      scale_fill_viridis('Transcripts')+
      ylim(axis_lims$lower_lim, axis_lims$upper_lim)+
      scale_x_log10(limits = c(50, 10000), breaks = c(100, 1000, 10000))+
      theme_bw()+
      theme(legend.position='none', 
            axis.text = element_text(size=14), 
            axis.title = element_text(size=18), 
            plot.title = element_text(hjust = 1, vjust = 0, size=12, face="bold"))+
      xlab(paste0(region, " length"))+
      ylab(expression(paste("mean ", Delta, " reactivity")))+
      ggtitle(r_label)
  }
  
  
  #plot mean delta vs GC_content
  r <- cor.test(x = delta_average$delta, delta_average$GC_content)
  r_label <- myR(x = r)
  axis_lims <- myAxisLims(delta_average$delta, 1)
  
  delta_GC_content_scatter <- ggplot(data = delta_average, aes(x = GC_content, y = delta))+
    stat_binhex(bins=40)+
    scale_fill_viridis('Transcripts')+
    ylim(axis_lims$lower_lim, axis_lims$upper_lim)+
    xlim(c(0.1,0.9))+
    theme_bw()+
    theme(legend.position='none', 
          axis.text = element_text(size=14), 
          axis.title = element_text(size=18), 
          plot.title = element_text(hjust = 1, vjust = 0, size=12, face="bold"))+
    xlab(paste0(region, " GC content"))+
    ylab(expression(paste("mean ", Delta, " reactivity")))+
    ggtitle(r_label)
  
  pdf(file = paste0(region, '/average_delta_vs_length_and_GC_scatters.pdf'), height = 4.5, width = 8)
  grid.arrange(delta_length_scatter, delta_GC_content_scatter, nrow = 1)
  dev.off()
}
