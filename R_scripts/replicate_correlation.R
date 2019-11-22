###This script was written by Joseph A.Waldron and produces panels for Figure S2 in Waldron et al. (2019) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accessions GSE134865 and GSE134888 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888

#Imports
library(gridExtra)
library(grid)
library(plyr)
library(tidyverse)
library(viridis)
library(corrplot)

#set home directory----
home <- '' #this needs to be set to the directory containing the data

#write functions----
#write a function to calculate replicate correlation by for every transcript
repCorFun <- function(xx)
{
  return(data.frame(COR = cor(xx$a, xx$b)))
}

#write a function to make a label from a cor test
myR <- function(x, y) {
  myCor <- cor.test(x, y)
  r <- round(as.numeric(myCor$estimate), digits = 3)
  r_label = paste('r =', r)
  return(r_label)
}

#load data----
#coverage data
coverage_data <- read_csv(file = file.path(home, 'plus_DMS_coverage.csv'), col_names = T) #download from GSE134865

#RT stop counts
#This file needs to be created using rtsc_correlation.py script which can be downloaded from from StructureFold2
#all 12 <.rtsc> files should be used for input.
#Note it will be a very large file
transcriptome_data <- read_csv(file.path(home, "all_rtsc.csv"))

#plot whole transcriptome scatters----
plot_list <- list()
for (condition in c("control", "hippuristanol")) {
  for (treatment in c("minus_DMS", "plus_DMS")) {
    for (comparison in list(c("1", "2"), c("1", "3"), c("2", "3"))) {
      a <- paste(condition, treatment, comparison[1], sep = "_")
      b <- paste(condition, treatment, comparison[2], sep = "_")
      
      df <- transcriptome_data[, c("transcript", "position", a, b)]
      names(df) <- c("transcript", "position", "a", "b")
      df$a <- log10(df$a + 1)
      df$b <- log10(df$b + 1)
      
      r_label <- myR(df$a, df$b)
      
      plot_list[[paste(condition, treatment, comparison[1], comparison[2], sep = "_")]] <- ggplot(data = df, aes(x = a, y = b))+
        stat_binhex(bins=100)+
        scale_fill_viridis('Transcripts', limits=c(0,25000), breaks = c(1, 5000, 10000, 15000, 20000))+
        ylim(0,5)+
        xlim(0,5)+
        theme_bw()+
        theme(legend.title = element_blank(),
              axis.text = element_text(size=14), 
              axis.title = element_text(size=18))+
        xlab(paste0("Log10(RT stops) ", comparison[1]))+
        ylab(paste0("Log10(RT stops) ", comparison[2]))+
        geom_abline(color = 'black', size = 0.5)+
        annotate("text", x = 1, y = 5, size=6, label = r_label)
    }
  }
}

#make headings
for (condition in c("Control", "Hippuristanol")) {
  for (treatment in c("DMS(-)", "DMS(+)")) {
    plot_list[[paste(condition, treatment, "lab", sep = "_")]] <- ggplot()+
      geom_point(aes(1,1), colour="white") +
      theme(plot.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
      annotate("text", x = 1, y = 1, size=10, label = paste(condition, treatment, sep = " / "))
  }
}

pdf(file = "panels/whole_transcriptome_replicate_correlation.pdf", width = 16, height = 20)  
grid.arrange(plot_list$`Control_DMS(-)_lab`,
             plot_list$control_minus_DMS_A_B, plot_list$control_minus_DMS_A_C, plot_list$control_minus_DMS_B_C,
             plot_list$`Control_DMS(+)_lab`,
             plot_list$control_plus_DMS_A_B, plot_list$control_plus_DMS_A_C, plot_list$control_plus_DMS_B_C,
             plot_list$`Hippuristanol_DMS(-)_lab`,
             plot_list$hippuristanol_minus_DMS_A_B, plot_list$hippuristanol_minus_DMS_A_C, plot_list$hippuristanol_minus_DMS_B_C,
             plot_list$`Hippuristanol_DMS(+)_lab`,
             plot_list$hippuristanol_plus_DMS_A_B, plot_list$hippuristanol_plus_DMS_A_C, plot_list$hippuristanol_plus_DMS_B_C,
             nrow = 8, layout_matrix = rbind(rep(1,3), 2:4, rep(5,3), 6:8, rep(9,3), 10:12, rep(13,3), 14:16), heights = c(1,4,1,4,1,4,1,4))
dev.off()

for (condition in c("control", "hippuristanol")) {
  for (treatment in c("minus_DMS", "plus_DMS")) {
    
    a <- paste(condition, treatment, '1', sep = "_")
    b <- paste(condition, treatment, '2', sep = "_")
    c <- paste(condition, treatment, '3', sep = "_")
    
    df <- read_csv(file = file.path(home, paste0("raw_data/replicate_correlation/coverage_0/", paste(a, b, c, sep = "_"), '_correlation.csv')), col_names = T)
    
    names(df) <- c("transcript", "position", "A", "B", "C")
    df$transcript <- str_replace(df$transcript, "\\|.+", "")
    transcriptome_list[[paste(condition, treatment, sep = "_")]] <- df
    
    #calculate correlation on a transcript level
    for (comparison in list(c("A", "B"), c("B", "C"), c("A", "C"))) {
      df %>%
        dplyr::select(transcript, position, comparison[1], comparison[2]) %>%
        dplyr::rename(a = comparison[1],
               b = comparison[2]) %>%
        ddply(.(transcript), repCorFun) %>%
        dplyr::mutate(comparison = rep(paste(comparison[1], "vs", comparison[2], sep = " ")),
               condition = rep(condition),
               treatment = rep(treatment)) -> transcript_list[[paste(condition, treatment, comparison[1], comparison[2], sep = "_")]]
    }
  }
}
#save(file = "R_objects/transcript_list_cov0.RData", transcript_list)
#save(file = "R_objects/transcriptome_list_cov0.RData", transcriptome_list)



transcript_data <- do.call("rbind", transcript_list)

#write out transcript level replicate correlation----
transcript_data %>%
  spread(key = comparison, value = COR) -> spread_transcript_data

spread_transcript_data[spread_transcript_data$condition == "control" & spread_transcript_data$treatment == "minus_DMS",] %>%
  rename(ctrl_minus_DMS_A_vs_B = `A vs B`,
         ctrl_minus_DMS_A_vs_C = `A vs C`,
         ctrl_minus_DMS_B_vs_C = `B vs C`) %>%
  select(-condition, -treatment) -> ctrl_minus_DMS_transcript_data

spread_transcript_data[spread_transcript_data$condition == "control" & spread_transcript_data$treatment == "plus_DMS",] %>%
  rename(ctrl_plus_DMS_A_vs_B = `A vs B`,
         ctrl_plus_DMS_A_vs_C = `A vs C`,
         ctrl_plus_DMS_B_vs_C = `B vs C`) %>%
  select(-condition, -treatment) -> ctrl_plus_DMS_transcript_data

spread_transcript_data[spread_transcript_data$condition == "hippuristanol" & spread_transcript_data$treatment == "minus_DMS",] %>%
  rename(hipp_minus_DMS_A_vs_B = `A vs B`,
         hipp_minus_DMS_A_vs_C = `A vs C`,
         hipp_minus_DMS_B_vs_C = `B vs C`) %>%
  select(-condition, -treatment) -> hipp_minus_DMS_transcript_data

spread_transcript_data[spread_transcript_data$condition == "hippuristanol" & spread_transcript_data$treatment == "plus_DMS",] %>%
  rename(hipp_plus_DMS_A_vs_B = `A vs B`,
         hipp_plus_DMS_A_vs_C = `A vs C`,
         hipp_plus_DMS_B_vs_C = `B vs C`) %>%
  select(-condition, -treatment) -> hipp_plus_DMS_transcript_data

inner_join(ctrl_minus_DMS_transcript_data, ctrl_plus_DMS_transcript_data, by = "transcript") %>%
  inner_join(hipp_minus_DMS_transcript_data, by = "transcript") %>%
  inner_join(hipp_plus_DMS_transcript_data, by = "transcript") -> merged_transcript_data

write.csv(file = file.path(home, "raw_data/replicate_correlation/coverage_0/all_transcripts_replicate_correlation.csv"), merged_transcript_data, row.names = F, quote = F)

#plot transcript level violins with different coverage thresholds----
transcript_data %>%
  filter(treatment == "plus_DMS") %>%
  inner_join(coverage_data, by = "transcript") -> merged_data

plot_list <- list()
for (coverage in c(0, 1, 10, 100)) {
  merged_data %>%
    filter(control_plus_DMS_1_coverage > coverage,
           control_plus_DMS_2_coverage > coverage,
           control_plus_DMS_3_coverage > coverage,
           hippuristanol_plus_DMS_1_coverage > coverage,
           hippuristanol_plus_DMS_2_coverage > coverage,
           hippuristanol_plus_DMS_3_coverage > coverage) -> df
  
  transcript_n <- n_distinct(df$transcript)
  
  plot_list[[paste("transcripts", coverage, sep = "_")]] <- ggplot(data = df, aes(x = condition, y = COR, fill = comparison))+
    #geom_violin(alpha = 0.5)+
    scale_fill_manual(values = c("#104BA9", "#FFE500", "#FF2800"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab('Correlation coefficient')+
    ylim(0,1)+
    ggtitle(paste("coverage >", coverage))+
    annotate("text", x = 2.2, y = 0.05, size = 5, label = paste("transcripts =", transcript_n))+
    theme_bw()+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  
}

pdf(file = "panels/transcript_replicate_correlation_coverage_thresholds.pdf", height = 12, width = 7)
grid.arrange(plot_list$transcripts_0, plot_list$transcripts_1, plot_list$transcripts_10, plot_list$transcripts_100,
             nrow = 4)
dev.off()



#plot correlation matrix----
#read in transcriptome data which has been pre-filtered to include only those transcripts with a coverage > 1 in all replicates
transcriptome_data_cov1 <- read_csv(file = file.path(home, "raw_data/replicate_correlation/coverage_1_all_replicates/all_rtsc_replicate_plus_DMS_cov1.csv"), col_names = T)

#calculate all comparison correlations
all_comparisons <- cor(transcriptome_data_cov1[,3:14], use = "complete.obs")

#change sample names
dimnames(all_comparisons)[[1]] <- c("Ctrl DMS(-) A", "Ctrl DMS(-) B", "Ctrl DMS(-) C", "Ctrl DMS(+) A", "Ctrl DMS(+) B", "Ctrl DMS(+) C",
                                    "Hipp DMS(-) A", "Hipp DMS(-) B", "Hipp DMS(-) C", "Hipp DMS(+) A", "Hipp DMS(+) B", "Hipp DMS(+) C")
dimnames(all_comparisons)[[2]] <- c("Ctrl DMS(-) A", "Ctrl DMS(-) B", "Ctrl DMS(-) C", "Ctrl DMS(+) A", "Ctrl DMS(+) B", "Ctrl DMS(+) C",
                                    "Hipp DMS(-) A", "Hipp DMS(-) B", "Hipp DMS(-) C", "Hipp DMS(+) A", "Hipp DMS(+) B", "Hipp DMS(+) C")

#plot
pdf(file = "panels/replicate_correlation_table_all_samples_DMS_plus_cov1_all_replicates.pdf", height = 10, width = 10)
corrplot(all_comparisons, type = "upper", method = "number")
dev.off()
