###This script was written by Joseph A. Waldron and produces Figure S2 from Waldron et al. (2020) Genome Biology
###Input data first needs to be generated using the the Shell scripts from this repository (see README file)

###Note that this script requires a lot of memory as the all_rtsc.csv is very large

#Imports----
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library(viridis)
library(corrplot)
library(ggplot2)

#write a function to calculate replicate correlation by transcript
repCorFun <- function(xx)
{
  return(data.frame(COR = cor(xx$a, xx$b)))
}

#write r label function
myR <- function(x, y) {
  myCor <- cor.test(x, y)
  r <- round(as.numeric(myCor$estimate), digits = 3)
  r_label = paste('r =', r)
  return(r_label)
}

#load data----
#coverage data
coverage_data <- read.csv(file = "plus_DMS_coverage_all_replicates.csv", header = T) #generate with SF2_pipeline_2_QC.sh

#rtsc stop counts
transcriptome_data <- read.csv(file = "all_rtsc.csv", header = T) #generate with SF2_pipeline_3a_QC.sh

#plot whole transcriptome scatters----
#the following for loop plots a scatter plot containing correlation coefficients for every comparison between replicates within each sample
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
      x_label <- paste0("Log10(RT stops) ", LETTERS[as.numeric(comparison[1])])
      y_label <- paste0("Log10(RT stops) ", LETTERS[as.numeric(comparison[2])])
      
      
      plot_list[[paste(condition, treatment, comparison[1], comparison[2], sep = "_")]] <- ggplot(data = df, aes(x = a, y = b))+
        stat_binhex(bins=100)+
        scale_fill_viridis('Transcripts', limits=c(0,2e6), breaks = c(5e5, 1e6, 1.5e6))+
        ylim(0,5)+
        xlim(0,5)+
        theme_bw()+
        theme(legend.title = element_blank(),
              axis.text = element_text(size=14), 
              axis.title = element_text(size=18))+
        xlab(x_label)+
        ylab(y_label)+
        geom_abline(color = 'black', size = 0.5)+
        annotate("text", x = 1, y = 5, size=6, label = r_label)
    }
  }
}

#make headings
heading_list <- list()
for (condition in c("Control", "Hippuristanol")) {
  for (treatment in c("DMS_minus", "DMS_plus")) {
    heading <- paste(condition, treatment, sep = " / ")
    heading <- gsub("DMS_minus", "DMS(-)", heading)
    heading <- gsub("DMS_plus", "DMS(+)", heading)
    
    heading_list[[paste(condition, treatment, "lab", sep = "_")]] <- ggplot()+
      geom_point(aes(1,1), colour="white") +
      theme(plot.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
      annotate("text", x = 1, y = 1, size=10, label = heading)
  }
}

pdf(file = "whole_transcriptome_replicate_correlation.pdf", width = 16, height = 20)
grid.arrange(heading_list$Control_DMS_minus_lab,
             plot_list$control_minus_DMS_1_2, plot_list$control_minus_DMS_1_3, plot_list$control_minus_DMS_2_3,
             heading_list$Control_DMS_plus_lab,
             plot_list$control_plus_DMS_1_2, plot_list$control_plus_DMS_1_3, plot_list$control_plus_DMS_2_3,
             heading_list$Hippuristanol_DMS_minus_lab,
             plot_list$hippuristanol_minus_DMS_1_2, plot_list$hippuristanol_minus_DMS_1_3, plot_list$hippuristanol_minus_DMS_2_3,
             heading_list$Hippuristanol_DMS_minus_lab,
             plot_list$hippuristanol_plus_DMS_1_2, plot_list$hippuristanol_plus_DMS_1_3, plot_list$hippuristanol_plus_DMS_2_3,
             nrow = 8, layout_matrix = rbind(rep(1,3), 2:4, rep(5,3), 6:8, rep(9,3), 10:12, rep(13,3), 14:16), heights = c(1,4,1,4,1,4,1,4))
dev.off()

#calculate replicate correlation within each transcript for all DMS(+) samples----
transcript_list <- list()
for (condition in c("control", "hippuristanol")) {
  for (comparison in list(c("1", "2"), c("1", "3"), c("2", "3"))) {
    a <- paste(condition, "plus_DMS", comparison[1], sep = "_")
    b <- paste(condition, "plus_DMS", comparison[2], sep = "_")
    
    df <- transcriptome_data[, c("transcript", "position", a, b)]
    names(df) <- c("transcript", "position", "a", "b")
  
    df %>%
      ddply(.(transcript), repCorFun) %>%
      dplyr::mutate(comparison = rep(paste(LETTERS[as.numeric(comparison[1])], "vs", LETTERS[as.numeric(comparison[2])], sep = " ")),
                    condition = rep(condition)) -> transcript_list[[paste(condition, treatment, comparison[1], comparison[2], sep = "_")]]
  }
}
transcript_data <- do.call("rbind", transcript_list)

#plot transcript level coverage data at increasing coverage thresholds----
transcript_data %>%
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
  
  plot_list[[paste("transcripts", coverage, sep = "_")]] <- ggplot(data = df, aes(x = factor(condition, levels = c("control", "hippuristanol"), labels = c("Ctrl", "Hipp"), ordered = T), y = COR, fill = comparison))+
    scale_fill_manual(values = c("#104BA9", "#FFE500", "#FF2800"))+
    geom_boxplot(width = 0.2, outlier.shape=NA)+
    ylab('Correlation coefficient')+
    ylim(0,1)+
    ggtitle(paste("coverage >", coverage))+
    annotate("text", x = 2.2, y = 0.05, size = 5, label = paste("transcripts =", transcript_n))+
    theme_bw()+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
}

pdf(file = "transcript_replicate_correlation_coverage_thresholds.pdf", height = 12, width = 7)
grid.arrange(plot_list$transcripts_0, plot_list$transcripts_1, plot_list$transcripts_10, plot_list$transcripts_100,
             nrow = 4)
dev.off()

#plot correlation matrix at a coverage threshold of 1----
coverage <- 1

#calculate all comparison correlations
transcriptome_data %>%
  inner_join(coverage_data, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage) -> filtered_data

all_comparisons <- cor(filtered_data[,3:14], use = "complete.obs")

#change sample names
dimnames(all_comparisons)[[1]] <- c("Ctrl DMS(-) A", "Ctrl DMS(-) B", "Ctrl DMS(-) C", "Ctrl DMS(+) A", "Ctrl DMS(+) B", "Ctrl DMS(+) C",
                                    "Hipp DMS(-) A", "Hipp DMS(-) B", "Hipp DMS(-) C", "Hipp DMS(+) A", "Hipp DMS(+) B", "Hipp DMS(+) C")
dimnames(all_comparisons)[[2]] <- c("Ctrl DMS(-) A", "Ctrl DMS(-) B", "Ctrl DMS(-) C", "Ctrl DMS(+) A", "Ctrl DMS(+) B", "Ctrl DMS(+) C",
                                    "Hipp DMS(-) A", "Hipp DMS(-) B", "Hipp DMS(-) C", "Hipp DMS(+) A", "Hipp DMS(+) B", "Hipp DMS(+) C")

#plot
pdf(file = "replicate_correlation_table_cov1.pdf", height = 10, width = 10)
corrplot(all_comparisons, type = "upper", method = "number")
dev.off()

