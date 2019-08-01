###This script uses a lot of memory (depending on the number of transcripts analysed). It would not run without generating a memory error (not enough RAM) 
#on my laptop so had to be run on the VM

#Imports
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library(ggplot2)
library("hexbin")


#set home and working directory
setwd('/mnt/data/JWALDRON/Papers/my_papers/Structure-seq/Figures/S3_structure_seq_QC/panels')
home <- '/mnt/data/JWALDRON/Structure_seq/MCF7_2015'

#set coverage
coverage <- 100

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
transcriptome_list <- list()
transcript_list <- list()

for (condition in c("control", "hippuristanol")) {
  for (treatment in c("minus_DMS", "plus_DMS")) {
    
    a <- paste(condition, treatment, '1', sep = "_")
    b <- paste(condition, treatment, '2', sep = "_")
    c <- paste(condition, treatment, '3', sep = "_")
    
    df <- read.csv(file = file.path(home, paste0('raw_data/replicate_correlation/coverage_', coverage, "/", paste(a, b, c, sep = "_"), '_correlation.csv')), header = T)
    
    names(df) <- c("transcript", "position", "A", "B", "C")
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

transcript_data <- do.call("rbind", transcript_list)

#plot violins----
plot_list <- list()
for (condition in c("control", "hippuristanol")) {
  for (treatment in c("minus_DMS", "plus_DMS")) {
    df <- transcript_data[transcript_data$condition == condition & transcript_data$treatment == treatment,]
    plot_list[[paste(condition, treatment, "transcript", sep = "_")]] <- ggplot(data = df, aes(x = factor(comparison), y = COR, fill = factor(comparison)))+
      geom_violin(alpha = 0.5)+
      scale_fill_manual(values = c("#104BA9", "#FFE500", "#FF2800"))+
      geom_boxplot(width = 0.2, outlier.shape=NA)+
      ylab('Correlation coefficient')+
      ylim(0,1)+
      theme_bw()+
      theme(legend.position='none',
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=14),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=18))
  }
}

#plot scatters----
for (condition in c("control", "hippuristanol")) {
  for (treatment in c("minus_DMS", "plus_DMS")) {
    for (comparison in list(c("A", "B"), c("B", "C"), c("A", "C"))) {
      transcriptome_list[[paste(condition, treatment, sep = "_")]] %>%
        dplyr::mutate(A = log10(A + 1),
               B = log10(B + 1),
               C = log10(C + 1)) %>%
        dplyr::select(transcript, position, comparison[1], comparison[2]) %>%
        dplyr::rename(a = comparison[1],
               b = comparison[2]) -> df
      
      r_label <- myR(df$a, df$b)
      
      plot_list[[paste(condition, treatment, comparison[1], comparison[2], sep = "_")]] <- ggplot(data = df, aes(x = a, y = b))+
        stat_binhex(bins=40)+
        scale_fill_distiller(palette='Spectral')+
        ylim(0,5)+
        xlim(0,5)+
        theme_bw()+
        theme(legend.position='none', 
              axis.text = element_text(size=14), 
              axis.title = element_text(size=18))+
        xlab(paste0("Log10(RT stops) ", comparison[1]))+
        ylab(paste0("Log10(RT stops) ", comparison[2]))+
        geom_abline(color = 'black', size = 0.5)+
        annotate("text", x = 1, y = 5, size=6, label = r_label)
    }
  }
}

#make headings----
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

#write figure----
pdf(file = paste0("replicate_correlation_cov", coverage, ".pdf"), width = 16, height = 20)  
grid.arrange(plot_list$`Control_DMS(-)_lab`,
             plot_list$control_minus_DMS_A_B, plot_list$control_minus_DMS_A_C, plot_list$control_minus_DMS_B_C, plot_list$control_minus_DMS_transcript,
             plot_list$`Control_DMS(+)_lab`,
             plot_list$control_plus_DMS_A_B, plot_list$control_plus_DMS_A_C, plot_list$control_plus_DMS_B_C, plot_list$control_plus_DMS_transcript,
             plot_list$`Hippuristanol_DMS(-)_lab`,
             plot_list$hippuristanol_minus_DMS_A_B, plot_list$hippuristanol_minus_DMS_A_C, plot_list$hippuristanol_minus_DMS_B_C, plot_list$hippuristanol_minus_DMS_transcript,
             plot_list$`Hippuristanol_DMS(+)_lab`,
             plot_list$hippuristanol_plus_DMS_A_B, plot_list$hippuristanol_plus_DMS_A_C, plot_list$hippuristanol_plus_DMS_B_C, plot_list$hippuristanol_plus_DMS_transcript,
             nrow = 8, layout_matrix = rbind(rep(1,4), 2:5, rep(6,4), 7:10, rep(11,4), 12:15, rep(16,4), 17:20), heights = c(1,4,1,4,1,4,1,4))
dev.off()

