###Imports
library(tidyverse)

#set home and working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Polysomes/panels')
home <- 'N:\\JWALDRON/Polysomes/library_prep/gradients/csv_files'

for (n in 1:3) {
  df <- read_csv(file = file.path(home, paste0("Polysomes_with_Hip_", n, ".csv")), col_names = T)
  
  ylim <- max(c(df$Ctrl, df$Hipp))
  
  df %>%
    gather(key = condition, value = absorbance, Ctrl, Hipp) %>%
    ggplot(aes(x = X1, y = absorbance, colour = factor(condition)))+
    geom_line(size = 1.5)+
    scale_x_continuous(limits = c(0.3, 10), breaks = 1:10)+
    ylim(c(0, ylim))+
    ylab("absorbance")+
    xlab("fraction")+
    coord_trans(limx = c(0.3, 10), limy = c(0, ylim))+
    theme_bw()+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 22),
          legend.title = element_blank()) -> trace
  
  pdf(file = paste0("polysome_trace_", n, ".pdf"), height = 4, width = 8)
  print(trace)
  dev.off()
}

