###Imports
library(tidyverse)

#set home and working directory----
setwd('N:\\JWALDRON/Structure_seq/Paper/Figures/R/Structure_seq_QC/panels')
home <- 'N:\\JWALDRON/Structure_seq'

#specificity----
#read data
specificity <- read_csv(file = file.path(home, "MCF7_2015/raw_data/specificity/specificity.csv"), col_names = T)

#reformat data
specificity %>%
  select(base, control_minus_DMS_specificity, control_plus_DMS_specificity, hippuristanol_minus_DMS_specificity, hippuristanol_plus_DMS_specificity) %>%
  gather(key = condition, value = percentage, control_minus_DMS_specificity, control_plus_DMS_specificity, hippuristanol_minus_DMS_specificity, hippuristanol_plus_DMS_specificity) %>%
  mutate(base = factor(base,
                       levels = c("T", "G", "C", "A"), 
                       labels = c("U", "G", "C", "A"),
                       ordered = T),
         condition = factor(condition, 
                            levels = c("hippuristanol_plus_DMS_specificity", "control_plus_DMS_specificity", "hippuristanol_minus_DMS_specificity", "control_minus_DMS_specificity"), 
                            labels = c("Hipp / DMS (+)", "Ctrl / DMS (+)", "Hipp / DMS (-)", "Ctrl / DMS (-)"),
                            ordered = T)) -> specificity

#plot   
specificity_plot <- ggplot(data = specificity, aes(y = percentage, x = condition, fill = base))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        panel.border = element_blank(),
        panel.grid = element_blank())+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=c(A="#086FA1", C="#FF3500", G="#00B64F", U="#FF8900"))+
  coord_flip()

pdf(file = 'specificity.pdf', height = 4, width = 6)
print(specificity_plot)
dev.off()

#ligation bias----
#read data
ligation_list <- list()
for (condition in c("control_minus_DMS", "control_plus_DMS", "hippuristanol_minus_DMS", "hippuristanol_plus_DMS")) {
  df <- read_csv(file = file.path(home, "ligation_bias", paste0(condition, "_ligation_bias.csv")), col_names = T)
  df %>%
    mutate(nt = factor(nt,
                       levels = c("T", "G", "C", "A"), 
                       labels = c("U", "G", "C", "A"),
                       ordered = T),
           condition = factor(condition, 
                              levels = c("hippuristanol_plus_DMS", "control_plus_DMS", "hippuristanol_minus_DMS", "control_minus_DMS"), 
                              labels = c("Hipp / DMS (+)", "Ctrl / DMS (+)", "Hipp / DMS (-)", "Ctrl / DMS (-)"),
                              ordered = T),
           position = factor(position),
           count = as.numeric(count)) -> df
  ligation_list[[condition]] <- df
}

ligation <- do.call("rbind", ligation_list)

#reformat data
ligation %>%
  filter(position == "all") %>%
  group_by(nt) %>%
  summarise(count = sum(count)) %>%
  mutate(condition = factor(rep("transcriptome"))) -> transcriptome_counts

transcriptome_counts %>%
  group_by(condition) %>%
  summarise(total_counts = sum(count)) %>%
  inner_join(transcriptome_counts, by = "condition") %>%
  mutate(percentage = count / total_counts) -> transcriptome_counts

ligation %>%
  filter(position == "ligated") %>%
  group_by(condition, nt) %>%
  summarise(count = sum(count)) -> sample_counts

sample_counts %>%
  group_by(condition) %>%
  summarise(total_counts = sum(count)) %>%
  inner_join(sample_counts, by = "condition") %>%
  mutate(percentage = count / total_counts) -> sample_counts

transcriptome_counts %>%
  bind_rows(sample_counts) %>%
  select(condition, nt, percentage) %>%
  mutate(condition = factor(condition, 
                            levels = c("transcriptome", "Hipp / DMS (+)", "Ctrl / DMS (+)", "Hipp / DMS (-)", "Ctrl / DMS (-)"),
                            ordered = T)) -> merged_data

#plot
ligation_plot <- ggplot(data = merged_data, aes(y = percentage, x = condition, fill = nt))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        panel.border = element_blank(),
        panel.grid = element_blank())+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=c(A="#086FA1", C="#FF3500", G="#00B64F", U="#FF8900"))+
  coord_flip()

pdf(file = 'ligation.pdf', height = 4, width = 6)
print(ligation_plot)
dev.off()
         