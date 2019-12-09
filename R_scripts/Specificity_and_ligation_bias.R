###This script was written by Joseph A. Waldron and produces panels S1F-G in Waldron et al. (2020) Genome Biology
###Input data first needs to be generated using the Shell scripts from this repository (see README file)

#Imports----
library(tidyverse)

#write theme----
my_theme <- theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        panel.border = element_blank(),
        panel.grid = element_blank())+

#load data----
#specificity
specificity <- read_csv(file = "specificity.csv", col_names = T) #generated with SF2_pipeline_3a_QC.sh

#ligation bias
ligation_list <- list()
for (condition in c("control_minus_DMS", "control_plus_DMS", "hippuristanol_minus_DMS", "hippuristanol_plus_DMS")) {
  for (replicate in 1:3) {
    df <- read_csv(file = "ligation_bias", paste(condition, replicate, "ligation_bias.csv", sep = "_"), col_names = T) #generated with SF2_pipeline_1.sh
    df %>%
      mutate(condition = rep(condition, nrow(df)),
             replicate = rep(replicate, nrow(df)),
             count = as.numeric(count)) -> df
    ligation_list[[paste(condition, replicate, sep = "_")]] <- df
  }
}
ligation <- do.call("rbind", ligation_list)

#reformat and plot data----
#specificity
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
  my_theme+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=c(A="#086FA1", C="#FF3500", G="#00B64F", U="#FF8900"))+
  coord_flip()

pdf(file = 'specificity.pdf', height = 4, width = 6)
print(specificity_plot)
dev.off()

#ligation bias
#Calculate the total counts
total_counts <- sum(ligation$count[ligation$position == "total_counts" & ligation$nt != "N"])

#for each nt, divide the total counts of that nt by the total counts to get a percentage
ligation %>%
  filter(position == "total_counts" & nt != "N") %>%
  group_by(nt) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  mutate(percentage = count / total_counts,
         condition = rep("transcriptome")) %>%
  select(condition, nt, percentage) -> transcriptome_counts

#calculate the total number counts at position_1 for each condition
ligation %>%
  filter(position == "position_1" & nt != "N") %>%
  group_by(condition) %>%
  summarise(total_counts = sum(count)) -> condition_counts

#for each nt within each condition, divide the total counts at position_1 of that nt by the total counts at psotion_1 in that condition to get a percentage
ligation %>%
  filter(position == "position_1" & nt != "N") %>%
  group_by(condition, nt) %>%
  summarise(count = sum(count)) %>%
  inner_join(condition_counts, by = "condition") %>%
  mutate(percentage = count / total_counts) %>%
  select(condition, nt, percentage) -> sample_counts

#merge data and mutate condition and nt into factors with appropriate labels for the figure
transcriptome_counts %>%
  bind_rows(sample_counts) %>%
  mutate(condition = factor(condition, levels = c("transcriptome", "hippuristanol_plus_DMS", "control_plus_DMS", "hippuristanol_minus_DMS", "control_minus_DMS"),
                            labels = c("transcriptome", "Hipp / DMS (+)", "Ctrl / DMS (+)", "Hipp / DMS (-)", "Ctrl / DMS (-)"), ordered = T),
         nt = factor(nt,
                     levels = c("T", "G", "C", "A"), 
                     labels = c("U", "G", "C", "A"),
                     ordered = T)) -> merged_data

#plot
ligation_plot <- ggplot(data = merged_data, aes(y = percentage, x = condition, fill = nt))+
  geom_bar(stat="identity")+
  my_theme+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=c(A="#086FA1", C="#FF3500", G="#00B64F", U="#FF8900"))+
  coord_flip()

pdf(file = 'ligation.pdf', height = 4, width = 6)
print(ligation_plot)
dev.off()
