###This script was written by Joseph A.Waldron and produces panels S1F-G in Waldron et al. (2019) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accessions GSE134865 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865

###Imports
library(tidyverse)

#set home directory----
home <- '' #this needs to be set to the directory containing the data

#write theme----
my_theme <- theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        panel.border = element_blank(),
        panel.grid = element_blank())

#specificity----
#read data
specificity <- read_csv(file = file.path(home, "specificity.csv"), col_names = T)

#reformat data
specificity %>%
  select(base, control_minus_DMS_specificity, control_plus_DMS_specificity, hippuristanol_minus_DMS_specificity, hippuristanol_plus_DMS_specificity) %>%
  gather(key = condition, value = percentage, control_minus_DMS_specificity, control_plus_DMS_specificity, hippuristanol_minus_DMS_specificity, hippuristanol_plus_DMS_specificity) %>%
  mutate(base = factor(base, levels = c("T", "G", "C", "A"), labels = c("U", "G", "C", "A"), ordered = T),
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
  my_theme+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=c(A="#086FA1", C="#FF3500", G="#00B64F", U="#FF8900"))+
  coord_flip()

pdf(file = 'ligation.pdf', height = 4, width = 6)
print(ligation_plot)
dev.off()
         
