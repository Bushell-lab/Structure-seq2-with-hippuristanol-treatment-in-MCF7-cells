###This script was written by Joseph A. Waldron and in combinatin with R10_analysis_1.R produces Figures 4e-g and Extended data Figures 6e-j in Schmidt et al. (in prep)
###Input data first needs to be generated using the R10_analysis_2.sh shell script from this repository and the variables from R10_analysis_1.R should be kept in the environment (see README file)

#windows MFE----
#load data
R10_structure_statistics <- read_csv(file = "hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp_filtered_IDs_structure_statistics.csv", col_names = T)
random_structure_statistics <- read_csv(file = "hippuristanol_fpUTR_50win_5step_filtered_random_IDs_structure_statistics.csv", col_names = T)

#merge data
R10_structure_statistics %>%
  mutate(step = str_replace(transcript, ".+_", ""),
         transcript = str_replace(transcript, "_.+", ""),
         motif = rep("AG")) -> R10_structure_statistics

random_structure_statistics %>%
  mutate(step = str_replace(transcript, ".+_", ""),
         transcript = str_replace(transcript, "_.+", ""),
         motif = rep("random")) -> random_structure_statistics

window_structure_statistics <- bind_rows(R10_structure_statistics, random_structure_statistics)

#plot
t <- wilcox.test(data = window_structure_statistics, CT_DeltaG ~ motif,
                 paired = F, alternative = "two.sided", var.equal = F, conf.int = T)
p_label <- myP(t)

window_structure_statistics %>%
  ggplot(aes(x = motif, y = CT_DeltaG, fill = motif))+
  geom_violin(alpha = 0.5, fill = "#00BFC4")+
  geom_boxplot(width = 0.2, outlier.shape=NA, fill = "#00BFC4")+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol)")))+
  ylim(c(-40,0))+
  violin_theme+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4)+
  ggtitle(p_label) -> windows_MFE_violin

#whole 5'UTR MFE----
#load structure statistics data for 50nt sliding windows with a step size of 10 for all filtered 5'UTRs with a length more than 100nt
fpUTR_structure_statistics <- read_csv(file = "hippuristanol_fpUTR_50win_10step_structure_statistics.csv", col_names = T)

#the following pipe calculates the mean MFE per transcript and creates a factor based on whether the transcript has an R10 motif in its 5'UTR
fpUTR_structure_statistics %>%
  mutate(step = as.numeric(str_replace(transcript, ".+_", "")),
         transcript = str_replace(transcript, "_.+", "")) %>%
  group_by(transcript) %>%
  summarise(avg_MFE = mean(CT_DeltaG)) %>%
  ungroup() %>%
  mutate(AG10 = factor(case_when(transcript %in% R10_IDs ~ "R10\n5'UTRs",
                                 !(transcript %in% R10_IDs) ~ "non-R10\n5'UTRs"),
                       levels = c("R10\n5'UTRs", "non-R10\n5'UTRs"), ordered = T)) -> fpUTR_MFE_data

#make matched size group
R10_n <- nrow(fpUTR_MFE_data[fpUTR_MFE_data$AG10 == "R10\n5'UTRs",])

set.seed(020588) #sets seed for random selection of size matched group

fpUTR_MFE_data %>%
  group_by(AG10) %>%
  sample_n(size = R10_n) %>%
  ungroup() -> fpUTR_MFE_size_matched

summary(fpUTR_MFE_size_matched)

#plot
t <- wilcox.test(data = fpUTR_MFE_size_matched, avg_MFE ~ AG10,
                 paired = F, alternative = "two.sided", var.equal = F, conf.int = T)
p_label <- myP(t)

#plot
fpUTR_MFE_size_matched %>%
  ggplot(aes(x = AG10, y = avg_MFE, fill = AG10))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  ylab(expression(paste("MFE (", Delta, "G kcal/mol)")))+
  ylim(c(-40,0))+
  violin_theme+
  stat_summary(fun.y=mean, geom='point', shape=16, size=4)+
  ggtitle(p_label) -> whole_fpUTR_avg_MFE_violin

pdf(file = "Fig_4f_GC_MFE_violins.pdf", width = 14, height = 4)
grid.arrange(windows_GC_content_violin, windows_MFE_violin, whole_fpUTR_GC_content_violin, whole_fpUTR_avg_MFE_violin, nrow = 1)
dev.off()

#fpUTR lengths----
#load all di_nt motifs in 5'UTR
for (di_nt in c("AG", "AC", "AT", "GC", "GT", "CT")) {
  df <- read_csv(file = paste0("control_fpUTR_hippuristanol_fpUTR_10_", di_nt, "_1.0_unique_50fp_50tp.csv"), col_names = T, col_types = col_types)
  df$di_nt <- rep(di_nt, nrow(df))
  data_list[[di_nt]] <- df
}

#filter
#the following pipe filters by 5' coverage
#it then removes any motifs that have less than 3 or either of the two nts composing the motif
#it then selects the most abundant transcript
fpUTR_lengths_list <- list()
for (di_nt in c("AG", "AC", "GT", "GC", "AT", "CT")) {
  data_list[[di_nt]] %>%
    inner_join(transcript_to_geneID, by = "transcript") %>%
    inner_join(abundance_data, by = "transcript") %>%
    inner_join(fp_coverage_data, by = "transcript") %>%
    mutate(motif = str_c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10),
           X = substr(di_nt, 1, 1), #identifies the first nt in the di_nt variable
           Xn = str_count(motif, X)) %>% #calculates the number of times that nt is present in the motif
    filter(control_minus_DMS_fp_10_coverage > fp_coverage,
           hippuristanol_minus_DMS_fp_10_coverage > fp_coverage,
           Xn > 2 & Xn < 8) %>% #this line filters any motifs that have less than 3 or more than 7 of the identified nt (ensuring at least 3 of each nt)
    group_by(gene) %>%
    top_n(n = 1, wt = abundance) %>%
    ungroup() %>%
    inner_join(FASTA_compositions_list$fpUTR, by = "transcript") %>%
    select(transcript, di_nt, length) -> df
  
  fpUTR_lengths_list[[di_nt]] <- df
}
fpUTR_lengths <- do.call("rbind", fpUTR_lengths_list)

#combine AC/GT/CT
fpUTR_lengths %>%
  mutate(di_nt_group = factor(case_when(di_nt == "AG" ~ "AG",
                                        di_nt == "GC" ~ "GC",
                                        di_nt == "AT" ~ "AU",
                                        di_nt == "AC" | di_nt == "GT" | di_nt == "CT" ~ "AC/GU/CU"),
                              levels = c("GC", "AG", "AC/GU/CU", "AU"), ordered = T)) %>%
  group_by(di_nt_group, transcript) %>%
  sample_n(size = 1) %>% #removes multiple transcripts per group
  ungroup() -> fpUTR_lengths_combined

#remove any transcripts in more than one group
fpUTR_lengths_combined %>%
  group_by(transcript) %>%
  count() %>%
  filter(n > 1) %>%
  pull(transcript) -> non_unique_IDs

fpUTR_lengths_combined %>%
  filter(!(transcript %in% non_unique_IDs)) -> filtered_fpUTR_lengths

#print summary
summary(filtered_fpUTR_lengths)

wilcox.test(filtered_fpUTR_lengths$length[filtered_fpUTR_lengths$di_nt_group == "AG"],
            filtered_fpUTR_lengths$length[filtered_fpUTR_lengths$di_nt_group == "AC/GU/CU"],
            paired = F, alternative = "two.sided", var.equal = F, conf.int = T)

#plot
density_plot <- ggplot(data = filtered_fpUTR_lengths, aes(x = length, colour = di_nt_group, fill = di_nt_group))+
  geom_density(size = 1, alpha = 0.2)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))+
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))+
  scale_alpha_continuous(range = c(0.1,1), guide = F)+
  xlab("5\'UTR length (nt)")+
  scale_x_log10(breaks=c(10,100,1000,10000),limits=c(10, 10000))+
  density_theme

pdf(file = "Extended_data_Fig__6h_fpUTR_length_density.pdf", width=6,height=4)
print(density_plot)
dev.off()

#translation analysis----
#the following pipe takes the translation data, selects only those transcripts that have an annotated 5'UTR
#then selects the most abundant transcript
#then creates a factor depending on whether there is an R10 motif in its 5'UTR
#then calculates the mean translation efficiency (TE) of each transcript by subtracting the
#log FPKM value for subs from the log FPKM value for polys for each replicate and taking the mean
translation_data %>%
  inner_join(transcript_to_geneID, by = "gene") %>%
  filter(transcript %in% FASTA_compositions_list$fpUTR$transcript) %>% #ensures only those transcripts whose 5'UTR has been annotated are included
  inner_join(abundance_data, by = "transcript") %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  mutate(AG10 = factor(case_when(transcript %in% R10_IDs ~ "R10\n5'UTRs",
                                 !(transcript %in% R10_IDs) ~ "non-R10\n5'UTRs"),
                       levels = c("R10\n5'UTRs", "non-R10\n5'UTRs"), ordered = T),
         TE_1 = mu_PD1_MCF7.gene - mu_SD1_MCF7.gene,
         TE_2 = mu_PD2_MCF7.gene - mu_SD2_MCF7.gene,
         TE_3 = mu_PD3_MCF7.gene - mu_SD3_MCF7.gene,
         mean_TE = rowMeans(cbind(TE_1, TE_2, TE_3))) %>%
  select(transcript, AG10, mean_TE) -> AG10_translation_data

#make a size matched group of non-R10 mRNAs
R10_n <- nrow(AG10_translation_data[AG10_translation_data$AG10 == "R10\n5'UTRs",])

set.seed(020588) #sets seed for random selection of size matched group

AG10_translation_data %>%
  group_by(AG10) %>%
  sample_n(size = R10_n) %>%
  ungroup() -> AG10_translation_size_matched

summary(AG10_translation_size_matched)

#plot
t <- wilcox.test(data = AG10_translation_size_matched, mean_TE ~ AG10,
                 paired = F, alternative = "two.sided", var.equal = F, conf.int = T)
p_label <- myP(t)

AG10_median <- median(AG10_translation_size_matched$mean_TE[AG10_translation_size_matched$AG10 == "R10\n5'UTRs"])
non_AG10_median <- median(AG10_translation_size_matched$mean_TE[AG10_translation_size_matched$AG10 == "non-R10\n5'UTRs"])

density_plot <- ggplot(data = AG10_translation_size_matched, aes(x = mean_TE, colour = AG10, fill = AG10))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab('TE')+
  xlim(c(-2.5,2.5))+
  geom_vline(xintercept = AG10_median, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = non_AG10_median, colour = "#fdae61", linetype="dashed", size = 2)+
  density_theme+
  ggtitle(p_label)

pdf(file = "Extended_data_Fig_6i_TE_density.pdf", width=6,height=4)
print(density_plot)
dev.off()

#AG10_position----
#the following pipe takes all 5'UTR AG motif positions, filters to include only those motifs with more than 2 adenines and guanines
#and only those 5'UTRs more than 100 nt
#it also filters by 5' end coverage and picks the most abundant transcript per gene
#it calculates the normalised position by dividing the mid-point of the 5'UTR by the 5'UTR length
data_list[["fpUTR"]] %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  inner_join(FASTA_compositions_list$fpUTR, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  mutate(motif = str_c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)) %>%
  filter(control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage,
         str_count(motif, "A") > 2 & str_count(motif, "G") > 2,
         length > 100) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>%
  ungroup() %>%
  select(transcript, start, length) %>%
  mutate(position = start + (size/2),
         normalised_position = position / length) -> df

med <- median(df$normalised_position)

print(paste("motif n =", nrow(df)))
print(paste("transcript n =", n_distinct(df$transcript)))

df %>%
  ggplot(aes(x = normalised_position))+
  geom_histogram(bins = 20, alpha = 0.5)+
  geom_vline(xintercept = med, linetype="dashed", size = 1)+
  density_theme+
  xlab("binned 5\'UTR length") -> R10_position_hist

pdf(file = "Extended_data_Fig_6j_fpUTR_positions.pdf", width = 4, height = 4)
print(R10_position_hist)
dev.off()
