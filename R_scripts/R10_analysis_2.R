###This script was written by Joseph A. Waldron and in combinatin with R10_analysis_1.R produces Figures 2H and Supplementary Figures 2H-J in Schmidt et al. NAR 2023
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

pdf("Supp_Fig_2J_MFE_windows.pdf", width = 4, height = 4)
print(windows_MFE_violin)
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

pdf(file = "Supp_Fig_2I_fpUTR_positions.pdf", width = 4, height = 4)
print(R10_position_hist)
dev.off()
