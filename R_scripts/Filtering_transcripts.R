#Imports
library(tidyverse)

#import variables----
source("Structure_seq_variables.R")

#load common data----
source("Structure_seq_common_data.R")

#filter----
#the following pipe will filter by the coverage and 5' end coverage thresholds and then pick the most abundant transcript per gene
FASTA_compositions_list$fpUTR %>%
  inner_join(coverage_data, by = "transcript") %>%
  inner_join(fp_coverage_data, by = "transcript") %>%
  inner_join(abundance_data, by = "transcript") %>%
  inner_join(transcript_to_geneID, by = "transcript") %>%
  filter(control_plus_DMS_1_coverage > coverage,
         control_plus_DMS_2_coverage > coverage,
         control_plus_DMS_3_coverage > coverage,
         hippuristanol_plus_DMS_1_coverage > coverage,
         hippuristanol_plus_DMS_2_coverage > coverage,
         hippuristanol_plus_DMS_3_coverage > coverage,
         control_minus_DMS_fp_10_coverage > fp_coverage,
         hippuristanol_minus_DMS_fp_10_coverage > fp_coverage) %>%
  group_by(gene) %>%
  top_n(n = 1, wt = abundance) %>% # selects most abundant transcript
  ungroup() -> filtered_transcripts

#the following two pipes splits the data into those transcripts with a 5'UTR length more than 100nt and less than or equal to 100nt
filtered_transcripts %>%
  filter(length > 100) %>%
  pull(transcript) -> filtered_100_plus_transcripts

filtered_transcripts %>%
  filter(length <= 100) %>%
  pull(transcript) -> filtered_minus_100_transcripts

#write out flat list of transcriptIDs----
write.table(filtered_100_plus_transcripts, file = "filtered_plus_100_transcripts.txt", row.names = F, col.names = F, quote = F)
write.table(filtered_minus_100_transcripts, file = "filtered_minus_100_transcripts.txt", row.names = F, col.names = F, quote = F)
