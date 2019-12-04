#coverage data----
coverage_data <- read_csv(file = 'plus_DMS_coverage_all_replicates.csv', col_names = T) #generated with SF2_pipeline_2.sh

#fp coverage data----
ctrl_fp_coverage_data <- read_csv(file = 'control_minus_DMS_FP_10.csv', col_names = T) #generated with SF2_pipeline_2.sh
hipp_fp_coverage_data <- read_csv(file = 'hippuristanol_minus_DMS_FP_10.csv', col_names = T) #generated with SF2_pipeline_2.sh

ctrl_fp_coverage_data %>%
  dplyr::rename(control_minus_DMS_fp_10_coverage = FP_coverage) %>%
  dplyr::inner_join(hipp_fp_coverage_data, by = "transcript") %>%
  dplyr::rename(hippuristanol_minus_DMS_fp_10_coverage = FP_coverage) -> fp_coverage_data

rm(ctrl_fp_coverage_data, hipp_fp_coverage_data)

#totals data
totals_data <- read_tsv(file = 'penn-DE.mmdiffMCF7', col_names = T, skip = 1) #download from GSE134888 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888
totals_data %>%
  dplyr::mutate(abundance = case_when(posterior_probability > positive_change ~ alpha1,
                               posterior_probability < positive_change ~ alpha0)) %>%
  dplyr::rename(transcript = feature_id) %>%
  dplyr::select(transcript, abundance) -> abundance_data
rm(totals_data)

#translation data
translation_data <- read_tsv(file = 'penn-DOD-gene.mmdiffMCF7', col_names = T, skip = 1) #download from GSE134888 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888
translation_data %>%
  dplyr::rename(gene = feature_id) %>%
  dplyr::mutate(DOD = eta1_1 - eta1_2,
         translation = factor(case_when(posterior_probability > positive_change & DOD < 0 ~ "4A-dep",
                                        posterior_probability < no_change ~ "4A-indep"), levels = c("4A-dep", "4A-indep"), ordered = T)) -> translation_data

#FASTA composition data
FASTA_compositions_list <- list()
for (region in c("fpUTR", "CDS", "tpUTR")) {
  df <- read_csv(file = 'MCF7_2015_', region, 's_composition.csv', col_names = T) #generate with SF2_pipeline_2.sh
  df$region <- rep(region, nrow(df))
  FASTA_compositions_list[[region]] <- df
}
FASTA_compositions <- do.call("rbind", FASTA_compositions_list)

#transcript to gene ID
transcript_to_geneID <- read_tsv(file = 'MCF7_2015_transcript_to_gene_map.txt', col_names = F) #download from the data folder of this repository
transcript_to_geneID %>%
  dplyr::rename(gene = X1, transcript = X2) -> transcript_to_geneID
