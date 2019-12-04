###This script was written by Joseph A.Waldron and produces panel S7A in Waldron et al. (2019) Genome Biology

#Imports----
library(tidyverse)
library(VennDiagram)

#import variables----
source("Structure_seq_variables.R")

#read in data----
#translation data
#download from GSE134888 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888
translation_data <- read_tsv(file = 'GSE134888_penn-DOD-gene.mmdiffMCF7.tsv', col_names = T, skip = 1)

#Iwasaki et al. data
#this is supplemental table S2A from Iwasaki et al. (2016) Nature
high_sensitivity_hipp <- read_csv(file = "nature17978-s2A.csv", col_names = T)

#Modelska et al. data
#this was created from Supplementary Table 5 from Modelska et al. (2015) Cell Death and Disease
fourAone_dep <- read_csv(file = "fourAone_dep_transcripts.csv", col_names = T)

#MCF7 gene names
#can be downloaded from the data folder of this repository
MCF7_IDs <- read_csv(file = "MCF7_2015_ensembl_IDs.csv", col_names = T)

#pull names of 4A-dep genes
translation_data %>%
  inner_join(transcript_to_geneID, by = c("feature_id" = "gene")) %>%
  inner_join(MCF7_IDs, by = "transcript") %>%
  filter(posterior_probability > positive_change & eta1_1 - eta1_2 < 0) %>%
  group_by(Ensembl_gene_symbol) %>%
  sample_n(size = 1) %>%
  pull(Ensembl_gene_symbol) -> fourAdep_gene_names

#pull unique gene names from Iwasaki table
high_sensitivity_hipp %>%
  group_by(Gene) %>%
  sample_n(size = 1) %>%
  pull(Gene) -> high_sensitivity_hipp_gene_names

#pull gene names from Modelska table
fourAone_dep %>%
  pull(Gene_name) -> fourAone_dep_gene_names

#plot venn----
#calculate all possible overlaps
overlap_123_IDs <- intersect(intersect(fourAdep_gene_names, high_sensitivity_hipp_gene_names), fourAone_dep_gene_names)
overlap_12_IDs <- intersect(fourAdep_gene_names, high_sensitivity_hipp_gene_names)
overlap_13_IDs <- intersect(fourAdep_gene_names, fourAone_dep_gene_names)
overlap_23_IDs <- intersect(high_sensitivity_hipp_gene_names, fourAone_dep_gene_names)

#plot Venn diagram
venn_plot <- draw.triple.venn(area1 = length(fourAdep_gene_names),
                 area2 = length(high_sensitivity_hipp_gene_names),
                 area3 = length(fourAone_dep_gene_names),
                 n123 = length(overlap_123_IDs),
                 n12 = length(overlap_12_IDs),
                 n13 = length(overlap_13_IDs),
                 n23 = length(overlap_23_IDs),
                 category = c("4A-dep", "Iwasaki et al.", "Modelska el at."),
                 lty = rep("blank",3),
                 fill = c("light blue", "red", "yellow"),
                 alpha = rep(0.5, 3), 
                 euler.d = T,
                 scaled = T,
                 overrideTriple = T)

pdf(file = "overlaps.pdf", width = 5, height = 4)
grid.draw(venn_plot)
dev.off()
