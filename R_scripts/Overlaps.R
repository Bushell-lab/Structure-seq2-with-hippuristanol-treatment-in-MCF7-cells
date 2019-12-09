###This script was written by Joseph A.Waldron and produces panel S7A in Waldron et al. (2019) Genome Biology
###Note that some of the data required for this script is generated following the mapping of the polysome sequencing reads to the Ensmebl version 90 transcriptome
###To avoid confusion and difficulties in naming the files appropriately, it was decided not to host this data within GSE134888
###Please contact us if you require this data as we would be happy to share. We also have a version of this script which uses the data available within GSE134888
###which prodcues a very similar Venn diagram, which we would also be happy to share on request.

#Imports----
library(tidyverse)
library(VennDiagram)

#import variables----
source("Structure_seq_variables.R")

#read in data----
#translation data
#contact us if you require this data
translation_data <- read_csv(file = "penn-DOD-gene.mmdiff90.csv", col_names = T)

#Iwasaki et al. data
#this is supplemental table S2A from Iwasaki et al. (2016) Nature
high_sensitivity_hipp <- read_csv(file = "nature17978-s2A.csv", col_names = T)

#Modelska et al. data
#this was created from Supplementary Table 5 from Modelska et al. (2015) Cell Death and Disease
fourAone_dep <- read_csv(file = "fourAone_dep_transcripts.csv", col_names = T)

#pull names of 4A-dep genes
translation_data %>%
  filter(posterior_probability > positive_change & eta1_1 - eta1_2 < 0) %>%
  pull(external_gene_name) -> fourAdep_gene_names

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
