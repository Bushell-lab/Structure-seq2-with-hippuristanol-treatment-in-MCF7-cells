###This script was written by Joseph A.Waldron and produces panel S7A in Waldron et al. (2019) Genome Biology
###Input data can be downloaded from the Gene Expression Omnibus (GEO) database accessions GSE134888 which can be found at 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888

#Imports
library(tidyverse)
library(VennDiagram)

#set home directory----
home <- '' #this needs to be set to the directory containing the data

#set variables----
#posterior probability threshold
positive_change <- 0.25

#read in translation data----
translation_data <- read_csv(file = file.path(home, "penn-DOD-gene.mmdiff90.csv"), col_names = T)

#pull names of 4A-dep genes
translation_data %>%
  filter(posterior_probability > positive_change & eta1_1 - eta1_2 < 0) %>%
  pull(external_gene_name) -> fourAdep_gene_names

#read in Iwasaki et al. data----
high_sensitivity_hipp <- read_csv(file = file.path(home, "nature17978-s2A.csv"), col_names = T) #this is supplemental table S2A from Iwasaki et al. (2016) Nature

#pull unique gene names from table
high_sensitivity_hipp %>%
  group_by(Gene) %>%
  sample_n(size = 1) %>%
  pull(Gene) -> high_sensitivity_hipp_gene_names

#read in Modelska et al. data----
fourAone_dep <- read_csv(file = file.path(home, "fourAone_dep_transcripts.csv"), col_names = T) #this was created from Supplementary Table 5 from Modelska et al. (2015) Cell Death and Disease

#pull gene names from table
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
