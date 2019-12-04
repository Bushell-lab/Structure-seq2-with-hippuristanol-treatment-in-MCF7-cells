# structure-seq

This repository contains all the Shell R and Python scripts used to create the figures within Waldron et al. (2020) Genome Biology.

Raw sequencing data are available in the Gene Expression Omnibus (GEO) database accessions GSE134865 and GSE134888 which can be found at 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865 and
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888.

# Download and splice FASTA file
First, download the MCF7 specific FASTA file named IsoSeq_MCF7_2015edition_polished.unimapped.fasta from http://datasets.pacb.com.s3.amazonaws.com/2015/IsoSeqHumanMCF7Transcriptome/list.html and then with this file in the current directory, run the reformatting_MCF7_2015_FASTA.py script. This generates a reformatted FASTA named MCF7_2015.fasta, that has a space between the transcript ID and the rest of the information on each header line.

Then splice this new FASTA into a 5'UTR, CDS and 3'UTR FASTA by running splicing_MCF7_2015_FASTA.py. This script uses the following three files, which are in the data folder of this repository;

All_R_V_HumanRefSeqs_Release85_CDSs_lengths.csv
NM_CDSs_to_MCF7_2015_blast_eris_start.tsv
NM_CDSs_to_MCF7_2015_blast_eris_end.tsv

so ensure that they and the FASTA are in the current directory when running the script.

# Processing Data
Then run the SF2_pipeline shell scripts to process the sequencing data and generate the files required by the R scripts to make the figures.
These SF2 scripts call python scripts from the StructureFold2 pipeline, which is available at https://github.com/StructureFold2/StructureFold2
Ensure all Software Dependencies for StructureFold2 are downloaded and that the StructureFold2 scripts are in the PATH.

SF2_pipeline_1.sh starts with the raw unprocessed <.fastq> files and generates <.rtsc> and <.react> files for each replicate from each sample. This script will take a long time to run and generates some very large intermediary files. It is not neccessary to run this script, as the <.rtsc> and <.react> files can also be downloaded from GSE134865 as supplementary files.

SF2_pipeline_2.sh needs to be run next as this generates files that are required for the SF2_pipeline_3.sh scripts

Filtering_transcripts.R then needs to be run to generate a flat text files of filtered transcripts with 5'UTR lengths of more or less than 100nt

TheSF2_pipeline_3.sh scripts can then be run to generate data that is required by the R scripts which produce the figures

# Generate figures
The R scripts can then be used to generate the figures. These scripts read in commonly used data and variables from the Structure_seq_common_data.R and Structure_seq_variables.R scripts respecitively, so ensure that these are in the current working directory when running these scripts.
