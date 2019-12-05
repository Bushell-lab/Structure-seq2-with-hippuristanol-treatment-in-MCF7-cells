# structure-seq

This repository contains all the Shell, R and Python scripts used to create the figures within Waldron et al. (2020) Genome Biology.

Raw sequencing data are available in the Gene Expression Omnibus (GEO) database accessions:
* [GSE134865](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865)
* [GSE134888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134888)

## Download and splice FASTA file
First, download the MCF7 specific FASTA file named IsoSeq_MCF7_2015edition_polished.unimapped.fasta 
[link](http://datasets.pacb.com.s3.amazonaws.com/2015/IsoSeqHumanMCF7Transcriptome/list.html). 
With this file in the current directory, first run the reformatting_MCF7_2015_FASTA.py script. 
This generates a reformatted FASTA named MCF7_2015.fasta, that has a space between the transcript ID 
and the rest of the information on each header line.

Then splice this new FASTA into a 5'UTR, CDS and 3'UTR FASTA, by running splicing_MCF7_2015_FASTA.py. 
It uses the output from a blast in which the MCF7_2015 sequences are the subjet and the NM_ RefSeq CDS 
sequences from release 85 are the query. The 5'UTR length is determined using the output of the blast 
filtered to the smallest value for the query start per subject transcript. Only when the query start equals 
position 1, i.e. it aligns with the start of the RefSeq CDS, is the 5'UTR length determined by this blast hit. 
The 3'UTR length is determined using the output of the blast filtered to the highest value for the query end per 
subject transcript. Only when the query end equals the length of the RefSeq CDS, i.e. it aligns with the end of the 
RefSeq CDS, is the 3'UTR length determined by this blast hit. Only transcripts for which both the 5'UTR and 
3'UTR lengths are determined, and when the resulting CDS is equally divisible by 3, are spliced.

This script uses the following three files, which are in the data folder of this repository;

All_R_V_HumanRefSeqs_Release85_CDSs_lengths.csv
NM_CDSs_to_MCF7_2015_blast_eris_start.tsv
NM_CDSs_to_MCF7_2015_blast_eris_end.tsv

so ensure that they and the FASTA file are in the current directory when running the script.

## Processing Data
Then run the SF2_pipeline shell scripts to process the sequencing data and generate the files required by the R scripts to make the figures.
These SF2 scripts call python scripts from the [StructureFold2](https://github.com/StructureFold2/StructureFold2) pipeline.
Ensure all Software Dependencies for StructureFold2 are downloaded and that the StructureFold2 scripts are in the PATH.

SF2_pipeline_1.sh starts with the raw unprocessed <.fastq> files and generates <.rtsc> and <.react> files for each replicate from each sample. 
This script will take a long time to run and generates some very large intermediary files. It is not neccessary to run this script, as the <.rtsc> 
and <.react> files can also be downloaded from GSE134865 as supplementary files.

SF2_pipeline_2.sh needs to be run next as this generates files that are required for the SF2_pipeline_3.sh scripts.

Filtering_transcripts.R then needs to be run to generate a flat text files of filtered transcripts with 5'UTR lengths of more or less than 100nt.

The SF2_pipeline_3.sh and Custom_scripts.sh scripts can then be run to generate the data required by the R scripts.

## Generate figures
The R scripts can then be used to generate the figures. These scripts read in commonly used data and variables from the 
Structure_seq_common_data.R and Structure_seq_variables.R scripts respecitively, so ensure that these scripts and all input data 
are in the current working directory when running these scripts.
