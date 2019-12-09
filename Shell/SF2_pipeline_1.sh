#!/bin/bash

#This shell script takes the unprocessed <.fastq> files, which are available in the Gene Expression Omnibus (GEO) database accession GSE134865
#at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134865 and processes them exactly as in Waldron et al. 2020 Genome Biology.
#It generates <.rtsc> and <.react> files for each replicate from each sample. These files can also be downloaded from GSE134865 as supplementary files.
#It also creates logs on the trimming/mapping/filtering of reads and a csv file that can be used to calculate ligation bias.
#It uses the StructureFold2 pipeline, which is available at https://github.com/StructureFold2/StructureFold2
#Ensure all Software Dependencies for StructureFold2 are downloaded and that the StructureFold2 scripts are in the PATH
#This script will take a long time to run and generates some very large intermediary files. It is up to the user whether they keep these intermediary files as they 
#will not be required for any analysis once the <.rtsc> and <.react> files have been generated. The <.rtsc> and <.react> will be much smaller in size.
#Run the script from within a directory containing the following sub-directories;
#fastq_files, trimmed_fastq_files, SAM_files, filtered_SAM_files, ligation_bias, rtsc_files, react_files, logs, FASTA
#Place all <.fastq.gz> files into the directory named fastq_files and the MCF7_2015.fasta file into the FASTA directory

#Unzip all <.fastq.gz> files
cd fastq_files
gunzip *

#Calculate ligation bias
check_ligation_bias.py
mv *.csv ../ligation_bias

#Trim fastq files of 5' and 3' adaptor sequences. Trim bases from the 3’ end with a NextSeq quality score below 30 and remove any reads that are less than 20nt after trimming
fastq_trimmer.py -log -nextseq

#Move trimmed fastq files to a seperate directory
mv *_trimmed.fastq ../trimmed_fastq_files
mv *.txt ../logs
gzip *.fastq
cd ../FASTA

#Build a Bowtie2 database of the MCF7_2015 transcriptome
mkdir bowtie2_database
bowtie2-build MCF7_2015.fasta bowtie2_database/MCF7_2015

#Map trimmed fastq files to MCF7_2015 transcriptome
cd ../trimmed_fastq_files
fastq_mapper.py ../FASTA/bowtie2_database/MCF7_2015 -log

#move mapped sam files to a seperate directory
mv *.sam ../SAM_files
mv *.csv ../logs

cd ../SAM_files
sam_filter.py -max_mismatch 4

#move filtered SAM files to a seperate directory
mv *_filtered.SAM ../filtered_SAM_files
mv *.csv ../logs

#Generate <.rtsc> files
cd ../filtered_SAM_files
sam_to_rtsc ../FASTA/MCF7_2015.fasta -trim _trimmed_mapped_filtered
mv *.rtsc ../rtsc_files
cd ../rtsc_files

#Generate <.react> files
rtsc_to_react.py control_minus_DMS_1.rtsc control_plus_DMS_1.rtsc ../FASTA/MCF7_2015.fasta
rtsc_to_react.py control_minus_DMS_2.rtsc control_plus_DMS_2.rtsc ../FASTA/MCF7_2015.fasta -scale control_minus_DMS_1_control_plus_DMS_1_ln_nrm.scale
rtsc_to_react.py control_minus_DMS_3.rtsc control_plus_DMS_3.rtsc ../FASTA/MCF7_2015.fasta -scale control_minus_DMS_1_control_plus_DMS_1_ln_nrm.scale

rtsc_to_react.py hippuristanol_minus_DMS_1.rtsc hippuristanol_plus_DMS_1.rtsc ../FASTA/MCF7_2015.fasta -scale control_minus_DMS_1_control_plus_DMS_1_ln_nrm.scale
rtsc_to_react.py hippuristanol_minus_DMS_2.rtsc hippuristanol_plus_DMS_2.rtsc ../FASTA/MCF7_2015.fasta -scale control_minus_DMS_1_control_plus_DMS_1_ln_nrm.scale
rtsc_to_react.py hippuristanol_minus_DMS_3.rtsc hippuristanol_plus_DMS_3.rtsc ../FASTA/MCF7_2015.fasta -scale control_minus_DMS_1_control_plus_DMS_1_ln_nrm.scale

mv *.react ../react_files
cd ..




