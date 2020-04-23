#!/bin/bash

#This script will produce the csv files required for the R10_analysis_2.R script
#Ensure that the input files are in the current directory before running

#use react_composition.py script to generate fasta and react files for downstream 65nt of AG motifs in 5'UTRs
react_composition.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta AG -size 10 -perc 1.0 -fp 0 -tp 65 -unique -fastaout -reactout

#use react_windows.py script to generate fasta and react files for 50nt windows with 5nt steps
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 50 -wstep 5 -fastaout -reactout

#trim hippuristanol 5'UTR AG motif react file and 5'UTR AG motif fasta file in order to fold the 31:50nt window
Trim_reacts.py hippuristanol_fpUTR_10_AG_1.0_unique_0fp_65tp.react 25
Trim_FASTA.py control_fpUTR_hippuristanol_fpUTR_10_AG_1.0_unique_0fp_65tp.fasta 25

#rename trimmed files
mv hippuristanol_fpUTR_10_AG_1.0_unique_0fp_65tp_trim_25_FP.react hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp.react
mv control_fpUTR_hippuristanol_fpUTR_10_AG_1.0_unique_0fp_65tp_trim_25_FP.fasta fpUTR_10_AG_1.0_unique_16-65tp.fasta

#rename sliding windows fasta
mv control_fpUTR_hippuristanol_fpUTR_50win_5step.fasta fpUTR_50win_5step.fasta

#using the batch_fold_rna.py script and the IDs outputed from R10_analysis_1.R, fold windows
batch_fold_rna.py filtered_10_AG_1.0_IDs.txt fpUTR_10_AG_1.0_unique_16-65tp.fasta 1 -r hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp.react -minl 50
batch_fold_rna.py filtered_random_IDs.txt fpUTR_50win_5step.fasta 1 -r hippuristanol_fpUTR_50win_5step.react -minl 50

#rename output folders
mv output_files_hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp.react_filtered_10_AG_1.0_IDs.txt_310.15_fpUTR_10_AG_1.0_unique_16-65tp.fasta_RNAstructure-mfe_sht_0_md_99999 \
output_files_hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp_filtered_IDs_310.15_RNAstructure-mfe_sht_0_md_99999
mv output_files_hippuristanol_fpUTR_50win_5step.react_filtered_random_IDs.txt_310.15_fpUTR_50win_5step.fasta_RNAstructure-mfe_sht_0_md_99999 \
output_files_hippuristanol_fpUTR_50win_5step_filtered_random_IDs_310.15_RNAstructure-mfe_sht_0_md_99999

#calculate structure statistics
structure_statistics.py -d output_files_hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp_filtered_IDs_310.15_RNAstructure-mfe_sht_0_md_99999/CT -mode F -offset 5 -na 0 -name hippuristanol_fpUTR_10_AG_1.0_unique_16-65tp_filtered_IDs_structure_statistics.csv
structure_statistics.py -d output_files_hippuristanol_fpUTR_50win_5step_filtered_random_IDs_310.15_RNAstructure-mfe_sht_0_md_99999/CT -mode F -offset 1 -na 0 -name hippuristanol_fpUTR_50win_5step_filtered_random_IDs_structure_statistics.csv

#generate a fasta and react files to fold all 50nt windows with the 5'UTRs of all filtered transcripts with a 5'UTR more than 100nt
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 50 -wstep 10 -restrict filtered_plus_100_transcripts.txt -fastaout -reactout

#rename fasta file
mv control_fpUTR_hippuristanol_fpUTR_50win_10step.fasta fpUTR_50win_10step.fasta

#extract all transcript_step IDs from FASTA
return_IDs.py fpUTR_50win_10step.fasta

#fold all windows
batch_fold_rna.py fpUTR_50win_10step_IDs.txt fpUTR_50win_10step.fasta 1 -r hippuristanol_fpUTR_50win_10step.react

#rename folders
mv output_files_hippuristanol_fpUTR_50win_10step.react_fpUTR_50win_10step_IDs.txt_310.15_fpUTR_50win_10step.fasta_RNAstructure-mfe_sht_0_md_99999 \
output_files_hippuristanol_fpUTR_50win_10step_filtered_plus_100_transcripts_310.15_RNAstructure-mfe_sht_0_md_99999

#extract MFEs and strandedness
structure_statistics.py -d output_files_hippuristanol_fpUTR_50win_10step_filtered_plus_100_transcripts_310.15_RNAstructure-mfe_sht_0_md_99999/CT -mode F -offset 1 -na 0 -name hippuristanol_fpUTR_50win_10step_structure_statistics.csv

echo All done