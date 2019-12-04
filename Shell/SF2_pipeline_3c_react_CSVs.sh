#!/bin/bash

#This script will produce the csv files required for the Positional_changes.R and dStruct.R scripts
#Ensure that the input files are in the current directory before running

#create csv files of averaged spliced react files
react_to_csv.py control_fpUTR.react MCF7_2015_fpUTRs.fasta -all_transcripts
react_to_csv.py control_CDS.react MCF7_2015_CDSs.fasta -all_transcripts
react_to_csv.py control_tpUTR.react MCF7_2015_tpUTRs.fasta -all_transcripts
react_to_csv.py hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -all_transcripts
react_to_csv.py hippuristanol_CDS.react MCF7_2015_CDSs.fasta -all_transcripts
react_to_csv.py hippuristanol_tpUTR.react MCF7_2015_tpUTRs.fasta -all_transcripts

#create csv files of replicate react files
react_to_csv.py control_1.react MCF7_2015.fasta -all_transcripts
react_to_csv.py control_2.react MCF7_2015.fasta -all_transcripts
react_to_csv.py control_3.react MCF7_2015.fasta -all_transcripts
react_to_csv.py hippuristanol_1.react MCF7_2015.fasta -all_transcripts
react_to_csv.py hippuristanol_2.react MCF7_2015.fasta -all_transcripts
react_to_csv.py hippuristanol_3.react MCF7_2015.fasta -all_transcripts

