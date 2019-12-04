#!/bin/bash

#This script needs to be run before any SF2_pipeline_3 script as these will depend on the output of this script
#Ensure all input files are in the current directory before running this script

#Calculate coverage
rtsc_coverage.py MCF7_2015.fasta -rtsc control_plus_DMS_1.rtsc control_plus_DMS_2.rtsc control_plus_DMS_3.rtsc \
hippuristanol_plus_DMS_1.rtsc hippuristanol_plus_DMS_2.rtsc hippuristanol_plus_DMS_3.rtsc -bases AC -name plus_DMS_coverage_all_replicates.csv

#Combine rtsc files
rtsc_combine.py control_minus_DMS_1.rtsc control_minus_DMS_2.rtsc control_minus_DMS_3.rtsc -name control_minus_DMS.rtsc
rtsc_combine.py control_plus_DMS_1.rtsc control_plus_DMS_2.rtsc control_plus_DMS_3.rtsc -name control_plus_DMS.rtsc
rtsc_combine.py hippuristanol_minus_DMS_1.rtsc hippuristanol_minus_DMS_2.rtsc hippuristanol_minus_DMS_3.rtsc -name hippuristanol_minus_DMS.rtsc
rtsc_combine.py hippuristanol_plus_DMS_1.rtsc hippuristanol_plus_DMS_2.rtsc hippuristanol_plus_DMS_3.rtsc -name hippuristanol_plus_DMS.rtsc

#Calculate 5'end coverage
rtsc_end_coverage.py control_minus_DMS.rtsc FP -length 10
rtsc_end_coverage.py hippuristanol_minus_DMS.rtsc FP -length 10

#Average reactivity data
react_average.py control_1.react control_2.react control_3.react -name control.react
react_average.py hippuristanol_1.react hippuristanol_2.react hippuristanol_3.react -name hippuristanol.react

#Splice reactivity files into 5'UTR, CDS, 3'UTR regions
splice_reacts_by_FASTA.py control.react MCF7_2015_fpUTRs.fasta MCF7_2015_CDSs.fasta MCF7_2015_tpUTRs.fasta
splice_reacts_by_FASTA.py hippuristanol.react MCF7_2015_fpUTRs.fasta MCF7_2015_CDSs.fasta MCF7_2015_tpUTRs.fasta

#Calculate composition data about each region
fasta_composition.py -single MCF7_2015_fpUTRs.fasta
fasta_composition.py -single MCF7_2015_CDSs.fasta
fasta_composition.py -single MCF7_2015_tpUTRs.fasta