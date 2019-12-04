#!/bin/bash

#This script will produce the data required for the replicate_correlation.R script and the Specificity_and_ligation_bias.R
#The check_ligation_bias.py script will also need to be run on all unzipped <.fastq> files
#Ensure that the input files are in the current directory before running

#Assess replicate correlation
rtsc_correlation.py control_minus_DMS_1.rtsc control_minus_DMS_2.rtsc control_minus_DMS_3.rtsc \
control_plus_DMS_1.rtsc control_plus_DMS_2.rtsc control_plus_DMS_3.rtsc hippuristanol_minus_DMS_1.rtsc \
hippuristanol_minus_DMS_2.rtsc hippuristanol_minus_DMS_3.rtsc hippuristanol_plus_DMS_1.rtsc \
hippuristanol_plus_DMS_2.rtsc hippuristanol_plus_DMS_3.rtsc -name all_rtsc.csv

#Calculate specificity
rtsc_specificity.py -index MCF7_2015.fasta -rtsc control_minus_DMS.rtsc control_plus_DMS.rtsc hippuristanol_minus_DMS.rtsc hippuristanol_plus_DMS.rtsc -name specificity.csv