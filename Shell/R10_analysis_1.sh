#!/bin/bash

#This script will produce the csv files required for the R10_analysis_1.R script
#Ensure that the input files are in the current directory before running

#use react_composition.py to generate csv files for AG motifs in 5'UTR CDS and 3'UTR
regions='fpUTR CDS tpUTR'
for region in $regions
do
react_composition.py 'control_'$region'.react' 'hippuristanol_'$region'.react' 'MCF7_2015_'$region's.fasta' AG -size 10 -perc 1.0 -fp 50 -tp 50 -unique
done

#use react_composition.py to generate csv files for all other di-nt motifs in 5'UTR
dints='AC AT GC GT CT'
for dint in $dints
do
react_composition.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta $dint -size 10 -perc 1.0 -fp 50 -tp 50 -unique
done

#use react_windows.py script to generate csv and fasta files of 20nt sliding windows for each region with 10nt steps
#rename fasta files and use fasta_composition.py to calculate GC content
regions='fpUTR CDS tpUTR'
for region in $regions
do
react_windows.py 'control_'$region'.react' 'hippuristanol_'$region'.react' 'MCF7_2015_'$region's.fasta' -wlen 20 -wstep 10 -fastaout
mv 'control_'$region'_hippuristanol_'$region'_20win_10step.fasta' $region'_20win_10step.fasta'
fasta_composition.py -single $region'_20win_10step.fasta'
done

echo All done

