#!/bin/bash

#Calculate gene uTIS scores
#download Dataset_S01 from Lee et al. 2012 PNAS at https://www.pnas.org/content/suppl/2012/08/27/1207846109.DCSupplemental and make into a <.csv> named "sd01.csv"
#and run the following script with "sd01.csv" in the current directory
uTIS_analysis_GTI.py

#Carry out G4 RNA screener analysis on 5'UTR sequences
#Ensure G4 RNA screener is first downloaded from http://gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener
screen.py MCF7_2015_fpUTRs.fasta -w 50 > fpUTR_G4_screener.tsv




