#!/bin/bash

#This script will produce data required for the GGC4.R script
#Ensure that the input files are in the current directory before running

#Calculate the reactivity of all (GGC)4 and (GCC)4 sequences in 5'UTRs
react_static_motif.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta GGCGGCGGCGGC -fp 0 -tp 0 -reactout
react_static_motif.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta GCCGCCGCCGCC -fp 0 -tp 0 -reactout

react_statistics.py -react control_fpUTR_GGCGGCGGCGGC_0fp_0tp.react -n 0 -m 0
react_statistics.py -react hippuristanol_fpUTR_GGCGGCGGCGGC_0fp_0tp.react -n 0 -m 0
react_statistics.py -react control_fpUTR_GCCGCCGCCGCC_0fp_0tp.react -n 0 -m 0
react_statistics.py -react hippuristanol_fpUTR_GCCGCCGCCGCC_0fp_0tp.react -n 0 -m 0

