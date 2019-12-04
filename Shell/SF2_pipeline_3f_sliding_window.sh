#!/bin/bash

#This script will produce data required for the Sliding_windows_reactivity.R script
#Ensure that the input files are in the current directory before running

#Carry out sliding window analysis
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 15 -wstep 3
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 20 -wstep 3
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 25 -wstep 3
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 30 -wstep 3
react_windows.py control_fpUTR.react hippuristanol_fpUTR.react MCF7_2015_fpUTRs.fasta -wlen 40 -wstep 3

