#!/bin/bash

#This script will produce data required for the Stats_plots.R script
#Ensure that the input files are in the current directory before running

#Calculate average reactivity and Gini coefficients
react_statistics.py -react control_fpUTR.react -n 0 -m 20
react_statistics.py -react hippuristanol_fpUTR.react -n 0 -m 20
react_statistics.py -react control_CDS.react -n 0 -m 20
react_statistics.py -react hippuristanol_CDS.react -n 0 -m 20
react_statistics.py -react control_tpUTR.react -n 125 -m 20
react_statistics.py -react hippuristanol_tpUTR.react -n 125 -m 20