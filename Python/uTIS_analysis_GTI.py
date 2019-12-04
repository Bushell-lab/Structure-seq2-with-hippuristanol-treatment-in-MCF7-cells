#!/usr/bin/env python

import sys

def read_in_data(csv_fyle):
	"""reads in the csv file and returns a dict with ID as key and dicts as values with 
	annotation as key and a list of LTM reads as values. Only returns aTIS and 5'UTR"""
	adict = {}
	with open(csv_fyle,'r') as f:
		for line in f:
			if line.startswith("NM_"):
				geneID = line.strip().split(",")[2]
				LTM_reads = int(line.strip().split(",")[7])
				pos = int(line.strip().split(",")[4])
				if pos <= 1:
					if geneID not in adict:
						adict[geneID] = {pos:LTM_reads}
					else:
						adict[geneID][pos] = LTM_reads
	return(adict)

def calculate_uTIS_score(adict):
	"""takes the output dict from read_in_data and returns a uTIS score, which is calculated by
	dividing 5'UTR counts by all counts"""
	aTIS_counts = {}
	fp_counts = {}
	for k,v in adict.items():
		atis_counter = 0
		fp_counter = 0
		for pos, LTM_reads in v.items():
			if pos == 1:
				atis_counter += LTM_reads
			if pos < 1:
				fp_counter += LTM_reads
		aTIS_counts[k] = atis_counter
		fp_counts[k] = fp_counter
					
	uTIS_scores = {}
	for k,v in fp_counts.items():
		uTIS_scores[k] = float(fp_counts[k]) / (fp_counts[k] + aTIS_counts[k])
	
	return(uTIS_scores)

def main():
	data = read_in_data("sd01.csv")
	scores = calculate_uTIS_score(data)
	with open("GTI_data_uTIS_scores.csv", 'w') as g:
		g.write(",".join(["Ensembl_gene_symbol", "uTIS_score"]) + "\n")
		for k,v in scores.items():
			g.write(",".join([k, str(v)]) + "\n")

if __name__ == '__main__': 
	main()