#!/usr/bin/env python

#Imports
from Bio import SeqIO
import argparse

#Functions
def read_in_fasta(afasta):
	'''Reads in a fasta file to a dictionary'''
	fasta_dict = {}
	fasta_sequences = SeqIO.parse(open(afasta),'fasta')
	for fasta in fasta_sequences:
		fasta_dict[fasta.id] = str(fasta.seq).upper()
	return fasta_dict

def write_fasta(dictionary, outfyle, LW=80):
	'''takes a dictionary and writes fasta'''
	with open(outfyle, 'w') as g:
		for k,v in dictionary.items():
			g.write('>'+k+'\n')
			for i in range(0,len(v),LW):
				g.write(v[i:i+LW]+'\n')

def main():
	parser = argparse.ArgumentParser(description = 'Trims FASTA files by the desired number of nts')
	parser.add_argument('fasta', type = str, help = 'FASTA file')
	parser.add_argument('n',type = int, help = '<int> number of nts to trim')
	parser.add_argument('-tp', action = 'store_true', help = 'Trim from the 3\' end (5\' default)')
	args = parser.parse_args()
	
	fasta = read_in_fasta(args.fasta)
	
	name_bit = 'FP' if args.tp == False else 'TP'
	out_name = '_'.join([args.fasta.replace('.fasta',''), 'trim', str(args.n), name_bit]) + '.fasta'
	trimmed_fasta = {}
	for k,v in fasta.items():
		if len(v)>args.n:
			trimmed_fasta[k] = v[args.n:] if args.tp == False else v[:-args.n]
			
	write_fasta(trimmed_fasta, out_name)


if __name__ == '__main__':
	main()