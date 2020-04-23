#!/usr/bin/env python

#Imports
import argparse
from itertools import islice

def read_in_react(areact):
	'''Reads in a reactivity file, returns a dictionary'''
	information = {}
	with open(areact,'r') as f:
		while True:
			next_n_lines = list(islice(f, 2))
			if not next_n_lines:
				break
			transcript,reactivities = [n.strip() for n in next_n_lines]
			information[transcript] = reactivities.split('\t')
	return information

def write_react(adict, outfyle):
	'''reads in a dictionary and writes a react file with tab seperated values'''
	with open(outfyle, 'w') as out:
		for k, v in adict.items():
			out.write(k+'\n'+'\t'.join(v)+'\n')

def main():
	parser = argparse.ArgumentParser(description='Trims n number of nts from either the 5\' or 3\' end of a react file')
	parser.add_argument('react_file', help = '<str> Operate on this react file')
	parser.add_argument('n', type=int, help='<int> number of nts from 5\' or 3\' end to remove')
	parser.add_argument('-tp',action='store_true', help='Start from the 3\' end (5\' default)')
	args = parser.parse_args()
	
	react_dict = read_in_react(args.react_file)
	name_bit = 'FP' if args.tp == False else 'TP'
	out_name = '_'.join([args.react_file.replace('.react',''), 'trim', str(args.n), name_bit]) + '.react'
	trimmed_react_dict = {}
	for k,v in react_dict.items():
		if len(v)>args.n:
			trimmed_react_dict[k] = v[args.n:] if args.tp == False else v[:-args.n]
			
	write_react(trimmed_react_dict, out_name)

if __name__ == '__main__': 
	main()
