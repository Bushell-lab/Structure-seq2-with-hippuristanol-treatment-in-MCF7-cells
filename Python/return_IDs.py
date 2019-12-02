#!/usr/bin/env python
import sys

def read_IDs(afastafile, outfyle):
	'''Read in Fasta, returns a list of IDs'''
	IDs = []
	with open(afastafile,'r') as f:
		with open(outfyle, 'w') as g:
			for line in f:
				if line.startswith('>'):
					g.write(line.strip()[1:])
					g.write('\n')

def main():
	fyle_name = sys.argv[1].split('.')[0]
	outname = fyle_name + "_IDs.txt"
	IDs = read_IDs(sys.argv[1], outname)

if __name__ == '__main__': 
	main()