#!/usr/bin/env python

from Bio import SeqIO

def reformat_fasta(afasta):
	'''Read in Fasta, returns dictionary with a space between transcript and info'''
	fasta_dict = {}
	fasta_sequences = SeqIO.parse(open(afasta),'fasta')
	for fasta in fasta_sequences:
		transcript = fasta.id.split('|')[0]
		info = '|'.join(fasta.id.split('|')[1:])
		fasta_dict[transcript + ' ' + info] = str(fasta.seq).upper()
	return fasta_dict

def write_fasta(dictionary, outfyle, LW=80):
	'''takes a dictionary and writes fasta'''
	with open(outfyle, 'w') as g:
		for k,v in dictionary.items():
			g.write('>'+k+'\n')
			for i in range(0,len(v),LW):
				g.write(v[i:i+LW]+'\n')

def main():
	input_dict=reformat_fasta('IsoSeq_MCF7_2015edition_polished.unimapped.fasta')
	write_fasta(input_dict,'MCF7_2015.fasta')

if __name__ == '__main__':
	main()