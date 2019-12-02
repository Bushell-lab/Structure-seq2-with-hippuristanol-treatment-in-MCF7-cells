#!/usr/bin/env python

'''
This script will splice the reformatted MCF7_2015.fasta file into a 5'UTR, CDS and 3'UTR fasta file. It uses the output from a blast in which the MCF7_2015
sequences are the subjet and the NM_ RefSeq CDS sequences from release 85 are the query. 
The 5'UTR length is determined using the output of the blast filtered to the smallest value for the query start per subject transcript. Only when the query start equals position 1, 
i.e. it aligns with the start of the RefSeq CDS, is the 5'UTR length determined by this blast hit.
The 3'UTR length is determined using the output of the blast filtered to the highest value for the query end per subject transcript. Only when the query end equals the length 
of the RefSeq CDS, i.e. it aligns with the end of the RefSeq CDS, is the 3'UTR length determined by this blast hit.
Only transcripts for which both the 5'UTR and 3'UTR lengths are determined are spliced, and only when the resulting CDS is equally divisible by 3. 
'''

from itertools import islice
from Bio import SeqIO

class Blast_Hit(object):
  '''makes a class for blast hits'''
  def __init__(self, query='',subject='',percent_id=0,align_len=0,mismatches=0,gaps=0,q_start=0,q_end=0,s_start=0,s_end=0,evalue=0,bitscore=0):
	'''Blast Hit'''
	self.query = query
	self.subject = subject
	self.percent_id = percent_id
	self.align_len = align_len
	self.mismatches = mismatches
	self.gaps = gaps
	self.q_start = q_start
	self.q_end = q_end
	self.s_start = s_start
	self.s_end = s_end
	self.evalue = evalue
	self.bitscore = bitscore

  @classmethod
  def from_tabular_line(klass, blast_line):
	'''Tabular format (BLAST+ -outfmt 6) -> Blast Hit, i.e.: my_hit = Blast_Hit.from_tabular_line(line)'''
	query, subject, percentid, alignmentlength, mismatches, gaps, qstart, qend, sstart, send, evalue, bitscore = blast_line.strip().split()
	return klass(query, subject, float(percentid), int(alignmentlength), int(mismatches), int(gaps), int(qstart), int(qend), int(sstart), int(send), float(evalue), float(bitscore))

def read_in_fasta(afasta):
	'''Reads in a fasta file to a dictionary'''
	fasta_dict = {}
	fasta_sequences = SeqIO.parse(open(afasta),'fasta')
	for fasta in fasta_sequences:
		fasta_dict[fasta.id] = str(fasta.seq).upper()
	return fasta_dict

def read_in_fasta_lengths_csv(csv_fyle):
	'''reads in a csv returns a dict'''
	adict = {}
	with open(csv_fyle, 'r') as f:
		line_count = 0
		for line in f:
			if line_count > 0:
				transcript = line.strip().split(',')[0]
				length = line.strip().split(',')[1]
				adict[transcript]=length
			line_count += 1
	return adict

def extract_fp_lengths(blast_results, PB_seqs):
	'''takes blast results filtered to smallest query start as tsv,
	returns a dict with MCF7_2015 ID as key and 5'UTR length as value'''
	adict = {}
	with open(blast_results, 'r') as f:
		for line in f:
			hit = Blast_Hit.from_tabular_line(line)
			if hit.q_start == 1:#checks that the blast hit starts at position 1
				fpUTR_length = hit.s_start - 1#calculates the length of the 5'UTR
				adict[hit.subject] = fpUTR_length
	return adict

def extract_tp_lengths(blast_results, rs_cds_lens, PB_seqs):
	'''takes blast results filtered to highest query end as tsv,
	returns a dict with MCF7_2015 ID as key and 3'UTR length as value'''
	adict = {}
	with open(blast_results, 'r') as f:
		for line in f:
			hit = Blast_Hit.from_tabular_line(line)
			PB_len = len(PB_seqs[hit.subject])#calculates the length of the MCF7_2015 transcript
			CDS_len = int(rs_cds_lens[hit.query])#extracts the length of the coding sequence
			if hit.q_end == CDS_len:#checks that the blast hit has extended to the end of the coding sequence
				tpUTR_length = PB_len - hit.s_end#calculates the length of the 3'UTR
				adict[hit.subject] = tpUTR_length#adds the 3'UTR length to the dictionary
	return adict

def annotate_transcripts(fp_lengths, tp_lengths, PB_seqs):
	'''takes the 5UTR lengths and 3UTR lengths dictionaries and for those transcripts which are in both,
	it calculates the CDS length by subtracting the 5'UTR and 3'UTR lengths from the total length. 
	If the calculated CDS length is equally divisible by 3 it stores the 5'UTR and 3'UTR lengths in a dictionary'''
	adict = {}
	for k,v in fp_lengths.items():
		fp_len = v
		if k in tp_lengths:
			tp_len = tp_lengths[k]
			ORF_length = (len(PB_seqs[k])) - fp_len - tp_len#calculates the coding sequence length
			if ORF_length > 0:
				if ORF_length % 3 == 0:#checks whether the coding sequence length is equally divisible by 3
					adict[k] = [fp_len, tp_len] #adds the 3'UTR and 5'UTR length to the dictionary
	return adict

def splice_features(MCF7_seqs, annotation_dict):
	'''Takes the FASTA dict and the feature length dicts and
	creates three new dicts with spliced feature sequences'''
	fp_seqs, cds_seqs, tp_seqs = {},{},{}
	for k,v in MCF7_seqs.items():
		if k in annotation_dict:
			fp_length = annotation_dict[k][0]
			tp_length = annotation_dict[k][1]
			if fp_length > 0:#some transcripts are 5'UTR less
				fp_seqs[k] = v[:fp_length]#splices 5'UTR sequences and adds to fp_seqs
			cds_seqs[k] = v[fp_length:-tp_length]#splices CDS sequences and adds to cds_seqs
			if tp_length > 0:#some transcripts are 3'UTR less
				tp_seqs[k] = v[-tp_length:]#splices 3'UTR sequences and adds to tp_seqs
	return fp_seqs, cds_seqs, tp_seqs

def write_fasta(dictionary, outfyle, LW=80):
	'''takes a dictionary and writes fasta'''
	with open(outfyle, 'w') as g:
		for k,v in dictionary.items():
			g.write('>' + k + '\n')
			for i in range(0, len(v), LW):
				g.write(v[i:i+LW] + '\n')

def main():
	MCF7_2015_seqs = read_in_fasta('MCF7_2015.fasta')
	refseq_CDS_lens = read_in_fasta_lengths_csv('All_R_V_HumanRefSeqs_Release85_CDSs_lengths.csv')
	fpUTR_lengths = extract_fp_lengths('NM_CDSs_to_MCF7_2015_blast_eris_start.tsv', MCF7_2015_seqs)
	tpUTR_lengths = extract_tp_lengths('NM_CDSs_to_MCF7_2015_blast_eris_end.tsv', refseq_CDS_lens, MCF7_2015_seqs)
	annotated_transcripts = annotate_transcripts(fpUTR_lengths, tpUTR_lengths, MCF7_2015_seqs)
	fp_dict, cds_dict, tp_dict = splice_features(MCF7_2015_seqs, annotated_transcripts)
	write_fasta(fp_dict, 'MCF7_2015_fpUTRs.fasta')
	write_fasta(cds_dict, 'MCF7_2015_CDSs.fasta')
	write_fasta(tp_dict, 'MCF7_2015_tpUTRs.fasta')

if __name__ == '__main__': 
	main()