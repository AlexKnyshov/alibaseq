from Bio import SeqIO
import sys

contigfilename = sys.argv[1]
fastafilename = sys.argv[2]

query = []
with open(contigfilename) as contigfile:
	for contig in contigfile:
		query.append(contig.strip())

outfile = open("extracted_contigs.fas", "w")
with open(fastafilename) as fastafile:
	seqs = SeqIO.parse(fastafile, "fasta")
	for seq in seqs:
		if seq.id in query:
			print >> outfile, ">"+seq.id+"\n"+seq.seq
			query.remove(seq.id)
		if len(query) == 0:
			break
outfile.close()