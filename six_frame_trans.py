from Bio import SeqIO
from Bio.Seq import Seq
import sys


if len(sys.argv) == 2:
	infname = sys.argv[1]

else:
	print "FORMAT: argument1 = file in fasta format"
	print "EXAMPLE: ./fasta"
	print "output is written to [file name]_trans.fas"
	sys.exit()

with open(infname) as infhandle:
	with open(infname+"_trans.fas", "w") as outhandle:
		seqs = SeqIO.parse(infhandle, "fasta")
		for seq in seqs:
			rem = len(seq.seq) % 3
			if rem != 0:
				seq.seq = seq.seq+Seq("N"*rem)
			for frame in range(3):
				print >> outhandle, ">"+seq.id+"_f"+str(frame+1)+"l"+str(len(seq.seq))
				print >> outhandle, (seq.seq[frame:]+Seq("N"*frame)).translate()
				print >> outhandle, ">"+seq.id+"_r"+str(frame+1)+"l"+str(len(seq.seq))
				print >> outhandle, (seq.seq.reverse_complement()[frame:]+Seq("N"*frame)).translate()
