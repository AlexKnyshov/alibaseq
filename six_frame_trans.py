from Bio import SeqIO
from Bio.Seq import Seq
import sys
import multiprocessing

if len(sys.argv) == 3:
	infname = sys.argv[1]
	threads = int(sys.argv[2])

else:
	print "FORMAT: argument1 = file in fasta format, threads"
	print "EXAMPLE: ./fasta 4"
	print "output is written to [file name]_trans.fas"
	sys.exit()

def six_frame_translate(seq):
	rem = len(seq.seq) % 3
	if rem != 0:
		seq.seq = seq.seq+Seq("N"*(3-rem))
	result = {}
	for frame in range(3):
		result[">"+seq.id+"_f"+str(frame+1)+"l"+str(len(seq.seq))] = (seq.seq[frame:]+Seq("N"*frame)).translate()
		result[">"+seq.id+"_r"+str(frame+1)+"l"+str(len(seq.seq))] = (seq.seq.reverse_complement()[frame:]+Seq("N"*frame)).translate()
	return result

with open(infname) as infhandle:
	with open(infname+"_trans.fas", "w") as outhandle:
		seqs = SeqIO.parse(infhandle, "fasta")
		p = multiprocessing.Pool(threads)
		results = p.imap(six_frame_translate,seqs)
		for result in results:
			for seq1 in result:
				print >> outhandle, seq1
				print >> outhandle, result[seq1]

print "done"