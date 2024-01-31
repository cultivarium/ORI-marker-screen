from Bio import SeqIO
import glob

for fn in glob.glob('./*.gb'):
	print(fn)
	f = open(fn.split(".gb")[0] + ".fasta", 'w+')
	with open(fn) as input_handle:
		for record in SeqIO.parse(input_handle, "genbank"):
			f.write(">" + record.id + "\n")
			f.write(str(record.seq) + "\n")
	f.close()
