#! /usr/bin/python3

from Bio import SeqIO
import re, os, sys

# Interpreting external parameters
nrparameters = len(sys.argv)
if nrparameters > 1:
	fastafile = sys.argv[1]
else:
	sys.exit("Program.py fastafile")

# Starting the files
input_file = open(fastafile, "r") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
output_file = open('Aminoacid_data.tsv','w')
output_file.write('Locus_tag\tProtein\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\tB\tZ\tX\tLength_(aa)\n')

#Data of interest
for cur_record in SeqIO.parse(input_file, "fasta") :
	protein_descriptor = cur_record.description
	if re.match("^\S+\s+\S+\s+\S+", protein_descriptor):
		target1 = re.search("^\S+\s+\[gene\=(\S+)\]\s+\[locus\_tag\=(\S+)\]", protein_descriptor)
		protein_name = target1.group(1)
		locus_tag = target1.group(2)
	A_count = cur_record.seq.count('A')
	C_count = cur_record.seq.count('C')
	D_count = cur_record.seq.count('D')
	E_count = cur_record.seq.count('E')
	F_count = cur_record.seq.count('F')
	G_count = cur_record.seq.count('G')
	H_count = cur_record.seq.count('H')
	I_count = cur_record.seq.count('I')
	K_count = cur_record.seq.count('K')
	L_count = cur_record.seq.count('L')
	M_count = cur_record.seq.count('M')
	N_count = cur_record.seq.count('N')
	P_count = cur_record.seq.count('P')
	Q_count = cur_record.seq.count('Q')
	R_count = cur_record.seq.count('R')
	S_count = cur_record.seq.count('S')
	T_count = cur_record.seq.count('T')
	V_count = cur_record.seq.count('V')
	W_count = cur_record.seq.count('W')
	Y_count = cur_record.seq.count('Y')
	B_count = cur_record.seq.count('B')
	Z_count = cur_record.seq.count('Z')
	X_count = cur_record.seq.count('X')
	length = len(cur_record.seq)
	output_line = '%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' % (locus_tag, protein_name, A_count, C_count, D_count, E_count, F_count, G_count, H_count, I_count, K_count, L_count, M_count, N_count, P_count, Q_count, R_count, S_count, T_count, V_count, W_count, Y_count, B_count, Z_count, X_count, length)
	output_file.write(output_line)
output_file.close()
input_file.close() 
