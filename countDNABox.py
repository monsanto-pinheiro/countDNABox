#!/usr/bin/env python
import sys, os, re
from optparse import OptionParser

vectAmbigous = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '*']
dt_ambigous = { 'R':'[AG]', 'Y':'[TC]', 'K':'[GT]', 'M':'[AC]', 'S':'[GC]', 
			'W':'[AT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]',
			'N':'[ACGT]', '*':'.' }

#complement
def complement(seq):  
	dict_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
		'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 
		'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 
		'N': 'N'}  
	complseq = [dict_complement[base] for base in seq]  
	return complseq


#reverse
def reverse_complement(seq):  
	seq = list(seq)  
	seq.reverse()   
	return ''.join(complement(seq))

def ambiguos_to_unambiguous(sequence):
	for ambig in vectAmbigous:
		sequence = sequence.replace(ambig, dt_ambigous[ambig])
	#print sequence
	return sequence


def FASTA_to_dict(file_, szType):
	"""Read FASTA file to a dictionary\n
	FASTA_to_dict(file)"""
	from Bio import SeqIO
	well = open(file_)
	spots_dict = SeqIO.to_dict(SeqIO.parse(well, szType))  # "fastq" "fasta"
	well.close()
	return spots_dict

## save output based on dictionary
def save_output(dt_data, handle_out, sz_header):
	
	handle_out.write("\n" + sz_header)
	vect_out = dt_data.keys()
	vect_out = sorted(vect_out)
	handle_out.write("\n\tgap\t#occurrence")
	n_total = 0
	for key in vect_out:
		handle_out.write("\n\t%d\t%s" % (key, dt_data[key]))
		n_total += int(dt_data[key])
	handle_out.write("\n\tTotal\t%d\n\n########################\n########################\n\n" % (n_total))

if __name__ == '__main__':

	parser = OptionParser(usage="%prog [-h] [-s] arguments", version="%prog 0.1", add_help_option=False,)
	parser.add_option("-f", "--format", type="string", dest="format", help="Format of in file, [fasta|fastq]", metavar="FORMAT")
	parser.add_option("-o", "--dir_output", type="string", dest="dir_output", help="Output directory", metavar="OUT_DIR")
	parser.add_option("-i", "--input", type="string", dest="input", help="Files separated by comma without space", metavar="IN_FILES")
	parser.add_option("-d", "--dna_box_1", type="string", dest="dna_box_1", help="DNA sequence to look for at first place", metavar="DNA")
	parser.add_option("-n", "--dna_box_2", type="string", dest="dna_box_2", help="Last DNA sequence to look for", metavar="DNA")
	parser.add_option("-s", action="store_true", dest="save_fail", help="Save sequences that doesn't have the DNA boxes in the sequence [False]", default=False)
	parser.add_option('-h', '--help', dest='help', action='store_true', help='show this help message and exit')

	(options, args) = parser.parse_args()
	
	if (options.help):
		parser.print_help()
		print "It can consume one or more files, separated them by comma, but without space"
		print "example: countDNABox -i file.fasta -d ATGTAGTAGTTGATGT -n TATGTAGATGATGAT -o out_dir -f fasta"
		print "example: countDNABox -i file.fasta,file1.fasta -d RYGTAGTAGTTGATGT -n TATGTAGATGATGAT -o out_dir -f fasta"
		print 
		print "\tprint ambiguous nucleotides:"
		count = 1
		for ambig in vectAmbigous:
			print "\t%d) %s -> %s" % (count, ambig, dt_ambigous[ambig])
			count += 1
		sys.exit(0)

	if (len(args) != 0):
		parser.error("incorrect number of arguments")
	
	if not options.format:   # 
		parser.error('Format not specified [fasta|fastq]')
	if not options.dir_output:   # 
		parser.error('Output directory not specified')
	if not options.input:   # 
		parser.error('File[s] not specified')
	if not options.dna_box_1:   # 
		parser.error('First DNA box not specified, could be degenerate')
	if not options.dna_box_2:   # 
		parser.error('Second DNA box not specified, could be degenerate')

	if (options.format <> "fasta" and options.format <> "fastq"):
		sys.exit("Error: input format file must be 'fasta' or 'fastq'")
	
	### need more dnaBox
	print "\nFind DNA boxes: " + options.dna_box_1 + "  " + options.dna_box_2
	
	pattern1 = re.compile(options.dna_box_1)
	pattern2 = re.compile(options.dna_box_2)
	
	patternReverse1 = re.compile(reverse_complement(options.dna_box_1))
	patternReverse2 = re.compile(reverse_complement(options.dna_box_2))
	
	### if dir doesn't exist make it
	if not os.path.exists(options.dir_output): os.makedirs(options.dir_output)

	handle_forward = open(options.dir_output + "/seq_forward.fasta", 'w')
	handle_reverse = open(options.dir_output + "/seq_reverse.fasta", 'w')
	handle_both = open(options.dir_output + "/seq_both.fasta", 'w')
	if (options.save_fail): handle_fail = open(options.dir_output + "/seq_fail.fasta", 'w')
	dt_data_forward = {}
	dt_data_reverse = {}
	dt_data_both = {}
	
	n_count = 0
	(n_count_forward, n_count_reverse, n_count_fail) = (0, 0, 0) 
	
	#### Testing files
	print "Testing files: "
	for sz_file in options.input.split(','):
		## start processing
		if (not os.path.exists(sz_file)):
			sys.exit("Error: FASTA file '%s' doesn't exist " % sz_file)
		print "\t" + sz_file + " -> exist"
	
	print	
	for sz_file in options.input.split(','):
		## start processing
		if (not os.path.exists(sz_file)):
			sys.exit("Error: FASTA file '%s' doesn't exist " % sz_file)

		dictRef = FASTA_to_dict(sz_file, options.format)
		if (len(dictRef) == 0):
			sys.exit("Error: FASTA file '%s' doesn't have any sequence " % sz_file)

		print "Processing %d sequences" % (len(dictRef))
		for key in dictRef:
			if ((n_count % 1000) == 0 and n_count > 1):
				print "%s -> %d sequences processed" % (sz_file, n_count)
			n_count += 1
			
			sz_seq = str(dictRef[key].seq).upper()
			n_start = -1
			n_end = -1
			for findRe in re.finditer(pattern1, sz_seq):
				n_start = findRe.start()
	
			for findRe in re.finditer(pattern2, sz_seq):
				n_end = findRe.start()
	
			## forward
			if (n_start < n_end and n_end <> -1 and n_start <> -1 and (n_start + len(options.dna_box_1)) < n_end):
				n_size = n_end - n_start - len(options.dna_box_1)
				n_count_forward += 1
				if (dt_data_forward.has_key(n_size)): dt_data_forward[n_size] += 1
				else: dt_data_forward[n_size] = 1
				if (dt_data_both.has_key(n_size)): dt_data_both[n_size] += 1
				else: dt_data_both[n_size] = 1
				handle_forward.write("%s\n" % (sz_seq[n_start+len(options.dna_box_1):n_end]))
				handle_both.write("%s\n" % (sz_seq[n_start+len(options.dna_box_1):n_end]))
				continue
			
			n_start = -1
			n_end = -1
			for findRe in re.finditer(patternReverse2, sz_seq):
				n_start = findRe.start()
			for findRe in re.finditer(patternReverse1, sz_seq):
				n_end = findRe.start()
	
			## forward
			if (n_start < n_end and n_end <> -1 and n_start <> -1 and (n_start + len(options.dna_box_2)) < n_end):
				n_size = n_end - n_start - len(options.dna_box_2)
				n_count_reverse += 1
				if (dt_data_reverse.has_key(n_size)): dt_data_reverse[n_size] += 1
				else: dt_data_reverse[n_size] = 1
				if (dt_data_both.has_key(n_size)): dt_data_both[n_size] += 1
				else: dt_data_both[n_size] = 1
				handle_reverse.write("%s\n" % (sz_seq[n_start+len(options.dna_box_2):n_end]))
				handle_both.write("%s\n" % (sz_seq[n_start+len(options.dna_box_2):n_end]))
				continue
	
			### fail
			if (options.save_fail): handle_fail.write(sz_seq + "\n")
			n_count_fail += 1
		
	print "%d sequences processed" % (n_count)
	print "writing reports"
	### save output
	handle_out = open(options.dir_output + "/report.txt", 'w')
	handle_out.write("File\t%s\nBox1\t%s\nBox2\t%s\nForward\t%d\nReverse\t%d\nForward + nReverse\t%d\nFail\t%d\nTotal seq\t%d\n" % (options.input,
						options.dna_box_1, options.dna_box_2, n_count_forward, n_count_reverse, n_count_reverse + n_count_forward, n_count_fail,
						n_count_reverse + n_count_forward + n_count_fail))
	save_output(dt_data_forward, handle_out, "Forward")
	save_output(dt_data_reverse, handle_out, "Reverse")
	save_output(dt_data_both, handle_out, "Both")
	
	handle_out.close()
	handle_forward.close()
	handle_reverse.close()
	if (options.save_fail): handle_fail.close()
	handle_both.close()
	print "finished"
