import sys
from Bio import SeqIO

fastq_location = str(sys.argv[1])

search_rnas = ["AGGGACGGGACGCGGTGCAGTG", "ACTGCCCCAGGTGCTGCTGG",	
"AGGGAATAGTTGCTGTGCTGTA", "TATCACAGTGGCTGTTCTTTTT",
"ATAAAGCTAGACAACCATTGA", "TGTGACTGGTTGACCAGAGGGG",	
"GGGGTTCCTGGGGATGGGATTT", "TTTGGTCCCCTTCAACCAGCTG",	
"CTGAAGCTCAGAGGGCTCTGAT", "TGGCTCAGTTCAGCAGGAACAG",	
"CTATACGGCCTCCTAGCTTTCC", "CCAATATTGGCTGTGCTGCTCC",
"TATTGCACTTGTCCCGGCCTGT", "AATCACTAACCACACGGCCAGG",	
"CTGGGAGAGGGTTGTTTACTCC", "AAAAGCTGGGTTGAGAGGGCGA",	
"ATGGTTCCGTCAAGCACCATGG", "TGGCAGTGTCTTAGCTGGTTGT",	
"TAAAGTGCTGACAGTGCAGAT", "TAACACTGTCTGGTAAAGATGG",
"TGAGGTAGTAGGTTGTATGGTT", "TAGCAGCACATAATGGTTTGTG",
"TAGCTTATCAGACTGATGTTGA", "TAGCACCATTTGAAATCAGTGTT",
"TGATATGTTTGATATATTAGGT"]

num_occurrences = 5
thresh = 3.0

for search_rna in search_rnas:
	barcodes = []
	barcode_counts = []
	found_counts = 0.0
	#ACTG
	freqs = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
	for y in SeqIO.parse(fastq_location, "fastq"):
		if search_rna.upper() in str(y.seq.upper()):
			read = str(y.seq.upper())
			fourN_1 = read[:4]
			fourN_2 = read[len(search_rna)+4:len(search_rna)+8]
			eightN = fourN_1 + fourN_2
			if len(eightN) == 8:
				found_counts += 1
				for count in range(8):
					if eightN[count] == 'A':
						freqs[count][0] += 1
					elif eightN[count] == 'C':
						freqs[count][1] += 1
					elif eightN[count] == 'T':
						freqs[count][2] += 1
					elif eightN[count] == 'G':
						freqs[count][3] += 1
				if eightN not in barcodes:
					barcodes.append(eightN)
					barcode_counts.append(1)
				else:
					barcode_counts[barcodes.index(eightN)] += 1	
			
	#print "PCR bias test"
	#print "[5' adapter-(N1)(N2)(N3)(N4)-miRNA-(N5)(N6)(N7)(N8)-3' adapter]"
	#print "8N seq    % occurrence"
	print search_rna
	for i in barcodes:
		if (barcode_counts[barcodes.index(i)]/found_counts)*100 >= thresh:
			print i, str("{:6.3f}".format((barcode_counts[barcodes.index(i)]/found_counts)*100))
	print "\n"
	#print found_counts
	# Prints the percent occurrence of A,C,T,G for each position in the 8N randomized sequence,
	# where the first four numbers represent A,C,T,G at the first position, the next four
	# represent A,C,T,G at the second position, and so on.
	# [5' adapter-(N1)(N2)(N3)(N4)-miRNA-(N5)(N6)(N7)(N8)-3' adapter]
	#Ncount = 1
	#for freq in freqs:
	#	print "% A,C,T,G @ position N"+str(Ncount)
	#	Ncount += 1
	#	for f in freq:
	#		print str("{:8.3f}".format((f/found_counts)*100))
