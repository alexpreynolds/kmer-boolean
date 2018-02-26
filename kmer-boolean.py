#!/usr/bin/env python

import sys
import os
import bitarray

# read FASTA
def read_sequences():
	global seqs
	seqs = []
	seq = ""
	for line in sys.stdin:
		if line.startswith('>'):
			if len(seq) > 0:
				seqs.append(seq)
				seq = ""
		else:
			seq += line.strip()

	seqs.append(seq)

# build and initialize bit array
def initialize_bitarray():
	global ba 
	ba = bitarray.bitarray(4**k)
	ba.setall(False)
	sys.stderr.write("Memory usage of bitarray.bitarray is [%ld] bytes\n" % (ba.buffer_info()[1]))

# process sequences
def process_sequences():
	global observed_kmers
	observed_kmers = {}
	for seq in seqs:
		for i in range(0, len(seq)):
			kmer = seq[i:i+k]
			if len(kmer) == k:
				observed_kmers[kmer] = None
				idx = 0
				for j in range(k-1, -1, -1):
					idx += 4**(k-j-1) * bm[kmer[j]]
				ba[idx] = True

def test_bitarray():
	test_idx = 0
	for j in range(k-1, -1, -1):
		test_idx += 4**(k-j-1) * bm[test_kmer[j]]
	test_result = ba[test_idx]
	if test_result:
		sys.stdout.write("%s found\n" % (test_kmer))
		sys.exit(os.EX_OK)
	else: 
		sys.stdout.write("%s not found\n" % (test_kmer))
		sys.exit(os.EX_DATAERR)

def main():
	global k
	k = int(sys.argv[1])

	global bm 
	bm = { 'A' : 0, 'C' : 1, 'T' : 2, 'G' : 3 }

	read_sequences()
	initialize_bitarray()
	process_sequences()

	global test_kmer
	try:
		global test_kmer
		test_kmer = sys.argv[2]
		if len(test_kmer) == k:
			test_bitarray()
		else:
			raise ValueError("test kmer (%s) should be of length k (%d)" % (test_kmer, k))
	except IndexError as err:
		keys = list(observed_kmers.keys())
		for i in range(0, len(keys)):
			sys.stdout.write("%s found\n" % (keys[i]))

	sys.exit(os.EX_OK)
  
if __name__== "__main__":
 	main()
