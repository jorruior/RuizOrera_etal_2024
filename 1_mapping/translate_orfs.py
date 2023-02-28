#! /usr/bin/python
__author__="jruiz"
__date__ ="$Aug 27, 2022 15:18:00 PM$"

import sys
import os
import Bio
from Bio import SeqIO
from Bio.Seq import Seq

##ARGUMENTS

fasta = sys.argv[1]
seqs = SeqIO.index(fasta, "fasta")

for seq in seqs:
	prot = seqs[seq].seq.translate()
	print(">" + seq + "\n" + str(prot) + "\n", end = "")

exit(0)
