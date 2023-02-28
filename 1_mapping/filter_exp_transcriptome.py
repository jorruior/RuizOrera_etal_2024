#! /usr/bin/python
__author__="jruizor"
__date__ ="$Dec 19, 2021 12:24:43 PM$"
'''FIlter transcripts with FPKM > cut-off
'''

import sys
import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

fpkm = 0.5
transcriptome = sys.argv[1]
exp = sys.argv[2]
tag = sys.argv[3]

exc = ("ENSMMUT00000062240","none")

expt = []
expg = []
for e in exp.split(","):
	for line in open(e):
		if line.startswith("#"):
			continue
		if not "\ttranscript\t" in line:
			continue
		t = line.split('transcript_id "')[1].split('"')[0]
		if t in exc:
			continue
		g = line.split('gene_id "')[1].split('"')[0]
		if float(line.split('FPKM "')[1].split('"')[0]) >= fpkm:
			expt.append(t)
			expg.append(g)

expt = list(set(expt))
expg = list(set(expg))
print(str(len(expt)) + " transcripts in transcriptome")
print(str(len(expg)) + " genes in transcriptome")

out = open(tag,"w+")
for line in open(transcriptome):
	if line.startswith("#"):
		continue
	if "\tgene\t" in line:
		g = line.split('gene_id "')[1].split('"')[0]
		if g in expg:
			out.write(line)
	else:
		t = line.split('transcript_id "')[1].split('"')[0]
		if t in expt:
			out.write(line)

out.close()
exit(0)








