import sys
import string
import subprocess
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

ens_gtf_file = sys.argv[1]
str_gtf_file = sys.argv[2]

exceptions = ("CABD030151935.1","CABD030151492.1","CABD030130636.1","CABD030130064.1","CABD030151826.1","CABD030151588.1","CABD030161071.1")

class trans_object:
	def __init__(self, chrm, gene, strand, start, end, biotype, biotype2, code, exonnum):
		self.chrm = chrm
		self.gene = gene
		self.strand = strand
		self.start = start
		self.end = end
		self.biotype = biotype
		self.biotype2 = biotype2
		self.code = code
		self.exonnum = exonnum

def parse_gtf_annot(gtf,field):
	#Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	total_genes = {}
	for line in open(gtf):
		if not "\t" + field + "\t" in line:
			continue
		t_name = line.split('transcript_id "')[1].split('"')[0]
		g_name = line.split('gene_id "')[1].split('"')[0]
		exonnum = int(line.split('exon_number "')[1].split('"')[0])

		if 'gene_biotype' in line:
			biot = line.split('gene_biotype "')[1].split('"')[0]
		else:
			biot = "unknown"

		if 'transcript_biotype' in line:
			biot2 = line.split('transcript_biotype "')[1].split('"')[0]
		else:
			biot2 = "unknown"	

		if not g_name in total_genes:
			total_genes[g_name] = [biot,[]]
		if not t_name in total_genes[g_name][1]:
			total_genes[g_name][1].append(t_name)
		
		trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[],biot,biot2,[], []))
		trans[t_name].code.append(line.split("\t")[7])
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))
		trans[t_name].exonnum.append(exonnum)

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]
	for x in trans:
		if trans[x].strand == "-":
			trans[x].exonnum = trans[x].exonnum[::-1]
			trans[x].code = trans[x].code[::-1]

	return trans,total_genes

def parse_gtf(gtf,field):
	#Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	codes = {}
	ids = {}
	gene_trans = {}
	for line in open(gtf):
		if "\ttranscript\t" in line:
			codes[line.split('transcript_id "')[1].split('"')[0].replace("_","-").replace(".","-")] = "="
			if "class_code" in line:
				codes[line.split('transcript_id "')[1].split('"')[0].replace("_","-").replace(".","-")] = line.split('class_code "')[1].split('"')[0]
			if "ref_gene_id" in line:
				ids[line.split('transcript_id "')[1].split('"')[0].replace("_","-").replace(".","-")] = line.split('ref_gene_id "')[1].split('"')[0]
				if not line.split('ref_gene_id "')[1].split('"')[0].replace("_","-").replace(".","-") in gene_trans:
					gene_trans[line.split('ref_gene_id "')[1].split('"')[0].replace("_","-").replace(".","-")] = []
	for line in open(gtf):
		if not "\t" + field + "\t" in line:
			continue
		t_name = line.split('transcript_id "')[1].split('"')[0].replace("_","-").replace(".","-")
		g_id = line.split('gene_id "')[1].split('"')[0].replace("_","-").replace(".","-")
		try:
			g_name = ids[t_name]
		except:
			g_name = g_id

		if 'gene_biotype' in line:
			biot = line.split('gene_biotype "')[1].split('"')[0]
		else:
			biot = "novel"

		if 'transcript_biotype' in line:
			biot2 = line.split('transcript_biotype "')[1].split('"')[0]
		else:
			biot2 = "novel"		


		trans.setdefault(t_name,trans_object(line.split("\t")[0],g_id + "/" + g_name,line.split("\t")[6],[],[],biot,biot2,codes[t_name],0))
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))		

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans


(ens_gtf_sc1,total_genes) = parse_gtf_annot(ens_gtf_file,"start_codon")
(ens_gtf_sc2,total_genes) = parse_gtf_annot(ens_gtf_file,"stop_codon")
(ens_gtf_CDS,total_genes) = parse_gtf_annot(ens_gtf_file,"CDS")
(ens_gtf_exon,total_genes) = parse_gtf_annot(ens_gtf_file,"exon")

str_gtf = parse_gtf(str_gtf_file,"exon")
lines = {}
boundaries = {}
for trans in str_gtf:
	gene = str_gtf[trans].gene.split("/")[1]
	code = str_gtf[trans].code
	gene_biotype = "processed_transcript"
	if ("AACZ0" in str_gtf[trans].chrm) or ("KV42" in str_gtf[trans].chrm) or (str_gtf[trans].chrm in exceptions):
		continue 
	if str_gtf[trans].strand == ".":
		continue
	if code in ("p","s","c"):
		continue
	
	if gene in total_genes:
		if total_genes[gene] == "protein_coding":
			if code != ("k","i","m","n"):
				continue
		gene_biotype = total_genes[gene][0]
	
	if not str_gtf[trans].gene.split("/")[0] in lines:
		lines[str_gtf[trans].gene.split("/")[0]] = [[],[]]

	if not str_gtf[trans].gene.split("/")[1] in lines[str_gtf[trans].gene.split("/")[0]][1]:
		lines[str_gtf[trans].gene.split("/")[0]][1].append(str_gtf[trans].gene.split("/")[1])	

	trans_biotype = "processed_transcript"
	trans_target = "none"
	if gene_biotype == "protein_coding":
		for trans2 in total_genes[gene][1]:
			if trans_biotype == "protein_coding":
				break
			if (len(str_gtf[trans].start) == 1) and (len(ens_gtf_exon[trans2].start) > 1):
				continue
			if (len(str_gtf[trans].start) > 1) and (len(ens_gtf_exon[trans2].start) == 1):
				continue			
			if (len(str_gtf[trans].start) == 1) or ((ens_gtf_exon[trans2].start[1:] == str_gtf[trans].start[1:]) and (ens_gtf_exon[trans2].end[:-1] == str_gtf[trans].end[:-1])):
				if not (trans2 in ens_gtf_sc1) or not (trans2 in ens_gtf_sc2):
					continue
				trans_biotype = ens_gtf_exon[trans2].biotype2					
				if (trans_biotype == "protein_coding") or (trans_biotype == "nonsense_mediated_decay"): 
					trans_target = trans2
					
		if trans_target != "none":
			if (str_gtf[trans].strand == "+") and ((ens_gtf_sc1[trans_target].start[0] < str_gtf[trans].start[0]) or (ens_gtf_sc2[trans_target].end[-1] > str_gtf[trans].end[-1])):
				trans_biotype = "processed_transcript"
			elif (str_gtf[trans].strand == "-") and ((ens_gtf_sc2[trans_target].start[0] < str_gtf[trans].start[0]) or (ens_gtf_sc1[trans_target].end[-1] > str_gtf[trans].end[-1])):
				trans_biotype = "processed_transcript"
			else:
				for n,st in enumerate(ens_gtf_sc1[trans_target].start):				
					lines[str_gtf[trans].gene.split("/")[0]][0].append(str_gtf[trans].chrm + "\tStringtie\tstart_codon\t" + str(ens_gtf_sc1[trans_target].start[n]) + "\t" + str(ens_gtf_sc1[trans_target].end[n]) + "\t.\t" + str_gtf[trans].strand + "\t.\tgene_id \"" + str_gtf[trans].gene.split("/")[0] + "\"; transcript_id \"" + trans + "\"; gene_name \"COMPLETE\"; transcript_name \"" + trans + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_biotype \"" + trans_biotype + "\"; exon_number \"" + str(ens_gtf_sc1[trans_target].exonnum[n]) + "\";\n")
				for n,st in enumerate(ens_gtf_sc2[trans_target].start):
					lines[str_gtf[trans].gene.split("/")[0]][0].append(str_gtf[trans].chrm + "\tStringtie\tstop_codon\t" + str(ens_gtf_sc2[trans_target].start[n]) + "\t" + str(ens_gtf_sc2[trans_target].end[n]) + "\t.\t" + str_gtf[trans].strand + "\t.\tgene_id \"" + str_gtf[trans].gene.split("/")[0] + "\"; transcript_id \"" + trans + "\"; gene_name \"COMPLETE\"; transcript_name \"" + trans + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_biotype \"" + trans_biotype + "\"; exon_number \"" + str(ens_gtf_sc2[trans_target].exonnum[n]) + "\";\n")			
				for n,st in enumerate(ens_gtf_CDS[trans_target].start):
					lines[str_gtf[trans].gene.split("/")[0]][0].append(str_gtf[trans].chrm + "\tStringtie\tCDS\t" + str(ens_gtf_CDS[trans_target].start[n]) + "\t" + str(ens_gtf_CDS[trans_target].end[n]) + "\t.\t" + str_gtf[trans].strand + "\t" + ens_gtf_CDS[trans_target].code[n] + "\tgene_id \"" + str_gtf[trans].gene.split("/")[0] + "\"; transcript_id \"" + trans + "\"; gene_name \"COMPLETE\"; transcript_name \"" + trans + "\"; protein_id \"CDS" + trans + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_biotype \"" + trans_biotype + "\"; exon_number \"" + str(ens_gtf_CDS[trans_target].exonnum[n]) + "\";\n")

	elif trans in ens_gtf_exon:
		trans_biotype = ens_gtf_exon[trans].biotype2

	boundaries[trans] = [str_gtf[trans].start[0],str_gtf[trans].end[-1]]
	if not str_gtf[trans].gene.split("/")[0] in boundaries:
		boundaries[str_gtf[trans].gene.split("/")[0]] = boundaries[trans]
	if boundaries[trans][0] < boundaries[str_gtf[trans].gene.split("/")[0]][0]:
		boundaries[str_gtf[trans].gene.split("/")[0]][0] = boundaries[trans][0]
	if boundaries[trans][1] > boundaries[str_gtf[trans].gene.split("/")[0]][1]:
		boundaries[str_gtf[trans].gene.split("/")[0]][1] = boundaries[trans][1]

	for n,st in enumerate(str_gtf[trans].start):
		if n == 0:
			lines[str_gtf[trans].gene.split("/")[0]][0].append(str_gtf[trans].chrm + "\tStringtie\ttranscript\t" + str(boundaries[trans][0]) + "\t" + str(boundaries[trans][1]) + "\t.\t" + str_gtf[trans].strand + "\t.\tgene_id \"" + str_gtf[trans].gene.split("/")[0] + "\"; transcript_id \"" + trans + "\"; gene_name \"COMPLETE\"; transcript_name \"" + trans + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_biotype \"" + trans_biotype + "\";\n")	
		ex = str(n+1)
		if str_gtf[trans].strand == "-":
			ex = str(len(str_gtf[trans].start) - n)
		lines[str_gtf[trans].gene.split("/")[0]][0].append(str_gtf[trans].chrm + "\tStringtie\texon\t" + str(str_gtf[trans].start[n]) + "\t" + str(str_gtf[trans].end[n]) + "\t.\t" + str_gtf[trans].strand + "\t.\tgene_id \"" + str_gtf[trans].gene.split("/")[0] + "\"; transcript_id \"" + trans + "\"; gene_name \"COMPLETE\"; transcript_name \"" + trans + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_biotype \"" + trans_biotype + "\"; exon_number \"" + ex + "\";\n")


for gene in lines:
	genename = []
	for line in lines[gene][1]:
		if not line in genename:	
			genename.append(line)
	genename = "-".join(genename)
	print(lines[gene][0][0].split("\t")[0] + "\tStringtie\tgene\t" + str(boundaries[gene][0]) + "\t" + str(boundaries[gene][1]) + "\t.\t" + lines[gene][0][0].split("\t")[6] + "\t.\tgene_id \"" + lines[gene][0][0].split('gene_id "')[1].split('"')[0] + "\"; gene_name \"" + genename + "\"; gene_biotype \"" + lines[gene][0][0].split('gene_biotype "')[1].split('"')[0] + "\";\n", end = "")
	for line in lines[gene][0]:
		print(line.replace("COMPLETE",genename),end="")



