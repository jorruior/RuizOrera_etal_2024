#! /usr/bin/python
__author__="jruizor"
__date__ ="$Jul 19, 2020 12:24:43 PM$"
'''Merge and quantify ORFs
'''

import sys
import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

orf_gtff = sys.argv[1]
trans_faf = sys.argv[2]
genome_faf = sys.argv[3]
cdsf = sys.argv[4]
gtff = sys.argv[5]
psites_cds = sys.argv[6]
try:
	expression = sys.argv[7] #stringtie.gtf
	expressiong = sys.argv[8] #abundance
except:
	expression = "none"
	expressiong = "none"
gencode = SeqIO.index("gencode.prot.fa", "fasta")
orf_gtff2 = orf_gtff.replace("_forprice","")

class trans_object:
	def __init__(self, chrm, gene, strand, start, end, desc):
		self.chrm = chrm
		self.gene = gene
		self.strand = strand
		self.start = start
		self.end = end
		self.desc = desc

def parse_gtf(gtf,field):
	#Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	trans_list = {}
	for line in open(gtf):
		if not "\t" + field + "\t" in line:
			continue
		g_name = line.split('gene_id "')[1].split('"')[0]
		
		if field == "CDS":
			t_name = line.split('ORF_id "')[1].split('"')[0]
			trans_list.setdefault(g_name,[])
			trans_list[g_name].append(t_name)
			trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[], line.split("\t")[8].rstrip("\n")))
		elif field == "transcript":
			t_name = line.split('transcript_id "')[1].split('"')[0]
			fpkm = float(line.split('FPKM "')[1].split('"')[0])
			trans_list.setdefault(g_name,[])
			trans_list[g_name].append([t_name,fpkm])
			trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[], str(fpkm)))
		else:
			t_name = line.split('transcript_id "')[1].split('"')[0] 
			trans_list.setdefault(g_name,[])
			trans_list[g_name].append(t_name)
			trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[], line.split("\t")[8].rstrip("\n")))
		
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans,trans_list

cds = SeqIO.index(cdsf, "fasta")
(trans_gtf,trans_list) = parse_gtf(gtff,"exon")
(orf_gtf,orf_list) = parse_gtf(orf_gtff,"CDS")
trans_fa = SeqIO.index(trans_faf, "fasta")

if expression != "none":
	exp_genes = {}
	(exp_gtf,exp_list) = parse_gtf(expression,"transcript")
	for line in open(expressiong):
		exp_genes[line.split("\t")[0]] = line.split("\t")[7]

os.system("grep -P \"\\tCDS\\t\" " + orf_gtff + " > " + orf_gtff2 + ".CDS.gtf")
os.system("sed -i \'s/transcript_id /transcript_name /\' " + orf_gtff2 + ".CDS.gtf")
os.system("sed -i \'s/ORF_id /transcript_id /\' " + orf_gtff2 + ".CDS.gtf")
os.system("gffread " + orf_gtff2 + ".CDS.gtf -g " + genome_faf + " -w " + orf_gtff2 + ".gffread.fa")
orfs_fa = SeqIO.index(orf_gtff2 + ".gffread.fa", "fasta")

new_names = {}
seqs = {}
for orf in orfs_fa:
	transcript = orf.split("_")[0]

	orf_seq = str(orfs_fa[orf].seq)
	if orf_seq[:3] == "ATG": #Remove if not near-cognate
		continue
	orf_seq = str(Seq(orf_seq).translate(cds=False))

	pf = orf_seq.find("*") 
	if (pf != -1) and (pf != len(orf_seq)-1): #Remove if internal stop codons were found
		continue

	if expression != "none":
		if not transcript in trans_gtf:
			c = []
			c.append([transcript,0])
		elif not trans_gtf[transcript].gene in exp_list:
			c = []
			c.append([transcript,0])
		else:
			c = []
			nc = []
			for elemento in exp_list[trans_gtf[transcript].gene]:
				c.append(elemento)

			c = sorted(c, key = lambda x: int(x[1]), reverse = True)

			if not transcript in c:
				c.append([transcript,0])
	else:
		c = []
		c.append([transcript,1])

	for t in c:
		if not t[0] in trans_fa:
			continue
		
		trans_seq = str(trans_fa[t[0]].seq)

		fa = str(Seq(trans_seq).translate(cds = False)).find(orf_seq)
		fb = str(Seq(trans_seq[1:]).translate(cds = False)).find(orf_seq)
		fc = str(Seq(trans_seq[2:]).translate(cds = False)).find(orf_seq)
		fi = [fa,fb,fc].index(max([fa,fb,fc]))

		if max(fa,fb,fc) == -1:
			continue

		f1 = max(fa,fb,fc)*3 + fi

		f2 = f1 + len(orf_seq)*3

		if f2 == f1:
			continue	

		if f2 > len(trans_seq):
			continue

		nucl = trans_seq[f1:f2]
		prot = str(Seq(nucl).translate(cds=False))

		#Check stop codons
		pf = prot.find("*")
		if (pf < len(prot)-1) or (pf == -1):
			continue

		name = t[0] + "_" + str(f1+1) + "_" + str(f2) + "_" + nucl[:3]
		new_names[orf] = name

		seqs[new_names[orf]] = (nucl,prot)
		break

psites = open(orf_gtff2 + ".psites.bed","w+")
tis = open(orf_gtff2 + ".TIS.bed","w+")
for orf in new_names:

	#Correct coordinates with stop codon
	all_pos = []
	for n,st in enumerate(trans_gtf[new_names[orf].split("_")[0]].start):
		for i in range(trans_gtf[new_names[orf].split("_")[0]].start[n],trans_gtf[new_names[orf].split("_")[0]].end[n]+1):
			all_pos.append(i)

	orf_gtf.setdefault(orf,trans_object(trans_gtf[new_names[orf].split("_")[0]].chrm,trans_gtf[new_names[orf].split("_")[0]].gene,trans_gtf[new_names[orf].split("_")[0]].strand,[],[], trans_gtf[new_names[orf].split("_")[0]].desc))

	if orf_gtf[orf].strand == "-":
		all_pos = all_pos[::-1]

	orf_gtf[orf].start = []
	orf_gtf[orf].end = []
	op = 0
	op2 = 0
	stop = 0
	stop2 = 0
	stop3 = 0
	ok = 0
	posp = 0
	un = int(new_names[orf].split("_")[1])-30
	un2 = int(new_names[orf].split("_")[1])+30+2	
	if un < 0:
		un = 0
	for n,pos in enumerate(all_pos):

		#TIS
		if n == int(new_names[orf].split("_")[1])-1:
			ok = 1
		if n == un:
			stop = 1
		elif n == un2:
			if orf_gtf[orf].strand == "+":
				tis.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tTIS\t" + orf_gtf[orf].strand + "\n")	
			stop = 0

		if stop == 1:
			if (orf_gtf[orf].strand == "-") and (stop2 % 3 == 2):
				if ok == 1:
					tis.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tTISa\t" + orf_gtf[orf].strand + "\n")
					ok = 0
				else:
					tis.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tTIS\t" + orf_gtf[orf].strand + "\n")
				stop3 += 1
			elif (orf_gtf[orf].strand == "+") and (stop2 % 3 == 2):
				if ok == 1:
					tis.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tTISa\t" + orf_gtf[orf].strand + "\n")
					ok = 0
				else:				
					tis.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tTIS\t" + orf_gtf[orf].strand + "\n")
				stop3 += 1	
			stop2 +=1
			
		#PSITES
		if (op != 0) and (abs(pos - posp) > 1):
			orf_gtf[orf].start.append(pos)
			orf_gtf[orf].end.append(posp)

		if n == int(new_names[orf].split("_")[1])-1:
			orf_gtf[orf].start.append(pos)
			op = 1
		elif n ==int(new_names[orf].split("_")[2])-1:
			orf_gtf[orf].end.append(pos)
			op = 0
			if orf_gtf[orf].strand == "+":
				psites.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tp2\t" + orf_gtf[orf].strand + "\n")	

		if op == 1:
			if (orf_gtf[orf].strand == "-") and (op2 % 3 == 0):
				psites.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tp2\t" + orf_gtf[orf].strand + "\n")
			elif (orf_gtf[orf].strand == "+") and (op2 % 3 == 2):
				psites.write(orf_gtf[orf].chrm + "\t" + str(pos) + "\t" + str(pos) + "\t" + new_names[orf] + "\tp2\t" + orf_gtf[orf].strand + "\n")	
			op2 += 1



		posp = pos

	if orf_gtf[orf].strand == "-":
		lista = orf_gtf[orf].start
		orf_gtf[orf].start = orf_gtf[orf].end
		orf_gtf[orf].end = lista
psites.close()
tis.close()

cds_ov = []
os.system("intersectBed -s -wo -a " + orf_gtff2 + ".psites.bed -b " + psites_cds + " | awk \'{if ($2 == $8) print $0}\' > " + orf_gtff2 + ".psites_vs_cds.ov")
for line in open(orf_gtff2 + ".psites_vs_cds.ov"):
	if (line.split("\t")[5] == "+") and (line.split("\t")[10] == "p0"):
		cds_ov.append(line.split("\t")[3])
	elif (line.split("\t")[5] == "-") and (line.split("\t")[10] == "p2"):
		cds_ov.append(line.split("\t")[3])	
cds_ov = list(set(cds_ov))

os.system("coverageBed -s -a " + orf_gtff2 + ".psites.bed -b " + orf_gtff2.replace(".price.gtf","_ORFquant_Detected_ORFs.gtf") + ".psites.bed > " + orf_gtff2 + "_PSITES.ov")
ovs = {}
for line in open(orf_gtff2 + "_PSITES.ov"):
	orf = line.split("\t")[3]
	c = float(line.split("\t")[-1].rstrip("\n"))
	if not orf in ovs:
		ovs[orf] = [0,0]
	ovs[orf][0] += 1
	if c >= 1:
		ovs[orf][1] += 1

#Writing outputs
outg = open(orf_gtff2.replace(".price.gtf","") + "_plus_NTG.gtf","w+")
outgn = open(orf_gtff2.replace(".price.gtf","") + "_plus_NTG.nucl.fa","w+")
outgp = open(orf_gtff2.replace(".price.gtf","") + "_plus_NTG.prot.fa","w+")
for line in open(orf_gtff2.replace(".price.gtf","_ORFquant_Detected_ORFs.gtf") + ".fixed.gtf"):
	outg.write(line)
for line in open(orf_gtff2.replace(".price.gtf","_ORFquant_Detected_ORFs.gtf") + ".fixed.nucl.fa"):
	outgn.write(line)
for line in open(orf_gtff2.replace(".price.gtf","_ORFquant_Detected_ORFs.gtf") + ".fixed.prot.fa"):
	outgp.write(line)

outov = open(orf_gtff2.replace(".price.gtf","") + "_NTG_ovs.gtf","w+")
for orf in new_names:
	if new_names[orf] in cds_ov:
		continue

	#ORF biotype
	o = seqs[new_names[orf]][1]
	t = trans_fa[new_names[orf].split("_")[0]].seq

	desc = trans_gtf[new_names[orf].split("_")[0]].desc.split('exon_number "')[0] + ' orf_id "' + new_names[orf]

	genes = []
	for g in orf_gtf[orf].gene.split("_"):
		genes.append(g)
	if new_names[orf].split("_")[0] in cds:
		c = str(cds[new_names[orf].split("_")[0]].seq).replace("*","")
		clase = "mRNA"
		if o == c:
			clase = "CDS"
		elif o in c:
			clase = "CDSa"
		elif c in o:
			clase = "CDSa"
		else:
			f1 = str(t.translate(cds = False)).find(str(c))
			f2 = str(t[1:].translate(cds = False)).find(str(c))
			f3 = str(t[2:].translate(cds = False)).find(str(c))
			fi = [f1,f2,f3].index(max([f1,f2,f3]))
			fs = max(f1,f2,f3)*3 + fi
			fe = fs + len(c)*3
			os = int(new_names[orf].split("_")[1])-1
			oe = int(new_names[orf].split("_")[2])
			if os < fs:
				if oe < fs:
					clase = "uORF"
				else:
					clase = "uoORF"
			elif oe > fe:
				if os > fe:
					clase = "dORF"
				else:
					clase = "doORF"
			else:
				clase = "intORF"

	else:
		if 'gene_biotype "protein_coding"' in desc:
			clase = "ncRNA-ORF"
			for g in genes:
				if g in trans_list:
					for cdsname in trans_list[g]:
						if cdsname in cds:
							c = str(cds[cdsname].seq).replace("*","")	
							if o == c:
								clase = "CDS"
								break
							elif o in c:
								clase = "CDSa"
								break
							elif c in o:
								clase = "CDSa"
							break

		elif "pseudo" in desc:
			clase = "pseudogene"
		else:
			clase = "lncRNA-ORF"

	#Gencode status
	gencode_status = "none"
	for gen in gencode:
		if "P1" in gen:
			if "rand" in gen:
				continue
			if (str(gencode[gen].seq).replace("*","") in o) or (o in str(gencode[gen].seq).replace("*","")):
				gencode_status = gen
		elif "ENST" in gen:
			if "rand" in gen:
				continue
			if str(gencode[gen].seq).replace("*","") == 0:
				gencode_status = gen	

	#Write output
	if ovs[new_names[orf]][1] > 0:
		if expression != "none":
			for n,st in enumerate(orf_gtf[orf].start):
				if not orf.split("_")[0] in exp_gtf:
					fpkm = "0"
				else:
					fpkm = exp_gtf[new_names[orf].split("_")[0]].desc

				fpkmg = "0"
				main = g[0]
				for g in genes:
					if g in exp_genes:
						if exp_genes[g] > fpkmg:
							fpkmg = exp_genes[g]
							main = g
				desc = desc.replace(orf_gtf[orf].gene,g) #Main gene is the one with the highest expression
				if gencode_status != "none":
					outov.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; FPKM "' + fpkm + '"; FPKMg "' + fpkmg + '"; orf_name "' + gencode_status + '"; psites_ov "' + str(ovs[new_names[orf]][1]/ovs[new_names[orf]][0]*100) + '";\n')
				else:
					outov.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; FPKM "' + fpkm + '"; FPKMg "' + fpkmg + '"; psites_ov "' + str(ovs[new_names[orf]][1]/ovs[new_names[orf]][0]*100) + '";\n')
		else:
			for n,st in enumerate(orf_gtf[orf].start):
				if gencode_status != "none":
					outov.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; orf_name "' + gencode_status + '"; psites_ov "' + str(ovs[new_names[orf]][1]/ovs[new_names[orf]][0]*100) + '";\n')
				else:
					outov.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; psites_ov "' + str(ovs[new_names[orf]][1]/ovs[new_names[orf]][0]*100) + '";\n')
	else:
		outgn.write(">" + new_names[orf] + " " + orf + "\n" + seqs[new_names[orf]][0] + "\n")
		outgp.write(">" + new_names[orf] + " " + orf + "\n" + seqs[new_names[orf]][1] + "\n")
		if expression != "none":
			for n,st in enumerate(orf_gtf[orf].start):
				if not orf.split("_")[0] in exp_gtf:
					fpkm = "0"
				else:
					fpkm = exp_gtf[new_names[orf].split("_")[0]].desc
				
				fpkmg = "0"
				main = g[0]
				for g in genes:
					if g in exp_genes:
						if exp_genes[g] > fpkmg:
							fpkmg = exp_genes[g]
							main = g
				desc = desc.replace(orf_gtf[orf].gene,g) #Main gene is the one with the highest expression
				if gencode_status != "none":
					outg.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; FPKM "' + fpkm + '"; FPKMg "' + fpkmg + '"; orf_name "' + gencode_status + '";\n')
				else:
					outg.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; FPKM "' + fpkm + '"; FPKMg "' + fpkmg + '";\n')
		else:
			for n,st in enumerate(orf_gtf[orf].start):
				if gencode_status != "none":
					outg.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '"; orf_name "' + gencode_status + '";\n')
				else:
					outg.write(orf_gtf[orf].chrm + "\tORFquant\tCDS\t" + str(orf_gtf[orf].start[n]) + "\t" + str(orf_gtf[orf].end[n]) + "\t.\t" + orf_gtf[orf].strand + "\t.\t" + desc + '"; orf_biotype "' + clase + '";\n')

outg.close()
outgn.close()
outgp.close()
outov.close()
exit(0)
