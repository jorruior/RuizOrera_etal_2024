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
from Bio.SeqRecord import SeqRecord

orf_gtff = sys.argv[1]
orf_faf = sys.argv[2]
trans_faf = sys.argv[3]
genome_faf = sys.argv[4]
cdsf = sys.argv[5] #FASTA with annotated CDS
gtff = sys.argv[6] #GTF annotaated and adapted to ORFquant
psites_cds = sys.argv[7] #BED with annotated psites
try:
	expression = sys.argv[8] #stringtie.gtf
	expressiong = sys.argv[9] #abundance
except:
	expression = "none"
	expressiong = "none"
gencode = SeqIO.index("gencode.prot.fa", "fasta")

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
			trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[], line.split("\t")[8].split(" exon_number")[0].rstrip("\n")))
		
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans,trans_list


cds = SeqIO.index(cdsf, "fasta")
cds2 = {}
(trans_gtf,trans_list) = parse_gtf(gtff,"exon")
(orf_gtf,orf_list) = parse_gtf(orf_gtff,"CDS")
trans_fa = SeqIO.index(trans_faf, "fasta")

if expression != "none":
	exp_genes = {}
	(exp_gtf,exp_list) = parse_gtf(expression,"transcript")
	for line in open(expressiong):
		exp_genes[line.split("\t")[0]] = line.split("\t")[7]


#Generate an adapted formated GTF with the CDSs
os.system("grep -P \"\\tCDS\\t\" " + orf_gtff + " > " + orf_gtff + ".CDS.gtf")
os.system("sed -i \'s/transcript_id /transcript_name /\' " + orf_gtff + ".CDS.gtf")
os.system("sed -i \'s/ORF_id /transcript_id /\' " + orf_gtff + ".CDS.gtf")

orfs_fa = SeqIO.index(orf_faf, "fasta")

#Check all the transcripts in the same gene and select the one with the highest expression that overlaps the gene.

outn = open(orf_gtff + ".fixed.nucl.fa","w+")
outp = open(orf_gtff + ".fixed.prot.fa","w+")
new_names = {}
excluded = 0

c = {}
for orf in orfs_fa:
	if "readthrough" in orf:
		continue
	transcript = orf.split("|")[0].split("_")[0]
	done = 0

	#Define sequence
	orf_seq = str(orfs_fa[orf].seq)

	if not transcript in trans_gtf:
		continue

	#If expression != none (default of 1), check all transcripts in the same locus and choose the one with the highest expression
	if expression != "none":
		if not trans_gtf[transcript].gene in exp_list:
			c[orf.split("|")[0]] = []
			c[orf.split("|")[0]].append([transcript,0])
		else:
			c[orf.split("|")[0]] = []
			for elemento in exp_list[trans_gtf[transcript].gene]: #For transcript in the host gene
				c[orf.split("|")[0]].append(elemento)

			c[orf.split("|")[0]] = sorted(c[orf.split("|")[0]], key = lambda x: int(x[1]), reverse = True) #sort list of transcripts by expression

			if not transcript in c[orf.split("|")[0]]: #Add expression of 0 in not in the list
				c[orf.split("|")[0]].append([transcript,0])
	else:
		c[orf.split("|")[0]] = []
		c[orf.split("|")[0]].append([transcript,1])

	for t in c[orf.split("|")[0]]: #t[0] is name, t[1] is expression
		if not t[0] in trans_fa:
			continue
		
		trans_seq = str(trans_fa[t[0]].seq) #Transcript sequence

		#Translate three frames and find the ORF
		fa = str(Seq(trans_seq).translate(cds = False)).find(orf_seq)
		fb = str(Seq(trans_seq[1:]).translate(cds = False)).find(orf_seq)
		fc = str(Seq(trans_seq[2:]).translate(cds = False)).find(orf_seq)
		fi = [fa,fb,fc].index(max([fa,fb,fc]))

		if max(fa,fb,fc) == -1: #The ORF is not present in the transcript
			continue

		f1 = max(fa,fb,fc)*3 + fi

		f2 = f1 + len(orf_seq)*3 + 3
		if f2 == f1:
			continue	

		if f2 > len(trans_seq): #If the length is not compatible with the transcript
			continue

		nucl = trans_seq[f1:f2]
		prot = str(Seq(nucl).translate(cds=False))

		#Check stop codons (some ORFquant structures do not contain stop codons or multiple ones, maybe fixed in following versions?)
		pf = prot.find("*")
		if (pf < len(prot)-1) or (pf == -1):
			continue

		name = t[0] + "_" + str(f1+1) + "_" + str(f2) #New name for the transcript

		new_names[orf.split("|")[0]] = name
		outn.write(">" + name + " " + orf.split("|")[0] + "\n" + nucl + "\n")
		outp.write(">" + name + " " + orf.split("|")[0] + "\n" + prot + "\n")

		#Check if the sequence is coding (partial or total)
		for cdsname in trans_list[trans_gtf[transcript].gene]:
			if cdsname in cds:
				cd = str(cds[cdsname].seq).replace("*","")
				if prot.replace("*","") == cd:
					cds2[transcript] = prot
				elif prot.replace("*","") in cd:
					cds2[transcript] = prot
				elif cd in prot.replace("*",""):
					cds2[transcript] = prot

		done = 1
		break

	if done == 0:
		excluded += 1
outn.close()
outp.close()

print(str(excluded) + " cases not considered")

#Append transcripts to the new GTF file, add FPKM information Stringtie
outg = open(orf_gtff + ".fixed.gtf","w+")
for trans in trans_gtf:
	if trans in cds2:
		trans_gtf[trans].desc = trans_gtf[trans].desc.replace('transcript_biotype "processed_transcript"','transcript_biotype "protein_coding"')
	for n,st in enumerate(trans_gtf[trans].start):
		if expression != "none":
			if not trans in exp_gtf:
				if trans_gtf[trans].gene in exp_genes:
					outg.write(trans_gtf[trans].chrm + "\tORFquant\texon\t" + str(trans_gtf[trans].start[n]) + "\t" + str(trans_gtf[trans].end[n]) + "\t.\t" + trans_gtf[trans].strand + "\t.\t" + trans_gtf[trans].desc + ' FPKM "0"; FPKMg "' + exp_genes[trans_gtf[trans].gene] + '";\n')
				else:
					outg.write(trans_gtf[trans].chrm + "\tORFquant\texon\t" + str(trans_gtf[trans].start[n]) + "\t" + str(trans_gtf[trans].end[n]) + "\t.\t" + trans_gtf[trans].strand + "\t.\t" + trans_gtf[trans].desc + ' FPKM "0"; FPKMg "0";\n')
			else:
				outg.write(trans_gtf[trans].chrm + "\tORFquant\texon\t" + str(trans_gtf[trans].start[n]) + "\t" + str(trans_gtf[trans].end[n]) + "\t.\t" + trans_gtf[trans].strand + "\t.\t" + trans_gtf[trans].desc + ' FPKM "' + exp_gtf[trans].desc + '"; FPKMg "' + exp_genes[trans_gtf[trans].gene] + '";\n')
		else:
			outg.write(trans_gtf[trans].chrm + "\tORFquant\texon\t" + str(trans_gtf[trans].start[n]) + "\t" + str(trans_gtf[trans].end[n]) + "\t.\t" + trans_gtf[trans].strand + "\t.\t" + trans_gtf[trans].desc + '\n')

#Write p-sites in a bed file
psites = open(orf_gtff + ".psites.bed","w+")
tis = open(orf_gtff + ".TIS.bed","w+")
for orft in orfs_fa:
	if "readthrough" in orft:
		continue

	orf = orft.split("|")[0]
	if not orf in new_names:
		print(orf + " not considered")
		continue

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

#Check cases that overlap with CDSs in-frame
cds_ov = []
os.system("intersectBed -s -wo -a " + orf_gtff + ".psites.bed -b " + psites_cds + " | awk \'{if ($2 == $8) print $0}\' > " + orf_gtff + ".psites_vs_cds.ov")
for line in open(orf_gtff + ".psites_vs_cds.ov"):
	if (line.split("\t")[5] == "+") and (line.split("\t")[10] == "p0"):
		cds_ov.append(line.split("\t")[3])
	elif (line.split("\t")[5] == "-") and (line.split("\t")[10] == "p2"):
		cds_ov.append(line.split("\t")[3])	
cds_ov = list(set(cds_ov))

#Define biotypes
tc = [] #Transcripts with protein-coding CDSs
ncrna_orfs = []
for orft in orfs_fa:
	if "readthrough" in orft:
		continue
	orf = orft.split("|")[0]
	if not orf in new_names:
		print(orf + " not considered")
		continue
	#ORF biotype
	o = str(orfs_fa[orft].seq).replace("*","")
	t = trans_fa[orf.split("_")[0]].seq

	desc = trans_gtf[new_names[orf].split("_")[0]].desc + ' orf_id "' + new_names[orf]

	if orf.split("_")[0] in cds2:
		desc = desc.replace('transcript_biotype "processed_transcript"','transcript_biotype "protein_coding"')

	gene = desc.split('gene_id "')[1].split('"')[0]	
	clase = "ncRNA-ORF"	
	if orf.split("_")[0] in cds:
		cs = str(cds[orf.split("_")[0]].seq).replace("*","")
		if o == cs:
			clase = "CDS"
		elif o in cs:
			clase = "CDSa"
		elif cs in o:
			clase = "CDSa"
	elif orf.split("_")[0] in cds2:
		cs = cds2[orf.split("_")[0]].replace("*","")
		if o == cs:
			clase = "CDSa"
		elif o in cs:
			clase = "CDSa"
		elif cs in o:
			clase = "CDSa"	

	if not "CDS" in clase:
		if (new_names[orf] in cds_ov) or (orf in cds_ov):
			clase = "CDSo"
		elif (new_names[orf].split("_")[0] in cds or new_names[orf].split("_")[0] in cds2):
			try:
				cs = str(cds[new_names[orf].split("_")[0]].seq).replace("*","")
			except:
				cs = cds2[new_names[orf].split("_")[0]].replace("*","")
			f1 = str(t.translate(cds = False)).find(str(cs))
			f2 = str(t[1:].translate(cds = False)).find(str(cs))
			f3 = str(t[2:].translate(cds = False)).find(str(cs))
			fi = [f1,f2,f3].index(max([f1,f2,f3]))
			fs = max(f1,f2,f3)*3 + fi
			fe = fs + len(cs)*3
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

	if clase == "ncRNA-ORF":
		if 'gene_biotype "protein_coding"' in desc:
			for cdsname in trans_list[gene]:
				if cdsname in cds:
					cs = str(cds[cdsname].seq).replace("*","")	
					if o == cs:
						clase = "CDS"
						break
					elif o in cs:
						clase = "CDSa"
						break
					elif cs in o:
						clase = "CDSa"
						break

		elif "pseudo" in desc:
			clase = "pseudogene"
		else:
			clase = "lncRNA-ORF"

	if "CDS" in clase:
		tc.append(orf.split("_")[0])
	elif clase == "ncRNA-ORF":
		ncrna_orfs.append(orf + "/" + o)

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
	if expression != "none":
		for n,st in enumerate(orf_gtf[orf].start):
			if not orf.split("_")[0] in exp_gtf:
				fpkm = "0"
			else:
				fpkm = exp_gtf[new_names[orf].split("_")[0]].desc
			if not gene in exp_genes:
				fpkmg = "0"
			else:
				fpkmg = exp_genes[gene]
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

#Check if ncRNA-ORFs might be uORFs instead
ncuorfs = open(orf_gtff + ".ncRNA_outORFs.tsv","w+")
for orf in ncrna_orfs:
	orf_seq = orf.split("/")[1]
	if not orf.split("/")[0] in c:
		continue
	for t in c[orf.split("/")[0]]: #t[0] is name, t[1] is expression
		if not t[0] in trans_fa:
			continue
		
		trans_seq = str(trans_fa[t[0]].seq) #Transcript sequence

		#Translate three frames and find the ORF
		fa = str(Seq(trans_seq).translate(cds = False)).find(orf_seq)
		fb = str(Seq(trans_seq[1:]).translate(cds = False)).find(orf_seq)
		fc = str(Seq(trans_seq[2:]).translate(cds = False)).find(orf_seq)
		fi = [fa,fb,fc].index(max([fa,fb,fc]))

		if max(fa,fb,fc) == -1: #The ORF is not present in the transcript
			continue

		f1 = max(fa,fb,fc)*3 + fi

		f2 = f1 + len(orf_seq)*3 + 3
		if f2 == f1:
			continue	

		if f2 > len(trans_seq): #If the length is not compatible with the transcript
			continue

		nucl = trans_seq[f1:f2]
		prot = str(Seq(nucl).translate(cds=False))

		#Check stop codons (some ORFquant structures do not contain stop codons or multiple ones, maybe fixed in following versions?)
		pf = prot.find("*")
		if (pf < len(prot)-1) or (pf == -1):
			continue

		if t[0] in tc:
			ncuorfs.write(orf.split("/")[0] + "\tpossible ouuORF\t" + t[0] + "\tFPKM:" + str(t[1]) + "\n")
			break
ncuorfs.close()

exit(0)
