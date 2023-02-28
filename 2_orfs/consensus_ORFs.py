#! /usr/bin/python
__author__="jruizor"
__date__ ="$Dec 19, 2020 12:24:43 PM$"
'''Merge and quantify ORFs
'''

import sys
import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

mergedgtf = sys.argv[1] #GTF with all merged ORFs (pooled reads or ORFs)
gtfs = sys.argv[2] #List with all individual BAM files, add class to second tab
orfsfile = sys.argv[3] #ORFs fasta
psitesbed = sys.argv[4]
psitesbedntg = sys.argv[5]
ensembl_gtf = sys.argv[6]
thr1 = int(sys.argv[7]) #Min samples with full ORF
thr2 = 15 #Min number of in-frame psites
thr3 = float(sys.argv[8]) #Min gene FPKM
tag = sys.argv[9]
expr = sys.argv[10] #yes,no,only_collapse
try:
	p_code = sys.argv[11]
except:
	p_code = "none"

class trans_object:
	def __init__(self, chrm, gene, strand, start, end, biotype,psites,desc,ns):
		self.chrm = chrm
		self.gene = gene
		self.strand = strand
		self.start = start
		self.end = end
		self.biotype = biotype
		self.psites = psites
		self.desc = desc
		self.ns = ns

def parse_gtf(gtf,field,thr,expr):
	#Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene
	c = []
	nv = []
	nc = []
	ps = []
	trans = {}
	for line in open(gtf):
		if 'gene_biotype' in line:
			biot = line.split('gene_biotype "')[1].split('"')[0]
		else:
			biot = "unknown"

		if not field in line:
			if thr != "none":
				if "\texon\t" in line:
					if expr == "yes":
						f = float(line.split('FPKMg "')[1].split('"')[0])
					else:
						f = thr
					if biot == "protein_coding":
						if f >= thr:
							c.append(line.split('gene_id "')[1].split('"')[0])
					elif "pseudogene" in biot:
						if f >= thr:
							ps.append(line.split('gene_id "')[1].split('"')[0])
					elif biot == "processed_transcript":
						if f >= thr:
							nv.append(line.split('gene_id "')[1].split('"')[0])	
					else:						
						if f >= thr:
							nc.append(line.split('gene_id "')[1].split('"')[0])					
			continue
		
		t_name = line.split('orf_id "')[1].split('"')[0]		
		g_name = line.split('gene_id "')[1].split('"')[0]
		desc = line.split("\t")[8]
	
		trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[],biot,0,desc,0))
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	if thr != "none":
		print("Expressed protein-coding genes:\t" + str(len(list(set(c)))))
		print("Expressed non-coding genes:\t" + str(len(list(set(nc)))))
		print("Expressed novel genes:\t" + str(len(list(set(nv)))))
		print("Expressed pseudogenes:\t" + str(len(list(set(ps)))))
	return trans


print("Generate collapsed fasta and bed for GENCODE script")
collapsed1 = open(mergedgtf.replace(".gtf",".collapsed.bed").replace(".bam.collapsed",".collapsed.bed.collapsed"),"w+")
collapsed2 = open(mergedgtf.replace(".gtf",".collapsed.fa").replace(".bam.collapsed",".collapsed.fa.collapsed"),"w+")
biotypes = {}
for gtf in open(gtfs):
	done = []
	sample = gtf.split("\t")[0].split("/")[-1].replace("_trimmed.fq.gzAligned.sortedByCoord.out.bam","")
	for line in open(gtf.split("\t")[0] + "_plus_NTG.gtf"):
		if "\tCDS\t" in line:
			if not line in done:
				collapsed1.write(line.split("\t")[0] + "\t" + line.split("\t")[3] + "\t" + line.split("\t")[4] + "\t" + line.split('orf_id "')[1].split('"')[0] + "\t" + tag + "_" + sample + "\t" + line.split("\t")[6] + "\n")
				biotypes[line.split('orf_id "')[1].split('"')[0]] = line.split('orf_biotype "')[1].split('"')[0]
				done.append(line)
	exc = 0
	for line in open(gtf.split("\t")[0] + "_plus_NTG.prot.fa"):
		if ">" in line:
			if line.split(" ")[0] in done:
				exc = 1
			else:
				exc = 0
				collapsed2.write(line.split(" ")[0] + "--" + tag + "_" + sample + " " + biotypes[line.split(" ")[0].replace(">","")] + "\n")
				done.append(line.split(" ")[0])
		else:
			if exc == 0:
				collapsed2.write(line)
collapsed1.close()
collapsed2.close()
if expr == "only_collapse":
	exit(0)

print("Preparing files")
merged = parse_gtf(mergedgtf,"orf_id",thr3,expr)
orfsf = {}
for line in open(orfsfile):
	if ">" in line:
		name = line.split(" ")[0].replace(">","").rstrip("\n")
		orfsf[name] = ""

	else:
		orfsf[name] = orfsf[name] + line.rstrip("\n")

os.system('grep -P "\texon\t" ' + ensembl_gtf + ' > ' + ensembl_gtf + '.exons.gtf')
os.system("coverageBed -s -a " + mergedgtf + " -b " + ensembl_gtf + ".exons.gtf > tmp/" + mergedgtf.split("/")[-1] + ".ensembl.ov")
ens_ov = {}
for line in open("tmp/" + mergedgtf.split("/")[-1] + ".ensembl.ov"):
	if "\texon\t" in line:
		continue
	orf = line.split('orf_id "')[1].split('"')[0]
	if not orf in ens_ov:
		ens_ov[orf] = [0,0]
	ens_ov[orf][0] = ens_ov[orf][0] + int(line.split("\t")[-3])
	ens_ov[orf][1] = ens_ov[orf][1] + int(line.split("\t")[-2])

ovs = {}
samples = []
clases = []
clasesu = []
for gtf in open(gtfs):
	#print(gtf)
	sample = gtf.split("\t")[0].split("/")[-1] + "_plus_NTG.gtf"
	ps1 = gtf.split("\t")[0] + "_ORFquant_Detected_ORFs.gtf.psites.bed"
	ps2 = gtf.split("\t")[0] + ".price.gtf.psites.bed"
	psites1 = gtf.split("\t")[0] + "_P_sites_uniq_plus.bedgraph"
	psites2 = gtf.split("\t")[0] + "_P_sites_uniq_minus.bedgraph"
	if not sample in samples:
		samples.append(sample)
	clases.append(gtf.split("\t")[1].rstrip("\n"))
	if not gtf.split("\t")[1].rstrip("\n") in clasesu:
		clasesu.append(gtf.split("\t")[1].rstrip("\n"))
	file = gtf.split("\t")[0] + "_plus_NTG.gtf"
	ind = parse_gtf(file,"orf_id","none",expr)
	done = []
	os.system("intersectBed -s -wao -a " + psitesbed + " -b " + ps1 + " > tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".ov")
	os.system("intersectBed -s -wao -a " + psitesbed + " -b " + ps2 + " >> tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".ov")
	os.system("intersectBed -wo -b " + psites1 + " -a " + psitesbed + " > tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.plus.ov")
	os.system("intersectBed -wo -b " + psites2 + " -a " + psitesbed + " > tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.minus.ov")
	if psitesbedntg != "none":
		os.system("intersectBed -s -wao -a " + psitesbedntg + " -b " + ps1 + " >> tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".ov")
		os.system("intersectBed -s -wao -a " + psitesbedntg + " -b " + ps2 + " >> tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".ov")
		os.system("intersectBed -wo -b " + psites1 + " -a " + psitesbedntg + " >> tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.plus.ov")
		os.system("intersectBed -wo -b " + psites2 + " -a " + psitesbedntg + " >> tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.minus.ov")
	for line in open("tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".ov"):
		orf1 = line.split('\t')[3]
		if orf1 in done:
			continue
		len1 = sum(merged[orf1].end) - sum(merged[orf1].start)
		if not orf1 in ovs:
			ovs[orf1] = {}
		orf2 = line.split('\t')[9]
		if not orf2 in ind:
			continue
		if p_code == "p0":
			if line.split('\t')[5] == "+":
				if int(line.split('\t')[1]) != int(line.split('\t')[7]):
					continue
			elif line.split('\t')[5] == "-":
				if (int(line.split('\t')[1])-1) != int(line.split('\t')[7]):
					continue	
		elif p_code == "none":
			if (line.split('\t')[4] != "p2") and (line.split('\t')[4] != "hq"):
				continue
			if int(line.split("\t")[1]) != int(line.split("\t")[7]):
				continue
		if orf1 == orf2:
			ovs[orf1][orf2 + "--" + sample] = [len1,len1,len1,orf1,orf2,ind[orf2].psites]
			done.append(orf1)
		else:
			len2 = sum(ind[orf2].end) - sum(ind[orf2].start)
			if not orf2 + "--" + sample in ovs[orf1]:
				ovs[orf1][orf2 + "--" + sample] = [3 , len1 , len2 , orf1 , orf2 , ind[orf2].psites]
			else:
				ovs[orf1][orf2 + "--" + sample][0] = ovs[orf1][orf2 + "--" + sample][0] + 3

out = {}
for orf1 in ovs:
	out[orf1] = []
	for sample in samples:
		orf22 = ["none",0,0,len(orf1)]
		for orf2 in ovs[orf1]:
			if sample in orf2:
				if float(ovs[orf1][orf2][0])/float(ovs[orf1][orf2][1]) > orf22[1]:
					orf22[1] = float(ovs[orf1][orf2][0])/float(ovs[orf1][orf2][1])
					orf22[0] = orf2.split("--")[0]
					orf22[2] = float(ovs[orf1][orf2][5])

		out[orf1].append(orf22)


print("Extracting P-sites - Unique reads")
ps = {}
for gtf in open(gtfs):
	sample = "tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.plus.ov"
	for line in open(sample):
		if not "\t+\t" in line:
			continue
		if p_code == "none":
			if not "p2" in line and not "hq" in line:
				continue
			if int(line.split("\t")[1]) != int(line.split("\t")[7]):
				continue
		else:
			if not p_code in line:
				continue
			if int(line.split("\t")[1]) != int(line.split("\t")[7]):
				continue			
		orf = line.split("\t")[3]
		if not orf in merged:
			continue
		if not orf in ps:
			ps[orf] = {}
		if not gtf.split("\t")[0] in ps[orf]:
			ps[orf][gtf.split("\t")[0]] = 0
		ps[orf][gtf.split("\t")[0]] = ps[orf][gtf.split("\t")[0]] + int(line.split("\t")[-2].rstrip("\n"))
		merged[orf].psites = merged[orf].psites + int(line.split("\t")[-2].rstrip("\n"))

for gtf in open(gtfs):
	sample = "tmp/" + mergedgtf.split("/")[-1] + "_" +  gtf.split("\t")[0].split("/")[-1] + ".psites.minus.ov"
	for line in open(sample):
		if not "\t-\t" in line:
			continue
		if p_code == "none":
			if not "p2" in line and not "hq" in line:
				continue
			if int(line.split("\t")[1]) != int(line.split("\t")[8]):
				continue
		else:
			if not p_code in line:
				continue
			if int(line.split("\t")[1]) != int(line.split("\t")[7]):
				continue
		orf = line.split("\t")[3]
		if not orf in merged:
			continue
		if not orf in ps:
			ps[orf] = {}
		if not gtf.split("\t")[0] in ps[orf]:
			ps[orf][gtf.split("\t")[0]] = 0
		ps[orf][gtf.split("\t")[0]] = ps[orf][gtf.split("\t")[0]] + int(line.split("\t")[-2].rstrip("\n"))
		merged[orf].psites = merged[orf].psites + int(line.split("\t")[-2].rstrip("\n"))



print("identifying candidate ORFs")
counts = open(mergedgtf + "_orf_replicability.txt","w+")
counts.write("orf\tclass\ttrans\tgene\tgene_bio\tT\tP\tN\tpsites\tFPKMg\talt_orfs\tlen")
hq = []
for clase in clasesu:
	counts.write("\tT--" + clase + "\tP--" + clase + "\tN--" + clase)
counts.write("\n")
for orf1 in out:
	t = {}
	p = {}
	nn = {}
	for clase in clasesu:
		t[clase] = 0
		p[clase] = 0
		nn[clase] = 0
	t["total"] = 0
	p["total"] = 0
	nn["total"] = 0
	mps = merged[orf1].psites
	names = []
	for n,elemento in enumerate(out[orf1]):
		if out[orf1][n][1] >= 0.9:
			t["total"] += 1
			t[clases[n]] += 1
			names.append(out[orf1][n][0])
		elif out[orf1][n][1] > 0:
			p["total"] += 1
			p[clases[n]] += 1
			names.append(out[orf1][n][0])
		else:
			nn[clases[n]] += 1
			nn["total"] += 1
		names = list(set(names))
		#print(orf1 + "\t" + samples[n] + "\t" + str(out[orf1][n]) + "\n", end = "")
	o = len(orfsf[orf1].replace("*",""))
	try:
		clase = merged[orf1].desc.split('orf_biotype "')[1].split('"')[0]
	except:
		clase = "unknown"
	try:
		tr = merged[orf1].desc.split('transcript_id "')[1].split('"')[0]
	except:
		tr = "unknown"
	try:
		gene_bio = merged[orf1].desc.split('gene_biotype "')[1].split('"')[0]
	except:
		gene_bio = "unknown"
	merged[orf1].ns = t["total"]
	if 'FPKMg "' in merged[orf1].desc:
		eg = merged[orf1].desc.split('FPKMg "')[1].split('"')[0]
	else:
		eg = str(thr3)
	counts.write(orf1 + "\t" + clase + "\t" + tr + "\t" + merged[orf1].gene + "\t" + gene_bio + "\t" + str(t["total"]) + "\t" + str(p["total"]) + "\t" + str(nn["total"]) + "\t" + str(mps) + "\t" + eg + "\t" + ";".join(names) + "\t" + str(o))
	if (t["total"] >= thr1) and (mps >= thr2) and (float(eg) >= thr3):
		hq.append(orf1)
	for clase in clasesu:
		counts.write("\t" + str(t[clase]) + "\t" + str(p[clase]) + "\t" + str(nn[clase]))
	counts.write("\n")
counts.close()



ps_file = open(mergedgtf + "_Psites.bed","w+")
for line in open(psitesbed):
	orf = line.split("\t")[3]
	if orf in out:
		if orf in hq:
			ps_file.write(line.replace("\tp2\t","\thq\t"))
		else:
			ps_file.write(line)
if psitesbedntg != "none":
	for line in open(psitesbedntg):
		orf = line.split("\t")[3]
		if orf in out:
			if orf in hq:
				ps_file.write(line.replace("\tp2\t","\thq\t"))
			else:
				ps_file.write(line)
ps_file.close()



print("Quantifying candidate ORFs")
c = []
nc = []
nv = []
pseudo = []
quant = open(mergedgtf + "_orf_quant_Psites.txt","w+")
quant.write("orf\tclass\ttrans\tgene\tgene_bio\tpsites\tFPKMg\tgencode\tlen\tseq")
for gtf in open(gtfs):
	quant.write("\t" + gtf.split("\t")[0].split("/")[-1].split("_R1_")[0] + "--" + gtf.split("\t")[1].rstrip("\n"))
quant.write("\n")
for orf in out:
	o = orfsf[orf].replace("*","")
	gene =  merged[orf].gene

	all_psites = []
	for gtf in open(gtfs):
		if orf in ps:
			if gtf.split("\t")[0] in ps[orf]:
				all_psites.append(str(ps[orf][gtf.split("\t")[0]]))
			else:
				all_psites.append("0")
		else:
			all_psites.append("0")

	try:
		clase = merged[orf].desc.split('orf_biotype "')[1].split('"')[0]
	except:
		clase = "unknown"
	try:
		gene_bio = merged[orf].desc.split('gene_biotype "')[1].split('"')[0]
	except:
		gene_bio = "unknown"
	try:
		tr = merged[orf].desc.split('transcript_id "')[1].split('"')[0]
	except:
		tr = "unknown"
	if gene_bio == "protein_coding":
		c.append(gene)
	elif "pseudogene" in gene_bio:
		pseudo.append(gene)
	elif gene_bio == "processed_transcript":
		nv.append(gene)
	else:
		nc.append(gene)
	if "CDS" in clase:
		clase = "CDS"

	if ens_ov[orf][0] == 0:
		clase = clase + "--novel"
	elif ens_ov[orf][0] < ens_ov[orf][1]:
		clase = clase + "--partial"	

	if "orf_name" in merged[orf].desc:
		gencode = merged[orf].desc.split('orf_name "')[1].split('"')[0]
	else:
		gencode = "none"

	if "phaseI_id" in merged[orf].desc:
		gencode = merged[orf].desc.split('phaseI_id "')[1].split('"')[0]
	else:
		gencode = "none"	

	if 'FPKMg "' in merged[orf1].desc:
		eg = merged[orf].desc.split('FPKMg "')[1].split('"')[0]
	else:
		eg = str(thr3)

	quant.write(orf + "\t" + clase + "\t" + tr + "\t" + gene + "\t" + gene_bio + "\t" + str(merged[orf].psites) + "\t" + eg + "\t" + gencode + "\t" + str(len(o)) + "\t" + o + "\t" + "\t".join(all_psites) + "\n")

quant.close()

if thr3 != "none":
	print("Translated protein-coding genes:\t" + str(len(list(set(c)))))
	print("Translated non-coding genes:\t" + str(len(list(set(nc)))))
	print("Translated novel genes:\t" + str(len(list(set(nv)))))
	print("Translated pseudogenes:\t" + str(len(list(set(pseudo)))))

hq_gtf = open(mergedgtf.replace(".gtf",".hq.gtf"),"w+")
hq_bed = open(mergedgtf.replace(".gtf",".hq.bed"),"w+")
hq_fa = open(orfsfile.replace(".fa",".hq.fa"),"w+")
for line in open(mergedgtf):
	if "\texon\t" in line:
		try:
			fpkm = float(line.split('FPKMg "')[1].split('"')[0])
		except:
			fpkm = thr3
		if fpkm >= thr3:
			hq_gtf.write(line)
	elif "\tCDS\t" in line:
		orf = line.split('orf_id "')[1].split('"')[0]
		if orf in hq:
			hq_gtf.write(line.rstrip("\n") + ' p-sites "' + str(merged[orf].psites) + '"; n_samples "' + str(merged[orf].ns) + '";\n')
			hq_bed.write(line.split("\t")[0] + "\t" + line.split("\t")[3] + "\t" + line.split("\t")[4] + "\t" + orf + "\t" + tag + "\t" + line.split("\t")[6] + "\n")

for orf in orfsf:
	if orf in hq:
		hq_fa.write(">" + orf + "--" + tag + "\n" + orfsf[orf] + "\n")

hq_gtf.close()
hq_bed.close()
hq_fa.close()
exit(0)


