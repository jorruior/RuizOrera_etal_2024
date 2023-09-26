#!/usr/bin/env python
__author__="jruizor"
__date__ ="$Jul 19, 2022 12:24:43 PM$"
'''Find the most likely homolog per each transcript
'''

##Activate conda -> Dependencies are pyliftover, bioPython, bedtools, gffread, BLAST

import sys
import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from pyliftover import LiftOver
import operator


gtf = sys.argv[1] #../annotation/custom/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf
annot = sys.argv[2] #../annotation/Homo_sapiens.GRCh38.98.gtf
fasta = sys.argv[3] #../annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.fa
sp1 = sys.argv[4] #human
make_blast = "no"
if sp1 == "human":
	species = ["chimp","gorilla","macaque","mouse","rat"] #the script is currently designed to work with these species
	species_short = ["hs","pt","gg","rm","mm","rn"]
elif sp1 == "chimp":
	species = ["human","macaque","mouse","rat"] #the script is currently designed to work with these species
	species_short = ["pt","hs","rm","mm","rn"]
elif sp1 == "macaque":
	species = ["human","chimp","mouse","rat"] #the script is currently designed to work with these species
	species_short = ["rm","hs","pt","mm","rn"]
elif sp1 == "gorilla":
	species = ["human"] #the script is currently designed to work with these species
	species_short = ["gg","hs"]	
elif sp1 == "mouse":
	species = ["human","chimp","macaque","rat"] #the script is currently designed to work with these species
	species_short = ["mm","hs","pt","rm","rn"]
elif sp1 == "rat":
	species = ["human","chimp","macaque","mouse"] #the script is currently designed to work with these species
	species_short = ["rn","hs","pt","rm","mm"]	

def get_coords(gtf,clase):
	'''Extract features from gtf'''
	coords_ex = {}
	# all_names = []
	for line in open(gtf):
		if clase == "exon":
			if '\t' + clase + '\t' in line:
				name = line.split('transcript_id "')[1].split('"')[0]
				if not name in coords_ex:
					coords_ex[name] = ["chr" + line.split("\t")[0],line.split("\t")[6],[],[],line.split('gene_name "')[1].split('"')[0],line.split('gene_biotype "')[1].split('"')[0],line.split('gene_id "')[1].split('"')[0]]
				coords_ex[name][2].append(line.split("\t")[3])
				coords_ex[name][3].append(line.split("\t")[4])
			elif '\tgene\t' in line:
				name = line.split('gene_id "')[1].split('"')[0]
				coords_ex[name + "_biotype"] = line.split('gene_biotype "')[1].split('"')[0]
				name = line.split('gene_name "')[1].split('"')[0]
				coords_ex[name + "_biotype"] = line.split('gene_biotype "')[1].split('"')[0]				

		elif clase == "gene":
			if '\t' + clase + '\t' in line:
				name = line.split('gene_id "')[1].split('"')[0]
				try:
					gname = line.split('gene_name "')[1].split('"')[0]
				except:
					gname = name
				if not name in coords_ex:
					coords_ex[name] = ["chr" + line.split("\t")[0],line.split("\t")[6],int(line.split("\t")[3]),int(line.split("\t")[4]),gname,line.split('gene_biotype "')[1].split('"')[0],line.split('gene_id "')[1].split('"')[0]]
			# elif '\ttranscript\t' in line:
			# 	name = line.split('gene_id "')[1].split('"')[0]
			# 	tname = line.split('transcript_id "')[1].split('"')[0]
			# 	coords_ex[t_name] = name

	return coords_ex


def load_chain(sp1,sp2):
	'''Load corresponding chain file'''
	eq = {}
	print("Running liftover")
	lo2 = "none"
	if sp1 == "human":
		if sp2 == "chimp":
			lo = LiftOver("chain/hg38ToPanTro5.over.chain") 
		elif sp2 == "gorilla":
			lo = LiftOver("chain/hg38ToGorGor4.over.chain")
		elif sp2 == "macaque":
			lo = LiftOver("chain/hg38ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif sp2 == "mouse":
			lo = LiftOver("chain/hg38ToMm10.over.chain") 
		elif sp2 == "rat":
			lo = LiftOver("chain/hg38ToRn6.over.chain") 
	elif sp1 == "chimp":
		if sp2 == "human":
			lo = LiftOver("chain/panTro5ToHg38.over.chain") 
		elif sp2 == "macaque":
			lo = LiftOver("chain/panTro5ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif sp2 == "mouse":
			lo = LiftOver("chain/panTro5ToMm10.over.chain") 
		elif sp2 == "rat":
			lo = LiftOver("chain/panTro5ToRn6.over.chain") 
	elif sp1 == "macaque":
		if sp2 == "human":
			lo = LiftOver("chain/rheMac10ToHg38.over.chain") 
		elif sp2 == "chimp":
			lo = LiftOver("chain/rheMac10ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToPanTro5.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif sp2 == "mouse":
			lo = LiftOver("chain/rheMac10ToMm10.over.chain")
		elif sp2 == "rat":
			lo = LiftOver("chain/rheMac10ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToRn6.over.chain")
	elif sp1 == "gorilla":
		if sp2 == "human":
			lo = LiftOver("chain/gorGor4ToHg38.over.chain") 
	elif sp1 == "mouse":		
		if sp2 == "human":
			lo = LiftOver("chain/mm10ToHg38.over.chain") 		
		if sp2 == "chimp":
			lo = LiftOver("chain/mm10ToPanTro5.over.chain") 
		elif sp2 == "macaque":
			lo = LiftOver("chain/mm10ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif sp2 == "rat":
			lo = LiftOver("chain/mm10ToRn6.over.chain") 
	elif sp1 == "rat":		
		if sp2 == "human":
			lo = LiftOver("chain/rn6ToHg38.over.chain") 		
		if sp2 == "chimp":
			lo = LiftOver("chain/rn6ToPanTro5.over.chain") 
		elif sp2 == "macaque":
			lo = LiftOver("chain/rn6ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif sp2 == "mouse":
			lo = LiftOver("chain/rn6ToMm10.over.chain") 
	return lo,lo2,eq	


def run_liftover(coords_ex,coords_annot,lo,lo2,eq,sp1,sp):
	'''Extract liftover coordinates'''
	status = {}
	g_to_t = {}
	liftex = open("liftover/" + sp1 + "_transcriptome_to_" + sp + ".exons.gtf","w+")
	statusex = open("liftover/" + sp1 + "_transcriptome_to_" + sp + ".exons.status.txt","w+")
	multiple = open("liftover/" + sp1 + "_transcriptome_to_" + sp + ".multiple.status.txt","w+")

	for name in coords_ex:
		if "_biotype" in name:
			continue
		#Define gene_name based on protein_coding gene, in case many are overlapping select the longest locus
		gene_name = coords_ex[name][4]
		if "-" in coords_ex[name][4].replace("-HSA-","_HSA_").replace("-PTR-","_PTR_").replace("-GGO-","_GGO_").replace("-MML-","_MML_").replace("-MMU-","_MMU_").replace("-RAT-","_RAT_"):
			l = 0
			for name2 in coords_ex[name][4].replace("-HSA-","_HSA_").replace("-PTR-","_PTR_").replace("-GGO-","_GGO_").replace("-MML-","_MML_").replace("-MMU-","_MMU_").replace("-RAT-","_RAT_").split("-"):
				if "STR" in name2:
					if l == 0:
						gene_name = name2
				else:
					l2 = coords_annot[name2][3] - coords_annot[name2][2]
					if l2 > l:
						l = l2
						gene_name = coords_annot[name2][4]
		elif "ENS" in coords_ex[name][4]:
			gene_name = coords_annot[coords_ex[name][4]][4]
		else:
			gene_name = coords_ex[name][4] #6?

		if not name in status:
			status[name] = ["aligned",gene_name]

		status["GENE-" + coords_ex[name][4]] = [gene_name,coords_ex[name][5],coords_ex[name][6]]
		if not gene_name in g_to_t:
			g_to_t[gene_name] = []

		#Define new transcriptome based on updated gene names

		for n,r in enumerate(coords_ex[name][2]):
			if lo2 == "none":
				try:
					st = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][2][n]), coords_ex[name][1])
					en = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][3][n]), coords_ex[name][1])
				except:
					st = []
					en = []
			else:
				try:
					st2 = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][2][n]), coords_ex[name][1])
					en2 = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][3][n]), coords_ex[name][1])
					st = lo2.convert_coordinate(st2[0][0], int(st2[0][1]), st2[0][2])
					en = lo2.convert_coordinate(en2[0][0], int(en2[0][1]), en2[0][2])
				except:
					st = []
					en = []				
			multiple.write(name + '\ttranscript\t' + str(st) + '\t' + str(en) + '\n')
			try:
				if (len(st) == 1) and (len(en) == 1):
					if st[0][0] != en[0][0]:
						statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tdiff_chrm\n')
						status[name] = ["rearranged",gene_name]
						g_to_t[gene_name].append(status[name][0])
					elif st[0][2] != en[0][2]:
						statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tdiff_strand\n')	
						status[name] = ["rearranged",gene_name]
						g_to_t[gene_name].append(status[name][0])
					elif abs(int(en[0][1])-int(st[0][1])) > abs(int(coords_ex[name][3][n])-int(coords_ex[name][2][n]))*1.5:
						statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tambiguous_synteny\n')
						status[name] = ["dubious_duplication",gene_name]
						g_to_t[gene_name].append(status[name][0])
					else:						
						chrm = st[0][0].replace("chr","").replace("Un_","").replace("v1_random","").replace("v1",".1")
						if 'NW' in chrm:
							try:
								chrm = eq["NW_" + chrm.split("NW_")[1]] 
							except:
								pass
						s1 = st[0][1]
						s2 = en[0][1]
						if st[0][2] != coords_ex[name][1]: #Fix to possible bug in pyliftover when converting one strand to another
							s1 = s1 - 2
						if en[0][2] != coords_ex[name][1]: #Fix to possible bug in pyliftover when converting one strand to another
							s2 = s2 - 2	
						if int(s2) >= int(s1):		
							liftex.write(chrm + "\tsynteny\texon\t" + str(s1+1) + "\t" + str(s2) + "\t.\t" + st[0][2] + '\t.\tgene_id "' + coords_ex[name][4] + '"; gene_name "' + gene_name + '"; transcript_id "' + name + '"; gene_biotype "' + coords_ex[name][5] + '";\n')
						else:
							liftex.write(chrm + "\tsynteny\texon\t" + str(s2) + "\t" + str(s1+1) + "\t.\t" + st[0][2] + '\t.\tgene_id "' + coords_ex[name][4] + '"; gene_name "' + gene_name + '"; transcript_id "' + name + '"; gene_biotype "' + coords_ex[name][5] + '";\n')
						statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tlifted\n')
						status[name] = ["aligned",gene_name]
						g_to_t[gene_name].append(status[name][0])
				elif (len(st) > 1) or (len(en) > 1):
					statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tmultiple_lifted\n')
					status[name] = ["dubious_duplication",gene_name]
					g_to_t[gene_name].append(status[name][0])
				else:
					statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tno_lifted\n')	
					status[name] = ["unaligned",gene_name]
					g_to_t[gene_name].append(status[name][0])
			except:
			 	statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + coords_ex[name][6] + '\t' + sp + '\texon\tno_lifted\n')
			 	status[name] = ["unaligned",gene_name]
			 	g_to_t[gene_name].append(status[name][0])
	
	liftex.close()
	statusex.close()
	multiple.close()
	return status,g_to_t


def identify_expression(sp1,sp,make_blast):
	#Check transcripts that overlap the syntenic regions
	ov = {}
	candidates = {}
	if sp == "human":
		trans_db = "../annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf"
	if sp == "chimp":
		trans_db = "../annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf"
	elif sp == "gorilla":
		trans_db = "../annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf"
	elif sp == "macaque":	
		trans_db = "../annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf"
	elif sp == "mouse":
		trans_db = "../annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf"
	elif sp == "rat":
		trans_db = "../annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf"	
	if make_blast == "yes":	
		os.system("intersectBed -wao -s -a liftover/" + sp1 + "_transcriptome_to_" + sp + ".exons.gtf -b " + trans_db + " > tmp/" + sp1 + "_transcriptome_to_" + sp + ".exons.ov")
	gene_ids = {}
	for line in open("tmp/" + sp1 + "_transcriptome_to_" + sp2 + ".exons.ov"):
		if line.count("\texon\t") < 2:
			continue
		t1 = line.split('transcript_id "')[1].split('"')[0]
		t2 = line.split('transcript_id "')[2].split('"')[0]
		gene_ids[t2] = line.split('gene_id "')[2].split('"')[0]
		if not t1 in ov:
			ov[t1] = {}
		if not t2 in ov[t1]:
			ov[t1][t2] = 0
		ov[t1][t2] = ov[t1][t2] + int(line.rstrip("\n").split("\t")[-1])

	for t1 in ov:
		candidates[t1] = ["none",0,"none"]
		for t2 in ov[t1]:
			if ov[t1][t2] > candidates[t1][1]:
				candidates[t1] = [t2,ov[t1][t2],gene_ids[t2]]

	return candidates


def blast_paralogs(fasta,status,sp1):
	#Run BLAST to find paralogs
	out = open("blast/" + sp1 + ".trans.paralogs.txt","w+")
	out.write("orf_id\tgene_id\tblast_sp\tquery\tquery_ov\tevalue\n")
	blasted = {}
	print("Running BLAST to same species")
	if sp1 == "human":		
	 	os.system("blastx -query " + fasta + " -db ../annotation/Homo_sapiens.GRCh38.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_trans_blastx.txt")
	elif sp1 == "chimp":
	 	os.system("blastx -query " + fasta + " -db ../annotation/Pan_troglodytes.Pan_tro_3.0.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_trans_blastx.txt")
	elif sp1 == "macaque":
	 	os.system("blastx -query " + fasta + " -db ../annotation/Macaca_mulatta.Mmul_10.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_trans_blastx.txt")
	elif sp1 == "mouse":
	 	os.system("blastx -query " + fasta + " -db ../annotation/Mus_musculus.GRCm38.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_trans_blastx.txt")
	elif sp1 == "rat":
	 	os.system("blastx -query " + fasta + " -db ../annotation/Rattus_norvegicus.Rnor_6.0.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_trans_blastx.txt")
			
	for line in open("tmp/" + sp1 + "_trans_blastx.txt"):
		name = line.split("\t")[0]
		if not name + "--" + sp1 in blasted:
			blasted[name + "--" + sp1] = []
		blasted[name + "--" + sp1].append(line.split("\t")[1] + "\t" + status[sp][name][1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10])

	for orf in blasted:
		for elemento in blasted[orf]:
				out.write(orf.split("--")[0] + "\t" + sp1 + "\t" + elemento + "\n")
	out.close()
	return blasted	


def blast_orthologs(fasta,status,sp1,species,coords_ex,make_blast):
	#Run BLAST to find orthologs
	out = open("blast/" + sp1 + ".trans.orthologs.txt","w+")
	out.write("trans_id\tmatch_id\tmatch_gene_id\tgene_id\tcoverage\tlength\n")
	blasted = {}
	for sp in species:
		if sp == "human":
			trans_db = "../annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed"
		elif sp == "chimp":
			trans_db = "../annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed"
		elif sp == "gorilla":
			trans_db = "../annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed"
		elif sp == "macaque":	
			trans_db = "../annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed"
		elif sp == "mouse":
			trans_db = "../annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed"
		elif sp == "rat":
			trans_db = "../annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed"	
		
		print("Finding genes")
		t_to_g = {}
		for line in open(trans_db + ".gtf"):
			if "\ttranscript\t" in line:
				t_to_g[line.split('transcript_id "')[1].split('"')[0]] = line.split('gene_id "')[1].split('"')[0]

		print("Running BLAST to " + sp)
		prog = "blastn"
		if sp in ("mouse","rat") and not sp1 in ("gorilla","mouse","rat","gorilla_cm","mouse_lv","rat_lv"):
			prog = "tblastx"
		##This step can be slow with tblastx, it can be also parallelized in a HPC
		if make_blast == "yes":
			os.system(prog + " -query " + fasta + " -db " + trans_db + " -strand plus -evalue 1e-4 -num_threads 4 -max_target_seqs 100 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq sframe\" > tmp/" + sp1 + "_" + sp + "_trans_" + prog + ".txt")
		intervals = {}
		for line in open("tmp/" + sp1 + "_" + sp + "_trans_" + prog + ".txt"):
			if line.rstrip("\n").split("\t")[-1] in ("-1","-2","-3"): #If match is to the opposite strand
				continue
			name = line.split("\t")[0]
			name2 = line.split("\t")[1]
			if not name in coords_ex:
				continue
			if not name in intervals:
				intervals[name] = {}
			if not name2 in intervals[name]:
				intervals[name][name2] = []

			for n in range(int(line.split("\t")[6]) , int(line.split("\t")[7])+1):
					intervals[name][name2].append(n)
			intervals[name][name2] = list(set(intervals[name][name2]))				

		for name in intervals:
			m = 0
			t = "none"
			for name2 in intervals[name]:
				l = len(intervals[name][name2])
				if l > m:
					t = name2
					m = l

			l_t = sum(map(int,coords_ex[name][3])) - sum(map(int,coords_ex[name][2])) + len(coords_ex[name][3])

			if t == "none":
				continue

			if not name in status[sp]:
				out.write(name + "\t" + t + "\t" + t_to_g[t] + "\t" + str(m) + "\t" + str(l_t) + "\n")
				continue

			if not status[sp][name][1] + "--" + sp in blasted:
				blasted[status[sp][name][1] + "--" + sp] = []

			blasted[status[sp][name][1] + "--" + sp].append(t_to_g[t])
			out.write(name + "\t" + t + "\t" + t_to_g[t] + "\t" + str(m) + "\t" + str(l_t) + "\n")

	return blasted


def get_exp(sp1, species):
	tpms_u = {}
	tpms = {}
	cats = []
	for file in glob.glob("../1_mapping/*sec/*29bp_5prime.fq.gzAligned.sortedByCoord.out.bam.unique.bam_abundance"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		for line in open(file):
			if line.startswith("Gene"):
				continue
			if not line.split("\t")[0] in tpms_u:
				tpms_u[line.split("\t")[0]] = {}
			tpms_u[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]
	
	for file in glob.glob("../1_mapping/*sec/*trimmed.fq.gzAligned.sortedByCoord.out.bam.unique.bam_abundance"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		for line in open(file):
			if line.startswith("Gene"):
				continue
			if not line.split("\t")[0] in tpms_u:
				tpms_u[line.split("\t")[0]] = {}
			tpms_u[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]		
	
	for file in glob.glob("../1_mapping/*sec/*29bp_5prime.fq.gz_abundance"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		cats.append(name)
		for line in open(file):
			if line.startswith("Gene"):
				continue
			if not line.split("\t")[0] in tpms:
				tpms[line.split("\t")[0]] = {}
			tpms[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]

	# if sp1 == "human":
	# 	for file in glob.glob("../all_lvs/*sec/*.29bp_5prime.fq.gz_abundance"):
	# 		name = file.split("/")[2]
	# 		cats.append(name)
	# 		for line in open(file):
	# 			if line.startswith("Gene"):
	# 				continue
	# 			if not line.split("\t")[0] in tpms:
	# 				tpms[line.split("\t")[0]] = {}
	# 			tpms[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]	
	
	for file in glob.glob("../1_mapping/*sec/*trimmed.fq.gz_abundance"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		cats.append(name)
		for line in open(file):
			if line.startswith("Gene"):
				continue
			if not line.split("\t")[0] in tpms:
				tpms[line.split("\t")[0]] = {}
			tpms[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]		
	
	return tpms,tpms_u,cats


def get_counts(sp1, species):
	counts = {}	
	for file in glob.glob("../1_mapping/*sec/*29bp_5prime.fq.gz.htseq.unique.counts"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		for line in open(file):
			if line.startswith("__"):
				continue
			if not line.split("\t")[0] in counts:
				counts[line.split("\t")[0]] = {}
			counts[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]

	# if sp1 == "human":
	# 	for file in glob.glob("../all_lvs/*sec/*trimmed.fq.gz.htseq.unique.counts"):
	# 		name = file.split("/")[2]
	# 		for line in open(file):
	# 			if line.startswith("__"):
	# 				continue
	# 			if not line.split("\t")[0] in counts:
	# 				counts[line.split("\t")[0]] = {}
	# 			counts[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]	
	
	for file in glob.glob("../1_mapping/*sec/*trimmed.fq.gz.htseq.unique.counts"):
		sp = file.split("/")[2].split("_")[0]
		name = file.split("/")[2]
		if not sp in species:
			continue
		for line in open(file):
			if line.startswith("__"):
				continue
			if not line.split("\t")[0] in counts:
				counts[line.split("\t")[0]] = {}
			counts[line.split("\t")[0]][name] = line.rstrip("\n").split("\t")[-1]		
	
	return counts


#1. Check synteny
coords_ex = get_coords(gtf,'exon')
coords_annot = get_coords(annot,'gene')
statuses = {}
g_to_t = {}
candidates = {}
for sp2 in species:
	(lo,lo2,eq) = load_chain(sp1,sp2)
	(statuses[sp2],g_to_t[sp2]) = run_liftover(coords_ex,coords_annot,lo,lo2,eq,sp1,sp2)
	candidates[sp2] = identify_expression(sp1,sp2,make_blast)
(tpms,tpms_u,cats) = get_exp(sp1,species_short)
counts = get_counts(sp1,species_short)


#2. Check homology, disabled for mouse, rat, gorilla
if sp1 in ("human","chimp","macaque"):
	# blastedp = blast_paralogs(fasta,statuses[species[0]],sp1)
	blasted = blast_orthologs(fasta,statuses,sp1,species,coords_ex,make_blast)
else:
	blasted = {}


#3. Find the most likely homolog gene and transcript
out1 = open("out/" + sp1 + ".trans.homologs.tsv","w+")
out2 = open("out/" + sp1 + ".gene.homologs.tsv","w+")
out3 = open("out/" + sp1 + ".gene.speciesexpression.tsv","w+")
out4 = open("out/" + sp1 + ".gene.unique.speciesexpression.tsv","w+")
out5 = open("out/" + sp1 + ".gene.unique.speciescounts.tsv","w+")
if sp1 == "human":
	out2.write("gene_id\tgtf_id\tgene_biotype\thom_status\tgtf_id_" + species[0] + "\tgtf_id_" + species[1] + "\tgtf_id_" + species[2] + "\tgtf_id_" + species[3]  + "\tgtf_id_" + species[4])
	out2.write("\thom_" + species[0] + "\tsyn_" + species[0] + "\thom_" + species[1] + "\tsyn_" + species[1])
	out2.write("\thom_" + species[2] + "\tsyn_" + species[2] + "\thom_" + species[3] + "\tsyn_" + species[3])
	out2.write("\thom_" + species[4] + "\tsyn_" + species[4] + "\n")
elif sp1 == "gorilla":
	out2.write("gene_id\tgtf_id\tgene_biotype\thom_status\tgtf_id_" + species[0])
	out2.write("\thom_" + species[0] + "\tsyn_" + species[0] + "\n")
else:
	out2.write("gene_id\tgtf_id\tgene_biotype\thom_status\tgtf_id_" + species[0] + "\tgtf_id_" + species[1] + "\tgtf_id_" + species[2] + "\tgtf_id_" + species[3])
	out2.write("\thom_" + species[0] + "\tsyn_" + species[0] + "\thom_" + species[1] + "\tsyn_" + species[1])
	out2.write("\thom_" + species[2] + "\tsyn_" + species[2] + "\thom_" + species[3] + "\tsyn_" + species[3] + "\n")	
out3.write("gene_id\tens_id\tgene_biotype")
out4.write("gene_id\tens_id\tgene_biotype")
out5.write("gene_id\tens_id\tgene_biotype")
for cat in cats:
	out3.write("\tTPM_" + cat)
	out4.write("\tTPMu_" + cat)
	out5.write("\tcounts_" + cat)
out3.write("\n")
out4.write("\n")
out5.write("\n")
all_genes = {}
for n,sp2 in enumerate(species):
	for trans in statuses[sp2]:
		if "GENE-" in trans:
			continue
		#Generate gene information
		if not statuses[sp2][trans][1] in all_genes:
			all_genes[statuses[sp2][trans][1]] = [[],statuses[sp2][trans][1],"none","none","none","none","none",0,0,0,0,0]
		if not trans in all_genes[statuses[sp2][trans][1]][0]:
			all_genes[statuses[sp2][trans][1]][0].append(trans)
		
		if trans in candidates[sp2]:
			out1.write(trans + "\t" + sp2 + "\t" + statuses[sp2][trans][1] + "\t" + candidates[sp2][trans][0] + "\t" + str(candidates[sp2][trans][1]) + "\t" + candidates[sp2][trans][2] + "\n")
			if candidates[sp2][trans][1] > all_genes[statuses[sp2][trans][1]][n+7]:
				all_genes[statuses[sp2][trans][1]][n+2] = candidates[sp2][trans][2]
				all_genes[statuses[sp2][trans][1]][n+7] = candidates[sp2][trans][1]
		else:
			out1.write(trans + "\t" + sp2 + "\t" + statuses[sp2][trans][1] + "\t" +  statuses[sp2][trans][0] + "\tnone\t0\n")

for gene in statuses[sp2]:
	if not "GENE-" in gene:
		continue
	if sp1 == "human":
		if statuses[sp2][gene][0] in all_genes:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tfull_synteny\t" + all_genes[statuses[sp2][gene][0]][2] + "\t" + all_genes[statuses[sp2][gene][0]][3] + "\t" + all_genes[statuses[sp2][gene][0]][4] + "\t" + all_genes[statuses[sp2][gene][0]][5] + "\t" + all_genes[statuses[sp2][gene][0]][6])
		else:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tno_full_synteny\tnone\tnone\tnone\tnone\tnone")
	elif sp1 == "gorilla":
		if statuses[sp2][gene][0] in all_genes:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tfull_synteny\t" + all_genes[statuses[sp2][gene][0]][2])
		else:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tno_full_synteny\tnone")	
	else:
		if statuses[sp2][gene][0] in all_genes:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tfull_synteny\t" + all_genes[statuses[sp2][gene][0]][2] + "\t" + all_genes[statuses[sp2][gene][0]][3] + "\t" + all_genes[statuses[sp2][gene][0]][4] + "\t" + all_genes[statuses[sp2][gene][0]][5])
		else:
			out2.write(statuses[sp2][gene][2] + "\t" + statuses[sp2][gene][0] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"] + "\tno_full_synteny\tnone\tnone\tnone\tnone")		
	for sp2 in species:
		non = 1
		if statuses[sp2][gene][0] + "--" + sp2 in blasted:
			for elemento in list(set(blasted[statuses[sp2][gene][0] + "--" + sp2])):
				if elemento in all_genes[statuses[sp2][gene][0]]: #Add the syntenic homolog
					out2.write("\t" + blasted[statuses[sp2][gene][0] + "--" + sp2][0] + "--syntenic")
					non = 0
					break
			if non == 1: #If it is not aligned, just add the first homolog
				out2.write("\t" + blasted[statuses[sp2][gene][0] + "--" + sp2][0])
		else:
			out2.write("\tno")
		if "aligned" in list(set(g_to_t[sp2][statuses[sp2][gene][0]])):
			out2.write("\taligned")
		elif "dubious_duplication" in list(set(g_to_t[sp2][statuses[sp2][gene][0]])):
			out2.write("\trearranged")
		elif "rearranged" in list(set(g_to_t[sp2][statuses[sp2][gene][0]])):
			out2.write("\trearranged")
		else:
			out2.write("\tunaligned")		
		#out2.write(";".join(list(set(g_to_t[sp2][statuses[sp2][gene][0]]))))
	out2.write("\n")

	out3.write(statuses[sp2][gene][0] + "\t" + statuses[sp2][gene][2] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"])
	out4.write(statuses[sp2][gene][0] + "\t" + statuses[sp2][gene][2] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"])
	out5.write(statuses[sp2][gene][0] + "\t" + statuses[sp2][gene][2] + "\t" + coords_ex[statuses[sp2][gene][2] + "_biotype"])

	for cat in cats:
		if species_short[0] + "_" in cat:
			#print(statuses[sp2][gene][2])
			try: 
				value = tpms[statuses[sp2][gene][2]][cat] 
			except:
				value = "0"	
		else:		
			for n,sp in enumerate(species_short[1:]):
				if sp + "_" in cat:
					try: 
						value = tpms[all_genes[statuses[sp2][gene][0]][n+2]][cat] 
					except:
						value = "0"
		out3.write("\t" + value)
	out3.write("\n")

	for cat in cats:
		if "SRR" in cat:
			continue
		if species_short[0] + "_" in cat:
			try: 
				value = tpms_u[statuses[sp2][gene][2]][cat] 
			except:
				value = "0"	
		else:		
			for n,sp in enumerate(species_short[1:]):
				if sp + "_" in cat:
					try: 
						value = tpms_u[all_genes[statuses[sp2][gene][0]][n+2]][cat] 
					except:
						value = "0"
		out4.write("\t" + value)
	out4.write("\n")

	for cat in cats:
		if species_short[0] + "_" in cat:
			try: 
				value = counts[statuses[sp2][gene][2]][cat] 
			except:
				value = "-1"	
		else:		
			for n,sp in enumerate(species_short[1:]):
				if sp + "_" in cat:
					try: 
						value = counts[all_genes[statuses[sp2][gene][0]][n+2]][cat] 
					except:
						value = "-1"
		out5.write("\t" + value)
	out5.write("\n")

out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
exit(0)
