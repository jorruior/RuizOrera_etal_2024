#!/usr/bin/env python
__author__="jruizor"
__date__ ="$Aug 5, 2022 9:43:43 PM$"
'''Find the most likely homolog per each transcript ORF and quantify
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


gtf = sys.argv[1] #../1_mapping/human_pooled/human.orfs.gtf_Psites.bed
fasta = sys.argv[2] #../1_mapping/human_pooled/human.orfs.fa
sp1 = sys.argv[3] #human
make_blast = "no"
if "human" in sp1:
	species = ["chimp","gorilla_cm","macaque","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["pt","gg","rm","mm","rn"]
elif "chimp" in sp1:
	species = ["human","macaque","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["hs","rm","mm","rn"]
elif "gorilla" in sp1:
	species = ["human"] #the script is currently designed to work with these species
	species_short = ["hs"]	
elif "macaque" in sp1:
	species = ["human","chimp","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["hs","pt","mm","rn"]
elif "mouse" in sp1:
	species = ["human","chimp","macaque","rat_lv"] #the script is currently designed to work with these species
	species_short = ["mm","hs","pt","rm","rn"]
elif "rat" in sp1:
	species = ["human","chimp","macaque","mouse_lv"] #the script is currently designed to work with these species
	species_short = ["rn","hs","pt","rm","mm"]	


def get_coords(gtf):
	'''Extract features from gtf'''
	coords_ex = {}
	for line in open(gtf):
		name = line.split("\t")[3]
		if line.split("\t")[5].rstrip("\n") == "+" and not "\thq\t" in line and not "\tp2\t" in line: #Positive: last site per codon
			continue
		elif line.split("\t")[5].rstrip("\n") == "-" and not "\tp0\t" in line: #Negative: first site per codon
			continue		
		if not name in coords_ex:
			coords_ex[name] = ["chr" + line.split("\t")[0],line.split("\t")[5].rstrip("\n"),[],[]]
		coords_ex[name][2].append(line.split("\t")[1])
		coords_ex[name][3].append(line.split("\t")[2])

	return coords_ex


def load_chain(sp1,sp2):
	'''Load corresponding chain file'''
	eq = {}
	lo2 = "none"
	print("Running liftover")
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


def run_liftover(coords_ex,lo,lo2,eq,sp1,sp):
	'''Extract liftover coordinates'''
	status = {}
	g_to_t = {}
	liftex = open("liftover/" + sp1 + "_orfeome_to_" + sp + ".psites.gtf","w+")
	statusex = open("liftover/" + sp1 + "_orfeome_to_" + sp + ".psites.status.txt","w+")

	for name in coords_ex:
		gene_name = name
		if not gene_name in g_to_t:
			g_to_t[gene_name] = []

		status[name] = ["unaligned",gene_name]

		for n,r in enumerate(coords_ex[name][2]):
			if lo2 == "none":
				try:
					st = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][2][n]), coords_ex[name][1])
				except:
					st = []
			else:
				try:
					st2 = lo.convert_coordinate(coords_ex[name][0], int(coords_ex[name][2][n]), coords_ex[name][1])
					st = lo2.convert_coordinate(st2[0][0], int(st2[0][1]), st2[0][2])
				except:
					st = []
			try:
				if (len(st) == 1):					
					chrm = st[0][0].replace("chr","").replace("Un_","").replace("v1_random","").replace("v1",".1")
					if 'NW' in chrm:
						try:
							chrm = eq["NW_" + chrm.split("NW_")[1]] 
						except:
							pass
			
					if int(st[0][1]) > 0:
						s1 = st[0][1]
						if st[0][2] != coords_ex[name][1]: #Possible bug in pyliftover when converting one strand to another
							s1 = s1 - 2
						if s1 > 0:
							liftex.write(chrm + "\tsynteny\texon\t" + str(s1) + "\t" + str(s1) + "\t.\t" + st[0][2] + '\t.\tgene_id "' + name + '"; gene_name "' + name + '"; orf_id "' + name + '"; initial_pos "' + coords_ex[name][0] + ":" + coords_ex[name][2][n] + ":" + coords_ex[name][1] + '";\n')
							statusex.write(name + '\torf\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tlifted\n')
							status[name] = ["aligned",gene_name]
						g_to_t[gene_name].append(status[name][0])
					else:
						statusex.write(name + '\torf\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tno_lifted\n')
						status[name] = ["unaligned",gene_name]
						g_to_t[gene_name].append(status[name][0])					
				elif len(st) > 1:
					statusex.write(name + '\torf\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tmultiple_lifted\n')
					status[name] = ["dubious_duplication",gene_name]
					g_to_t[gene_name].append(status[name][0])
				else:
					statusex.write(name + '\torf\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tno_lifted\n')	
					status[name] = ["unaligned",gene_name]
					g_to_t[gene_name].append(status[name][0])
			except:
				statusex.write(name + '\torf\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tno_lifted\n')
				status[name] = ["unaligned",gene_name]
				
	liftex.close()
	statusex.close()
	return status,g_to_t


def identify_expression(sp1,sp,make_blast):
	#Check ORFs that overlap the syntenic regions
	out = open("liftover/" + sp1 + "_" + sp + ".orfs_to_transcripts_heart.synteny.txt","w+")
	#out.write("orf_id\ttrans_id\ttissue\n")
	ov = {}
	all_ov = {}
	candidates = {}
	done = []
	trans_db = "../1_mapping/" + sp + "_pooled/" + sp + ".orfs.allframes.bed"
	if "human" in sp:
		trans_db2 = "../annotation/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf"
	elif "chimp" in sp:
		trans_db2 = "../annotation/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf"
	elif "gorilla" in sp:
		trans_db2 = "../annotation/custom/Gorilla_gorilla.gorGor4.98.stringtie2021.heart.gtf"
	elif "macaque" in sp:
		trans_db2 = "../annotation/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf "
	elif "mouse" in sp:
		trans_db2 = "../annotation/custom/Mus_musculus.GRCm38.98.stringtie2021.heart.gtf"
	elif "rat" in sp:
		trans_db2 = "../annotation/custom/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.heart.gtf"

	if make_blast == "yes":
		os.system("intersectBed -wao -s -a liftover/" + sp1 + "_orfeome_to_" + sp + ".psites.gtf -b " + trans_db + " > tmp/" + sp1 + "_orfeome_to_" + sp + ".psites.ov")
		os.system("intersectBed -wao -s -a liftover/" + sp1 + "_orfeome_to_" + sp + ".psites.gtf -b " + trans_db2 + " > tmp/" + sp1 + "_orfeome_to_" + sp + "_heart_transcripts.ov")
	gene_ids = {}
	for line in open("tmp/" + sp1 + "_orfeome_to_" + sp2 + ".psites.ov"):
		t1 = line.split('orf_id "')[1].split('"')[0]
		if not t1 in ov:
			ov[t1] = {}
			all_ov[t1] = {}

		if ((line.split("\t")[6] == "+") and (line.split("\t")[3] == line.split("\t")[10]) and (':+"' in line) and (line.split("\t")[13] == "p2")) or ((line.split("\t")[6] == "+") and (line.split("\t")[3] == line.split("\t")[10]) and (':-"' in line) and (line.split("\t")[13] == "p1")) or ((line.split("\t")[6] == "-") and (line.split("\t")[3] == line.split("\t")[10]) and (':-"' in line) and (line.split("\t")[13] == "p0")) or ((line.split("\t")[6] == "-") and (line.split("\t")[3] == line.split("\t")[10]) and (':+"' in line) and (line.split("\t")[13] == "p2")):
			t2 = line.split("\t")[-4]
			if not t2 in ov[t1]:
				ov[t1][t2] = 0
			ov[t1][t2] = ov[t1][t2] + 1

		if line.split("\t")[3] == line.split("\t")[10]:
			t2 = line.split("\t")[-4]
			if not t2 in all_ov[t1]:
				all_ov[t1][t2] = 0
			all_ov[t1][t2] = all_ov[t1][t2] + 1			


	for t1 in ov:
		candidates[t1] = ["none",0]
		for t2 in ov[t1]:
			if ov[t1][t2] > candidates[t1][1]:
				candidates[t1] = [t2,ov[t1][t2]]

	for line in open("tmp/" + sp1 + "_orfeome_to_" + sp2 + "_heart_transcripts.ov"):
		t1 = line.split('orf_id "')[1].split('"')[0]
		if line.count("\texon\t") == 2:
			t2 = line.split('transcript_id "')[1].split('"')[0]
			if not t1 + "-" + t2 in done:
				done.append(t1 + "-" + t2)
				out.write(t1 + "\t" + t2 + "\t" + sp + "\texon\theart\n")
		elif line.count("\tCDS\t") > 0:
			t2 = line.split('transcript_id "')[1].split('"')[0]
			if not t1 + "-" + t2 in done:
				done.append(t1 + "-" + t2)
				out.write(t1 + "\t" + t2 + "\t" + sp + "\tCDS\theart\n")

	for t1 in all_ov:
		for t2 in all_ov[t1]:
			if t2 in ov[t1]:
				out.write(t1 + "\t" + t2 + "\t" + sp + "\tif_ORF\theart\n")
			else:
				out.write(t1 + "\t" + t2 + "\t" + sp + "\tof_ORF\theart\n")

	out.close()
	return candidates


def blast_paralogs(fasta,sp1,make_blast):
	#Run BLAST to find paralogs
	out = open("blast/" + sp1 + ".orfs_to_trans.paralogs.txt","w+")
	out.write("orf_id\tblast_sp\tquery\tquery_ov\tevalue\n")
	blasted = {}
	print("Running BLAST to same species")
	if make_blast == "yes":
		if sp1 == "human":		
		 	os.system("blastp -query " + fasta + " -db ../annotation/Homo_sapiens.GRCh38.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe sframe\" > tmp/" + sp1 + "_orf_blastp.txt")
		elif sp1 == "chimp":
		 	os.system("blastp -query " + fasta + " -db ../annotation/Pan_troglodytes.Pan_tro_3.0.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe sframe\" > tmp/" + sp1 + "_orf_blastp.txt")
		elif sp1 == "macaque":
		 	os.system("blastp -query " + fasta + " -db ../annotation/Macaca_mulatta.Mmul_10.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe sframe\" > tmp/" + sp1 + "_orf_blastp.txt")
		elif sp1 == "mouse":
		 	os.system("blastp -query " + fasta + " -db ../annotation/Mus_musculus.GRCm38.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe sframe\" > tmp/" + sp1 + "_orf_blastp.txt")
		elif sp1 == "rat":
		 	os.system("blastp -query " + fasta + " -db ../annotation/Rattus_norvegicus.Rnor_6.0.pep.all -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe sframe\" > tmp/" + sp1 + "_orf_blastp.txt")
			
	for line in open("tmp/" + sp1 + "_orf_blastp.txt"):
		name = line.split("\t")[0]
		if not name + "--" + sp1 in blasted:
			blasted[name + "--" + sp1] = []
		blasted[name + "--" + sp1].append(line.split("\t")[1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10])
		out.write(line.split("\t")[0] + "\t" + sp1 + "\t" + line.split("\t")[1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10] + "\n")
	out.close()
	return blasted


def blast_orthologs(fasta,sp1,species,make_blast):
	'''Run BLAST to find orthologs'''
	out = open("blast/" + sp1 + ".orfs_to_trans.orthologs.txt","w+")
	out.write("orf_id\tblast_sp\tquery\tquery_ov\tevalue\n")
	blasted = {}

	for sp in species:
		trans_db1 = "../1_mapping/" + sp + "_pooled/" + sp + ".orfs.fa.db"		
		print("Running BLAST to " + sp)
		if make_blast == "yes":
			if not sp1 in ("human_gencode","human_paper"):
				os.system("blastp -query " + fasta + " -db " + trans_db1 + " -evalue 1e-2 -num_threads 1 -max_target_seqs 100 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_" + sp + "_orf_blastp.txt")
		for line in open("tmp/" + sp1 + "_" + sp + "_orf_blastp.txt"):
			name = line.split("\t")[0]
			if not name + "--" + sp in blasted:
				blasted[name + "--" + sp] = []
			blasted[name + "--" + sp].append(line.split("\t")[1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10])
			out.write(line.split("\t")[0] + "\t" + sp + "\t" + line.split("\t")[1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10] + "\n")

	return blasted


def evidence_othertissues(fasta,sp1,species,make_blast):
	'''Run BLAST to find orthologs'''
	out = open("blast/" + sp1 + ".orfs_to_trans_othertissues.orthologs.txt","w+")
	out.write("orf_id\tsp\ttis\torf_ov\torf_ov_len\tquery\tquery_ov\tevalue\n")

	for sp in species:
		for trans_db1 in glob.glob("../../other_tissues/" + sp.split("_")[0] + "*_pooled/*.orfs.fa"):
			blasted = {}
			tis = trans_db1.split(sp.split("_")[0] + "_")[1].split("_pooled")[0]
			print("Running BLAST to " + trans_db1)
			if make_blast == "yes":
				if not sp1 in ("human_gencode","human_paper"):
					os.system("blastp -query " + fasta + " -db " + trans_db1 + ".db -evalue 1e-2 -num_threads 1 -max_target_seqs 10 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" > tmp/" + sp1 + "_" + sp.split("_")[0] + "_" + tis + "_orf_blastp.txt")
					os.system("intersectBed -wao -s -a liftover/" + sp1 + "_orfeome_to_" + sp + ".psites.gtf -b " + trans_db1.replace("orfs.fa","orfs.allframes.bed") + " > tmp/" + sp1 + "_orfeome_to_" + sp.split("_")[0] + "_" + tis + ".psites.ov")
			for line in open("tmp/" + sp1 + "_" + sp.split("_")[0] + "_" + tis + "_orf_blastp.txt"):
				name = line.split("\t")[0]
				if not name in blasted:
					blasted[name] = []
				blasted[name].append(line.split("\t")[1] + "\t" + str(int(line.split("\t")[7]) - int(line.split("\t")[6])) + "\t" + line.split("\t")[10])
			
			ov = {}
			for line in open("tmp/" + sp1 + "_orfeome_to_" + sp.split("_")[0] + "_" + tis + ".psites.ov"):
				t1 = line.split('orf_id "')[1].split('"')[0]
				if not t1 in ov:
					ov[t1] = {}

				if ((line.split("\t")[6] == "+") and (line.split("\t")[3] == line.split("\t")[10]) and (':+"' in line) and (line.split("\t")[13] == "p2")) or ((line.split("\t")[6] == "+") and (line.split("\t")[3] == line.split("\t")[10]) and (':-"' in line) and (line.split("\t")[13] == "p1")) or ((line.split("\t")[6] == "-") and (line.split("\t")[3] == line.split("\t")[10]) and (':-"' in line) and (line.split("\t")[13] == "p0")) or ((line.split("\t")[6] == "-") and (line.split("\t")[3] == line.split("\t")[10]) and (':+"' in line) and (line.split("\t")[13] == "p2")):
					t2 = line.split("\t")[-4]
					if not t2 in ov[t1]:
						ov[t1][t2] = 0
					ov[t1][t2] = ov[t1][t2] + 1

			candidates = {}
			for t1 in ov:
				candidates[t1] = ["none",0]
				for t2 in ov[t1]:
					if ov[t1][t2] > candidates[t1][1]:
						candidates[t1] = [t2,ov[t1][t2]]

			for t1 in candidates:
				if t1 in blasted:
					for b in blasted[t1]:
						out.write(t1 + "\t" + sp.split("_")[0] + "\t" + tis + "\t" + candidates[t1][0] + "\t" + str(candidates[t1][1]) + "\t" + b + "\n")
				else:
					out.write(t1 + "\t" + sp.split("_")[0] + "\t" + tis + "\t" + candidates[t1][0] + "\t" + str(candidates[t1][1]) + "\tnone\t0\t1000\n")

	out.close()


def get_exp(sp1, species):
	'''Extract values to calculate TPMs'''
	tpms = {}
	counts = {}
	mapped = {}
	headers = {}
	for line in open("../../primates_lv/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000
	for line in open("../../ipsc_cm/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000	
	for line in open("../../rodents_lv/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000	
	for line in open("../../cardioids/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000

	#Ribo-seq
	for sp2 in [sp1] + species:
		if sp2 in ("mouse","rat"):
			sp2 = sp2 + "_lv"
		if sp2 == "gorilla":
			sp2 = sp2 + "_cm"			
		tpms[sp2] = {}
		counts[sp2] = {}
		for line in open("../1_mapping/" + sp2 + "_pooled/" + sp2 + ".orfs.gtf_orf_quant_Psites.txt"):
			if line.startswith("orf"):
				tpms[sp2]["main"] = []
				counts[sp2]["main"] = []
				if (sp1 == sp2):
					for l in line.rstrip("\n").split("\t"):
						if "_Ri_" in l:
							tpms[sp2]["main"].append(l.split("_Ri_")[0] + "_Ri--" + sp2)
							counts[sp2]["main"].append(l.split("_Ri_")[0] + "_Ri--" + sp2)
						else:
							tpms[sp2]["main"].append(l)
							counts[sp2]["main"].append(l)
					
					tpms[sp2]["main"] = "\t".join(tpms[sp2]["main"])
					counts[sp2]["main"] = "\t".join(counts[sp2]["main"])
				else:
					for l in line.rstrip("\n").split("\t")[10:]:
						if "_Ri_" in l:
							tpms[sp2]["main"].append(l.split("_Ri_")[0] + "_Ri--" + sp2)
							counts[sp2]["main"].append(l.split("_Ri_")[0] + "_Ri--" + sp2)
					
					tpms[sp2]["main"] = "\t".join(tpms[sp2]["main"])
					counts[sp2]["main"] = "\t".join(counts[sp2]["main"])

			elif not line.split("\t")[0] in tpms[sp2]:
				#print(sp2 + "\t" + line.split("\t")[0])
				tpms[sp2][line.split("\t")[0]] = []
				counts[sp2][line.split("\t")[0]] = []
				if (sp1 == sp2):
					for y,l in enumerate(line.rstrip("\n").split("\t")[:10]): #Non-numeric ORF info
						tpms[sp2][line.split("\t")[0]].append(l)
						counts[sp2][line.split("\t")[0]].append(l)
					for l in line.rstrip("\n").split("\t")[10:]: #Counts
						y +=1
						if "ENST" in line.split("\t")[0] or "merged" in line.split("\t")[0]:
							n = float(l) / ((float(line.split("\t")[0].replace("-rand","").split("_")[-1]) - float(line.split("\t")[0].replace("-rand","").split("_")[-2])) / 1000)
						else:
							n = float(l) / (float(line.split("\t")[0].split(":")[-1].replace("__v",""))/1000)
							tpms[sp2][line.split("\t")[0]].append(str(n / mapped[tpms[sp2]["main"].split("\t")[y].split("--")[0]]))
						counts[sp2][line.split("\t")[0]].append(str(l))

					tpms[sp2][line.split("\t")[0]] = "\t".join(tpms[sp2][line.split("\t")[0]])
					counts[sp2][line.split("\t")[0]] = "\t".join(counts[sp2][line.split("\t")[0]])
				else:
					for y,l in enumerate(line.rstrip("\n").split("\t")[10:]): #Counts						
						if "ENST" in line.split("\t")[0] or "merged" in line.split("\t")[0]:
							n = float(l) / ((float(line.split("\t")[0].replace("-rand","").split("_")[-1]) - float(line.split("\t")[0].replace("-rand","").split("_")[-2])) / 1000)
						else:
							n = float(l) / (float(line.split("\t")[0].split(":")[-1].replace("__v",""))/1000)
						tpms[sp2][line.split("\t")[0]].append(str(n / mapped[tpms[sp2]["main"].split("\t")[y].split("--")[0]]))
						counts[sp2][line.split("\t")[0]].append(str(l))

					tpms[sp2][line.split("\t")[0]] = "\t".join(tpms[sp2][line.split("\t")[0]])
					counts[sp2][line.split("\t")[0]] = "\t".join(counts[sp2][line.split("\t")[0]])

	#RNA-seq
	for sp2 in [sp1] + species:
		if sp2 in ("mouse","rat"):
			sp2 = sp2 + "_lv"
		if sp2 == "gorilla":
			sp2 = sp2 + "_cm"
		headers[sp2] = []
		for file in glob.glob("../1_mapping/*/*out.bam." + sp2 + ".htseq.rna.counts"):
			if "_iPS_" in file:
				continue
			if "_mR_" in file:
				headers[sp2].append(file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR--" + sp2)
			for line in open(file):
				if line.startswith("__"):
					continue
				if not line.split("\t")[0] + "-RNA" in tpms[sp2]:
					tpms[sp2][line.split("\t")[0] + "-RNA"] = []
					counts[sp2][line.split("\t")[0] + "-RNA"] = []
				n = float(line.split("\t")[1].rstrip("\n")) / (float(line.split("\t")[0].split(":")[-1].replace("__v",""))/1000)
				if "_mR_" in file:
					n = n / mapped[file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR"]
				tpms[sp2][line.split("\t")[0] + "-RNA"].append(str(n))	
				counts[sp2][line.split("\t")[0] + "-RNA"].append(str(line.split("\t")[1].rstrip("\n")))	

	return (tpms,counts,headers)


#1. Check synteny
coords_ex = get_coords(gtf)
statuses = {}
g_to_t = {}
candidates = {}
for sp2 in species:
	(lo,lo2,eq) = load_chain(sp1.split("_")[0],sp2.split("_")[0])
	(statuses[sp2],g_to_t[sp2]) = run_liftover(coords_ex,lo,lo2,eq,sp1,sp2)
	candidates[sp2] = identify_expression(sp1,sp2,make_blast)

seqs = {}
for line in open(fasta):
	if ">" in line:
		name = line.rstrip("\n").replace(">","")
		seqs[name] = ""
	else:##
		seqs[name] = line.rstrip("\n").replace(">","")

lens = {}
for s in seqs:
	lens[s.split("--")[0]] = len(seqs[s])*3
(tpms,counts,headers) = get_exp(sp1,species)


#2. Check homology
if sp1 in ("human","chimp","macaque"):
	blastedp = blast_paralogs(fasta,sp1,make_blast)
	blasted = blast_orthologs(fasta,sp1,species,make_blast)
	evidence_othertissues(fasta,sp1,species,make_blast)
else:
	blastedp = {}
	blasted = {}


#3. Find the most likely homolog gene and transcript
out2 = open("out/" + sp1 + ".orf.homologs.tsv","w+")
out3 = open("out/" + sp1 + ".orf.speciesexpression.tsv","w+")
out4 = open("out/" + sp1 + ".orf.speciescounts.tsv","w+")
out2.write("orf_id\tlen")
for sp2 in species:
	out2.write("\tsyn_" + sp2 + "\tgtf_id_" + sp2 + "\torf_aligned_" + sp2 + "\thom_" + sp2 + "\tsyn_aligned_" + sp2)
out2.write("\n")

if (sp1 == "human_gencode") or (sp1 == "human_paper"):
	out3.write("orf_id")

for cat in [sp1] + species:
	if cat == sp1:
		out3.write(tpms[cat]["main"])
		out4.write(tpms[cat]["main"])
	else:
		out3.write("\t" + tpms[cat]["main"])
		out4.write("\t" + tpms[cat]["main"])
	out3.write("\t" + "\t".join(headers[cat]))
	out4.write("\t" + "\t".join(headers[cat]))
out3.write("\n")
out4.write("\n")

for gene in coords_ex:
	if not gene in lens:
		continue
	out2.write(gene + "\t" + str(lens[gene]))
	for sp2 in species:
		if gene in candidates[sp2]:
			out2.write("\tfull_synteny\t" + candidates[sp2][gene][0] + "\t" + str(candidates[sp2][gene][1]))
		else:
			out2.write("\tno_full_synteny\tnone\t0")

		hom = "none"
		if gene + "--" + sp2 in blasted:
			hom = blasted[gene + "--" + sp2][0].split("\t")[0]
			out2.write("\t" + hom)
		else:
			out2.write("\tno")

		if gene in g_to_t[sp2]:
			try:
				out2.write("\t" + str(round(g_to_t[sp2][gene].count("aligned")/len(g_to_t[sp2][gene])*100,2)))
			except:
				out2.write("\t0")	
		else:
			out2.write("\t0")
	out2.write("\n")

	if not sp1 in tpms:
		out3.write(gene)
		out4.write(gene)
	else:
		if not gene in tpms[sp1]:
			continue
	
		#out3.write(gene + "\t" + coords_ex[gene][6] + "\t" + coords_ex[gene][5])
		out3.write(str(tpms[sp1][gene]) + "\t" + "\t".join(tpms[sp1][gene + "-RNA"]))
		out4.write(str(counts[sp1][gene]) + "\t" + "\t".join(counts[sp1][gene + "-RNA"]))

	for cat in species:

		#Ribo-seq	
		if cat.split("_")[0] == sp1.split("_")[0]:
			try:
				value = str(tpms[cat][gene])
			except:
				value = "\t".join(['00' for i in range(len(tpms[cat]["main"].split("\t")))])
		elif not gene in candidates[cat]:
			value = "\t".join(['0' for i in range(len(tpms[cat]["main"].split("\t")))])
		elif candidates[cat][gene][0] in tpms[cat]:
			value = tpms[cat][candidates[cat][gene][0]]
		elif hom in tpms[cat]:
			value = tpms[cat][hom]		
		else:
			value = "\t".join(['0' for i in range(len(tpms[cat]["main"].split("\t")))])
		
		#RNA-seq
		if cat.split("_")[0] == sp1.split("_")[0]:
			try:
				value2 = "\t".join(tpms[cat][gene + "-RNA"])
			except:
				value2 = "\t".join(['000' for i in range(len(headers[cat]))])	
		elif not gene in candidates[cat]:		
			value2 = "\t".join(['0' for i in range(len(headers[cat]))])	
		elif candidates[cat][gene][0] + "-RNA" in tpms[cat]:
			value2 = "\t".join(tpms[cat][candidates[cat][gene][0] + "-RNA"])
		elif hom + "-RNA" in tpms[cat]:
			value2 = "\t".join(tpms[cat][hom + "-RNA"])		
		else:	
			value2 = "\t".join(['0' for i in range(len(headers[cat]))])
		
		out3.write("\t" + str(value) + "\t" + str(value2))


		#Ribo-seq	
		if cat.split("_")[0] == sp1.split("_")[0]:
			try:
				value = str(counts[cat][gene])
			except:
				value = "\t".join(['0' for i in range(len(counts[cat]["main"].split("\t")))])
		elif not gene in candidates[cat]:
			value = "\t".join(['0' for i in range(len(counts[cat]["main"].split("\t")))])
		elif candidates[cat][gene][0] in counts[cat]:
			value = counts[cat][candidates[cat][gene][0]]
		elif hom in counts[cat]:
			value = counts[cat][hom]		
		else:
			value = "\t".join(['0' for i in range(len(counts[cat]["main"].split("\t")))])
		
		#RNA-seq
		if cat.split("_")[0] == sp1.split("_")[0]:
			try:
				value2 = "\t".join(counts[cat][gene + "-RNA"])
			except:
				value2 = "\t".join(['0' for i in range(len(headers[cat]))])	
		elif not gene in candidates[cat]:		
			value2 = "\t".join(['0' for i in range(len(headers[cat]))])	
		elif candidates[cat][gene][0] + "-RNA" in counts[cat]:
			value2 = "\t".join(counts[cat][candidates[cat][gene][0] + "-RNA"])
		elif hom + "-RNA" in counts[cat]:
			value2 = "\t".join(counts[cat][hom + "-RNA"])	
		else:	
			value2 = "\t".join(['0' for i in range(len(headers[cat]))])
		
		out4.write("\t" + str(value) + "\t" + str(value2))	
	out3.write("\n")
	out4.write("\n")

out2.close()
out3.close()
out4.close()

exit(0)
