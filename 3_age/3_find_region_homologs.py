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


gtf = sys.argv[1] #../1_mapping/human_cm_pooled/human.orfs.hq.bed
psites = sys.argv[2] #../1_mapping/human_cm_pooled/human.orfs.gtf_Psites.bed
sp1 = sys.argv[3] #human
if "human" in sp1:
	species = ["chimp","gorilla_cm","macaque","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["pt","gg","rm","mm","rn"]
	sp_short = "hs"
elif "chimp" in sp1:
	species = ["human","macaque","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["hs","rm","mm","rn"]
	sp_short = "pt"
elif "macaque" in sp1:
	species = ["human","chimp","mouse_lv","rat_lv"] #the script is currently designed to work with these species
	species_short = ["hs","pt","mm","rn"]
	sp_short = "rm"

def get_coords(gtf):
	'''Extract features from gtf'''
	coords_ex = {}
	for line in open(gtf):
		name = line.split("\t")[3]
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
		elif "gorilla" in sp2:
			lo = LiftOver("chain/hg38ToGorGor4.over.chain") 	
		elif sp2 == "macaque":
			lo = LiftOver("chain/hg38ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif "mouse" in sp2:
			lo = LiftOver("chain/hg38ToMm10.over.chain") 
		elif "rat" in sp2:
			lo = LiftOver("chain/hg38ToRn6.over.chain") 
	elif sp1 == "chimp":
		if sp2 == "human":
			lo = LiftOver("chain/panTro5ToHg38.over.chain") 
		elif sp2 == "macaque":
			lo = LiftOver("chain/panTro5ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToRheMac10.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif "mouse" in sp2:
			lo = LiftOver("chain/panTro5ToMm10.over.chain") 
		elif "rat" in sp2:
			lo = LiftOver("chain/panTro5ToRn6.over.chain") 
	elif sp1 == "macaque":
		if sp2 == "human":
			lo = LiftOver("chain/rheMac10ToHg38.over.chain") 
		elif sp2 == "chimp":
			lo = LiftOver("chain/rheMac10ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToPanTro5.over.chain")
			for line in open("/fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.equivalents"):
				eq[line.split("\t")[6]] = line.split("\t")[4]
		elif "mouse" in sp2:
			lo = LiftOver("chain/rheMac10ToMm10.over.chain")
		elif "rat" in sp2:
			lo = LiftOver("chain/rheMac10ToRheMac8.over.chain")
			lo2 = LiftOver("chain/rheMac8ToRn6.over.chain")
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
	liftex = open("liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf","w+")
	statusex = open("liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.status.txt","w+")
	multiple = open("liftover/" + sp1 + "_orfeome_to_" + sp + ".multiple.status.txt","w+")

	for name in coords_ex:
		gene_name = name
		if not gene_name in g_to_t:
			g_to_t[gene_name] = []

		status[name] = ["unaligned",gene_name]

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

			multiple.write(name + '\tCDS\t' + str(st) + '\t' + str(en) + '\n')
			try:
				if (len(st) == 1) and (len(en) == 1):
					if st[0][0] != en[0][0]:
						statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tdiff_chrm\n')
						status[name] = ["rearranged",gene_name]
						g_to_t[gene_name].append(status[name][0])
					elif st[0][2] != en[0][2]:
						statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tdiff_strand\n')	
						status[name] = ["rearranged",gene_name]
						g_to_t[gene_name].append(status[name][0])
					elif abs(int(en[0][1])-int(st[0][1])) > abs(int(coords_ex[name][3][n])-int(coords_ex[name][2][n]))*1.5:
						statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tambiguous_synteny\n')
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
							liftex.write(chrm + "\tsynteny\texon\t" + str(s1) + "\t" + str(s2) + "\t.\t" + st[0][2] + '\t.\tgene_id "' + name + '"; transcript_id "' + name + '";\n')
						else:
							liftex.write(chrm + "\tsynteny\texon\t" + str(s2) + "\t" + str(s1) + "\t.\t" + st[0][2] + '\t.\tgene_id "' + name + '"; transcript_id "' + name + '";\n')
						statusex.write(name + '\ttranscript\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tlifted\n')
						status[name] = ["aligned",gene_name]
						g_to_t[gene_name].append(status[name][0])
				elif (len(st) > 1) or (len(en) > 1):
					statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tmultiple_lifted\n')
					status[name] = ["dubious_duplication",gene_name]
					g_to_t[gene_name].append(status[name][0])
				else:
					statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tno_lifted\n')	
					status[name] = ["unaligned",gene_name]
					g_to_t[gene_name].append(status[name][0])
			except:
			 	statusex.write(name + '\tCDS\t' + coords_ex[name][0] + '\t' + coords_ex[name][2][n] + '-' + coords_ex[name][3][n] + '\t' + coords_ex[name][1] + '\t' + sp + '\texon\tno_lifted\n')
			 	status[name] = ["unaligned",gene_name]
			 	g_to_t[gene_name].append(status[name][0])
	
	
	liftex.close()
	statusex.close()
	multiple.close()

	if "human" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	elif "chimp" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	elif "gorilla" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Gorilla_gorilla.gorGor4.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	elif "macaque" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	elif "mouse" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	elif "rat" in sp:
		os.system("gffread liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.gtf -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa -w liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa")
	return status,g_to_t


def analyze_sequence(sp1,sp):
	'''Analyze ORF structures'''
	seqs = {}
	ss = {}
	for line in open("liftover/" + sp1 + "_orfeome_to_" + sp + ".exons.fa"):
		if ">" in line:
			name = line.rstrip("\n").replace(">","")
			seqs[name] = ""
		else:
			seqs[name] = seqs[name] + line.rstrip("\n")

	for name in seqs:
		translated = str(Seq(seqs[name]).translate(cds = False))
		tis = seqs[name][:3]
		try:
			if translated.count("*") > 0:
				stop = round(float(translated.find("*")+1) / len(translated)*100,3)
			else:
				stop = 101 #No stop
			ss[name] = [tis,str(stop),translated]
		except:
			pass

	return ss


def get_exp(sp1, species, species_short, sp_short, psites, gtf):
	'''Extract values to calculate TPMs'''
	tpms = {}
	counts = {}
	mapped = {}
	headers = {}
	for line in open("../../primates_lv/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000
	for line in open("../../ipsc_cm/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000	
	for line in open("../../rodents_lv/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000
		elif "_Ri_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_Ri_")[0] + "_Ri"] = float(line.split("\t")[3]) / 1000000	
	for line in open("../../cardioids/notebooks/total_stats.txt"):
		if "_mR_" in line:
			mapped[line.split("\t")[0].split("/")[-2].replace("_trimmed.29bp_5prime.fq.gz_STARpass1/","").replace("_trimmed.fq.gz_STARpass1/","").split("__")[0].split("_mR_")[0] + "_mR"] = float(line.split("\t")[3]) / 1000000

	#RNA-seq
	headers[sp1.split("_")[0]] = []
	tpms[sp1.split("_")[0]] = {}
	counts[sp1.split("_")[0]] = {}
	for file in glob.glob("../1_mapping/*29bp*_sec/" + sp_short + "*Aligned.sortedByCoord.out.bam"):
		if "pooled" in file:
			continue
		if "iPS" in file:
			continue			

		os.system("coverageBed -split -s -a " + gtf + " -b " + file + " > tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + ".exons.ov.tmp")
		
		print(file)
		headers[sp1.split("_")[0]].append(file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR--" + sp1.split("_")[0])
		all_orfs = {}
		for line in open("tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + ".exons.ov.tmp"):
				orf = line.split("\t")[3]				
				if not orf in tpms[sp1.split("_")[0]]:
					tpms[sp1.split("_")[0]][orf] = []
					counts[sp1.split("_")[0]][orf] = []
				if not orf in all_orfs:
					all_orfs[orf] = 0
				
				all_orfs[orf] = all_orfs[orf] + int(line.split("\t")[-4])
		#os.remove("tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + ".exons.ov.tmp")

		for orf in all_orfs:
			if "ENST" in orf or "merged" in orf:
				le =  ((float(orf.replace("-rand","").split("_")[-1]) - float(orf.replace("-rand","").split("_")[-2])) / 1000)
				tpms[sp1.split("_")[0]][orf].append( ( float(all_orfs[orf]) / le ) / mapped[file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR"] )
				counts[sp1.split("_")[0]][orf].append(float(all_orfs[orf]))
			elif "P1_" in orf:
				tpms[sp1.split("_")[0]][orf].append( ( float(all_orfs[orf]) / (float(orf.split(":")[-1].replace("__v",""))/1000) ) / mapped[file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR"] )
				counts[sp1.split("_")[0]][orf].append(float(all_orfs[orf]))

	for n,sp2 in enumerate(species):
		headers[sp2] = []
		tpms[sp2] = {}
		counts[sp2] = {}
		for file in glob.glob("../1_mapping/*29bp*_sec/" + species_short[n] + "*Aligned.sortedByCoord.out.bam"):
			if "pooled" in file:
				continue
			if "iPS" in file:
				continue
			os.system("coverageBed -split -s -a liftover/" + sp1 + "_orfeome_to_" + sp2 + ".exons.gtf -b " + file + " > tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + "_" + sp2 + ".exons.ov.tmp")

			print(file)
			headers[sp2].append(file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR--" + sp2)
			all_orfs = {}
			for line in open("tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + "_" + sp2 + ".exons.ov.tmp"):
					orf = line.split('transcript_id "')[1].split('"')[0]				
					if not orf in tpms[sp2]:
						tpms[sp2][orf] = []
						counts[sp2][orf] = []
					if not orf in all_orfs:
						all_orfs[orf] = 0
					
					all_orfs[orf] = all_orfs[orf] + int(line.split("\t")[-4])
			#os.remove("tmp/" + file.split("/")[-1].split("__")[0].split("_mR_")[0] + sp1 + "_" + sp2 + ".exons.ov.tmp")

			for orf in all_orfs:
				if "ENST" in orf or "merged" in orf:
					le =  ((float(orf.replace("-rand","").split("_")[-1]) - float(orf.replace("-rand","").split("_")[-2])) / 1000)
					tpms[sp2][orf].append( ( float(all_orfs[orf]) / le ) / mapped[file.split("/")[-1].split("_mR_")[0] + "_mR"] )
					counts[sp2][orf].append(float(all_orfs[orf]))
				elif "P1_" in orf:
					#print(orf)
					tpms[sp2][orf].append( ( float(all_orfs[orf]) / (float(orf.split(":")[-1].replace("__v",""))/1000) ) / mapped[file.split("/")[-1].split("__")[0].split("_mR_")[0] + "_mR"] )
					counts[sp2][orf].append(float(all_orfs[orf]))

	#Ribo-seq
	psites_cov = {}
	for file in glob.glob("../1_mapping/*_sec/" + sp_short + "*_P_sites_uniq_plus.bedgraph"):
		if "pooled" in file:
			continue			

		psites_cov[file] = {}
		os.system("intersectBed -wao -a " + psites + " -b " + file + " | grep -P \"\\t\\+\\t\" > tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + ".psites.ov.tmp")
		os.system("intersectBed -wao -a " + psites + " -b " + file.replace("_plus","_minus") + " | grep -P \"\\t\\-\\t\" >> tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + ".psites.ov.tmp")
		
		print(file)
		headers[sp1.split("_")[0]].append(file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri--" + sp1.split("_")[0])
		all_orfs = {}
		for line in open("tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + ".psites.ov.tmp"):
			if "hq\t+" in line or "p0\t-" in line:
				orf = line.split("\t")[3]
				if not orf in tpms[sp1.split("_")[0]]:
					continue
				if not orf in all_orfs:
					all_orfs[orf] = 0
					psites_cov[file][orf] = {}
				if not line.split("\t")[1] in psites_cov[file][orf]:
					psites_cov[file][orf][line.split("\t")[1]] = 0
				if line.split("\t")[1] == line.split("\t")[7]:
					all_orfs[orf] = all_orfs[orf] + int(line.split("\t")[-2].replace(".","0"))
					psites_cov[file][orf][line.split("\t")[1]] = psites_cov[file][orf][line.split("\t")[1]] + int(line.split("\t")[-2].replace(".","0"))
		#os.remove("tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1.split("_")[0] + ".psites.ov.tmp")

		for orf in tpms[sp1.split("_")[0]]:
			if not orf in all_orfs:
				tpms[sp1.split("_")[0]][orf].append(0)
				counts[sp1.split("_")[0]][orf].append(0)
				continue
			if "ENST" in orf or "merged" in orf:
				le =  ((float(orf.replace("-rand","").split("_")[-1]) - float(orf.replace("-rand","").split("_")[-2])) / 1000)
				tpms[sp1.split("_")[0]][orf].append( ( float(all_orfs[orf]) / le ) / mapped[file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri"] )
				counts[sp1.split("_")[0]][orf].append(float(all_orfs[orf]))
			elif "P1_" in orf:
				tpms[sp1.split("_")[0]][orf].append( ( float(all_orfs[orf]) / (float(orf.split(":")[-1].replace("__v",""))/1000) ) / mapped[file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri"] )
				counts[sp1.split("_")[0]][orf].append(float(all_orfs[orf]))	

	for n,sp2 in enumerate(species):
		for file in glob.glob("../1_mapping/*_sec/" + species_short[n] + "*_P_sites_uniq_*plus*.bedgraph"):
			if "pooled" in file:
				continue
			if "AssumboCG06" in file:
				continue
			os.system("intersectBed -wao -b liftover/" + sp1 + "_orfeome_to_" + sp2 + ".psites.gtf -a " + file + " | grep -P \"\\t\\+\\t\" > tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + "_" + sp2 + ".psites.ov.tmp") #P-sites are generated in script 2
			os.system("intersectBed -wao -b liftover/" + sp1 + "_orfeome_to_" + sp2 + ".psites.gtf -a " + file.replace("_plus","_minus") + " | grep -P \"\\t\\-\\t\" >> tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + "_" + sp2 + ".psites.ov.tmp")

			print(file)
			headers[sp2].append(file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri--" + sp2)
			all_orfs = {}
			for line in open("tmp/" + file.split("/")[-1].split("__")[0].split("_Ri_")[0] + sp1 + "_" + sp2 + ".psites.ov.tmp"):
					orf = line.split('orf_id "')[1].split('"')[0]
					if not orf in tpms[sp2]:
						continue
					if not orf in all_orfs:
						all_orfs[orf] = 0
					if line.split("\t")[2] == line.split("\t")[7]:
						all_orfs[orf] = all_orfs[orf] + int(line.split("\t")[3].replace(".","0"))
			#os.remove("tmp/" + file.split("/")[-1].split("_Ri_")[0] + sp2 + ".psites.ov.tmp")

			for orf in tpms[sp2]:
				if not orf in all_orfs:
					tpms[sp2][orf].append(0)
					counts[sp2][orf].append(0)
					continue
				if "ENST" in orf or "merged" in orf:
					le =  ((float(orf.replace("-rand","").split("_")[-1]) - float(orf.replace("-rand","").split("_")[-2])) / 1000)
					tpms[sp2][orf].append( ( float(all_orfs[orf]) / le ) / mapped[file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri"] )
					counts[sp2][orf].append(float(all_orfs[orf]))
				elif "P1_" in orf:
					tpms[sp2][orf].append( ( float(all_orfs[orf]) / (float(orf.split(":")[-1].replace("__v",""))/1000) ) / mapped[file.split("/")[-1].split("__")[0].split("_Ri_")[0] + "_Ri"] )
					counts[sp2][orf].append(float(all_orfs[orf]))

	out = open("out/" + sp1 + ".orf_psites.tsv","w+")
	out.write("orf_id\tfile\tcat\tpos\tcov\n")
	for file in psites_cov:
		for orf in psites_cov[file]:
			for pos in psites_cov[file][orf]:
				out.write(orf + "\t" + file + "\tunique\t" + str(pos) + "\t" + str(psites_cov[file][orf][pos]) + "\n")
	
	out.close()
	return (tpms,counts,headers)


coords_ex = get_coords(gtf)
statuses = {}
g_to_t = {}
strs = {}
candidates = {}
for sp2 in species:
	(lo,lo2,eq) = load_chain(sp1.split("_")[0],sp2)
	(statuses[sp2],g_to_t[sp2]) = run_liftover(coords_ex,lo,lo2,eq,sp1,sp2)
	strs[sp2] = analyze_sequence(sp1,sp2)

(tpms,counts,headers) = get_exp(sp1,species,species_short,sp_short,psites,gtf)


out3 = open("out/" + sp1 + ".orf_region.speciesexpression.tsv","w+")
out3.write("orf_id\tsp2\tsample\tTIS\ttranslated\tstr\tstatus\tTPM\traw_counts\n")
for cat in [sp1.split("_")[0]] + species:
	for orf in tpms[cat]:
		for n,header in enumerate(headers[cat]):
			if cat != sp1.split("_")[0]:
				if not orf in strs[cat]:
					strs[cat][orf] = ["XXX","-1","XXX"]
				if not orf in statuses[cat]:
					statuses[cat][orf] = ["unaligned","unknown"]
				try:
					out3.write(orf + "\t" + cat + "\t" + header + "\t" + "\t".join(strs[cat][orf]) + "\t" + statuses[cat][orf][0] + "\t" + str(tpms[cat][orf][n]) + "\t" + str(counts[cat][orf][n]) + "\n")
				except:
					out3.write(orf + "\t" + cat + "\t" + header + "\t" + "\t".join(strs[cat][orf]) + "\t" + statuses[cat][orf][0] + "\t0\t0\n")
			else:
				try:
					out3.write(orf + "\t" + cat + "\t" + header + "\tnone\tnone\tnone\tnone\t" + str(tpms[cat][orf][n]) + "\t" + str(counts[cat][orf][n]) + "\n")
				except:
					out3.write(orf + "\t" + cat + "\t" + header + "\tnone\tnone\tnone\tnone\t0\t0\t\n")				

out3.close()

exit(0)
