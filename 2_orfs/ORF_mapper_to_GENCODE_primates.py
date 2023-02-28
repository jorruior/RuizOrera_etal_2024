#!/usr/bin/env python
import sys
import string
import subprocess
import os
import random
import string
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from optparse import OptionParser

import functions

__author__ = "Jorge Ruiz-Orera"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com"


def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print("Error: " + str(s) + " argument not given\n")
		exit()

def check_file (file_str):
	'''Check if input really exists'''
	try:
		open("%s" %file_str)
	except:
		print("Error: " + file_str + " input not found\n")
		exit()

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-d","--input_dir",action="store",dest="folder",help="(Required and modified) Name of the tag")
	parser.add_option("-f","--input_fasta",action="store",dest="orfs_fa_file",help="(Required) File with all translated candidate ORFs. A FASTA can be generated from a BED file using the script 'bed1_to_fasta.sh'")
	parser.add_option("-b","--bed",action="store",dest="orfs_bed_file",help="(Required) File with 1-based BED coordinates of all translated candidate ORFs. ORF names should match in both fasta and bed files. If -a is activated, the BED coordinates will be written into this file.")
	parser.add_option("-o","--output",action="store",dest="out_name",default="input_bed",help="Tag name for the generated output files. Default: BED filename")	
	parser.add_option("-l","--min_len_cutoff",action="store",type=int,dest="len_cutoff",default="16",help="Minimum ORF length threshold in amino acids (without stop codon). default=16")
	parser.add_option("-L","--max_len_cutoff",action="store",type=int,dest="max_len_cutoff",default="999999999999",help="Maximum ORF length threshold in amino acids (without stop codon). default=none")
	parser.add_option("-c","--collapse_cutoff",action="store",type=float,dest="col_thr",default="0.9",help="Minimum required fraction to collapse ORFs with similar stretches of overlapping amino-acid sequences. default=0.9")
	parser.add_option("-m","--collapse_method",action="store",dest="method",default="longest_string",help="Method to cluster ORF variants. 'longest_string' will collapse ORFs if the longest shared string is above -c threshold (default). 'psite_overlap' is a slower method that collapses ORFs if the fraction of shared psites is above -c threshold")
	parser.add_option("-a","--make_annot_bed",action="store",dest="calculate_coordinates",default="no",help="If ORF BED file is not available, generate it from the FASTA file. WARNING: BED file will be writen to the filename given by -b/--bed. (ATG/NTG/XTG/no, default = no)")
	parser.add_option("--multiple",action="store",dest="mult",default="yes",help="If -a option is activated, include ORFs that map to multiple genomic coordinates. (yes/no, default = yes)")
	parser.add_option("-g","--genomic",action="store",dest="genomic",default="none",help="If -a option is activated, this optional argument uses a BED file with ORF genomic coordinates to HELP mapping ORF sequences to the correct exon-intron positions. ORFs that cannot be mapped to these genomic regions will be mapped to alternative positions if possible. Use 'none' if no file is given or -a is not activated. (default = none)")
	parser.add_option("-G","--force_genomic",action="store",dest="fgenomic",default="none",help="If -a option is activated, this optional argument uses a BED file with ORF genomic coordinates to FORCE mapping ORF sequences to the correct exon-intron positions. ORFs that cannot be mapped to these genomic regions will be discarded. Use 'none' if no file is given or -a is not activated. (default = none)")
	parser.add_option("-C","--add_cds",action="store",dest="cds_cases",default="no",help="If 'yes', include CDS and pseudogenes in the output GTF. (default = 'no')")
	parser.add_option("-s","--stringtie",action="store",dest="stringtie",help="(Required and modified) GTF with transcript expression from Stringtie")


	(opt,args)=parser.parse_args()

	check_arg(opt.folder,"--folder")
	check_arg(opt.orfs_fa_file,"--input_fasta")
	check_arg(opt.orfs_bed_file,"--input_bed")
	check_arg(opt.stringtie,"--stringtie")
	check_file(opt.orfs_fa_file)
	if opt.calculate_coordinates != "no":
		print("A new BED file will be generated in " + opt.orfs_bed_file)
	else:
		check_file(opt.orfs_bed_file)
	folder = opt.folder
	orfs_fa_file = opt.orfs_fa_file
	orfs_bed_file = opt.orfs_bed_file
	len_cutoff = opt.len_cutoff
	max_len_cutoff = opt.max_len_cutoff
	col_thr = opt.col_thr
	method = opt.method
	calculate_coordinates = opt.calculate_coordinates
	mult = opt.mult
	genomic = opt.genomic
	fgenomic = opt.fgenomic
	cds_cases = opt.cds_cases
	if opt.out_name == "input_bed":
		out_name = orfs_bed_file
	else:
		out_name = opt.out_name
	if not calculate_coordinates in ("ATG","NTG","XTG","no"):
		print("Error: " + method + " is not a valid -a argument\n")	
		exit()
	if not method in ("longest_string","psite_overlap"):
		print("Error: " + calculate_coordinates + " is not a valid -m argument\n")		
		exit()
	if not opt.mult in ("yes","no"):
		print("Error: " + opt.mult + " is not a valid --multiple argument\n")
		exit()
	if (calculate_coordinates == "no") and (mult != "yes"):
		print("Error: " + opt.mult + " is not a valid --multiple argument when -a is not enabled\n")
		exit()
	if (calculate_coordinates == "no") and (genomic != "none"):
		print("Error: " + opt.genomic + " is not a valid --multiple argument when -a is not enabled\n")
		exit()
	if (calculate_coordinates == "no") and (fgenomic != "none"):
		print("Error: " + opt.fgenomic + " is not a valid --multiple argument when -a is not enabled\n")
		exit()
	elif (genomic != "none") and (fgenomic != "none"):
		print("Error: -g and -G are mutually exclusive parameters, please use only one of them\n")
		exit()		
	elif genomic != "none":
		check_file(genomic)
	elif fgenomic != "none":
		check_file(fgenomic)	

	#Annotation files
	transcriptome_fa_file = "../1_mapping/annotation/" + folder + ".fixed.fa" #include both mRNA and ncRNA, only transcript ID in header, e.g. cat Ens103/Homo_sapiens.GRCh38.cdna.all.fa Ens103/Homo_sapiens.GRCh38.ncrna.fa | cut -d"." -f1,1 > prueba; mv prueba Ens103/Homo_sapiens.GRCh38.trans.fa
	transcriptome_gtf_file =  "../1_mapping/annotation/" + folder + ".fixed.sorted.gtf" #sort -k1,1 -k4,4n -k5,5n Ens103/Homo_sapiens.GRCh38.gtf > Ens103/Homo_sapiens.GRCh38.sorted.gtf

	check_file(transcriptome_fa_file)
	check_file(transcriptome_gtf_file)
	check_file(opt.stringtie)

	try:
   		os.mkdir(folder + "/tmp/")
	except:
		pass 

	(orfs_fa,transcriptome_fa) = functions.load_fasta(orfs_fa_file,transcriptome_fa_file)
	gtf = functions.parse_gtf(transcriptome_gtf_file,"exon")
	if calculate_coordinates != "no":
		functions.make_bed(orfs_fa,transcriptome_fa,gtf,len_cutoff,max_len_cutoff,calculate_coordinates,orfs_bed_file,out_name,mult,genomic,fgenomic)
	supp = functions.read_support(opt.stringtie)
	(overlaps,overlaps_cds,other_overlaps,total_studies,seed) = functions.insersect_orf_gtf(orfs_bed_file,transcriptome_gtf_file,folder)
	other_overlaps = functions.pseudo_or_cds_ov(orfs_bed_file,transcriptome_gtf_file,other_overlaps,folder,seed)
	(candidates,trans_orfs,coord_psites) = functions.orf_tags(overlaps,overlaps_cds,orfs_fa,transcriptome_fa,gtf,len_cutoff,max_len_cutoff,folder,seed)
	(exc,variants,variants_names,datasets) = functions.exclude_variants(orfs_fa,trans_orfs,col_thr,candidates,method,coord_psites,folder,seed)
	functions.write_output(orfs_fa_file,orfs_bed_file,candidates,exc,variants,variants_names,datasets,supp,gtf,transcriptome_fa,len_cutoff,max_len_cutoff,col_thr,total_studies,other_overlaps,folder,out_name,method,seed,genomic,fgenomic,cds_cases)

if __name__ == '__main__':
	main()

exit(0)

#Add maximum length
#Upload again


