import sys
import os

class trans_object:
		def __init__(self, chrm, gene, strand, start, end, code):
				self.chrm = chrm
				self.gene = gene
				self.strand = strand
				self.start = start
				self.end = end
				self.code = code

def parse_gtf(gtf,field):
	#Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	for line in open(gtf):
		if not "\t" + field + "\t" in line:
			continue
		t_name = line.split('transcript_id "')[1].split('"')[0]
		g_name = line.split('gene_id "')[1].split('"')[0]

		if 'gene_biotype' in line:
			biot = line.split("\t")[1]
		else:
			biot = "unknown"
		
		trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[],biot))
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans

def parse_bed(bed):
	#Read a bed and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	for line in open(bed):
		t_name = line.split('\t')[4]
		
		trans.setdefault(t_name,trans_object(line.split("\t")[0],line.split("\t")[4],line.split("\t")[5].rstrip('\n'),[],[],line.split("\t")[3]))
		trans[t_name].start.append(int(line.split("\t")[1]))
		trans[t_name].end.append(int(line.split("\t")[2]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans


gtff = sys.argv[1]
annot = sys.argv[2]
out = open(gtff + "_psites_plus_partial.bed","w+")

if gtff.endswith("bed"):
	gtf = parse_bed(gtff)
elif gtff.endswith("gtf"):
	gtf = parse_gtf(gtff,"CDS")

status = {}
if (annot == "yes") and (gtff.endswith("gtf")): #Remove uncomplete proteins 
		for line in open(gtff):
			if "\tCDS\t" in line:
				if not line.split('transcript_id "')[1].split('"')[0] in status:
					status[line.split('transcript_id "')[1].split('"')[0]] =  0				
			if "\tstart_codon\t" in line:
				if not line.split('transcript_id "')[1].split('"')[0] in status:
					status[line.split('transcript_id "')[1].split('"')[0]] =  1
				elif status[line.split('transcript_id "')[1].split('"')[0]] == 2:
					status[line.split('transcript_id "')[1].split('"')[0]] = 3
				else:
					status[line.split('transcript_id "')[1].split('"')[0]] =  1
			elif "\tstop_codon\t" in line:
				if not line.split('transcript_id "')[1].split('"')[0] in status:
					status[line.split('transcript_id "')[1].split('"')[0]] = 2
				elif status[line.split('transcript_id "')[1].split('"')[0]] == 1:
					status[line.split('transcript_id "')[1].split('"')[0]] = 3
				else:
					status[line.split('transcript_id "')[1].split('"')[0]] = 2

for orf in gtf:
	if orf in status:
		if status[orf] == 0:
			print(orf + " excluded")
			continue
		if status[orf] == 1:
			stat = 1
		elif status[orf] == 2:
			stat = 2
		elif status[orf] == 3:
			stat = 3
	else:
		stat = 3
	i = 0
	x1 = 0
	x2 = 0
	x3 = 0
	y1 = 0
	y2 = 0
	y3 = 0
	for n,ex in enumerate(gtf[orf].start):
		for co in range(int(gtf[orf].start[n]),int(gtf[orf].end[n])+1):
			if (gtf[orf].strand == "+") and ((stat == 1) or (stat == 3)):
				if (i % 3) == 2:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp0\t" + gtf[orf].strand + "\n")
					y3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST2\t" + gtf[orf].strand + "\n"
					if x3 == 0:
						x3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG2\t" + gtf[orf].strand + "\n"
				elif (i % 3) == 0:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp1\t" + gtf[orf].strand + "\n")
					y1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST0\t" + gtf[orf].strand + "\n"
					if x1 == 0:
						x1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG0\t" + gtf[orf].strand + "\n"
				else:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp2\t" + gtf[orf].strand + "\n")
					y2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST1\t" + gtf[orf].strand + "\n"
					if x2 == 0:
						x2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG1\t" + gtf[orf].strand + "\n"
			

			elif (gtf[orf].strand == "+") and (stat == 2):
				l = (sum(gtf[orf].end) - sum(gtf[orf].start) + len(gtf[orf].end)) % 3
				if (i % 3) == (l+2) % 3:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp0\t" + gtf[orf].strand + "\n")
					y3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST2\t" + gtf[orf].strand + "\n"
					if x3 == 0:
						x3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG2\t" + gtf[orf].strand + "\n"
				elif (i % 3) == l % 3:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp1\t" + gtf[orf].strand + "\n")
					y1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST0\t" + gtf[orf].strand + "\n"
					if x1 == 0:
						x1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG0\t" + gtf[orf].strand + "\n"
				else:
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp2\t" + gtf[orf].strand + "\n")
					y2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST1\t" + gtf[orf].strand + "\n"
					if x2 == 0:
						x2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG1\t" + gtf[orf].strand + "\n"


			elif (gtf[orf].strand == "-") and ((stat == 2) or (stat == 3)):
				if (i % 3) == 0:
					x1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG2\t" + gtf[orf].strand + "\n"
					if y1 == 0:
						y1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST2\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp0\t" + gtf[orf].strand + "\n")
				elif (i % 3) == 1:
					x2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG1\t" + gtf[orf].strand + "\n"
					if y2 == 0:
						y2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST1\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp1\t" + gtf[orf].strand + "\n")
				else:
					x3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG0\t" + gtf[orf].strand + "\n"
					if y3 == 0:
						y3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST0\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp2\t" + gtf[orf].strand + "\n")


			elif (gtf[orf].strand == "-") and (stat == 1):
				l = (sum(gtf[orf].end) - sum(gtf[orf].start) + len(gtf[orf].end)) % 3
				if (i % 3) == l % 3:
					x1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG2\t" + gtf[orf].strand + "\n"
					if y1 == 0:
						y1 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST2\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp0\t" + gtf[orf].strand + "\n")
				elif (i % 3) == (l+1) % 3:
					x2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG1\t" + gtf[orf].strand + "\n"
					if y2 == 0:
						y2 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST1\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp1\t" + gtf[orf].strand + "\n")
				else:
					x3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpATG0\t" + gtf[orf].strand + "\n"
					if y3 == 0:
						y3 = gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tpST0\t" + gtf[orf].strand + "\n"
					out.write(gtf[orf].chrm + "\t" + str(co) + "\t" + str(co) + "\t" + orf + "\tp2\t" + gtf[orf].strand + "\n")

			i += 1
	out.write(x1)
	out.write(x2)
	out.write(x3)
	out.write(y1)
	out.write(y2)
	out.write(y3)
out.close()
