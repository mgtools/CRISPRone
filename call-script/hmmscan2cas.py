from __future__ import print_function
#!/usr/bin/env python

#infer putative cas genes from hmmscan results
#latest update Nov 2020 YY, compatible with python 3

import sys
import operator
import os
import re

def get_HMGI_gff(srsid):
	gff = "tmp/" + srsid + "/" + srsid + ".gff"
	return gff 

#add any (for local version) Ye, Nov 2016
def get_any_gff(fdir, srsid): 
        gff = fdir + "/" + srsid + ".gff"
        return gff

def get_bacteria_faa(folder="/home/data1/groups/quzhang/Streptococcus-pyogenes-thermophilus-db-Dec16/Bacteria/"):
	bacteria_faa = {}
	for dir in os.listdir(folder):
		thisdir = folder + dir
		if os.path.isdir(thisdir):
			for afile in os.listdir(thisdir):
				if ".faa" in afile:
					subs = afile.split(".")
					bacteria_faa[subs[0]] = thisdir

	print ("total faa from " + folder + str(len(bacteria_faa)))
	return bacteria_faa

def get_complete_gff(faaid, bacteria_faa):
	#print "faaid", faaid
	file = bacteria_faa[faaid] + "/" + faaid + ".gff"
	return file

def get_draft_gff(faaid, bacteria_faa):
	#print "faaid", faaid
	file = bacteria_faa[faaid] + "/" + faaid + ".gff"
	return file

def get_hmminfo(gafile):
	hmminfo = {}
	infile = open(gafile, "r")
	print("get " + gafile)
	for aline in infile:
		if aline[0] == '#':
			continue
		subs = aline.strip().split("\t")
		
		hmminfo[subs[0]] = [float(subs[1]), subs[2]]
		print("hmm " + subs[0] + " " + subs[1] + " " + subs[2])
	infile.close()
	return hmminfo

def getcrispr(crisprhmmfile):
	infile = open(crisprhmmfile, "r")
	#print "get", crisprhmmfile
	crisprid, crisprdes = [], []
	for aline in infile:
		if aline[0] == '#':
			continue
		subs = aline.strip().split("\t")
		crisprid.append(subs[0])
		crisprdes.append(subs[2])
	infile.close()
	return crisprid, crisprdes
		
def cas_contig(hmmscanfile, crisprhmm, hmminfo, gfffile, ifhmmsearch=False):
	infile = open(hmmscanfile, "r")
	all, valid = 0, 0
	cascontig = []
	genelist = []
	gene2contig = []
	#get significant cas genes
	for aline in infile:
		if aline[0] == '#':
			continue	
		subs = aline.split()
		if ifhmmsearch:
			seqid, famacc = subs[0], subs[3] 
		else:
			seqid, famacc = subs[3], subs[1]
		if famacc not in crisprhmm:
			continue
		ga = hmminfo[famacc][0]
		score = float(subs[13])
		#print "hmm", subs[0], "ga-cut", ga, "score", score
		#raw_input("type enter to continue..")
		if score >= ga:
			#subs2 = seqid
			#subs2 = seqid.split(";")
			#if len(subs2) > 1:
			#	contigid = subs2[0]
			#	geneid = subs2[1]
			#	if contigid not in cascontig:
			#		cascontig.append(contigid)
			#	if geneid not in genelist:
			#		genelist.append(geneid)
			#		gene2contig.append(cascontig.index(contigid))
			#else:
			geneid = seqid 
			if geneid not in genelist:
				genelist.append(geneid)
				gene2contig.append(-1)
	infile.close()
	#print hmmscanfile, "total casgene", len(genelist), "total contig", len(cascontig)

	#get contigs that contain the significant cas genes
	if genelist and gfffile and (not cascontig):
		infile = open(gfffile, "r")
		for aline in infile:
			if "##FASTA" in aline:
				break
			if aline[0] == '#':
				continue 
			subs = aline.split("\t")
			if subs[2] == 'CDS':
				g = re.findall(r'ID=([^\;]+);Name=([^\;]+)',subs[-1]) 
				if g:
					id = g[0][1] #use Name
				else:
					subs2 = subs[-1][3:].split(";")
					id = subs2[0]
				if (id in genelist):
					contigid = subs[0]
					if contigid not in cascontig:
						contigidx = len(cascontig)
						cascontig.append(contigid)
					else:
						contigidx = cascontig.index(contigid) 
					gene2contig[genelist.index(id)] = contigidx 
		infile.close()
		
	#print "total cas-contig candidates", len(cascontig)
	#for gidx in range(len(genelist)):
	#	print "cas-gene", genelist[gidx], "contig", cascontig[gene2contig[gidx]]

	return (cascontig, genelist, gene2contig)

def contig_gene(cascontig, gfffile):
	idx = -1 
	lastcontig = ""
	contiggene = []

	infile = open(gfffile, "r")
	for aline in infile:
		if "##FASTA" in aline:
			break
		if aline[0] == '#':
			continue
		subs = aline.split("\t")
		if subs[2] != "CDS":
			continue
		if subs[0] != lastcontig:
			valid = 0
			if (not contiggene) or contiggene[-1]:
				idx += 1
				if idx >= len(cascontig):
					break
				currcontig = cascontig[idx]	
				contiggene.append([])
			if subs[0] == currcontig:
				valid = 1
		if valid:
			g = re.findall(r'ID=([^\;]+);Name=([^\;]+)',subs[-1]) 
			if g:
				id = g[0][1] #use Name
			else:
				subs2 = subs[-1][3:].split(";")
				id = subs2[0]
			if (not contiggene[-1]) or (id != contiggene[-1][-1][0]):
				contiggene[-1].append([id, int(subs[3]), int(subs[4]), subs[6], "", "", ""])
		lastcontig = subs[0]	
	infile.close()

	return contiggene
		

def cas_genes(hmmscanfile, crisprhmm, hmminfo, cascontig, contiggene, flank=3, ifhmmsearch=False, clusteronly=True, outfile="", filtered=""):
	#get hits for all genes 
	infile = open(hmmscanfile, "r")
	for aline in infile:
		if aline[0] == '#':
			continue
		subs = aline.split()
		if ifhmmsearch:
			hmmid, geneid, score = subs[4], subs[0], float(subs[13])
		else:
			#hmmid, geneid, score = subs[1], subs[3], float(subs[13])
						      #acc
			hmmid, geneid, score = subs[0], subs[3], float(subs[13])
						      #name
		#g = re.findall(r'zzzzz\|([^\|]+)', geneid) 
		#if g:
		#	geneid = g[0]
		if hmmid not in hmminfo:
			sys.exit("hmm not found '" + hmmid + "' line " + aline)
		ga = hmminfo[hmmid][0]
		if score >= ga:
			#print "score",score," ga",ga
			sig = 2 
		else:
			sig = 0 
		for cidx in range(len(cascontig)):
			for gidx in range(len(contiggene[cidx])):
				if geneid == contiggene[cidx][gidx][0]:
					if not contiggene[cidx][gidx][4]: #only keep the best one
						contiggene[cidx][gidx][4] = hmmid
						contiggene[cidx][gidx][5] = score
						contiggene[cidx][gidx][6] = sig
					else:
						if sig and (not contiggene[cidx][gidx][6]):
							contiggene[cidx][gidx][4] = hmmid
							contiggene[cidx][gidx][5] = score
							contiggene[cidx][gidx][6] = sig
						elif sig == contiggene[cidx][gidx][6]:
							#don't replace the sig one with insig but higher score
							if score > contiggene[cidx][gidx][5]: 
								contiggene[cidx][gidx][4] = hmmid
								contiggene[cidx][gidx][5] = score
								contiggene[cidx][gidx][6] = sig
	infile.close()

	#recruit cas genes with weak similarities & remove isolated cas genes
	casgene_all = []
	casunk_all = []
	tagnum = 0
	for cidx in range(len(cascontig)):
		allgene_sorted = sorted(contiggene[cidx], key=operator.itemgetter(1))
		isolate = [0] * len(allgene_sorted)
		#recruit weak cas in the neighborhood
		oldcas = 0
		casgene = []
		casunk = []
		#print "contig: ", cascontig[cidx]
		#print "contig-gene", contiggene[cidx]
		for gidx in range(len(allgene_sorted)):
			if (allgene_sorted[gidx][-1] == 2) and (allgene_sorted[gidx][4] in crisprhmm):
				casgene.append(gidx)
				#print "seed", allgene_sorted[gidx]
		#print "seed cas", len(casgene)
		if len(casgene) == 0:
			print("Error: cascontig " + cascontig[cidx] + " has NO casgene!!")
			#raw_input("type enter to continue..")
			sys.exit()
		iter = 0
		while (len(casgene) > oldcas):
			oldcas = len(casgene)
			oldgene = casgene[:]
			for bait in oldgene:
				#for gidx in range(max(0, bait - flank), min(bait + flank, len(allgene_sorted))):
				#bug, Ye, May 21, 2013
				for gidx in range(max(0, bait - flank), min(bait + flank + 1, len(allgene_sorted))):
					if (allgene_sorted[gidx][-1] == 0) and (allgene_sorted[gidx][4] in crisprhmm):	
						casgene.append(gidx)
						allgene_sorted[gidx][-1] = 1
			iter += 1
			#print "iter", iter, "casgene", len(casgene)
		#print "total casgene", len(casgene)
		casgene_all.extend(casgene)
		#raw_input("type enter to continue..")

		unkstart = -1
		for gidx in range(1, len(allgene_sorted) - 1):
			if (not allgene_sorted[gidx][-1]) and (allgene_sorted[gidx-1][-1]):
				unkstart = gidx
			if (not allgene_sorted[gidx][-1]) and (allgene_sorted[gidx+1][-1]):	
				unkend = gidx
				if (unkstart != -1) and (unkend - unkstart + 1 <= flank):
					for idx2 in range(unkstart, unkend + 1):
						casunk.append(allgene_sorted[idx2][0])
						allgene_sorted[idx2][-1] = 1
						#print "CRISPR-associated UNK:", casunk[-1]
						#raw_input("found new CAS protein..")
				unkstart = -1
		if casunk:
			casunk_all.extend(casunk)

		#tag isolated cas genes
		if clusteronly:
			for gidx in range(len(allgene_sorted)):
				hmmid = allgene_sorted[gidx][4]
				if hmmid:
					clustersize = 0
					for b in range(max(0, gidx - 10), min(len(allgene_sorted), gidx + 10)):
						if allgene_sorted[b][4]:
							clustersize += 1
					if clustersize < 2:
						isolate[gidx] = 1
						if allgene_sorted[gidx][6]:
							tagnum += 1
		
		if outfile:
			outfile.write("##" + cascontig[cidx] + " " + str(len(allgene_sorted)) + " derived-fm " + hmmscanfile + "\n")
		if filtered and casgene:
			filtered.write("##" + cascontig[cidx] + " " + str(len(casgene) + len(casunk)) + " derived-fm " + hmmscanfile + "\n")
		subs = hmmscanfile.split("/")
		subs2 = subs[-1].split("-")
		sampleid = subs2[0]
		for gidx in range(len(allgene_sorted)):
			hmmid = allgene_sorted[gidx][4]
			if hmmid:
				hmmacc = allgene_sorted[gidx][4]
				hmmdes = hmminfo[hmmacc][1]
				subs = hmmdes.split(":")
				shortdes = hmmacc + ":" + subs[0]
			else:
				shortdes = "unk"
			if outfile:
				outfile.write(" ".join([allgene_sorted[gidx][0], sampleid, cascontig[cidx], str(allgene_sorted[gidx][1]), str(allgene_sorted[gidx][2]), str(allgene_sorted[gidx][3]), shortdes, str(allgene_sorted[gidx][6])]) + "\n")
			if filtered and allgene_sorted[gidx][6] and (not isolate[gidx]):
				filtered.write(" ".join([str(allgene_sorted[gidx][0]), sampleid, cascontig[cidx], str(allgene_sorted[gidx][1]), str(allgene_sorted[gidx][2]), str(allgene_sorted[gidx][3]), shortdes, str(allgene_sorted[gidx][6])]) + "\n")
	print(hmmscanfile + " all known-cas " + str(len(casgene_all)) + " unknown-cas " + str(len(casunk_all)) + " isolate-cas-tagged " + str(tagnum))
	return casgene_all, casunk_all

def main():
	gadir = "/u/yye/CRISPRpk/CRISPRone/scripts/"
	gafile0 = gadir + "local/cas-db/cas-CRISPR.GA"
	crisprfile0 = gadir + "local/cas-db/cas-CRISPR.GA"
	#both are GA files
	gafile, crisprhmmfile = gafile0, crisprfile0
	hmmsearchfile, hmmscanfile, gfffile = "", "", ""
	collection = "HMGI"
	outfile, outfiltered = "summary.txt", "summary-filtered.txt"
	clusteronly = True
	fdir = "tmp";
	for idx in range(len(sys.argv)):
		if (sys.argv[idx] == "-gafile") and (len(sys.argv) > idx + 1):
			crisprhmmfile = gafile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-crisprhmm") and (len(sys.argv) > idx + 1):
			crisprhmmfile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-collection") and (len(sys.argv) > idx + 1):
			collection = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-gff") and (len(sys.argv) > idx + 1):
			gfffile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-hmmscan") and (len(sys.argv) > idx + 1):
			hmmscanfile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-hmmsearch") and (len(sys.argv) > idx + 1):
			hmmsearchfile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-outfile") and (len(sys.argv) > idx + 1):
			outfile = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-outfiltered") and (len(sys.argv) > idx + 1):
			outfiltered = sys.argv[idx + 1]
		elif (sys.argv[idx] == "-isolate"):
			clusteronly = False	
	if not (gafile and (hmmscanfile or hmmsearchfile) and crisprhmmfile):
		print(sys.argv[0] + " -gafile filename -crisprhmm <-collection HMGI/complete/draft> <-hmmscan filename-or-folder>/<-hmmsearch filename-or-folder> -newfile filename")
		sys.exit()

	if hmmsearchfile:
		file0 = hmmsearchfile
		ifhmmsearch = True
	else:
		file0 = hmmscanfile
		ifhmmsearch=False
	allfiles = []
	if os.path.isdir(file0):
		fdir = file0
		files = os.listdir(file0)
		for afile in files:
			if ".domtblout" in afile:
				allfiles.append(file0 + "/" + afile)
	else:
		allfiles.append(file0)

	out = open(outfile, "w")
	out2 = open(outfiltered, "w")
	out.write("#" + " ".join(sys.argv) + "\n")
	out2.write("#" + " ".join(sys.argv) + "\n")
	hmminfo = get_hmminfo(gafile)
	crisprid, crisprdes = getcrispr(crisprhmmfile)
	casgene_add, casunk_add = 0, 0

	print("collection " + collection)
	if collection == 'complete':
		bacteria_faa = get_bacteria_faa()
	elif collection == 'draft':
		bacteria_faa = get_bacteria_faa(folder="/home/data1/groups/quzhang/Streptococcus-pyogenes-thermophilus-db-Dec16/Bacteria_draft/")

	for afile in allfiles:
		subs0 = afile.split("/")
		subs = subs0[-1].split("-")
		if collection == 'HMGI':
			gfffile = get_HMGI_gff(subs[0])
		elif collection == 'any':
			gfffile = get_any_gff(fdir, subs[0])
		elif collection == 'complete':
			gfffile = get_complete_gff(subs[0], bacteria_faa)
		elif collection == 'draft':
			gfffile = get_draft_gff(subs[0], bacteria_faa)
		elif not gfffile:
			print("Error: gfffile not found for " + subs[0])
			sys.exit()
		if not (os.path.exists(gfffile) and os.path.exists(afile)):
			continue
		cascontig, casgene, gene2contig = cas_contig(afile, crisprid, hmminfo, gfffile, ifhmmsearch=ifhmmsearch)
		contig_allgene = contig_gene(cascontig, gfffile)
		casgene_all, casunk_all = cas_genes(afile, crisprid, hmminfo, cascontig, contig_allgene, ifhmmsearch=ifhmmsearch, outfile=out, filtered=out2, clusteronly=clusteronly)
		casgene_add += len(casgene_all)
		casunk_add += len(casunk_all)
	out.close()
	out2.close()

	print("total casgene " + str(casgene_add) + " unknown-cas " + str(casunk_add))
	print("results saved to " + outfile + " and " + outfiltered)

if __name__ == '__main__':
	main()
