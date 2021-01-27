#!/usr/bin/env python

# verify types of CRISPR-Cas systems / verify using predictions
import sys
import os
import operator
import subprocess

cmd = sys.argv[0]
resultdir, resultfile = "", ""
genome = True
for idx in range(len(sys.argv)):
	if (sys.argv[idx] == "-d") and (len(sys.argv) > idx + 1):
		resultdir = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-f") and (len(sys.argv) > idx + 1):
		resultfile = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-g") and (len(sys.argv) > idx + 1):
		genome = sys.argv[idx + 1]
cmddir = os.path.dirname(cmd)
if not cmddir:
	modulefile = "../local/cas-db/interference-module.txt"
	cdhitest = "../bin/cd-hit-v4.6.1-2012-08-27/cd-hit-est"
else:
	modulefile = cmddir + "/../local/cas-db/interference-module.txt"
	cdhitest = cmddir + "/../bin/cd-hit-v4.6.1-2012-08-27/cd-hit-est"

plotfiles = []
if resultfile:
	plotfiles.append(resultfile)
if resultdir:
	for afile in os.listdir(resultdir):
	    if afile.endswith("plot.txt"):
		plotfiles.append(resultdir + "/" + afile)
if len(plotfiles) < 1:
	print("Usage: " + sys.argv[0] + " <-d folder> <-f plot.txt> <-s False/True>")
	print("  if  -d folder is given, all .plot.txt files under that folder will be checked\n")
	print("  use -f filename if only use a single .plot.txt file\n")
	sys.exit("Error: no plot.txt files provided/found\n")

imptfam = []
impthmm = []
inf = open(modulefile, "r")
for aline in inf:
	subs = aline.split()	
	if subs[0] not in impthmm:
		impthmm.append(subs[0])
	if subs[2] not in imptfam:
		imptfam.append(subs[2][:-1])
inf.close()
	
imptcas = 0
crispr = 0
cas9 = 0
cpf1 = 0
repeat, typeI, typeII, typeIII, typeIV, typeV = 0, 0, 0, 0, 0, 0
typelist = ["Type-I", "Type-II", "Type-III", "Type-IV", "Type-V"]
typegene = [0] * 5
typegeneimp = {}
repeatseq = []
castot = 0
subtypespt, subtypesptimp = {}, {}
for afile in plotfiles:
	#print("plotfile " + afile)
	inf = open(afile, "r")
	for aline in inf:
		subs = aline.strip().split()
		namedes = subs[4]
		typedes = subs[5]
		if "repeat" in namedes:
			repeat += 1
			repeatseq.append(namedes[7:])
		elif (namedes != "unk") and (namedes != "antiRepeat"):
			castot += 1
			(famid, famname) = namedes.split(":")
			if famid in impthmm:
				imptcas += 1
			if famname == "cas9":
				cas9 += 1
			elif famname == "cpf1":
				cpf1 += 1
			for idx in [4, 3, 2, 1, 0]:
				if typelist[idx][1:] in typedes:
					typegene[idx] += 1
					if famid in impthmm:
						if typelist[idx] in typegeneimp:
							typegeneimp[typelist[idx]] = typegeneimp[typelist[idx]] + 1
						else:
							typegeneimp[typelist[idx]] = 1
					break
			if "Subtype" in typedes:
				if typedes not in subtypespt:
					subtypespt[typedes] = 1
					subtypesptimp[typedes] = 0
				else:
					subtypespt[typedes] = subtypespt[typedes] + 1
				if famid in impthmm:
					subtypesptimp[typedes] = subtypesptimp[typedes] + 1
	inf.close()

typedit = {typelist[0]:typegene[0], typelist[1]:typegene[1], typelist[2]:typegene[2], typelist[3]:typegene[3], typelist[4]:typegene[4]}
sorted_typedit = sorted(typedit.items(), key=operator.itemgetter(1), reverse=True)
sorted_subtypespt = sorted(subtypespt.items(), key=operator.itemgetter(1), reverse=True)

repeatnr = 0
if len(repeatseq) > 1:
	out = open("repeat-temp-0000.fa", "w")
	for idx in range(len(repeatseq)):
		out.write(">repeat" + str(idx+1)+"\n")
		out.write(repeatseq[idx]+"\n")
	out.close()
	cmd = cdhitest + " -c 0.9 -i repeat-temp-0000.fa -o repeat-temp-0000-nr0.9.fa >/dev/null 2>&1"
	os.system(cmd)
	inf = open("repeat-temp-0000-nr0.9.fa", "r")
	for aline in inf:
		if aline[0] == '>':
			repeatnr += 1
	inf.close()
	os.remove("repeat-temp-0000.fa")
	os.remove("repeat-temp-0000-nr0.9.fa")
	os.remove("repeat-temp-0000-nr0.9.fa.clstr")
else:
	repeatnr = repeat

typefound = []
typeunflt = []
if (not imptcas) and (not repeat) and (genome == True):
	print("NO CRISPR-Cas found")
elif (imptcas or repeat):
	print("CRISPR-Cas found")
	print("#repeat/cas all nr/important")
	print("Repeat " + str(repeat) + " " + str(repeatnr))
	print("Cas " + str(castot) + " " + str(imptcas))
	for idx in [0, 1, 2, 3, 4]:
		#if sorted_typedit[idx][1] > 0:
		if (sorted_typedit[idx][1] > 0) and (idx == 0):
			if sorted_typedit[idx][0] in typegeneimp:
				print(sorted_typedit[idx][0] + " " + str(sorted_typedit[idx][1]) + " " + str(typegeneimp[sorted_typedit[idx][0]]))
			else:
				print(sorted_typedit[idx][0] + " " + str(sorted_typedit[idx][1]) + " 0")
		elif (sorted_typedit[idx][1] > 0) and (idx > 0): 
			if ((sorted_typedit[idx][0] in typegeneimp) and (typegeneimp[sorted_typedit[idx][0]] > 1)) or ((cas9>0) and (sorted_typedit[idx][0] == "Type-II")) or ((cpf1>0) and (sorted_typedit[idx][0] == "Type-V")):
				if sorted_typedit[idx][0] in typegeneimp:
					print(sorted_typedit[idx][0] + " " + str(sorted_typedit[idx][1]) + " " + str(typegeneimp[sorted_typedit[idx][0]]))
				else:
					print(sorted_typedit[idx][0] + " " + str(sorted_typedit[idx][1]) + " 0")
	for idx in range(len(subtypespt)):
		if sorted_subtypespt[idx][1] > 0:
			if sorted_subtypespt[idx][0] in subtypesptimp:
				print(sorted_subtypespt[idx][0] + " " + str(sorted_subtypespt[idx][1]) + " " + str(subtypesptimp[sorted_subtypespt[idx][0]]))
			else:
				print(sorted_subtypespt[idx][0] + " " + str(sorted_subtypespt[idx][1]) + " 0")
