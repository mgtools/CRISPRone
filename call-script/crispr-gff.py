from __future__ import print_function
#!/usr/bin/env python
#latest update Nov 2020 YY, compatible with python 3
import sys
import os
from os import listdir
from os.path import isfile, join
from os import path
import subprocess

if "check_output" not in dir( subprocess ): #old python doesn't support check_output
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

if len(sys.argv) < 4:
    sys.exit(sys.argv[0] + " crisprone-infor-file plot-dir crisprone-gff-outputfile <old-gff-file>")
inforfile, plotdir, gffoutfile = sys.argv[1:4]
dname = os.path.dirname(sys.argv[0])
colorf = dname + "/../local/cas-db/subtype.color"
typepg = dname + "/crispr-type.pl " + colorf
oldname, oldgene = {}, {}
if len(sys.argv) > 4:
    inf = open(sys.argv[4], "r")
    for aline in inf:
        subs = aline.strip().split("\t")
        subs2 = subs[-1].split(";")
        if len(subs2) > 2:
            #print subs2 
            oldname[subs2[0][3:]] = subs2[2][4:]
            oldgene[subs2[0][3:]] = subs2[3][5:]
    inf.close()
	
inf = open(inforfile, "r")
seqlen = {}
for aline in inf:
    if aline[0] == '#':
        break
    tmp = aline.strip().split()
    seqlen[tmp[0]] = tmp[1]
inf.close()

allfiles = [f for f in listdir(plotdir) if isfile(join(plotdir, f))]

for suffix in [".plot.txt", ".suspicious"]: 
   plotfiles = [f for f in allfiles if (suffix in f)]
   if suffix == ".suspicious": 
        out = open(gffoutfile + ".sus", "w")	
   else:
        out = open(gffoutfile, "w")
   out.write("##gff-version 3\n")
   out.write("#!processor CRISPRone\n")
   for plotfile in plotfiles:
       inf = open(plotdir + "/" + plotfile, "r")
       lines = inf.readlines()
       inf.close()
       if len(lines) < 1:
           continue
       tmps = lines[0].split()
       sid, spos, epos, strand, des = tmps[0:5]
       if len(tmps) > 5:
           stype = tmps[5]
       else:
           stype = "sus"
       slen = seqlen[sid]
       items = [sid, "Genbank", "region", "1", slen, ".", "+", "."]
       if suffix == ".plot.txt":
           cmd="perl " + typepg + " " + plotdir + "/" + plotfile
           if sys.version_info.major == 2:
               output = subprocess.check_output(cmd, shell=True)
           else:
               output = subprocess.check_output(cmd, shell=True, encoding='UTF-8')
           print ("type: " + cmd)
           print ("# " + output)
           tmps = output.split("\n")
           items.append("ID=" + sid + ";array=" + tmps[0] + ";gene=" + tmps[1] + ";type=" + tmps[2])
       else:
           items.append("ID=" + sid)
       out.write("\t".join(items) + "\n")
       for aline in lines:
           (sid, spos, epos, strand, des, stype) = aline.split()
           if "repeat" in des:
               tmp = stype.split(":")
               note = "what=CRISPR;des=repeat:" + des[7:] + ";" + "note=copy:" + str(len(tmp))
               items = [sid, "metaCRT", "direct_repeat", spos, epos, ".", strand, ".", note]
               out.write("\t".join(items) + "\n")
           elif "antiRepeat" in des:
               note = "what=antiCRISPR;des=antiRepeat:" + des[7:] + ";" + "note=" + stype
               items = [sid, "CRISPRone", "partial_repeat", spos, epos, ".", strand, ".", note] 
               out.write("\t".join(items) + "\n")
           else:
               thisid = "_".join([sid, spos, epos, strand])
               note = "what=cas;des=" + des + ";" + "note=" + stype
               if thisid in oldname:
                   note = note + ";Name=" + oldname[thisid]
               if thisid in oldgene:
                   note = note + ";gene=" + oldgene[thisid]
               items = [sid, "FGS-hmm", "CDS", spos, epos, ".", strand, ".", note]
               out.write("\t".join(items) + "\n")
   out.close()
