from __future__ import print_function
#!/usr/bin/env python
import sys
import os
import subprocess
import time
from time import gmtime, strftime

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

def file_put_contents(filename, content):
	subprocess.check_output("echo \"" + content + "\" > " + filename, shell=True)  	

cleanup=True

if len(sys.argv) < 4:
   sys.exit("usage: crisprone-local.py output-base output-prefix input-fna <input-gff>\n")


cmd=sys.argv[0]
crisprone=os.path.dirname(os.path.realpath(__file__))
print("crisprone " + crisprone)
if os.path.exists(crisprone + "/local") == False:
   sys.exit("CRISPRone directory not set correctly")

output_base=sys.argv[1]
output_prefix=sys.argv[2]

output_dir=output_base + "/" + output_prefix
orign_input_name=output_prefix
output = output_base + "/" + output_prefix + "/" + output_prefix
old_input_file=sys.argv[3]
old_gff_file="doesnotexist"
old_gff_valid=False
if len(sys.argv) >= 5:
   old_gff_file=sys.argv[4]
   old_gff_valid=True 

#initial setup (folders, file names)
seq_infor=output + ".infor"
seq_header=output + ".des"

plot_dir=output_dir + "/plot/"
spacer_dir=output_dir + "/spacer/"
if os.path.exists(output_dir) == False:
	subprocess.check_output("mkdir -p " + output_dir, shell=True)

new_input_file=output + ".fna"
crt_output=output + ".crt"
gff_file=output + ".gff"
gff_final=output + "-sm.gff"
repeat_spacer=output + ".repeat-spacer"

cas_hmmer_domt=output + "-vs-cas.domtblout"
fraggenescan_out_faa=output + ".faa"
fraggenescan_out_faa_new_tmp=output + ".faa.tmp"

cas_hmmer_summary=output + ".summary"
cas_hmmer_summary_filter=output + ".summary-filter"
plot_order_with_spacer_seq=output + ".plot.spacer"

if os.path.exists(gff_final):
	sys.exit("results found in folder " + output_dir + "--if you need to redo the annotation, delete the result folder then rerun this script\n")

#check input file
if sys.version_info.major == 2:
	check_result = subprocess.check_output("perl " + crisprone + "/call-script/check_file.pl " + old_input_file + " " + new_input_file + " " + seq_infor, shell=True)
else:
	check_result = subprocess.check_output("perl " + crisprone + "/call-script/check_file.pl " + old_input_file + " " + new_input_file + " " + seq_infor, shell=True, encoding='UTF-8')

if "error" in check_result:
  	sys.exit("error found" + check_result)

subprocess.check_output("cp " + seq_infor + " " + seq_header, shell=True)

#call metaCRT to predict CRISPRs
if not os.path.exists(crt_output):
	print ("run metaCRT...")
	check_result=subprocess.check_output("java -cp " + crisprone + "/bin/metaCRT.jar crt " + new_input_file + " " + crt_output, shell=True)
else:
	print ("metaCRT was done. skip.")
pg = "perl " + crisprone + "/plot-script/build-repeat-spacer-table.pl "
check_result=subprocess.check_output(" ".join([pg, crt_output, repeat_spacer]), shell=True)

#prepare protein file when gff file is given
if os.path.exists(old_gff_file):
	if old_gff_file != gff_file:
		subprocess.check_output("cp " + old_gff_file + " " + gff_file, shell=True)
	print ("prepare protein sequences according to the gff file...")
	pg = "perl " + crisprone + "/call-script/generate-faa-given-gff.pl"
	#print " ".join([pg, gff_file, new_input_file, fraggenescan_out_faa])
	result=subprocess.check_output(" ".join([pg, gff_file, new_input_file, fraggenescan_out_faa]), shell=True)
	if os.path.getsize(fraggenescan_out_faa) <= 10:
		old_gff_valid=False

#predict proteins using FragGeneScan if gff was not provided, or the given gff did not contain information about proteins
#if not old_gff_valid:
#to be checked Oct 21, 2020
if not old_gff_valid or ( not os.path.exists(fraggenescan_out_faa) or os.stat(fraggenescan_out_faa).st_size == 0):
	print ("run FragGeneScan...")
	fraggenescan=crisprone + "/bin/FragGeneScan1.31/FragGeneScan -s " + new_input_file + " -o " + output + " -w 1 -t complete"
	try:
		subprocess.check_output(fraggenescan, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)
	subprocess.check_output("rm -f " + output + ".ffn " + output, shell=True)

	subprocess.check_output("perl " + crisprone + "/call-script/change-fraggenescan-seq-head.pl " + fraggenescan_out_faa + " " + fraggenescan_out_faa_new_tmp, shell=True)
	subprocess.check_output("perl " + crisprone + "/call-script/generate-gff.pl " + fraggenescan_out_faa + " " + gff_file, shell=True)

if not os.path.exists(fraggenescan_out_faa):
    	sys.exit("Generation of protein sequences failed")

#cas prediction
if not os.path.exists(cas_hmmer_domt):
	run_hmmer=crisprone + "/bin/hmmer-3.1b1/bin/hmmscan --cpu 4 --noali --domtblout " + cas_hmmer_domt + " -E 0.01 " + crisprone + "/local/cas-db/cas-CRISPR.LIB " + fraggenescan_out_faa
	subprocess.check_output(run_hmmer, shell=True)

create_cas_order_table="python " + crisprone + "/call-script/hmmscan2cas.py -collection unk -hmmscan " + output + "-vs-cas.domtblout -gff " + gff_file + " -isolate -outfile " + cas_hmmer_summary + " -outfiltered " + cas_hmmer_summary_filter + " -gafile " + crisprone + "/local/cas-db/cas-CRISPR.GA"
subprocess.check_output(create_cas_order_table, shell=True)
subprocess.check_output("rm -f " + cas_hmmer_summary, shell=True)
 
print ("final touches...")
#prepare to plot the CAS locus & CRISPR array
plot_cmd="perl " + crisprone + "/plot-script/insert-repeat-with-spacer-v4.pl " + repeat_spacer + " " + cas_hmmer_summary_filter + " " + plot_order_with_spacer_seq
subprocess.check_output(plot_cmd, shell=True)
subprocess.check_output("mkdir -p " + plot_dir + " " + spacer_dir, shell=True)
subprocess.check_output("perl " + crisprone + "/plot-script/ready-for-plot-v2.pl " + plot_order_with_spacer_seq + " " + plot_dir + " " + spacer_dir, shell=True)

#predict type II tracrRNA
if os.path.getsize(repeat_spacer) > 10:
	adjust_typeII="perl " + crisprone + "/plot-script/adjust-typeII-v2.pl " + plot_dir
	gen_crRNA="perl " + crisprone + "/call-script/identify-tracrRNA-v2.pl " + plot_dir + " " + repeat_spacer + " " + new_input_file
	adjust_crRNA="perl " + crisprone + "/call-script/adjust-tracrRNA-result.pl " + plot_dir
	subprocess.check_output(adjust_typeII, shell=True)
	subprocess.check_output(gen_crRNA, shell=True)
	subprocess.check_output(adjust_crRNA, shell=True)
	#modify_plot="perl " + crisprone + "/plot-script/modify-plot-file.pl " + plot_dir
	#subprocess.check_output(modify_plot, shell=True) #not needed anymore

#check for mock CRISPRs
if os.path.getsize(repeat_spacer) > 10: 
	mock_crispr="perl " + crisprone + "/call-script/check-mock-CRISPR.pl " + crt_output
	subprocess.check_output(mock_crispr, shell=True)

#check for suspicious cas
cas_check="perl " + crisprone + "/call-script/check-suspicious-cas.pl " + plot_dir + " " + crt_output
subprocess.check_output(cas_check, shell=True)
regenerate_plot="perl " + crisprone + "/call-script/regenerate-plot.pl " + seq_infor + " " + plot_dir
subprocess.check_output(regenerate_plot, shell=True)

#summarize in one standard gff file (YY June 2018)
generate_gff="python " + crisprone + "/call-script/crispr-gff.py " + seq_infor + " " + plot_dir + " " + gff_final
if old_gff_valid:
	generate_gff = generate_gff + " " + gff_file
subprocess.check_output(generate_gff, shell=True)

#get cas gene sequences (in one file, YY June 2018)
generate_seq="perl " + crisprone + "/call-script/get_cas_seq-v2.pl " + new_input_file + " " + gff_final + " " + output + "-sm"
subprocess.check_output(generate_seq, shell=True)

#cleaning up immediate files
if cleanup:
	subprocess.check_output("rm -r -f " + plot_dir + " " + spacer_dir, shell=True)
	subprocess.check_output("rm -f " + " ".join([output + ".fna", output + ".repeat-spacer", output + ".gff", output + ".gff.old", output + ".faa", output + "summary-filter", output + ".plot.spacer", output + ".fna.n*", output + ".crt.*", output + ".infor"]), shell=True)

#complete php file
#os.environ['TZ'] = 'EST+05EDT,M4.1.0,M10.5.0'
#time.tzset()
#timestamp = time.strftime('%X %x %Z')
#job_show_result=" http://omics.informatics.indiana.edu/CRISPRone/show-col.php?col=tmp&id=" + output_prefix
