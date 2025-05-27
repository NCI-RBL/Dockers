

"""
@package docstring
Documentation for this module.
More details.

:Program Name: post_process_eCLIP.py
:Developer: Carl McIntosh <mcintoshc@nih.gov>
:Date: May 13, 2025
:Version: May 13, 2025

Usage Example
=============
``post_process_eCLIP.py``

Methods
=======

"""

import argparse
import yaml
import json
import sys
import glob
import subprocess
from tabulate import tabulate
import os
import pandas as pd


# pip install dash
# pip install dash-bio==1.0.1

## Test
# python /home/mcintoshc/eCLIP_TestZone/process_eclip.py \
# --yaml_file /mnt/gridftp/guibletwm/CCBRRBL16/manifests/ZAP.yaml \
# --directory /mnt/gridftp/guibletwm/CCBRRBL16


## Get program arguements ######################################################
def get_args():
	"""*get_args* - parses program's arg values.

	:returns: (*dict*) Contains user provided variables.
	"""
	parser = argparse.ArgumentParser()

	## Required
	group_required = parser.add_argument_group('required arguments')
	group_required.add_argument("-d", "--directory", help="Provide path to data directory. ", required=True)
	group_required.add_argument("-y", "--yaml_file", help="Provide path to data directory. ", required=True)

	## options
	group_optional = parser.add_argument_group('optional arguments')
	group_optional.add_argument("--verbose", help=". ", action="store_true", required=False)

	return parser.parse_args()

def summarize_cutadapt(in_dict, in_args):

	in_args.SUMMARY_HTML_FILE.write("<h2>Cutadapt Information</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p><b><i>Cutadapt</i></b> information is provided below. Documentation can be found at <a href="https://cutadapt.readthedocs.io/en/stable/#">Cutadapt</a>.</p>\n')

	summary=[]
	headers=["Total reads processed","Reads with adapters","Reads that were too short","Reads written (passing filters)","Total basepairs processed","Quality-trimmed","Total written (filtered)"]
	for key in in_dict.keys():
		if "_metric_" in key:
			metric_row={"Cutadapt Metrics File":os.path.basename(in_dict[key])}
			with open(in_dict[key],'r') as METRIC_FILE:
				for index in range(1, 20):
					if index > 20: break
					line = METRIC_FILE.readline()
					line = line.split(':')
					if line[0] in headers:
						value = line[1].replace(' ','')
						value = value.replace('bp',' bp')
						value = value.replace('(' ,' (')
						metric_row[line[0]]=value.rstrip('\n')
				summary.append(metric_row)
	df = pd.DataFrame(summary)
# 	pd.set_option('display.max_rows', None)  # Show all rows
# 	pd.set_option('display.max_columns', None)  # Show all columns
# 	pd.set_option('display.width', 1000)  # Adjust the width as needed
#
# 	print(df.to_string(index=False))
# 	print("")

	print_nice_table("Cutadapt Metrics", df, in_args)

	return df

def print_nice_table(title, in_df, in_args):
	print("\n#### " + title + " ####")
	table = in_df.to_markdown(index=False, tablefmt='simple')
	print(table)

	in_args.SUMMARY_HTML_FILE.write("<h3>" + title + "</h3>\n")




	in_args.SUMMARY_HTML_FILE.write(in_df.to_html(index=False, border=3, table_id='table_fmt'))


def summarize_STAR_FINAL_LOGS(star_log_files, in_args):

	in_args.SUMMARY_HTML_FILE.write("<h2>STAR Aligner Information</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p><b><i>STAR Aligner</i></b> information is provided below. Documentation can be found at <a href="https://github.com/alexdobin/STAR">STAR Aligner</a>.</p>\n')

	star_final_logs_info = []
	general_info = []
	uniq_reads_info = []
	multi_mapped_reads_info = []
	unmapped_reads_info = []
	chimeric_reads_info = []

	# print("Number of Log Files: " + str(len(star_log_files)))

	for star_log in star_log_files:
		short_file_name = os.path.basename(star_log)
		short_file_name = short_file_name.replace(".sorted.STARLog.final.out","")
		log_header = "STAR Log File (.sorted.STARLog.final.out)"

		star_final_logs_info.append({log_header: short_file_name, "Full Path": star_log})

		tmp_dict = {"STAR Log File":short_file_name}
		with open(star_log,'r') as STAR_LOG_FILE:
			file_lines = STAR_LOG_FILE.readlines()
			for line in file_lines:
				if "UNIQUE READS:" in line:
					general_info.append(tmp_dict)
					tmp_dict = {"STAR Log File":short_file_name}
					# tmp_dict = {}
				elif "MULTI-MAPPING READS:" in line:
					uniq_reads_info.append(tmp_dict)
					tmp_dict = {"STAR Log File":short_file_name}
					# tmp_dict = {}
				elif "UNMAPPED READS:" in line:
					multi_mapped_reads_info.append(tmp_dict)
					tmp_dict = {"STAR Log File":short_file_name}
					# tmp_dict = {}
				elif "CHIMERIC READS:" in line:
					unmapped_reads_info.append(tmp_dict)
					tmp_dict = {"STAR Log File":short_file_name}
					# tmp_dict = {}
				elif "|" in line:
					line = line.lstrip(" ")
					line = line.rstrip("\n")
					(header,value) = line.split('|')
					header = header.rstrip(" ")
					header = header.replace("ddd","ddd")
					header = header.replace("Number of input reads","Input Read Cts")
					header = header.replace("Average input read length","Ave Input Read Len")
					header = header.replace("Mapping speed, Million of reads per hour","Map Speed (Million Reads/hour")
					header = header.replace("Average mapped length","Ave. Map Len")
					header = header.replace("Number of splices","Splice Number")
					header = header.replace("Deletion rate per base","Deletion Rate/Base")
					header = header.replace("Insertion rate per base","Insertion Rate/Base")
					header = header.replace("Number of splices","Splice Cts")
					header = header.replace("Number of ","# of ")
					header = header.replace("Average","Ave.")
					header = header.replace("Count","Cts.")
					header = header.replace("Count","Cts.")


					value = value.replace("%","")
					value = value.replace("\t","")

					tmp_dict[header] = value

			chimeric_reads_info.append(tmp_dict)
			tmp_dict = {}

	star_final_logs_info_df = pd.DataFrame(star_final_logs_info)
	print_nice_table("STAR Log Files - Path Info Table", star_final_logs_info_df, in_args)

	general_info_df = pd.DataFrame(general_info)
	print_nice_table("STAR Final Log Files - General Info Table", general_info_df, in_args)

	uniq_reads_info_df = pd.DataFrame(uniq_reads_info)
	print_nice_table("STAR Final Log Files - Uniq Reads Table", uniq_reads_info_df, in_args)

	multi_mapped_reads_info_df = pd.DataFrame(multi_mapped_reads_info)
	print_nice_table("STAR Final Log Files - Uniq Reads Table", multi_mapped_reads_info_df, in_args)

	unmapped_reads_info_df = pd.DataFrame(unmapped_reads_info)
	print_nice_table("STAR Final Log Files - Unmapped Reads Table", unmapped_reads_info_df, in_args)

	chimeric_reads_info_df = pd.DataFrame(chimeric_reads_info)
	print_nice_table("STAR Final Log Files - Chimeric Reads Table", chimeric_reads_info_df, in_args)

def extract_fastqc_data_files(fastqc_data_file):
	out_data_dict = {"Filename":""}
	with open(fastqc_data_file,'r') as FASTQC_DATA_FILE:
		file_lines = FASTQC_DATA_FILE.readlines()
		for line in file_lines:
			line = line.rstrip("\n")
			if ">>Basic Statistics" in line:                  out_data_dict["Basic Stats"] = line.split('\t')[1]
			elif "Filename" in line:                          out_data_dict["Filename"] = line.split('\t')[1]
			# elif "File type" in line:                         out_data_dict["File type"] = line.split('\t')[1]
			# elif "Encoding" in line:                          out_data_dict["Encoding"] = line.split('\t')[1]
			elif "Total Sequences" in line:                   out_data_dict["Total Seqs"] = line.split('\t')[1]
			elif "Sequences flagged as poor quality" in line: out_data_dict["Seqs Flagged as Poor Quality"] = line.split('\t')[1]
			elif "Sequence length" in line:                   out_data_dict["Seqs Length Range"] = line.split('\t')[1]
			elif "%GC" in line:                               out_data_dict["%GC"] = line.split('\t')[1]
			elif ">>Per base sequence quality" in line:       out_data_dict["Per Base Seq Quality"] = line.split('\t')[1]
			elif ">>Per tile sequence quality" in line:       out_data_dict["Per Tile Seq Quality"] = line.split('\t')[1]
			elif ">>Per sequence quality scores" in line:     out_data_dict["Per Seq Quality scores"] = line.split('\t')[1]
			elif ">>Per base sequence content" in line:       out_data_dict["Per base Seq Content"] = line.split('\t')[1]
			elif ">>Per sequence GC content" in line:         out_data_dict["Per Seq GC Content"] = line.split('\t')[1]
			elif ">>Per base N content" in line:              out_data_dict["Per Base N Content"] = line.split('\t')[1]
			elif ">>Sequence Length Distribution" in line:    out_data_dict["Seq Length Distribution"] = line.split('\t')[1]
			elif ">>Sequence Duplication Levels" in line:     out_data_dict["Seq Duplication Levels"] = line.split('\t')[1]
			elif ">>Overrepresented sequences" in line:       out_data_dict["Overrepresented Seq"] = line.split('\t')[1]
			elif ">>Adapter Content" in line:                 out_data_dict["Adapter Content"] = line.split('\t')[1]
	return out_data_dict

def summarize_fastqc(fastqc_data_files, fastqc_html_files, in_args):

	in_args.SUMMARY_HTML_FILE.write("<h2>FastQC Information</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>Individual <b><i>FastQC</i></b> HTML reports are provided for each of the input FASTQ files. Documentation can be found at <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">Babraham Bioinformatics</a>.</p>\n')

	results_dir = in_args.directory

	results_dir = results_dir.rstrip('/')
	print("Files located: " + results_dir)

	summary_fastq_info = {
		"FastQC HTML Files":[fastqc_html_file.replace(results_dir + "/","") for fastqc_html_file in fastqc_html_files],
		"FastQC Data Files":[fastqc_data_file.replace(results_dir + "/","") for fastqc_data_file in fastqc_data_files]
	}

	general_info_fastqc_df = pd.DataFrame(summary_fastq_info)
	print_nice_table("FastQC - File Summary Table", general_info_fastqc_df, in_args)



	fastq_summary_data=[]
	for fastqc_data_file in fastqc_data_files:
		fastqc_data = extract_fastqc_data_files(fastqc_data_file)
		fastq_summary_data.append(fastqc_data)

	fastq_summary_data_df = pd.DataFrame(fastq_summary_data)
	print_nice_table("FastQC - Basic Summary Table", fastq_summary_data_df, in_args)

def set_up_readme_summary_html(in_args):
	eCLIP_ReadMeSummary_html = in_args.directory.rstrip('/') + "/eCLIP_ReadMeSummary.html"
	in_args.SUMMARY_HTML_FILE = open (eCLIP_ReadMeSummary_html, "w")
	in_args.SUMMARY_HTML_FILE.write('<!DOCTYPE html>\n<html lang="en">\n<head>\n')
	in_args.SUMMARY_HTML_FILE.write('\t<style>\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t h1 {text-align: center;}\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t h2 {text-decoration: underline;}\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t h3 {text-decoration: underline;}\n')

	in_args.SUMMARY_HTML_FILE.write('\t\t p {font-size: 20px;}\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t table {border-collapse: collapse;}\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t th {font-size: 18px; padding: 8px; text-align: center;}\n')
	in_args.SUMMARY_HTML_FILE.write('\t\t td {font-size: 18px; padding: 8px}\n')



	in_args.SUMMARY_HTML_FILE.write('\t\t tr:nth-child(even) {background-color: #AED6F1;}\n')

	in_args.SUMMARY_HTML_FILE.write('\t</style>\n')



	# in_args.SUMMARY_HTML_FILE.write('\t\tth, td {padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }\n')

	# in_args.SUMMARY_HTML_FILE.write('\t\tddd\n')



	in_args.SUMMARY_HTML_FILE.write('\t<meta charset="utf-8" />\n\t<title>Summary eCLIP Workflow Results</title>\n</head>\n<body>\n')
	in_args.SUMMARY_HTML_FILE.write("<h1><i>eCLIP</i> Workflow Report</h1>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>A full description Yeo Lab <b><i>eCLIP</i></b> software can be found at this <a href="https://github.com/YeoLab/eCLIP?tab=readme-ov-file">GitHub link</a>.</p>\n')


	in_args.SUMMARY_HTML_FILE.write("<h2>MultiQC Information</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>An aggregate HTML summary report was generated using <b><i>MultiQC</i></b>. Documentation can be found at <a href="https://seqera.io/multiqc/">Seqera MultiQC</a>. The file, <b><i>eCLIP_multiqc_report.html</i></b>, can be located in the run directory.</p>\n')

	in_args.SUMMARY_HTML_FILE.write('<p></p>\n')


	return eCLIP_ReadMeSummary_html


def close_readme_summary_html(in_args):
	## Close
	in_args.SUMMARY_HTML_FILE.write("\n</body>\n</html>\n")

	in_args.SUMMARY_HTML_FILE.close()


######## main - start ##############################################################################
def main():
	"""*main* - main function.

	"""
	######## Parse arguements ########
	args = get_args()

	YAML_FILE = args.yaml_file
	RESULTS_DIR = args.directory.rstrip('/')

	args.eCLIP_ReadMeSummary_html = set_up_readme_summary_html(args)

	yaml_dict=None

	with open(YAML_FILE, 'r') as file:
		yaml_dict = yaml.safe_load(file)

	# print(json.dumps(yaml_dict,sort_keys=True, indent=4))

	dataset_name = yaml_dict["dataset"]

	number_of_repliacates = len(yaml_dict["samples"])

	print("Dateset Name: " + dataset_name)
	print("Replicate Count: " + str(number_of_repliacates))

	replicates = []
	replicate_dict = {}
	for replicate in yaml_dict["samples"]:
		# print(json.dumps(replicate[0],sort_keys=True, indent=4))
		replicates.append([replicate[0][ "name"],replicate[1][ "name"]])

	print(replicates)

	## Bam files
	for replicate in replicates:

		## IP Files
		# deduplicate reads based on Unique Molecular Identifiers(UMIs) remove PCR Duplicates
		replicate_dict["IP_tr_metric_" + replicate[0]]    = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + ".*.fqTr.metrics")[0]
		replicate_dict["IP_tr_tr_metric_" + replicate[0]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + ".*.fqTrTr.metrics")[0]
		replicate_dict["IP_bam_" + replicate[0]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + ".*.rmDupSo.bam")[0]
		replicate_dict["IP_rpm_norm_pos_bw_" + replicate[0]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + "*.rmDupSo.norm.pos.bw")[0]
		replicate_dict["IP_rpm_norm_neg_bw_" + replicate[0]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + "*.rmDupSo.norm.neg.bw")[0]

		## Input Files
		# deduplicate reads based on Unique Molecular Identifiers(UMIs) remove PCR Duplicates
		replicate_dict["Input_tr_metric_" + replicate[0]]    = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[1] + ".*.fqTr.metrics")[0]
		replicate_dict["Input_tr_tr_metric_" + replicate[0]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[1] + ".*.fqTrTr.metrics")[0]
		replicate_dict["Input_bam_" + replicate[1]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[1] + ".*.rmDupSo.bam")[0]
		replicate_dict["Input_rpm_norm_pos_bw_" + replicate[1]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[1] + "*.rmDupSo.norm.pos.bw")[0]
		replicate_dict["Input_rpm_norm_neg_bw_" + replicate[1]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[1] + "*.rmDupSo.norm.neg.bw")[0]

		## Combined files
		replicate_dict["CLIPperPeaksBed" + replicate[1]] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + ".*.peakClusters.bed")[0]
		replicate_dict["input_norm_peaks"] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + ".*.peakClusters.normed.compressed.bed")[0]
		replicate_dict["Blacklist-filtered peaks"] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + "*.blacklist-removed.bed")[0]
		replicate_dict["Blacklist-filtered bigBeds"] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + "*blacklist-removed.fx.bb")[0]
		replicate_dict["Blacklist-filtered narrowPeaks"] = glob.glob(RESULTS_DIR + "/" + dataset_name + "." + replicate[0] + "*narrowPeak")[0]

		## STAR Log files
		replicate_dict["STAR_FINAL_LOGS"] = glob.glob(RESULTS_DIR + "/" + dataset_name + "*.STARLog.final.out")


	# print(json.dumps(replicate_dict, indent=4))
	print("\n##### MultiQC Report #####")
	replicate_dict["MultiQC Report"] = glob.glob(RESULTS_DIR + "/" + "*multiqc_report.html")[0]
	print (replicate_dict["MultiQC Report"])


	print("\n##### FastQC Report #####")
	replicate_dict["FastQC Data Files"] = glob.glob(RESULTS_DIR + "/" + "*.fastqc_data.txt")
	replicate_dict["FastQC HTMLs"] = glob.glob(RESULTS_DIR + "/" + "*.fastqc_report.html")
	summarize_fastqc(replicate_dict["FastQC Data Files"], replicate_dict["FastQC HTMLs"], args)

	print("\n######## Cutadapt ########")
	summarize_cutadapt(replicate_dict, args)
	print("\n########## STAR ##########")
	summarize_STAR_FINAL_LOGS(replicate_dict["STAR_FINAL_LOGS"], args)




	# Close
	close_readme_summary_html(args)
	print(args.eCLIP_ReadMeSummary_html)

	## Peak Cluster Cluster

#
#
# 	ZAP.rep1_INPUT.*
#
# 	"dataset": "ZAP"
# 		[rep1_IP,rep1_INPUT]
# 		[rep2_IP,rep2_INPUT]

# [mcintoshc@fsitgl-head01p CCBRRBL16]$ ls -l ZAP.*.rmDupSo.bam
# -rw-rw-r-- 1 guibletwm guibletwm 109788 Apr 24 14:09 ZAP.rep1_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
# -rw-rw-r-- 1 guibletwm guibletwm  34160 Apr 24 19:38 ZAP.rep1_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
# -rw-rw-r-- 1 guibletwm guibletwm 132729 Apr 24 14:00 ZAP.rep2_INPUT.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
# -rw-rw-r-- 1 guibletwm guibletwm  27047 Apr 24 14:40 ZAP.rep2_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
#
#
# [mcintoshc@fsitgl-head01p CCBRRBL16]$ ls -l *.rmDupSo.peakClusters.bed
# -rw-rw-r-- 1 guibletwm guibletwm 7033 Apr 24 20:38 ZAP.rep1_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed
# -rw-rw-r-- 1 guibletwm guibletwm 3719 Apr 24 15:42 ZAP.rep2_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed

if __name__ == "__main__":
	main()
######## main - end ################################################################################
