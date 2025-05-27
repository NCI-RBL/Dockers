

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
import numpy as np

# pip install dash
# pip install dash-bio==1.0.1

## Test
# module load mamba
# mamba activate cwl_env
#
# python /home/mcintoshc/eCLIP_WF/src/summarize_merge_peaks_wf.py \
# --yaml_file /mnt/gridftp/guibletwm/CCBRRBL16/manifests/ZAP.merge.yaml \
# --directory /mnt/gridftp/guibletwm/CCBRRBL16/ZAP
#
# ln -s /mnt/gridftp/guibletwm/CCBRRBL16/ZAP/eCLIP_MergePeaks_ReadMeSummary.html reports/eCLIP_MergePeaks_ReadMeSummary.html



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

def set_up_readme_summary_html(in_args):
	eCLIP_MergePeaks_ReadMeSummary_html = ""
	try:
		eCLIP_MergePeaks_ReadMeSummary_html = in_args.directory.rstrip('/') + "/eCLIP_MergePeaks_ReadMeSummary.html"
		in_args.SUMMARY_HTML_FILE = open (eCLIP_MergePeaks_ReadMeSummary_html, "w")
	except:
		print("Warning(PermissionError): File will be created in current working directory.")
		eCLIP_MergePeaks_ReadMeSummary_html = "./eCLIP_MergePeaks_ReadMeSummary.html"
		in_args.SUMMARY_HTML_FILE = open (eCLIP_MergePeaks_ReadMeSummary_html, "w")

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

	in_args.SUMMARY_HTML_FILE.write('\t<meta charset="utf-8" />\n\t<title>Summary eCLIP Merge Peaks Workflow Results</title>\n</head>\n<body>\n')
	in_args.SUMMARY_HTML_FILE.write("<h1><i>eCLIP</i> Merge Peaks Workflow Report</h1>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>A full description Yeo Lab <b><i>eCLIP Merge Peaks</i></b> software can be found at this <a href="https://github.com/YeoLab/merge_peaks">GitHub link</a>.</p>\n')

	in_args.SUMMARY_HTML_FILE.write("<p>NOTE: For each table below, only up to " + str(in_args.max_row_count) + " rows are displayed.</p1>\n")


	in_args.SUMMARY_HTML_FILE.write('<p></p>\n')

def close_readme_summary_html(in_args):
	## Close
	in_args.SUMMARY_HTML_FILE.write("\n</body>\n</html>\n")

	in_args.SUMMARY_HTML_FILE.close()

	return eCLIP_MergePeaks_ReadMeSummary_html

def get_yaml_data(in_args):
	with open(in_args.yaml_file, 'r') as YAML_FILE:
		yaml_data = yaml.safe_load(YAML_FILE)

	in_args.merged_peaks_bed       = in_args.directory + "/" + yaml_data['merged_peaks_bed']
	in_args.merged_peaks_custombed = in_args.directory + "/" + yaml_data['merged_peaks_custombed']

	print("\nLocation of Merge Peak BED file:\n\t" + in_args.merged_peaks_bed)
	print("Location of \n\t" + in_args.merged_peaks_custombed)

	## Merged Peaks BED
	colnames = ["Peak Chr","Peak Start","Peak End","Conservative p-value","geomean of the l2FC","Strand"]
	merged_peaks_bed_data_df = pd.read_csv(in_args.merged_peaks_bed, sep="\t", names=colnames, header=None, nrows=in_args.max_row_count)
	# merged_peaks_bed_data_df = merged_peaks_bed_data_df.sort_values(by=['Chr'])

	# try:
	# merged_peaks_bed_data_df['Coordinate'] = pd.concat([merged_peaks_bed_data_df['Chr'], ":", merged_peaks_bed_data_df['Start'], "-",merged_peaks_bed_data_df['End']])
	#
	# merged_peaks_bed_data_df['Coordinate'] = [''.join(i) for i in zip(df["Year"].map(str),df["quarter"])]
	# data = np.array(['g', 'e', 'e', 'k', 's'])
	#
	# y = np.arange(6, dtype=np.double)
	# np.full_like(y, 0.1)
	# np.full_like(x, 1)
	#
	# ser = pd.Series(data)
	# print(ser)
	#
	# except:
	# 	pass

	## Custom BED
	custom_bed_colnames = ["IDR Region","Reproducible Reak Region","geomean of the l2FC","rep1 log2FC","rep2 log2FC","rep1 -log10 pvalue","rep2 -log10 pvalue"]
	custombed_df = pd.read_csv(in_args.merged_peaks_custombed, sep="\t", names=custom_bed_colnames, header=None, nrows=in_args.max_row_count)
	custombed_df['Conservative p-value'] = custombed_df[['rep1 -log10 pvalue', 'rep2 -log10 pvalue']].min(axis=1)

	custombed_df = custombed_df.sort_values(by=['Conservative p-value'], ascending=False)

	return (merged_peaks_bed_data_df, custombed_df)

def print_nice_table(title, in_df, in_args):
	print("\n#### " + title + " ####")
	table = in_df.to_markdown(index=False, tablefmt='simple')
	print(table)

	in_args.SUMMARY_HTML_FILE.write("<h3>" + title + "</h3>\n")
	in_args.SUMMARY_HTML_FILE.write(in_df.to_html(index=False, border=3, table_id='table_fmt'))

def write_merged_peaks_bed_data_df(merged_peaks_bed_data_df, in_args):
	in_args.SUMMARY_HTML_FILE.write("<h2>Reproducible Peaks BED</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>Yeo Lab "<i>merged_peaks_bed: this is the BED6 file containing reproducible peaks as determined by entropy-ordered peaks between two replicates.</i>"</p>\n')
	in_args.SUMMARY_HTML_FILE.write('<p>Column(<b><i>Conservative p-value</i></b>) is the more conservative p-value or greater of the replicates i.e "<b>rep1 -log10 pvalue</b>" vs "<b>rep2 -log10 pvalue</b>" of <b><u>Custom BED</u></b> table above.</p>\n')


	merged_peaks_bed_data_df = merged_peaks_bed_data_df.sort_values(by="Conservative p-value", ascending=False)

	print_nice_table("Reproducible Peaks BED", merged_peaks_bed_data_df, in_args)

def write_custom_bed_data_df(custombed_df, in_args):
	in_args.SUMMARY_HTML_FILE.write("<h2>Replicate Information</h2>\n")
	in_args.SUMMARY_HTML_FILE.write('<p>Yeo Lab - "<i>*.custombed: contains individual replicate information.</i>"</p>\n')
	print_nice_table("Custom BED", custombed_df, in_args)

def write_full_files(in_args):
	full_files = glob.glob(in_args.directory + "/*.full")

	header = "Chr","Start","End","Name (colon separated region)","Reads in CLIP","Reads in Input","p-value","chi value or (F)isher","(F)isher or (C)hi square test","Change","-log10 p-value (400 if > threshold)","log2FC","Entropy"
	in_args.SUMMARY_HTML_FILE.write("<h2>Full Files</h2>\n")

	for full_file in full_files:
		print(full_file)
		full_file_name = os.path.basename(full_file)
		df = pd.read_csv(full_file, sep="\t", names=header, header=None, nrows=in_args.max_row_count)
		print_nice_table(full_file_name, df, in_args)

# def	IDR(in_args):
# 	header=["Chr","Start","End","Name","Score","Strand","signalValue","p-value",
# 		"q-value","summit","Local IDR","Global IDR","ddd","ddd","ddd","ddd"]
# 	full_files = glob.glob(in_args.directory + "/*.idr.out")[0]
# 	print(full_files)



######## main - start ##############################################################################
def main():
	"""*main* - main function.

	"""
	######## Parse arguements ########
	args = get_args()

	args.max_row_count = 50

	YAML_FILE = args.yaml_file
	RESULTS_DIR = args.directory.rstrip('/')

	args.eCLIP_MergePeaks_ReadMeSummary_html = set_up_readme_summary_html(args)
	(merged_peaks_bed_data_df, custombed_df) = get_yaml_data(args)

	## Write merged_peaks_bed_data_df
	write_custom_bed_data_df(custombed_df, args)
	write_merged_peaks_bed_data_df(merged_peaks_bed_data_df, args)
	# write_full_files(args)

	# IDR(args)


if __name__ == "__main__":
		main()
######## main - end ################################################################################
