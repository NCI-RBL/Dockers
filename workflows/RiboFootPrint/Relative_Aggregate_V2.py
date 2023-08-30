# Relative_Aggregate_V2.py

# Relative Aggregates with split into 5utr, cds, 3utr

__author__ = 'Wilfried Guiblet'

import sys
import argparse



import numpy as np


def IndexSubparts(SubPartsFile):

	SubpartsIndex = {}

	for line in open(SubPartsFile, 'rt'):
		chrom, start_5utr, start_cds, start_3utr, end_3utr = line.strip().split('\t')

		SubpartsIndex[chrom] = [start_5utr, start_cds, start_3utr, end_3utr]

	return SubpartsIndex


def Scale_and_Aggregate(Infile, Outfile, SubpartsIndex):

	scaled_positions = {}

	for line in open(Infile, 'rt'):

		chrom, start, end, pos, depth = line.strip().split('\t')
		#length = int(end) - int(start)

		try:
			start_5utr, start_cds, start_3utr, end_3utr = SubpartsIndex[chrom]

			if int(start_5utr) < int(pos) <= int(start_cds):
				subpart = '5utr'
				length = int(start_cds) - int(start_5utr)
				relative_pos = int(pos) - int(start_5utr)

			elif int(start_cds) < int(pos) <= int(start_3utr):
				subpart = 'cds'
				length = int(start_3utr) - int(start_cds)
				relative_pos = int(pos) - int(start_cds)

			elif int(start_3utr) < int(pos) <= int(end_3utr):
				subpart = '3utr'
				length = int(end_3utr) - int(start_3utr)
				relative_pos = int(pos) - int(start_3utr)
			else:
				print('Subpart Error')


			if subpart not in scaled_positions:
				scaled_positions[subpart] = {}
			
			relative_start = round((float(relative_pos) - 1 )/ length, 5)
			relative_end   = round((float(relative_pos)    )/ length, 5)
			for relative_pos in np.arange(relative_start, relative_end, 0.00001):
				relative_pos = round(relative_pos, 5) # getting read of error in float
				if relative_pos not in scaled_positions[subpart]:
					scaled_positions[subpart][relative_pos] = 0
				scaled_positions[subpart][relative_pos] += float(depth)

		except:
			#print('Missing transcript ID in index? : '+chrom)
			continue

	output = open(Outfile, 'w+')

	for subpart in scaled_positions:
		for relative_pos in scaled_positions[subpart]:
			output.write(subpart+'\t'+str(relative_pos)+'\t'+str(scaled_positions[subpart][relative_pos])+'\n')



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--Infile', type=str, required=True)
	parser.add_argument('--SubPartsFile', type=str, required=True)
	parser.add_argument('--Outfile', type=str, required=True)

	args = parser.parse_args()
	Scale_and_Aggregate(args.Infile, args.Outfile, IndexSubparts(args.SubPartsFile))


if __name__ == '__main__':
	main()









