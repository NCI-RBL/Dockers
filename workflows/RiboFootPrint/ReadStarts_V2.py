# ReadStarts_V2.py

__author__ = 'Wilfried Guiblet'

import sys
import argparse

def IndexTranscriptome(TranscriptomeFile):

	TranscriptomeIndex = {}

	for line in open(TranscriptomeFile, 'rt'):
		chrom, transcript_start, transcript_end = line.strip().split('\t')

		TranscriptomeIndex[chrom] = [transcript_start, transcript_end]

	return TranscriptomeIndex

def ReadStarts(Infile, Outfile, TranscriptomeIndex):

	infile = open(Infile, 'rt')
	outfile = open(Outfile, 'w+')
	starts = {}

	for line in infile:

		array = line.strip().split('\t')
		chrom, start, end = array[0:3]
		#subpart_start, subpart_end = array[7], array[8]

		try:
			transcript_start, transcript_end = TranscriptomeIndex[chrom]


			if int(start) >= int(transcript_start):

				key = chrom+'|'+transcript_start+'|'+transcript_end+'|'+str(int(start)+1)

				if key not in starts:
					starts[key] = 0
				starts[key] += 1

		except:
			#print('Missing transcript ID in index? : '+chrom)
			continue

	for key in starts:
		chrom, transcript_start, transcript_end, pos = key.split('|')
		outfile.write(chrom+'\t'+transcript_start+'\t'+transcript_end+'\t'+pos+'\t'+str(starts[key])+'\n')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--Infile', type=str, required=True)
	parser.add_argument('--Outfile', type=str, required=True)
	parser.add_argument('--TranscriptomeFile', type=str, required=True)


	args = parser.parse_args()
	ReadStarts(args.Infile, args.Outfile, IndexTranscriptome(args.TranscriptomeFile))


if __name__ == '__main__':
	main()
