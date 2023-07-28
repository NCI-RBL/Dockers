#!/usr/bin/env python

__author__ = 'Wilfried Guiblet'

import sys
import pandas as pd
import argparse



def mergefunction(uniqueCounts, FracMMCounts, totalCounts):

    merged = uniqueCounts.merge(FracMMCounts, how='outer', on=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'])
    merged = merged.merge(totalCounts, how='outer', on=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'])
    merged.columns = ['ID', 'chr', 'start', 'end', 'strand', 'Length', 'Counts_unique', 'Counts_fracMM', 'Counts_total']
    return(merged)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--uniqueCountsFile', type=str, required=True)
    parser.add_argument('--FracMMCountsFile', type=str, required=True)
    parser.add_argument('--totalCountsFile', type=str, required=True)
    parser.add_argument('--outName', type=str, required=True)


    args = parser.parse_args()

    uniqueCountsFile = pd.read_csv(args.uniqueCountsFile, sep='\t', header=1, skip_blank_lines=True)
    FracMMCountsFile = pd.read_csv(args.FracMMCountsFile, sep='\t', header=1, skip_blank_lines=True)
    totalCountsFile = pd.read_csv(args.totalCountsFile, sep='\t', header=1, skip_blank_lines=True)

    merged_table = mergefunction(uniqueCountsFile, FracMMCountsFile, totalCountsFile)
    merged_table.to_csv(args.outName, sep='\t', header=True, index=False)

if __name__ == '__main__':
	main()