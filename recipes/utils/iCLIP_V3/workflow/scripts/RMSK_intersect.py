# RMSK_intersect.py

__author__ = 'Wilfried Guiblet'

import pandas as pd
import itertools
import sys
import re
import argparse
import numpy as np



def ProcessRMSK(InputFile):

    RMSK_df = pd.read_csv(InputFile, sep = '\t', header = None,\
                names = ['chrom', 'start', 'end', 'ID', 'score', 'strand', 'length', 'counts_unique', 'counts_fracMM',\
                     'counts_total', 'repChrom', 'repStart', 'repEnd', 'repName', 'repClass', 'repStrand', 'repFamily'])

    RMSK_df['repStart'] = RMSK_df['repStart'].astype('str') # change type necessary for aggregate
    RMSK_df['repEnd']   = RMSK_df['repEnd'].astype('str') # change type necessary for aggregate

    # https://stackoverflow.com/questions/36271413/pandas-merge-nearly-duplicate-rows-based-on-column-value
    RMSK_df = RMSK_df.groupby(['chrom', 'start', 'end', 'ID', 'score', 'strand', 'length', 'counts_unique', 'counts_fracMM',\
                     'counts_total']).agg({'repStart': ', '.join,
                                           'repEnd': ', '.join,
                                           'repName' : ', '.join,
                                           'repClass': ', '.join,
                                           'repStrand': ', '.join,
                                           'repFamily': ', '.join,
                                                                 }).reset_index()



    print(RMSK_df)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--Input', type=str, required=True)
    #parser.add_argument('--Output', type=str, required=True)

    args = parser.parse_args()
    ProcessRMSK(args.Input)


if __name__ == '__main__':
    main()
