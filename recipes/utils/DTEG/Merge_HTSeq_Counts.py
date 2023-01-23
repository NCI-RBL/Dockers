#!/usr/bin/env python

import pandas as pd
import glob
import sys
from os import listdir
from os.path import isfile, join

# Get a list of all files that end in "counts.csv"
countfiles = sys.argv[3:] 
outputfile = sys.argv[2]
path = sys.argv[1]

# Create an empty list to hold the DataFrames
dfs = []

# Loop through the files and read them into a DataFrame
for file in countfiles:
    df = pd.read_csv(path+'/'+file, sep = '\t')
    df.columns = ['gene', str(file).replace('_counts.txt', '')]
    dfs.append(df)

# Merge the DataFrames on the "common_column"
merged_df = dfs[0]

for df in dfs[1:]:
        merged_df = merged_df.merge(df, on='gene', how="outer")

merged_df.to_csv(path+'/'+outputfile, sep = '\t', index = False)
