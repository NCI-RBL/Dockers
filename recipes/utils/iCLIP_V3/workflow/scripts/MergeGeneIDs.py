# MergeGeneIDs.py

import pandas as pd
import sys

PeakFile = pd.read_csv(sys.argv[1], sep = '\t', header=None).fillna('NA')
PeakFile.columns = ['transcript_id', 'Peak_Start', 'Peak_End', 'Peak_ID', 'Score', 'Strand', 'Length', 'Counts_Unique', 'Counts_fracMM', 'Counts_Total', 'Feature', 'transcript_id_dup', 'Feature_Start', 'Feature_End']


GenCodeFile = pd.read_csv(sys.argv[2], sep = '\t', header=0, low_memory=False)

GenCodeFile = GenCodeFile[['gene_name', 'transcript_id']].drop_duplicates()


Merged_df = PeakFile.merge(GenCodeFile, how='inner', on=['transcript_id'])
Merged_df = Merged_df[['gene_name', 'transcript_id','Peak_Start', 'Peak_End', 'Peak_ID', 'Length', 'Counts_Unique', 'Counts_fracMM', 'Counts_Total', 'Feature', 'Feature_Start', 'Feature_End']]

Merged_df.loc[Merged_df['Feature'] == '.', 'Feature'] = ''
Merged_df.loc[Merged_df['Feature_Start'] == -1, 'Feature_Start'] = ''
Merged_df.loc[Merged_df['Feature_End'] == -1, 'Feature_End'] = ''

Merged_df.to_csv(sys.argv[3], index=False)
