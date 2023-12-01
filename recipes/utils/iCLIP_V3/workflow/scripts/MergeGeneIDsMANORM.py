# MergeGeneIDsMANORM.py

import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ManormFile', type=str, required=True)
    parser.add_argument('--PeakFile', type=str, required=True)
    parser.add_argument('--BackgroundPeaks', type=str, required=True)
    parser.add_argument('--GenCodeFile', type=str, required=True)
    parser.add_argument('--outName', type=str, required=True)

    args = parser.parse_args()

    PeakFile = pd.read_csv(args.PeakFile, sep = '\t')
    PeakFile.columns = ['transcript_id', 'Peak_Start', 'Peak_End', 'Peak_ID', 'Score', 'Strand', 'Length', 'Counts_Unique_Sample', 'Counts_fracMM_Sample',\
                        'Counts_Total_Sample', 'Feature', 'transcript_id_dup', 'Feature_Start', 'Feature_End']

    BackgroundPeaks =  pd.read_csv(args.BackgroundPeaks, sep = '\t')
    BackgroundPeaks.columns = ['transcript_id', 'Peak_Start', 'Peak_End', 'Peak_ID', 'Score', 'Strand', 'Length', 'Counts_Unique_Background',\
                               'Counts_fracMM_Background', 'Counts_Total_Background']

    PeakFile = PeakFile.merge(BackgroundPeaks, how='inner', on=['transcript_id', 'Peak_Start', 'Peak_End', 'Peak_ID', 'Score', 'Strand', 'Length'])

    ManormFile = pd.read_csv(args.ManormFile, sep = '\t')
    ManormFile.columns = ['transcript_id', 'Peak_Start', 'Peak_End', 'MANORM_PEAK', 'Score', 'Strand', 'm_value', 'p_value', 'peak_group',\
                          'normalized_read_density_in_sample', 'normalized_read_density_in_background']

    PeakFile = PeakFile.merge(ManormFile.drop(columns=['MANORM_PEAK']), how='inner', on=['transcript_id', 'Peak_Start', 'Peak_End', 'Score', 'Strand'])

    GenCodeFile = pd.read_csv(args.GenCodeFile, sep = '\t', header=0, low_memory=False)
    GenCodeFile = GenCodeFile[['gene_name', 'transcript_id']].drop_duplicates()

    Merged_df = PeakFile.merge(GenCodeFile, how='inner', on=['transcript_id'])
    Merged_df = Merged_df[['gene_name', 'transcript_id','Peak_Start', 'Peak_End', 'Peak_ID', 'Length', 'Counts_Unique_Sample', 'Counts_fracMM_Sample',\
                           'Counts_Total_Sample', 'Counts_Unique_Background', 'Counts_fracMM_Background', 'Counts_Total_Background', 'Feature',\
                           'Feature_Start', 'Feature_End', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_in_sample',\
                           'normalized_read_density_in_background']]

    Merged_df.loc[Merged_df['Feature'] == '.', 'Feature'] = ''
    Merged_df.loc[Merged_df['Feature_Start'] == -1, 'Feature_Start'] = ''
    Merged_df.loc[Merged_df['Feature_End'] == -1, 'Feature_End'] = ''

    Merged_df.to_csv(args.outName, index=False)

if __name__ == '__main__':
        main()
