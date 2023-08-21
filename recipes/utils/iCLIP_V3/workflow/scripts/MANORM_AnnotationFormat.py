# MANORM_AnnotationFormat.py

__author__ = 'Wilfried Guiblet'

# Inputs:  Outputs of bedtools intersect for GenCode, RepeatMasker, and Introns (custom)
# Goals: - Aggregates multiple rows with same read
#        - Combine the three annotations in a curated table

import pandas as pd
import itertools
import sys
import re
import argparse
import numpy as np


def CombineTypes(row):

    # Prioritizing features in a summarized column 

    for strand in ['Same', 'Oppo']:

        if row[strand+'_gene_type'] == '.' and row[strand+'_Repeat'] == '.' and str(row[strand+'_intron_number']) in ['.', ' ', ''] and row[strand+'_ncRNA'] == '.':
            row[strand+'_Comb_type'] = 'no Feature'

        if row[strand+'_gene_type'] == 'lncRNA':
            if str(row[strand+'_intron_number']) not in ['.', ' ', '']:
                row[strand+'_Comb_type'] = 'lncRNA: Intronic'

        if row[strand+'_gene_type'] == 'protein_coding':
            if str(row[strand+'_intron_number']) not in ['.', ' ', '']:
                row[strand+'_Comb_type'] = 'protein_coding: Intronic'

        if row[strand+'_gene_type'] == 'lncRNA':
            if str(row[strand+'_exon_number']) not in ['.', ' ', '']:
                row[strand+'_Comb_type'] = 'lncRNA: Exonic'

        if row[strand+'_gene_type'] == 'processed_pseudogene':
            row[strand+'_Comb_type'] = 'pseudogene'

        if row[strand+'_Repeat'] not in ['.', ' ', '', 'lncRNA']:
            row[strand+'_Comb_type'] = 'Repeat Element'

        if row[strand+'_gene_type'] == 'protein_coding':
            if str(row[strand+'_exon_number']) not in ['.', ' ', '']:
                row[strand+'_Comb_type'] = 'protein_coding: Exonic'

        if row[strand+'_ncRNA'] != '.':
            row[strand+'_Comb_type'] = 'ncRNA' 

    return row



def MergeDataframes(SameStrandGen_df, SameStrandIntrons_df, SameStrandRMSK_df, SameStrandNCRNA_df, OppoStrandGen_df, OppoStrandIntrons_df, OppoStrandRMSK_df, OppoStrandNCRNA_df):

    # Create empty dataframe
    Out_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background',\
                                 'Same_ensembl_gene_id',   'Same_external_gene_name',    'Same_gene_type',\
                                 'Same_transcript_type',   'Same_transcript_name',       'Same_feature',\
                                 'Same_exon_number',       'Same_intron_number',         'Same_intron_Overlap',\
                                 'Same_exon_Overlap',      'Same_Repeat',                'Same_ncRNA',           'Same_Comb_type',\
                                 'Oppo_ensembl_gene_id',   'Oppo_external_gene_name',    'Oppo_gene_type',\
                                 'Oppo_transcript_type',   'Oppo_transcript_name',       'Oppo_feature',\
                                 'Oppo_exon_number',       'Oppo_intron_number',         'Oppo_intron_Overlap',\
                                 'Oppo_exon_Overlap',      'Oppo_Repeat',                'Oppo_ncRNA',           'Oppo_Comb_type'])


    SameStrandMerged_df = SameStrandGen_df.merge(SameStrandIntrons_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])
    SameStrandMerged_df = SameStrandMerged_df.merge(SameStrandRMSK_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])
    SameStrandMerged_df = SameStrandMerged_df.merge(SameStrandNCRNA_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])

    OppoStrandMerged_df = OppoStrandGen_df.merge(OppoStrandIntrons_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])
    OppoStrandMerged_df = OppoStrandMerged_df.merge(OppoStrandRMSK_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])
    OppoStrandMerged_df = OppoStrandMerged_df.merge(OppoStrandNCRNA_df, how='outer', on=['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background'])
    
    Out_df[['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background']] = SameStrandMerged_df[['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background']]
    Out_df['Same_ensembl_gene_id'] = SameStrandMerged_df['gene_id']
    Out_df['Same_external_gene_name'] = SameStrandMerged_df['gene_name']
    Out_df['Same_gene_type'] = SameStrandMerged_df['gene_type']
    Out_df['Same_transcript_type'] = SameStrandMerged_df['transcript_type']
    Out_df['Same_transcript_name'] = SameStrandMerged_df['transcript_name']
    Out_df['Same_feature'] = SameStrandMerged_df['feature']
    Out_df['Same_exon_number'] = SameStrandMerged_df['exon_number']
    Out_df['Same_intron_number'] = SameStrandMerged_df['intron_number']
    #Out_df['Same_intron_LargeOL'] = SameStrandMerged_df['???'] FROM BEDTOOLS INTERSECT?
    #Out_df['Same_exon_LargeOL'] = SameStrandMerged_df['???'] FROM BEDTOOLS INTERSECT?
    Out_df['Same_Repeat'] = SameStrandMerged_df['repClass']
    Out_df['Same_ncRNA'] = SameStrandMerged_df['ncrnaName']

    Out_df['Oppo_ensembl_gene_id'] = OppoStrandMerged_df['gene_id']
    Out_df['Oppo_external_gene_name'] = OppoStrandMerged_df['gene_name']
    Out_df['Oppo_gene_type'] = OppoStrandMerged_df['gene_type']
    Out_df['Oppo_transcript_type'] = OppoStrandMerged_df['transcript_type']
    Out_df['Oppo_transcript_name'] = OppoStrandMerged_df['transcript_name']
    Out_df['Oppo_feature'] = OppoStrandMerged_df['feature']
    Out_df['Oppo_exon_number'] = OppoStrandMerged_df['exon_number']
    Out_df['Oppo_intron_number'] = OppoStrandMerged_df['intron_number']
    #Out_df['Oppo_intron_LargeOL'] = OppoStrandMerged_df['???'] FROM BEDTOOLS INTERSECT?
    #Out_df['Oppo_exon_LargeOL'] = OppoStrandMerged_df['???'] FROM BEDTOOLS INTERSECT?    
    Out_df['Oppo_Repeat'] = OppoStrandMerged_df['repClass']
    Out_df['Oppo_ncRNA'] = OppoStrandMerged_df['ncrnaName']

    Out_df = Out_df.apply(CombineTypes, axis=1)

    return(Out_df)

def ProcessGencode(InputFile):
    Gen_df = pd.read_csv(InputFile, sep = '\t', header = None,\
                names = ['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background',\
                         'genChrom', 'genStart', 'genEnd', 'Source', 'genScore', 'genStrand', 'feature', 'frame',\
                         'gene_id', 'gene_type', 'gene_name', 'level', 'mgi_id', 'havana_gene', 'transcript_id',\
                         'transcript_type', 'transcript_name', 'transcript_support_level', 'tag', 'havana_transcript',\
                         'exon_number', 'exon_id', 'protein_id', 'ccdsid', 'ont', 'genOverlap'])

    Gen_df.iloc[:, 10:] = Gen_df.iloc[:, 10:].astype('str') # change type necessary for aggregate
    #print(Gen_df.dtypes)

    # Aggregate reads found in multiple entries
    # https://stackoverflow.com/questions/36271413/pandas-merge-nearly-duplicate-rows-based-on-column-value
    Gen_df = Gen_df.groupby(['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group','normalized_read_density_sample',\
                             'normalized_read_density_background']).agg({   'genChrom'          : ', '.join, 'genStart'          : ', '.join, 'genEnd'                   : ', '.join,
                                                                            'Source'            : ', '.join, 'genScore'          : ', '.join, 'genStrand'                : ', '.join,
                                                                            'feature'           : ', '.join, 'frame'             : ', '.join, 'gene_id'                  : ', '.join,
                                                                            'gene_type'         : ', '.join, 'gene_name'         : ', '.join, 'level'                    : ', '.join,
                                                                            'mgi_id'            : ', '.join, 'havana_gene'       : ', '.join, 'transcript_id'            : ', '.join,
                                                                            'transcript_type'   : ', '.join, 'transcript_name'   : ', '.join, 'transcript_support_level' : ', '.join,
                                                                            'tag'               : ', '.join, 'havana_transcript' : ', '.join, 'exon_number'              : ', '.join,
                                                                            'exon_id'           : ', '.join, 'protein_id'        : ', '.join, 'ccdsid'                   : ', '.join,
                                                                            'ont'               : ', '.join, 'genOverlap'        : '; '.join 
                                                                            }).reset_index()

    # Remove duplicate terms in cells, remove 'nan'
    Gen_df.iloc[:, 10:] = Gen_df.iloc[:, 10:].applymap(lambda x: ', '.join(set(x.strip().replace(' ', '').replace('nan,', '').replace('nan', '').split(','))))
    return(Gen_df)


def ProcessIntrons(InputFile):
    Introns_df = pd.read_csv(InputFile, sep = '\t', header = None,\
                names = ['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background',\
                         'intronChrom', 'intronStart', 'intronEnd', 'intronID', 'intronScore', 'intronStrand', 'intronOverlap', 'intron_number'])

    #Introns_df = Introns_df.astype('str') # change type necessary for aggregate
    Introns_df.iloc[:, 10:] = Introns_df.iloc[:, 10:].astype('str') # change type necessary for aggregate
   
    # https://stackoverflow.com/questions/36271413/pandas-merge-nearly-duplicate-rows-based-on-column-value
    Introns_df = Introns_df.groupby(['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group','normalized_read_density_sample',\
                                     'normalized_read_density_background']).agg({'intronChrom'   : ', '.join, 'intronStart'   : ', '.join, 'intronEnd'    : ', '.join,
                                                                                 'intronID'      : ', '.join, 'intronScore'   : ', '.join, 'intronStrand' : ', '.join,
                                                                                 'intronOverlap' : ', '.join, 'intron_number' : ', '.join }).reset_index()

    # Remove duplicate terms in cells, remove 'nan'
    Introns_df.iloc[:, 10:] = Introns_df.iloc[:, 10:].applymap(lambda x: ', '.join(set(x.strip().replace(' ', '').replace('nan,', '').replace('nan', '').split(','))))
    return(Introns_df)


def ProcessRMSK(InputFile):
    RMSK_df = pd.read_csv(InputFile, sep = '\t', header = None, \
                names = ['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample', 'normalized_read_density_background',\
                         'repChrom', 'repStart', 'repEnd', 'repName', 'repClass', 'repStrand', 'repFamily', 'repOverlap'])

    #RMSK_df = RMSK_df.astype('str') # change type necessary for aggregate
    RMSK_df.iloc[:, 10:] = RMSK_df.iloc[:, 10:].astype('str') # change type necessary for aggregate

    # https://stackoverflow.com/questions/36271413/pandas-merge-nearly-duplicate-rows-based-on-column-value
    RMSK_df = RMSK_df.groupby(['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group','normalized_read_density_sample',\
                               'normalized_read_density_background']).agg({ 'repStart'  : ', '.join, 'repEnd'   : ', '.join, 'repName' : ', '.join,
                                                                            'repClass'  : ', '.join, 'repStrand': ', '.join, 'repFamily': ', '.join,
                                                                            'repOverlap': ', '.join,}).reset_index()

    # Remove duplicate terms in cells, remove 'nan'
    RMSK_df.iloc[:, 10:] = RMSK_df.iloc[:, 10:].applymap(lambda x: ', '.join(set(x.strip().replace(' ', '').replace('nan,', '').replace('nan', '').split(','))))
    return(RMSK_df)


def ProcessNCRNA(InputFile):
    NCRNA_df = pd.read_csv(InputFile, sep = '\t', header = None, \
                names = ['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group', 'normalized_read_density_sample',\
                         'normalized_read_density_background', 'ncrnaChrom', 'ncrnaStart', 'ncrnaEnd', 'ncrnaName', 'ncrnaClass', 'ncrnaStrand', 'ncrnaOverlap'])

    #RMSK_df = RMSK_df.astype('str') # change type necessary for aggregate
    NCRNA_df.iloc[:, 10:] = NCRNA_df.iloc[:, 10:].astype('str') # change type necessary for aggregate

    # https://stackoverflow.com/questions/36271413/pandas-merge-nearly-duplicate-rows-based-on-column-value
    NCRNA_df = NCRNA_df.groupby(['chrom', 'start', 'end', 'type', 'score', 'strand', 'm_value', 'p_value', 'peak_group','normalized_read_density_sample',\
                                 'normalized_read_density_background']).agg({  'ncrnaChrom'  : ', '.join, 'ncrnaStart': ', '.join, 'ncrnaEnd'   : ', '.join,
                                                                                'ncrnaName'   : ', '.join, 'ncrnaClass': ', '.join, 'ncrnaStrand': ', '.join,
                                                                                'ncrnaOverlap': ', '.join,}).reset_index()

    # Remove duplicate terms in cells, remove 'nan'
    NCRNA_df.iloc[:, 10:] = NCRNA_df.iloc[:, 10:].applymap(lambda x: ', '.join(set(x.strip().replace(' ', '').replace('nan,', '').replace('nan', '').split(','))))
    return(NCRNA_df)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--SameStrandRMSK', type=str, required=True)
    parser.add_argument('--SameStrandGenCode', type=str, required=True)
    parser.add_argument('--SameStrandIntrons', type=str, required=True)
    parser.add_argument('--SameStrandncRNA', type=str, required=True)
    parser.add_argument('--OppoStrandRMSK', type=str, required=True)
    parser.add_argument('--OppoStrandGenCode', type=str, required=True)
    parser.add_argument('--OppoStrandIntrons', type=str, required=True)
    parser.add_argument('--OppoStrandncRNA', type=str, required=True)
    parser.add_argument('--Output', type=str, required=True)


    args = parser.parse_args()
    SameStrandRMSK_df = ProcessRMSK(args.SameStrandRMSK)
    SameStrandGen_df = ProcessGencode(args.SameStrandGenCode)
    SameStrandIntrons_df = ProcessIntrons(args.SameStrandIntrons)
    SameStrandNCRNA_df = ProcessNCRNA(args.SameStrandncRNA)

    OppoStrandRMSK_df = ProcessRMSK(args.OppoStrandRMSK)
    OppoStrandGen_df = ProcessGencode(args.OppoStrandGenCode)
    OppoStrandIntrons_df = ProcessIntrons(args.OppoStrandIntrons)
    OppoStrandNCRNA_df = ProcessNCRNA(args.OppoStrandncRNA)

    Out_df = MergeDataframes(SameStrandGen_df, SameStrandIntrons_df, SameStrandRMSK_df, SameStrandNCRNA_df, OppoStrandGen_df, OppoStrandIntrons_df, OppoStrandRMSK_df, OppoStrandNCRNA_df)
    Out_df.to_csv(args.Output, sep = '\t', index=False)

if __name__ == '__main__':
    main()
