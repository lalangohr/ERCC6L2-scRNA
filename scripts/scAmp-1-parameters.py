"""
Part 1 of TP53 mutation status identification from single-cell amplicon sequencing (scAmp-seq) data.

Parse and specify input parameters and write them into a file to be used in scAmp-seq data preprocessing.
"""


import argparse
import os

import pandas as pd
import regex as re

from time import gmtime, strftime


# Directories

dir_metadata = './metadata'
dir_results = './results'


# Functions

def parse_tp53_reference():
    """
    Parse TP53 reference sequence
    """
    
    file_name_cds = 'TP53_CDS_NM_000546-5.fasta'
    file_cds = os.path.join(dir_metadata, file_name_cds)

    f1=open(file_cds,'r')
    for line in f1:
        tp53_ref = line
    f1.close()

    return tp53_ref


def parse_primer_seq(tp53_ref, scAmp, read_len):
    """
    Parse used primer sequence and its length
    """

    # Get length of current primer
    file_name_primers = 'TP53_primers.csv'
    file_primers = os.path.join(dir_metadata, file_name_primers)
    df_primers = pd.read_csv(file_primers)
    primer_infos = df_primers[(df_primers['primer name'] == scAmp)].values.tolist()
    primer_seq_used = primer_infos[0][2]
    primer_len = len(primer_seq_used)

    # Extract substring of TP53 CDS reference after primer: the subsequence of TP53 CDS that is sequenced
    start_pos = int(primer_infos[0][1])+primer_len-1
    tp53_ref_sub = tp53_ref[start_pos:start_pos+read_len-primer_len]

    # Primer start position in TP53 CDS
    primer_pos = int(primer_infos[0][1])

    return primer_seq_used, primer_len, tp53_ref_sub, primer_pos


def get_variant_type(variant):
    """
    Identify whether variant is a substition, insertion or deletion
    """
    
    variant_type = ''
    if bool(re.search('>', variant)):
        variant_type = 'sub'
    elif bool(re.search('del', variant)):
        variant_type = 'del'
    elif bool(re.search('ins', variant)):
        variant_type = 'ins'

    return variant_type


def calculate_variant_position(variant, variant_type, primer_len, primer_pos):

    if variant_type=='sub' or variant_type=='del':

        variant_pos_from_start = int(re.findall('\d+', variant)[0]) # Variant position in TP53 CDS
        variant_pos = str(variant_pos_from_start - primer_pos - primer_len) # Variant position in read

    elif variant_type=='ins':
         # Variant position in TP53 CDS
        variant_start_pos_from_start = int(re.findall('\d+', variant.split('_')[:1][0])[0])
        variant_end_pos_from_start = int(re.findall('\d+', variant.split('_')[1:][0])[0])
    
        # Variant position in read
        variant_start_pos = variant_start_pos_from_start - primer_pos - primer_len
        variant_end_pos = variant_end_pos_from_start - primer_pos - primer_len
        variant_pos = str(variant_start_pos) + '_' + str(variant_end_pos)

    return variant_pos


# Main

def main():

    # Command line parameters
    parser = argparse.ArgumentParser(
        prog = 'python scAmp-1-parameters.py',
        description = 'Part 1 of TP53 mutation status identification from scAmp-seq data: \
                       Parse metadata and specify input parameters for TP53 mutation status identification'
    )
    parser.add_argument('-s', '--sample',  type=str, required=True, help='sample id')
    parser.add_argument('-v', '--variant', type=str, required=True, help='variant')
    parser.add_argument('-l', '--length',  type=int, default=250,   help='read length (default: 250)')

    args = parser.parse_args()

    sample = args.sample
    variant = args.variant
    read_length = args.length


    tp53_ref = parse_tp53_reference()

    # Parse TP53 variant,primer ID -pair informations
    file_name_mut_blocks = 'samples_variant_primer_ID.csv'
    file_mut_blocks = os.path.join(dir_metadata, file_name_mut_blocks)
    df_tp53_variants = pd.read_csv(file_mut_blocks)

    # Get TP53 variant position infos of a specific variant and sample
    df_variant_infos = df_tp53_variants[(df_tp53_variants['sample'] == sample) & \
                                        (df_tp53_variants['TP53 c. mutation'] == variant)]

    for _, row in df_variant_infos.iterrows():

        primer_id = row['primer ID']

        variant_in_file_name = variant.replace('.','').replace('>','')
        time_in_file_name = strftime('%Y%m%d-%H%M%S', gmtime())
        file_name = 'scAmp-parameters-' + primer_id + '-' + sample + '-' + variant_in_file_name + '-' + time_in_file_name + '.csv'
        file_out = os.path.join(dir_results, file_name)

        primer_seq_used, primer_len, tp53_ref_sub, primer_pos = parse_primer_seq(tp53_ref, primer_id, read_length)

        variant_type = get_variant_type(variant)

        variant_pos = calculate_variant_position(variant, variant_type, primer_len, primer_pos)

        # Parse R1 file locations
        file_name_scAmp_paths = 'samples_primer_ID_scAmp-seq_R1_files.csv'
        file_scAmp_paths = os.path.join(dir_metadata, file_name_scAmp_paths)
        df_files = pd.read_csv(file_scAmp_paths)

        # Specify R1 and R2 files
        df_R1 = df_files[((df_files['sample'] == sample) & (df_files['primer ID'] == primer_id))]
        file_R1 = df_R1['scAmp-seq R1 file'].iloc[0]
        file_R2 = file_R1.replace('R1','R2')

        # Combine all parameters needed
        col = ['sample', 'variant', 'variant position in substring',
              'primer seq', 'substring of TP53 CDS after primer', 'R1 file', 'R2 file']
        data = [[sample, variant, variant_pos,
                primer_seq_used, tp53_ref_sub, file_R1, file_R2]]

        df_out = pd.DataFrame(data, columns=col)
        df_out.to_csv(file_out, index=False)

        print('Result written to: ' + file_out)


if __name__ == "__main__":
    main()
