"""
Part 8 of TP53 mutation status identification from single-cell amplicon sequencing (scAmp-seq) data.

Define cells (barcodes) as WT or MUT based on genotype likelihoods (GL) or mpileup majority (MM)
"""


import argparse
import csv
import os

import pandas as pd
import numpy as np
import regex as re


# Directories

dir_metadata = './metadata'


# Functions

def parse_parameters(file_param):
    """
    Parse TP53 variant parameters
    """

    des = []
    data = []

    with open(file_param, newline='') as f:
        reader = csv.reader(f)
        data_in = list(reader)

    for i in data_in:
        if des == []:
            des = i
        else:
            data = i

    sample = data[0]
    variant = data[1]

    # Variant type [sub, del, ins]
    variant_type = ''
    if bool(re.search('>', variant)):
        variant_type = 'sub'
    elif bool(re.search('del', variant)):
        variant_type = 'del'
    elif bool(re.search('ins', variant)):
        variant_type = 'ins'

    # Bases [WT,MUT]
    if variant_type=='sub':
        bases = [variant.split('>')[:1][0][-1], variant.split('>')[1:][0]]
    elif variant_type=='del':
        bases = [variant.split('del')[1:][0], '']
    elif variant_type=='ins':
        bases = ['', variant.split('ins')[1:][0]]

    print('sample: ' + sample)
    print('variant: ' + variant)
    print('variant_type: ' + variant_type)
    print('WT base:  ' + bases[0])
    print('MUT base: ' + bases[1])

    return variant, variant_type, bases


def parse_pos_GRCh38(variant):
    """
    Parse TP53 variant position in reference genome
    """
    
    # Parse TP53 variant,primer ID -pair informations
    file_name_mut_blocks = 'samples_variant_primer_ID.csv'
    file_mut_blocks = os.path.join(dir_metadata, file_name_mut_blocks)
    df_tp53_variants = pd.read_csv(file_mut_blocks)

    variant_infos = df_tp53_variants[(df_tp53_variants['TP53 c. mutation'] == variant)]

    pos_GRCh38 = variant_infos['GRCh38'].values.tolist()[0]

    print(pos_GRCh38)

    return pos_GRCh38


def get_lines_that_contain_position(string, fp):
    """
    Get all lines in file that contain the given position
    """

    return [line for line in fp if string in line]


def add_wtmut_to_dict(bc_wtmut, bc, mut, DP):
    """
    Add current WT/MUT status to barcode (bc) in dictionary
    """

    if bc in bc_wtmut:
        bc_wtmut[bc][mut] = bc_wtmut[bc][mut] + 1
        bc_wtmut[bc][2] = bc_wtmut[bc][2] + int(DP)

    else: 
        wtmut = [0,0,0]
        wtmut[mut] = wtmut[mut] + 1
        wtmut[2] = wtmut[2] + int(DP)
        bc_wtmut.update({bc : wtmut})

    return bc_wtmut


def get_genotypes(folder_bc, input_suffix, pos_GRCh38, bases, variant_type):
    """
    Get genotype (GT) of highest genotype likelihood (GL)
    """

    bc_wtmut = dict()

    GL_of_WT = []
    GL_of_MUT = []

    n_reads_wt = 0
    n_reads_mut = 0

    for filename in os.listdir(folder_bc):
        file = os.path.join(folder_bc, filename)

        # Check if it is a file and file ends with .haploid.vcf
        if os.path.isfile(file) and filename.endswith(input_suffix):

            prefix=filename.split('.')[0]
            bc = prefix.split('_')[0]
            umi = prefix.split('_')[1]

            print(bc)

            # Get GT (of highest GL) at variant position

            with open(file, "r") as fp:
                for line in get_lines_that_contain_position(pos_GRCh38.split(':')[1], fp):

                    print(line)

                    line_as_list = re.split(r'\t+', line.rstrip())

                    if input_suffix == '.haploid.vcf':
        
                        GT_PL = line_as_list[len(line_as_list)-1]
                        
                        if ":" in GT_PL:
                            GT = GT_PL.split(':')[0]
                            PL = GT_PL.split(':')[1]
                            PL_splitted = PL.split(',')
                            GL = pow(10, (-int(PL_splitted[int(GT)])/10) )
                        else:
                            GT='0'
                            GL=1

                        DP = line_as_list[7].split(';')[0].split('=')[1]

                        if int(DP)>1: # Check once again in case BCFtools mpileup filters more reads
                            if GT == '0' or GT == '1':
                                GL_of_MUT.append(GL)
                                bc_wtmut = add_wtmut_to_dict(bc_wtmut, bc, int(GT), DP)
                            else:
                                print('GT was either 0 nor 1 in ' + filename)
                            
                            DP4 = [s for s in line_as_list[7].split(';') if "DP4" in s]
                            DP4_list = DP4[0].split('=')[1].split(',')
                            n_wt = int(DP4_list[0]) + int(DP4_list[1])
                            n_mut = int(DP4_list[2]) + int(DP4_list[3])

                            n_reads_wt = n_reads_wt + n_wt
                            n_reads_mut = n_reads_mut + n_mut

                    elif input_suffix == '.mpileup':
                        
                        bases_dict = {'A':'t', 'T':'a', 'G':'c', 'C':'g', '':'*', 'CCTCTTGCT':'c+9agcaagagg'}

                        bases_reads = line_as_list[len(line_as_list)-2]

                        n_wt = bases_reads.count(bases_dict[bases[0]])

                        if bases[1] == 'X': # benchmark (only substitutions)
                            n_mut = len(bases_reads) - n_wt
                        else:
                            n_mut = bases_reads.count(bases_dict[bases[1]])

                        if variant_type=='ins':
                            n_wt = n_wt - (n_mut * bases_dict[bases[1]].count(bases_dict[bases[0]]))

                        DP = n_mut+n_wt
                        if  DP>1: # Check once again in case SAMtools mpileup filters more reads
                            if n_mut==0 and n_wt==0:
                                print('Neither WT nor MUT bases detected at variant position in ' + filename)
                            elif n_mut > n_wt:
                                bc_wtmut = add_wtmut_to_dict(bc_wtmut, bc, 1, n_wt+n_mut)
                            else: # i.e. n_mut <= n_wt:
                                bc_wtmut = add_wtmut_to_dict(bc_wtmut, bc, 0, n_wt+n_mut)

                            GL_of_MUT.append(0) # NO GL in case of mpileup majority pipeline. For now, just append zeros.
                                
                            n_reads_wt = n_reads_wt + n_wt
                            n_reads_mut = n_reads_mut + n_mut

    return bc_wtmut, GL_of_WT, GL_of_MUT, n_reads_wt, n_reads_mut


def write_genotypes_to_file(bc_transcripts_wtmut, folder_results, output_suffix):
    """
    Write genotypes to file
    """

    columns=['barcode','%mut','#wt','#mut','#transcripts','#reads']
    df_bc = pd.DataFrame.from_dict(bc_transcripts_wtmut, orient='index', columns=['#wt','#mut','#reads'])
    df_bc['#transcripts'] = df_bc['#wt'] + df_bc['#mut']
    df_bc['%mut'] = df_bc['#mut'] / df_bc['#transcripts']
    df_bc['barcode'] = df_bc.index
    df_bc = df_bc.reindex(columns, axis=1)

    file_bc = folder_results + 'barcodes_'+output_suffix+'.csv'
    df_bc.to_csv(file_bc, index=False, line_terminator='\n')

    return df_bc


def write_genotype_threshold_to_file(df_bc, folder_results, output_suffix):
    """
    Generate barcodes below/above thresholds (WT/MUT and WT/HET/HOM) and write to file
    """

    df_bc['barcode'] = df_bc['barcode'] + '-1'
    df_bc['%mut'] = df_bc['%mut'].astype(float)

    mut = df_bc['%mut']

    cond_list = [mut <= 0.1, mut >= 0.3, (mut>0.1) & (mut<0.3)]
    choice_list = ['0', '1', 'NA'] # WT, MUT, NA
    df_bc['TP53_gt_WTMUT'] = np.select(cond_list, choice_list)

    cond_list = [mut <= 0.1, (mut >= 0.4) & (mut <= 0.6), mut >= 0.9, ((mut>0.1) & (mut<0.4)) | ((mut>0.6) & (mut<0.9))]
    choice_list = ['0', '1', '2', 'NA'] # WT, HET, HOM, NA
    df_bc['TP53_gt_WTHETHOM'] = np.select(cond_list, choice_list)

    df_bc_wtmut = df_bc[df_bc['TP53_gt_WTMUT'] != 'NA'][['barcode', 'TP53_gt_WTMUT']]
    df_bc_wthethom = df_bc[df_bc['TP53_gt_WTHETHOM'] != 'NA'][['barcode', 'TP53_gt_WTHETHOM']]

    file_bc_wtmut = folder_results + 'barcodes_'+output_suffix+'_wtmut.csv'
    df_bc_wtmut[['barcode', 'TP53_gt_WTMUT']].to_csv(file_bc_wtmut, index=False, line_terminator='\n')

    file_bc_wthethom = folder_results + 'barcodes_'+output_suffix+'_wthethom.csv'
    df_bc_wthethom[['barcode', 'TP53_gt_WTHETHOM']].to_csv(file_bc_wthethom, index=False, line_terminator='\n')

    return df_bc_wtmut, df_bc_wthethom


def write_stats_gt(df_bc, n_reads_wt, n_reads_mut, df_bc_wtmut, df_bc_wthethom, folder_results, output_suffix):
    """
    Write genotype statistics to file
    """

    data=[]
    data.append([
        'Reads', 
        float(df_bc['#reads'].sum()), 
        n_reads_wt,
        n_reads_mut,
        0.0, 
        0.0
    ])
    data.append([
        'Transcripts', 
        float(df_bc['#transcripts'].sum()), 
        float(df_bc['#wt'].sum()), 
        float(df_bc['#mut'].sum()), 
        0.0, 
        0.0,
    ])
    data.append([
        'Barcodes', 
        float(len(df_bc)), 
        float(len(df_bc_wtmut[df_bc_wtmut['TP53_gt_WTMUT'] == '0'])), 
        float(len(df_bc_wtmut[df_bc_wtmut['TP53_gt_WTMUT'] == '1'])), 
        float(len(df_bc_wthethom[df_bc_wthethom['TP53_gt_WTHETHOM'] == '1'])), 
        float(len(df_bc_wthethom[df_bc_wthethom['TP53_gt_WTHETHOM'] == '2']))
    ])

    columns=['Type','#genotyped','#genotyped-WT','#genotyped-MUT','#genotyped-HET','#genotyped-HOM']
    df_stats_gt = pd.DataFrame(data, columns=columns)

    # Calculate fractions of genotyping stats
    for col in ['genotyped-WT', 'genotyped-MUT', 'genotyped-HET', 'genotyped-HOM']:
        df_stats_gt['%'+ col] = df_stats_gt['#' + col] / df_stats_gt['#genotyped']

    # Write genotyping stats to file
    file_stats_gt_all = folder_results + 'stats_'+output_suffix+'_wtmut_wthethom.csv'
    df_stats_gt.to_csv(file_stats_gt_all, index=True, line_terminator='\n')


def write_mean_std(folder_results, output_suffix, GL_of_WT, GL_of_MUT):
    """
    # Write mean and standard deviation of GL_of_WT, GL_of_MUT to file
    """

    file_name = folder_results + 'likelihoods_'+output_suffix+'.csv'
    with open(file_name, 'w') as file_handle:
        writer = csv.writer(file_handle)
        writer.writerow(['mean(GL_of_WT)', 'std(GL_of_WT)', 'mean(GL_of_MUT)', 'std(GL_of_MUT)'])
    try:
        file_handle = open(file_name, 'a')
        writer = csv.writer(file_handle)
        writer.writerow([np.mean(GL_of_WT), np.std(GL_of_WT), np.mean(GL_of_MUT), np.std(GL_of_MUT)])

    finally:
        file_handle.close()


# Main

def main():

    # Command line parameters

    parser = argparse.ArgumentParser(
        prog = 'python scAmp-8-define-cells-as-WT-MUT.py',
        description = 'Part 8 of TP53 mutation status identification from scAmp-seq data: \
                       Define cells (barcodes) as WT or MUT'
    )

    parser.add_argument('-amp', '--scampparam', type=str, required=True, help='scAmp parameter file')
    parser.add_argument('-b', '--bcs', type=str, required=True, help='bc folder (with haploid.vcf files per transcript)')
    parser.add_argument('-r', '--results', type=str, required=True, help='(final) results folder')
    parser.add_argument('-i', '--input',  type=str, required=True, help='files of which type should be used as input (.mpileup, .haploid.vcf)')
    parser.add_argument('-s', '--suffix',  type=str, required=True, help='which suffix to use for output files')
    parser.add_argument('-g', '--grch38',  type=str, required=False, help='GRCh38 position of variant')

    args = parser.parse_args()

    file_param = args.scampparam
    folder_results = args.results
    folder_bc = args.bcs

    input_suffix = args.input
    output_suffix = args.suffix

    if args.grch38 is not None:
        pos_GRCh38 = args.grch38
    else:
        pos_GRCh38 = ''


    # Identify variant position for given variant
    variant, variant_type, bases = parse_parameters(file_param)

    if args.grch38 is None:
        pos_GRCh38 = parse_pos_GRCh38(variant)


    # Iterate over *haploid.vcf -files in given bc folder
    bc_transcripts_wtmut, GL_of_WT, GL_of_MUT, n_reads_wt, n_reads_mut = get_genotypes(folder_bc, input_suffix, pos_GRCh38, bases, variant_type)


    # Define bc as wt/mut and wt/het/mut
    df_bc = write_genotypes_to_file(bc_transcripts_wtmut, folder_results, output_suffix)
    df_bc_wtmut, df_bc_wthethom = write_genotype_threshold_to_file(df_bc, folder_results, output_suffix)


    # Write mean and standard deviation of GL_of_WT, GL_of_MUT to file
    write_mean_std(folder_results, output_suffix, GL_of_WT, GL_of_MUT)

    # Write statistics to file
    write_stats_gt(df_bc, n_reads_wt, n_reads_mut, df_bc_wtmut, df_bc_wthethom, folder_results, output_suffix)

    print('Finished')


if __name__ == "__main__":
    main()
