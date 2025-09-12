"""
Part 2 of TP53 mutation status identification from single-cell amplicon sequencing (scAmp-seq) data.

Perform quality control (QC) on scAmp-seq data, consisting of the following steps:
1) Parse and count raw data
2) Barcode QC
3) Sequencing QC
4) Primer QC
5) Read QC
 
Requires: 
1) Paramaters file written by scAm-part1-parameters-NovaSeq.py and 
2) Barcode whitelist file.
Produces:
1) Reads splitted into k files such that indentical barcodes are found in same file
2) QC statistics
"""


import sys
import csv
import gzip
import math
import argparse

import pandas as pd
import regex as re

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.QualityIO import FastqPhredIterator
from statistics import mean


# Classes

class ScAmp_counts: 
    """
    Class of scAmp-seq QC counts
    """

    num_all = 0
    num_after_bccorr = 0
    num_after_seqqc = 0
    num_after_primerqc = 0
    num_after_readqc = 0


class ScAmp_stats:
    """
    Class of scAmp-seq QC statistics
    """
    
    # Create instances of classes
    bc = ScAmp_counts()
    transcript = ScAmp_counts()
    read = ScAmp_counts()

    def write_to_file(self, file):
        """
        Write statistics to file
        """

        lines = []
        lines.append(['Step', '#bc', '#transcript', '#read'])
        lines.append(['1) Raw data', self.bc.num_all, self.transcript.num_all, self.read.num_all])
        lines.append(['2) Barcode QC', self.bc.num_after_bccorr, self.transcript.num_after_bccorr, self.read.num_after_bccorr])
        lines.append(['3) Sequencing QC', self.bc.num_after_seqqc, self.transcript.num_after_seqqc, self.read.num_after_seqqc])
        lines.append(['4) Primer QC', self.bc.num_after_primerqc, self.transcript.num_after_primerqc, self.read.num_after_primerqc])
        lines.append(['5) Read QC', self.bc.num_after_readqc, self.transcript.num_after_readqc, self.read.num_after_readqc])
        
        with open(file, 'w') as f:
            writer = csv.writer(f)
            for line in lines:
                writer.writerow(line)
    

# Functions

def add_transcript_to_dict(transcripts, bc, umi):
    """
    Add unique (bc, umi) -pairs to transcripts -dict
    """

    if bc in transcripts:
        if umi not in transcripts[bc]:
            transcripts[bc] = transcripts[bc] + [umi]
    else:
        transcripts[bc] = [umi]

    return transcripts


def parse_bcs_perform_qc(file, split, local_stats, whitelist, folder_tmp):
    """
    Parse barcodes and perform barcode QC
    """
    
    bc_occ = dict()
    transcripts = dict()

    num_bc_all = 0
    num_bc_after_bccorr = 0
    num_reads_all = 0

    # Parse whitelist
    whitelist_bc=[]
    with gzip.open(whitelist,'rt') as f:
        for bc in f:
            whitelist_bc.append(bc.replace('\n',''))

    # Parse barcodes (bcs)
    with gzip.open(file, "rt") as handle:
        for i, record in enumerate(FastqPhredIterator(handle)):
            
            if (i % 1000) == 0:
                sys.stderr.write('Line {}, list: {}\n'.format(i, len(bc_occ))); sys.stdout.flush()
            
            bc = str(record.seq[:split])
            umi = str((record.seq[split:]))

            if bc in bc_occ:
                bc_occ[bc] = bc_occ[bc] + 1
            else:
                bc_occ[bc] = 1
            
            transcripts = add_transcript_to_dict(transcripts, bc, umi)
    

    # Barcode QC

    set_bc = set(bc_occ.keys())
    set_exact_bc = set_bc.intersection(whitelist_bc)
    list_nonexact_bc = list(set_bc - set_exact_bc)

    # Count number of bcs
    num_bc_all = len(set_bc)
    num_bc_after_bccorr = len(set_exact_bc)

    # Count number of transcript
    num_transcripts_all = 0
    for transcripts_of_one_bc in transcripts.values():
        num_transcripts_all = num_transcripts_all + len(transcripts_of_one_bc)

    [transcripts.pop(key) for key in list_nonexact_bc]

    num_transcripts_after_bccorr = 0
    for transcripts_of_one_bc in transcripts.values():
        num_transcripts_after_bccorr = num_transcripts_after_bccorr + len(transcripts_of_one_bc)

    # Count number of reads
    num_reads_all = sum(bc_occ.values())

    # Save statistics
    write_bc_read_distribution_to_file(bc_occ, 'all', folder_tmp)

    # Remove non-eaxct barcodes from dict bc_occ
    [bc_occ.pop(key) for key in list_nonexact_bc]

    num_reads_after_bccorr = sum(bc_occ.values())

    # Save statistics
    write_bc_read_distribution_to_file(bc_occ, 'after_bcqc', folder_tmp)

    # Statistics
    local_stats.bc.num_all = num_bc_all
    local_stats.bc.num_after_bccorr = num_bc_after_bccorr

    local_stats.transcript.num_all = num_transcripts_all
    local_stats.transcript.num_after_bccorr = num_transcripts_after_bccorr

    local_stats.read.num_all = num_reads_all
    local_stats.read.num_after_bccorr = num_reads_after_bccorr


    return bc_occ, local_stats


def write_bc_read_distribution_to_file(bc_occ, postfix, folder_tmp):
    """
    Write barcode read distribution to file
    """
    
    df_bc = pd.DataFrame.from_dict(bc_occ, orient='index', columns=['n_reads'])
    df_bc = df_bc.sort_values(by=['n_reads'], ascending=False)

    df_bc.to_csv(folder_tmp+"bc_reads_"+postfix+".csv", index_label='barcode')


def sort_bcs_into_buckets(bc_occurrences, k):
    """
    Sort barcodes into k buckets
    """
    
    df_bc = pd.DataFrame.from_dict(bc_occurrences, orient='index', columns=['n_reads'])
    df_bc['n_cumsum'] = df_bc.sort_values(by=['n_reads'], ascending=True).n_reads.cumsum()

    k_all = k
    k_size_all = df_bc.n_cumsum.max() / k_all
    excessive_bcs = df_bc[df_bc['n_reads'] > k_size_all].index
    n_excessive_bcs = excessive_bcs.size

    k_rest = k_all - n_excessive_bcs
    k_size = df_bc[df_bc['n_reads'] <= k_size_all].n_cumsum.max() / k_rest

    # Assign buckets to bcs with amount of reads of maximum k% of all reads
    # Round by round(x*1000000)/1000000 to omit rounding problems causing bucket to be k_rest+1 due to ceil(k_rest.0000000001)=k_rest+1
    df_bc['bucket'] = (round((df_bc[df_bc['n_reads'] <= k_size_all].n_cumsum / k_size)*1000000)/1000000).apply(math.ceil)
    df_bc.sort_values(by=['n_reads'], ascending=True)

    # Assign buckets to bcs with amount of reads exceeding k% of all reads
    i=0
    excessive_bcs_buckets = []
    for bc in excessive_bcs:
        df_bc.loc[bc, 'bucket'] = k_all - i
        excessive_bcs_buckets.append(k_all - i)
        i=i+1

    return df_bc


def check_whether_primer_passes_primer_qc(primer, primer_used):
    """
    Check whether primer passes primer QC or not
    """
    
    cnt = 0
    for i in range(len(primer)):
        if primer[i] != primer_used[i]:
            cnt = cnt + 1
            if cnt > 1:
                break

    if cnt < 2:
        # Keep row as primer sequence differs max 1bp
        return True

    return False


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
    variant_pos_string = data[2] # variant position in substring
    primer_used = data[3]        # primer sequence (reference)
    tp53_ref_sub = data[4]       # substring of TP53 CDS reference
    file_R1 = data[5]
    file_R2 = data[6]

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

    primer_len = len(primer_used)

    print('sample: ' + sample)
    print('variant: ' + variant)
    print('variant_type: ' + variant_type)
    print('WT base:  ' + bases[0])
    print('MUT base: ' + bases[1])
    print('primer seq: ' + primer_used)
    print('primer has length: ' + str(primer_len))
    print('substring of TP53 CDS after primer: ' + tp53_ref_sub)
    print('substring has length: ' + str(len(tp53_ref_sub)))
    print('variant position in substring: ' + variant_pos_string + ' (counting starting from 0)')

    return variant, variant_type, bases, primer_used, tp53_ref_sub, variant_pos_string, file_R1, file_R2


def calculate_parameters(variant, variant_pos_string, tp53_ref_sub):
    """
    Calculate additional TP53 variant parameters
    """
    
    # Variant type [sub, del, ins]
    variant_type = ''
    if bool(re.search('>', variant)):
        variant_type = 'sub'
    elif bool(re.search('del', variant)):
        variant_type = 'del'
    elif bool(re.search('ins', variant)):
        variant_type = 'ins'

    if variant_type=='sub' or variant_type=='del':
        variant_pos = int(variant_pos_string)
    elif variant_type=='ins':
        variant_pos = int(variant_pos_string.split('_')[:1][0])

    if variant_type=='sub' or variant_type=='del':
        pos_start = variant_pos - 10
    elif variant_type=='ins':
        pos_start = variant_pos - 9

    if pos_start < 0:
        pos_start = 0
        
    pos_end = variant_pos + 11

    if pos_end >= len(tp53_ref_sub):
        pos_end = len(tp53_ref_sub)
        
    return variant_type, variant_pos, pos_start, pos_end


def check_whether_read_passes_read_qc(read, variant, variant_pos_string, tp53_ref_sub, bases):
    """
    Check whether read passes read QC or not
    Retains only reads which have
    * the WT or variant nucleotide(s) at the variant position(s), i.e.
      - the WT or MUT nucleotide in case variant is substitution (e.g. c.659A>G)
      - the WT or insertion nucleotides in case variant is insertion (e.g. c.795_796insCCTCTTGCT)
      - the WT or next WT nucleotide in case variant is deletion (e.g. c.723delC), and
    * 10bps up- and downstream of the variant position a sequence that matches the TP53 CDS reference
    """

    variant_type, variant_pos, pos_start, pos_end = calculate_parameters(variant, variant_pos_string, tp53_ref_sub) # Could to be calculated only onces, not for each line

    ref_wt = tp53_ref_sub[pos_start:pos_end]

    if variant_type=='sub':
        ref_mut = tp53_ref_sub[pos_start:variant_pos] + bases[1] + tp53_ref_sub[variant_pos+1:pos_end]
    elif variant_type=='del':
        ref_mut = tp53_ref_sub[pos_start:variant_pos] + tp53_ref_sub[variant_pos+1:pos_end+1]
    elif variant_type=='ins':
        pos_end_mut_ins = pos_end
        if pos_end + len(bases[1]) > len(tp53_ref_sub):
            pos_end_mut_ins = len(tp53_ref_sub) - len(bases[1])
        ref_mut = tp53_ref_sub[pos_start:variant_pos+1] + bases[1] + tp53_ref_sub[variant_pos+1:pos_end_mut_ins]

    # Will be used for subsitutions only when checking whether +/-10bp max 1 mismatch occurs
    seq_pre = tp53_ref_sub[pos_start:variant_pos]
    seq_post = tp53_ref_sub[variant_pos+1:pos_end]
    
    read_sub = read[pos_start:pos_end]
    
    if variant_type=='sub' or variant_type=='del':
        read_sub_mut = read_sub
    if variant_type=='ins':
        read_sub_mut = read[pos_start:pos_end+len(bases[1])]

    cnt = 0

    if read_sub == ref_wt or read_sub_mut == ref_mut:
        # Substring of read matches wt reference sequence, either exactly or with variant nucleotide
        return True
    
    if variant_type=='sub' and read[variant_pos] != bases[0] and read[variant_pos] != bases[1]:
        # Neither wt nor mut base at variant position, but other base (only for subsitutions)
        return False
    
    if variant_type=='sub':
        # Check whether substring differ in max 1bp (only for subsitutions)
        if variant_type=='sub':
            for i in range(len(seq_pre)):
                if read_sub[i] != seq_pre[i]:
                    cnt = cnt + 1
                if cnt > 1:
                    # Substring differs in more than max 1bp
                    return False
            for i in range(len(seq_post)):
                if read_sub[len(seq_pre)+1+i] != seq_post[i]:
                    cnt = cnt + 1
                if cnt > 1:
                    # Substring differs in more than max 1bp
                    return False            
        if cnt < 2:
            # Substring differ in max 1bp
            return True

    # An insertion or deletion that was not an exact match
    return False


def write_to_fastq_file(f_handle, title, seq, qual):
    """
    Write read to FASTQ file
    """
    
    f_handle.write('@' + title + '\n')
    f_handle.write(seq + '\n')
    f_handle.write('+' + '\n')
    f_handle.write(qual + '\n')


def parse_fastq_files_perform_qc(file1, split1, file2, split2, bc_occ, df_bc_buckets, local_stats, primer_used, variant, variant_pos_string, tp53_ref_sub, bases, k, folder_tmp):
    """
    Parse input data line by line and write into k files according to barcode buckets
    Perform QC steps 3-5: 
    3) Sequencing QC
    4) Primer QC
    5) Read QC
    Save statistics (how many reads, mRNA transcripts, and cells we kept)
    """
        
    bcs_after_seqqc = set()
    bcs_after_primerqc = set()
    bcs_after_readqc = set()

    transcripts_after_seqqc = dict()
    transcripts_after_primerqc = dict()
    transcripts_after_readqc = dict()

    num_transcripts_after_seqqc = 0
    num_transcripts_after_primerqc = 0
    num_transcripts_after_readqc = 0

    num_reads_after_seqqc = 0
    num_reads_after_primerqc = 0
    num_reads_after_readqc = 0
    
    fastq_R1_file_names = []
    fastq_R2_file_names = []
    for i in range(k):
        fastq_R1_file_names.append(folder_tmp + 'R1_after_QC_part'+str(i+1)+'.fastq')
        fastq_R2_file_names.append(folder_tmp + 'R2_after_QC_part'+str(i+1)+'.fastq')
    
    for file_name in fastq_R1_file_names:
        with open(file_name, 'w') as file_handle:
            file_handle.truncate(0)
    for file_name in fastq_R2_file_names:
        with open(file_name, 'w') as file_handle:
            file_handle.truncate(0)

    # Read data line by line and write to bucket (and subsample for bcs with excessive amount of reads)
    fastq_R1_file_handles = []
    fastq_R2_file_handles = []
    try:
        for file_name in fastq_R1_file_names:
            fastq_R1_file_handles.append(open(file_name, 'a'))
        for file_name in fastq_R2_file_names:
            fastq_R2_file_handles.append(open(file_name, 'a'))
        
        with gzip.open(file1, "rt") as handle1, gzip.open(file2, "rt") as handle2:
            for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):

                bc = str(seq1[:split1])
                umi = str(seq1[split1:])
                qual_bc_umi = mean([ord(letter)-33 for letter in qual1])
                primer = str(seq2[:split2])
                read = str(seq2[split2:])
                qual_primer_read = mean([ord(letter)-33 for letter in qual2])
            
                if (qual_bc_umi > 30) & (qual_primer_read > 30):

                    if bc in bc_occ:
                        
                        curr_bucket = int(df_bc_buckets.loc[bc]['bucket'])

                        # 6) Primer QC
                        if check_whether_primer_passes_primer_qc(primer, primer_used):

                            # 7) Read QC
                            if check_whether_read_passes_read_qc(read, variant, variant_pos_string, tp53_ref_sub, bases):

                                # Write filtered data into files
                                write_to_fastq_file(fastq_R1_file_handles[curr_bucket-1], title1, seq1, qual1)
                                write_to_fastq_file(fastq_R2_file_handles[curr_bucket-1], title2, seq2, qual2)
    
                                # Stats: Add unique bcs to bc -set; Add unique (bc, umi) -pairs to transcripts -set; Count reads
                                if bc not in bcs_after_readqc: bcs_after_readqc.add(bc)
                                transcripts_after_readqc = add_transcript_to_dict(transcripts_after_readqc, bc, umi)
                                num_reads_after_readqc = num_reads_after_readqc +1
                        
                            # Stats: Add unique bcs to bc -set; Add unique (bc, umi) -pairs to transcripts -set; Count reads
                            if bc not in bcs_after_primerqc: bcs_after_primerqc.add(bc)    
                            transcripts_after_primerqc  = add_transcript_to_dict(transcripts_after_primerqc, bc, umi)
                            num_reads_after_primerqc = num_reads_after_primerqc +1
                                                
                        # Stats: Add unique bcs to bc -set; Add unique (bc, umi) -pairs to transcripts -set; Count reads
                        if bc not in bcs_after_seqqc: bcs_after_seqqc.add(bc)
                        transcripts_after_seqqc = add_transcript_to_dict(transcripts_after_seqqc, bc, umi)
                        num_reads_after_seqqc = num_reads_after_seqqc +1


        for transcripts in transcripts_after_seqqc.values():
            num_transcripts_after_seqqc = num_transcripts_after_seqqc + len(transcripts)

        for transcripts in transcripts_after_primerqc.values():
            num_transcripts_after_primerqc = num_transcripts_after_primerqc + len(transcripts)

        for transcripts in transcripts_after_readqc.values():
            num_transcripts_after_readqc = num_transcripts_after_readqc + len(transcripts)

    finally:
        for file_handle in fastq_R1_file_handles:
            file_handle.close()
        for file_handle in fastq_R2_file_handles:
            file_handle.close()

    local_stats.bc.num_after_seqqc       = len(bcs_after_seqqc)
    local_stats.bc.num_after_primerqc    = len(bcs_after_primerqc)
    local_stats.bc.num_after_readqc      = len(bcs_after_readqc)

    local_stats.transcript.num_after_seqqc       = num_transcripts_after_seqqc
    local_stats.transcript.num_after_primerqc    = num_transcripts_after_primerqc
    local_stats.transcript.num_after_readqc      = num_transcripts_after_readqc

    local_stats.read.num_after_seqqc       = num_reads_after_seqqc
    local_stats.read.num_after_primerqc    = num_reads_after_primerqc
    local_stats.read.num_after_readqc      = num_reads_after_readqc

    return local_stats


# Main

def main():

    # Command line parameters

    parser = argparse.ArgumentParser(
        prog = 'python scAmp-2-QC.py',
        description = 'Part 2 of TP53 mutation status identification from scAmp-seq data: \
                    Perform quality control (QC)'
    )

    parser.add_argument('-amp', '--scampparam', type=str, required=True, help='scAmp parameter file')
    parser.add_argument('-w', '--whitelist', type=str, required=True, help='barcode whitelist file')
    parser.add_argument('-k', '--part', type=int, required=True, help='specify k, i.e. how many parts to split data into')
    parser.add_argument('-t', '--tmp', type=str, required=True, help='tmp (intermediate results) folder')
    
    args = parser.parse_args()

    file_param = args.scampparam
    whitelist = args.whitelist
    k = args.part
    folder_tmp = args.tmp
    

    # Create instance of class
    my_stats = ScAmp_stats()


    # Parse input parameters
    variant, variant_type, bases, primer_used, tp53_ref_sub, variant_pos_string, file_R1, file_R2 = parse_parameters(file_param)
    primer_len = len(primer_used)

    print('Input parameters:')
    print('primer seq: ' + primer_used)
    print('primer has length: ' + str(primer_len))
    print('R1 file: ' + file_R1)
    print('R2 file: ' + file_R2)
    sys.stdout.flush()


    # Perform QC steps 1-2: 1) Parse raw data and 2) perform Barcode QC
    print('Parsing barcodes and barcode QC...')
    sys.stdout.flush()
    
    bc_occurrences, my_stats = parse_bcs_perform_qc(file_R1, 16, my_stats, whitelist, folder_tmp)


    # In order to be able to split the data into k buckets, define which barcodes goes into which bucket
    print('Sort barcodes into buckets...')
    sys.stdout.flush()

    df_bc = sort_bcs_into_buckets(bc_occurrences, k)


    # Parse input data line by line and write into k files according to barcode buckets
    # Perform QC steps 3-5: 3) Sequencing QC, 4) Primer QC, and 5) Read QC
    # Count how many reads, mRNA transcripts, and cells we kept
    print('Parse input data and write into files according to buckets...')
    sys.stdout.flush()

    my_stats = parse_fastq_files_perform_qc(file_R1, 16, file_R2, primer_len, bc_occurrences, df_bc, my_stats, primer_used, variant, variant_pos_string, tp53_ref_sub, bases, k, folder_tmp)

    # Write statistics to file
    file_stats = folder_tmp + 'scAmp-2-QC_stats.csv'
    my_stats.write_to_file(file_stats)

    print('Finished')


if __name__ == "__main__":
    main()
