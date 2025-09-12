"""
Part 3 of TP53 mutation status identification from single-cell amplicon sequencing (scAmp-seq) data.

Generate separate FASTQ files for each mRNA transcript
"""


import argparse

from Bio.SeqIO.QualityIO import FastqGeneralIterator


# Functions

def write_to_fastq_file(f_handle, title, seq, qual):
    """
    Write read to FASTQ file
    """

    f_handle.write('@' + title + '\n')
    f_handle.write(seq + '\n')
    f_handle.write('+' + '\n')
    f_handle.write(qual + '\n')


def write_fastq_for_each_transcript(folder_tmp, part, folder_bc):
    """
    Generate separate FASTQ files for each mRNA transcript
    """

    bc_occ = dict()

    file_all_bcs_in_k_R1 = folder_tmp + 'R1_after_QC_'+part+'.fastq'
    file_all_bcs_in_k_R2 = folder_tmp + 'R2_after_QC_'+part+'.fastq'


    with open(file_all_bcs_in_k_R1, "rt") as handle1, open(file_all_bcs_in_k_R2, "rt") as handle2:
        for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):

            split1 = 16
            bc = str(seq1[:split1])
            umi = str(seq1[split1:])

            file_name_bc_R1 = folder_bc + bc + '_' + umi + '_R1.fastq'
            file_name_bc_R2 = folder_bc + bc + '_' + umi + '_R2.fastq'

            how = 'w' # if first time, use write and not append
            if bc in bc_occ:
                if umi in bc_occ[bc]:
                    how = 'a' # if not first time, use append
                else:
                    bc_occ[bc] = bc_occ[bc] + [umi]
            else:
                bc_occ[bc] = [umi]

            file_handle_bc_R1 = open(file_name_bc_R1, how)
            file_handle_bc_R2 = open(file_name_bc_R2, how)

            write_to_fastq_file(file_handle_bc_R1, title1, seq1, qual1)
            write_to_fastq_file(file_handle_bc_R2, title2, seq2, qual2)


# Main

def main():

    # Command line parameters

    parser = argparse.ArgumentParser(
        prog = 'python scAmp-3-FASTQ-per-mRNA-transcript.py',
        description = 'Part 3 of TP53 mutation status identification from scAmp-seq data: \
                       Generate separate FASTQ files for each mRNA transcript'
    )
    parser.add_argument('-t', '--tmp', type=str, required=True, help='tmp (intermediate results) folder')
    parser.add_argument('-k', '--part', type=int, required=True, help='specify k, i.e. part of intermediate results to further process')
    parser.add_argument('-b', '--bcs', type=str, required=True, help='bc folder (for fastq and bam files per transcript)')
    args = parser.parse_args()

    folder_tmp = args.tmp
    part = 'part' + str(args.part)
    folder_bc = args.bcs


    write_fastq_for_each_transcript(folder_tmp, part, folder_bc)

    print('Finished')


if __name__ == "__main__":
    main()
