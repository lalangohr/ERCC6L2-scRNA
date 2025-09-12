#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-5-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify directory, where BAM files per mRNA transcript are stored and BAI files will be stored
bc_dir=./results/E1a_c490AG/bc/

# Create a BAI file for each BAM file
for R2_file in $(cat ${bc_dir}fastq_to_bam.txt); do 
    bamfile=${R2_file}.Aligned.sortedByCoord.out.bam

    /usr/bin/time -v \
    srun samtools index $bamfile

done
