#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-4-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify which STAR to use and where the genome indices are stored
star=./software/STAR
genome_indices_dir=./STAR-genome-indices/refdata-gex-GRCh38-2020-A_star-20211012-200507/

# Specify directory, where FASTQ files per mRNA transcript are stored and BAM files will be stored
bc_dir=./results/E1a_c490AG/bc/

# Remove FASTQ files with only 1 read (i.e. 4 lines in FASTQ file)
find ${bc_dir} -name "*R1.fastq" -type f > ${bc_dir}fastq_files.txt
find ${bc_dir} -name "*R2.fastq" -type f >> ${bc_dir}fastq_files.txt
for file in $(cat ${bc_dir}fastq_files.txt); do n=$(cat $file | wc -l); if (($n == 4)); then rm $file; fi; done

# Create list of remaining FASTQ files (i.e. with at least 2 reads)
find ${bc_dir} -name "*R2.fastq" -type f > ${bc_dir}fastq_to_bam.txt

# Create a BAM file for each remaining FASTQ files (i.e. one BAM file for each barcode-UMI-pair)
for R2_file in $(cat ${bc_dir}fastq_to_bam.txt); do 
    output_prefix=${R2_file}.

    /usr/bin/time -v \
    srun $star \
    --genomeDir $genome_indices_dir \
    --readFilesIn $R2_file \
    --readNameSeparator "space" \
    --outSAMmultNmax "1" \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $output_prefix \
    
done
