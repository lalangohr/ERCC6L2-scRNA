#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-6-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify which BCFtools to use, where the genome reference is stored, and the start and end position of TP53
bcftools=./software/bcftools
hg38ref=./hg38/Homo_sapiens_assembly38.fasta
tp53position=chr17:7668421-7687490

# Specify directory, where BAM and BAI files per mRNA transcript are stored, and VCF files will be stored
bc_dir=./results/E1a_c490AG/bc/

# Create VCF files per mRNA transcript (needed before haploid VCF files can be created)
for R2_file in $(cat ${bc_dir}fastq_to_bam.txt); do 
    bamfile=${R2_file}.Aligned.sortedByCoord.out.bam
    outfile=${R2_file}.Aligned.sortedByCoord.out.bam.vcf

    /usr/bin/time -v \
    srun $bcftools mpileup \
    --count-orphans \
    --max-depth 99999999 \
    --min-BQ 30 \
    --ignore-overlaps \
    --region $tp53position \
    --fasta-ref $hg38ref \
    --output-type v \
    --output $outfile $bamfile

done
