#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-7-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify which BCFtools to use, and the start and end position of TP53
bcftools=./software/bcftools
tp53position=chr17:7668421-7687490

# Specify directory, where VCF files per mRNA transcript are stored
bc_dir=./results/E1a_c490AG/bc/

# Create haploid VCF files per mRNA transcript
for R2_file in $(cat ${bc_dir}fastq_to_bam.txt); do 
    bcffile=${R2_file}.Aligned.sortedByCoord.out.bam.vcf
    outfile=${R2_file}.Aligned.sortedByCoord.out.bam.haploid.vcf

    /usr/bin/time -v \
    srun cat $bcffile | $bcftools call --consensus-caller --ploidy 1 > $outfile

done
