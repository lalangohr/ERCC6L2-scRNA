#!/bin/bash
#SBATCH --time 7-00:00:00
#SBATCH --mem 20G
#SBATCH --partition cpu
#SBATCH --array=1-10
#SBATCH --output ./log/scAmp-3-FASTQ-%A-%a.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Set array numbers in the SBATCH parameter to the number of parts data was splitted into

# Specify directory, where FASTQ files after QC are stored, and directory, where FASTQ files per mRNA transcript will be stored
tmp_dir=./results/E1a_c490AG/tmp/
bc_dir=./results/E1a_c490AG/bc/

mkdir $bc_dir

# Generate separate FASTQ files for each mRNA transcript
/usr/bin/time -v \
srun python src/scAmp-3-FASTQ-per-mRNA-transcript.py \
-k $SLURM_ARRAY_TASK_ID \
-t $tmp_dir \
-b $bc_dir
