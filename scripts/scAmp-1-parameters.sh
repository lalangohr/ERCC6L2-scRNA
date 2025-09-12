#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-1-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify sample, variant and read length
sample='E1a'
variant='c.818G>A'
read_len=270

# Create parameter file for sample and variant
/usr/bin/time -v \
srun python src/scAmp-1-parameters.py \
-s $sample \
-v $variant \
-l $read_len
