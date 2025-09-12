#!/bin/bash
#SBATCH --time 7-00:00:00
#SBATCH --mem 20G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-8-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify sample and variant file, directory, where final results will be stored, and directory, where VCF files are stored
param_file=./results/scAmp-parameters-scAmp3-E1a-c490AG-DATE-TIME.csv
results_dir=./results/E1a_c490AG/
bc_dir=./results/E1a_c490AG/bc/

# Specify input VCF file and output file suffices
input_suffix=.haploid.vcf
output_suffix=GL

# Define cells (barcodes) as WT or MUT
/usr/bin/time -v \
srun python src/scAmp-8-define-cells-as-WT-MUT.py \
-amp $param_file \
-r $results_dir \
-b $bc_dir \
-i $input_suffix \
-s $output_suffix
