#!/bin/bash
#SBATCH --time 1-00:00:00
#SBATCH --mem 10G
#SBATCH --partition cpu
#SBATCH --output ./log/scAmp-2-%j.out

echo "node: $HOSTNAME"
echo "time: $(date)"


# Specify where 10X barcode whitelist is stored
wl_file=./10X-barcode-whitelist/3M-february-2018.txt.gz

# Specify sample and variant file, create directory, where final results will be stored, and directory, where FASTQ files after QC will be stored
param_file=./results/scAmp-parameters-scAmp3-E1a-c490AG-DATE-TIME.csv
results_dir=./results/E1a_c490AG/
tmp_dir=./results/E1a_c490AG/tmp/

mkdir $results_dir
mkdir $tmp_dir

# Specify k, the number of parts data should be splitted into
k=10

# Run quality control (QC)
/usr/bin/time -v \
srun python src/scAmp-2-QC.py \
-k $k \
-amp $param_file \
-w $wl_file \
-t $tmp_dir
