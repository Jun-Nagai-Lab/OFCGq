#!/bin/bash

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioenv

# Set directories
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script directory: $SCRIPT_DIR"

RAW_DIR="/mnt/g/60-819312648/N2301043_60-819312648_Lib/230209-A00151B_L001"
OUT_DIR="/mnt/g/60-819312648/N2301043_60-819312648_Lib/trimmed"
SUFFIX=".fq.gz"

# Create output folder if it doesn't exist
mkdir -p "$OUT_DIR"

# List samples
cd "$RAW_DATA"
samples=$(ls *_L1_1$SUFFIX | sed "s/_L1_1${SUFFIX}//g")

echo "Samples to process:"
echo "$samples"

# Run Trim Galore
for sample in $samples; do
    echo "Processing: $sample"
    trim_galore --paired --illumina \
        -o "$OUT_DIR" \
        "${sample}_L1_1${SUFFIX}" "${sample}_L1_2${SUFFIX}"
done

# Rename output files for clarity
cd "$OUT_DIR"
for f in *_val_1.fq.gz; do
    base=$(basename "$f" _val_1.fq.gz)
    mv "$f" "${base}_R1.cor_trim.fq.gz"
done

for f in *_val_2.fq.gz; do
    base=$(basename "$f" _val_2.fq.gz)
    mv "$f" "${base}_R2.cor_trim.fq.gz"
done

echo "✅ Trimming complete. Results saved in: $OUT_DIR"
