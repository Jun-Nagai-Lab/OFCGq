#!/bin/bash

# === Activate environment (adjust as needed) ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioenv

# === Setup ===
DIR=~/bio_project
READ_DIR="/mnt/g/60-819312648/N2301043_60-819312648_Lib/trimmed"
OUT_STAR="/mnt/g/60-819312648/N2301043_60-819312648_Lib/STAR_mapped"
TMP_DIR_BASE="$HOME/bio_project/tmp_star"
OUT_COUNTS="/mnt/g/60-819312648/N2301043_60-819312648_Lib/featureCounts"
GENOME="/mnt/g/60-819312648/N2301043_60-819312648_Lib/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"   # Path to genome
GTF="/mnt/g/60-819312648/N2301043_60-819312648_Lib/genome/Mus_musculus.GRCm39.112.gtf"
SUFFIX=".fq.gz"

mkdir -p "$OUT_STAR"
mkdir -p "$OUT_COUNTS"

# === Detect Samples ===
samples=""
for f in "$READ_DIR"/*_L1_1_R1.cor_trim$SUFFIX; do
    [ -e "$f" ] || continue
    sample=$(basename "$f" _L1_1_R1.cor_trim$SUFFIX)
    samples+="$sample "
done

echo "Detected samples: $samples"

mkdir -p ~/bio_project/genome_index

STAR --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir ~/bio_project/genome_index \
  --genomeFastaFiles $GENOME \
  --sjdbGTFfile $GTF 

GENOME_INDEX="$HOME/bio_project/genome_index"

mkdir -p "$TMP_DIR_BASE"

# === STAR Mapping ===
for sample in $samples; do
    echo "Mapping: $sample"
	rm -rf "$TMP_DIR_BASE/${sample}_tmp" 
    STAR --runMode alignReads \
        --genomeDir "$GENOME_INDEX" \
        --readFilesIn "$READ_DIR/${sample}_L1_1_R1.cor_trim$SUFFIX" "$READ_DIR/${sample}_L1_2_R2.cor_trim$SUFFIX" \
        --readFilesCommand zcat \
        --runThreadN 2 \
        --quantMode GeneCounts \
        --outFileNamePrefix "$OUT_STAR/${sample}_" \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped Fastx \
	--outTmpDir "$TMP_DIR_BASE/${sample}_tmp"
done

# === Collect BAM Files for FeatureCounts ===
bam_files=()
for sample in $samples; do
    bam_files+=("$OUT_STAR/${sample}_Aligned.out.bam")
done

# === Run FeatureCounts ===
echo "Running featureCounts..."
featureCounts -p --countReadPairs -t exon -g gene_id \
    -a "$GTF" \
    -o "$OUT_COUNTS/featurecounts_output.txt" \
    -T 2 \
    "${bam_files[@]}"

rm -rf ~/bio_project/genome_index
rm -rf "$TMP_DIR_BASE"

echo "✅ All done. Counts saved to: $OUT_COUNTS/featurecounts_output.txt"
