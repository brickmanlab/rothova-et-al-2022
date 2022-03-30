#!/bin/bash

set -e

if [ "$#" -ne 3 ]; then
	echo usage : run_velocity.sh [path/to/fastq] [batch] [threads]
	exit
fi

FASTQ_PATH=$1
BATCH=$2
THREADS=$3
REF_DIR="/home/projects/dan_bri/data/DataBase/mm9-m1/star-index/"
OUTPUT=output

echo "[run_velocity.sh]> running velocity"

check if _temp exists
if [ ! -d "_temp" ]; then
    mkdir -p _temp
    mkdir -p _temp/converted
fi

echo "[run_velocity.sh]> splitting FastQ files"
fastp \
  -i $FASTQ_PATH/Undetermined_S0_R1_001.fastq.gz \
  -I $FASTQ_PATH/Undetermined_S0_R2_001.fastq.gz \
  -o _temp/Undetermined_S0_R1_001.fastq.gz \
  -O _temp/Undetermined_S0_R2_001.fastq.gz \
  --split_by_lines 4000000 \
  --thread $THREADS \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g

echo "[run_velocity.sh]> converting reads to 10X format"
python parse_parallel.py \
  --input _temp/ \
  --threads $THREADS \
  --output _temp/converted/

echo "[run_velocity.sh]> concatenate converted reads"
for f in _temp/converted/*R1*.fastq.gz; do cat $f >> Undetermined_S0_R1_001.fastq.gz; done
for f in _temp/converted/*R2*.fastq.gz; do cat $f >> Undetermined_S0_R2_001.fastq.gz; done

echo "[run_velocity.sh]> trimming poly-T and low quality reads"
cutadapt \
  --cores $THREADS \
  -m 20 \
  -A "T{68}" \
  --pair-filter=both \
  -o R1.trim.fastq.gz \
  -p R2.trim.fastq.gz \
  Undetermined_S0_R1_001.fastq.gz \
  Undetermined_S0_R2_001.fastq.gz

echo "[run_velocity.sh]> construct whitelist for StarSolo"
python build_whitelist.py \
  --batch $BATCH \
  --input $FASTQ_PATH

echo "[run_velocity.sh]> running StarSolo"
STAR \
  --genomeDir $REF_DIR \
  --readFilesIn R2.trim.fastq.gz R1.trim.fastq.gz \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 11 \
  --soloUMIstart 12 \
  --soloUMIlen 8 \
  --runThreadN 35 \
  --soloCBwhitelist ./whitelist.txt \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --readFilesCommand zcat

echo "[run_velocity.sh]> cleaning temporary files"
rm -rf Aligned.out.sam SJ.out.tab _STARtmp/ _temp/ *trim.fastq.gz

# move files to output directory
echo "[run_velocity.sh]> moving files into $OUTPUT/$BATCH"
mkdir -p $OUTPUT/$BATCH
mv fastp* $OUTPUT/$BATCH
mv Log* $OUTPUT/$BATCH
mv whitelist.txt $OUTPUT/$BATCH
mv Solo.out $OUTPUT/$BATCH

# genozip new raw files
genozip Undetermined_S0_R1_001.fastq.gz -@ $THREADS -o $OUTPUT/$BATCH/Undetermined_S0_R1_001.fastq.genozip
genozip Undetermined_S0_R2_001.fastq.gz -@ $THREADS -o $OUTPUT/$BATCH/Undetermined_S0_R2_001.fastq.genozip
rm Undetermined*.gz

echo "[run_velocity.sh]> Done"