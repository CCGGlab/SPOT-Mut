#!/bin/bash

set -euo pipefail

# Variables
INPUT_DIR="" # raw fastq files directory
OUTPUT_DIR="" # output directory for aligned BAM files
REF_GENOME="hg38.fa" # indexed UCSC hg38 reference genome file
GERMLINE="somatic_hg38/af-only-gnomad.hg38.vcf.gz" # Common/known germline variants file from gnomad
## Download from:
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
##
TARGETS="data/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed" # WES target regions file
THREADS=32 ## CPU threads used during alignment and post-processing

mkdir -p $OUTPUT_DIR

# Loop through paired FASTQ files
# R1=$1
for R1 in $(find $INPUT_DIR -type f -name "*_good_1.fq.gz"); do
	# Get sample name
	R2="${R1:0:${#R1}-7}2${R1:${#R1}-6}"
	filename=$(basename $R1)
	Sample=${filename:0:10}

	echo "Processing sample: $Sample"
	OUTDIR="${OUTPUT_DIR}/$Sample"
	mkdir -p $OUTDIR
	# Output files
	ALIGNED="$OUTDIR/$Sample"	
	echo 
	# Alignment, including read group information definition
	printf "@RG\tID:${Sample}\tLB:1\tPL:Illumina_NovaSeq_X_plus\tSM:${Sample}\tPU:LH00330\n"
	bwa-mem2 mem -v 1 -t $THREADS -R "@RG\tID:${Sample}\tLB:1\tPL:Illumina_NovaSeq_X_plus\tSM:${Sample}\tPU:LH00330" "$REF_GENOME" "${R1}" "$R2" > "${ALIGNED}.sam"
			
	# Convert to sorted BAM
	echo "Converting to BAM and sorting."
	samtools view -@ $THREADS -Sb "${ALIGNED}.sam" > ${ALIGNED}_unsorted.bam
    samtools sort -@ $THREADS -o "${ALIGNED}.bam" ${ALIGNED}_unsorted.bam
	# remove sam
	rm ${ALIGNED}.sam ${ALIGNED}_unsorted.bam

	# Index BAM
	echo "Indexing"
	samtools index -@ $THREADS "${ALIGNED}.bam"
			
	gatk MarkDuplicates -I "${ALIGNED}.bam" -O "${ALIGNED}_nodup.bam" -M "${ALIGNED}_mkdup_metrics.txt" --TMP_DIR $OUTDIR

	## remove unmarked
	rm ${ALIGNED}.bam

	gatk BaseRecalibrator -I ${ALIGNED}_nodup.bam -R $REF_GENOME --known-sites $GERMLINE -O ${ALIGNED}_recal_data.table

	gatk ApplyBQSR -I ${ALIGNED}_nodup.bam -R $REF_GENOME --bqsr-recal-file ${ALIGNED}_recal_data.table -O ${ALIGNED}_nodup_recal.bam	

	rm ${ALIGNED}_nodup.bam

	echo "Finished sample: $Sample"
done