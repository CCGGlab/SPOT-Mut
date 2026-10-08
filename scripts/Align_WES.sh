#!/bin/bash

# Usage: ./Align_WES.sh sample_number

set -euo pipefail

# Variables
OUTPUT_DIR="" # output directory for aligned BAM files
REF_GENOME="hg38.fa" # indexed UCSC hg38 reference genome file
GERMLINE="somatic_hg38/af-only-gnomad.hg38.vcf.gz" # Common/known germline variants file from gnomad
## Download GERMLINE from:
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
##
TARGETS="data/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed" # WES target regions file
THREADS=24 # CPU threads used during alignment and post-processing

ZIP_DIR="raw/" # directory where the zip files containing fastq files are located
OUT_DIR="${1}" # directory where the zip files are extracted

mkdir -p $OUT_DIR

for zip_arch in "$ZIP_DIR"/*.zip; do
    zip_base=$(basename "$zip_arch" .zip)
    target_dir="$OUT_DIR/$zip_base/"
    
    # Read all matching files into an array
    readarray -t matching_files < <(zipinfo -1 "$zip_arch" | grep -E ".*/${1}_R[12]\.fastq\.gz$")

    # If we have any matches, extract them
    if (( ${#matching_files[@]} > 0 )); then
        mkdir -p "$target_dir"
        echo "Extracting from: $zip_arch"
        for file in "${matching_files[@]}"; do
            echo "  -> $file  ==>  $target_dir"
            unzip -j "$zip_arch" "$file" -d "$target_dir"
        done
    fi
done


# # Loop through paired FASTQ files
for SD in "$OUT_DIR"/*/; do
	echo $SD
	for R1 in "${SD}/*_R1.fastq.gz"; do

		# Get sample name
		R2="$SD/${1}_R2.fastq.gz"
		RUN=$(basename $SD)
		echo "Processing sample: $1 in sequencing run $RUN"

		mkdir -p "${OUTPUT_DIR}/${1}"

		# Output files
		ALIGNED="$OUTPUT_DIR/${1}/${RUN}_${1}"

		## unzip
		TMP_R1="${OUTPUT_DIR}/${1}/${1}_R1.fastq"
		TMP_R2="${OUTPUT_DIR}/${1}/${1}_R2.fastq"
		TMP_R1_L="${OUTPUT_DIR}/${1}/${1}_R1_L"
		TMP_R2_L="${OUTPUT_DIR}/${1}/${1}_R2_L"
		#if [[ ! -f "${TMP_R1}" || ! -f "${TMP_R2}" ]]; then
		echo "Decompressing FASTQ files..."
		pigz -p $THREADS -d -c $R1 > "${TMP_R1}"
		pigz -p $THREADS -d -c $R2 > "${TMP_R2}"
		#fi

		## remove gzipped fastq files
		rm -r $SD

		echo "splitting fastq based on sequencing lane reported in fastq header :x:"	
		# Split fastq files according to sequencing lane
		# Assume both R1 and R2 reads have the same number of sequencing lanes included in the SAM file
		# Assume the lane is readable from the fastq record header as a single digit between two colons: ':lane_digit:'
		# Assume there aren't more than 10 lanes with numbers within [0-9]
			awk -v prefix="$TMP_R1_L" '{
				header = $0;
				getline seq; getline plus; getline qual;
			if (match(header, /:[0-9]:/)) {
					digit = substr(header, RSTART + 1, 1);  # Extract the digit between colons
					filename = prefix digit ".fastq";
					print header "\n" seq "\n" plus "\n" qual >> filename;
			} else {
					print "Warning: No digit found in header: " header > "/dev/stderr"
			}
			}' "$TMP_R1"
			
			awk -v prefix="$TMP_R2_L" '{
			            header = $0;
			            getline seq; getline plus; getline qual;
			    if (match(header, /:[0-9]:/)) {
			            digit = substr(header, RSTART + 1, 1);  # Extract the digit between colons
			            filename = prefix digit ".fastq";
			            print header "\n" seq "\n" plus "\n" qual >> filename;
			    } else {
			            print "Warning: No digit found in header: " header > "/dev/stderr"
			    }
				}' "$TMP_R2"

		# Assume the library barcode is the last element of the fastq header line. An element is split by ':'
		# The library barcode can differ 1 base from the reference. Thus, there are many possible strings that refer to the same barcode.
		# We count all the library barcodes in the file and pick the most occurring one as the true library barcode.
		# Extract barcode from fastq files and compare if they are the same
		# Count barcodes:
		echo "counting barcodes"
		barcode_counts_1=$(head -n 4000000 "$TMP_R1" | awk -F'[: ]' 'NR % 4 == 1 { gsub(/\+/, "", $NF); count[$NF]++ } 
			END {
				for (b in count) print count[b], b
			}' | sort -nr)

		# Log the barcodes
		echo $barcode_counts_1 > "${OUT_DIR}/barcode_counts_${RUN}_1.txt"

		# Get the most common barcode
		barcode1=$(echo "$barcode_counts_1" | head -n 1 | awk '{print $2}')
		TMP_R2=$OUT_DIR/27_R2.fastq
		echo "counting barcodes R2"
		barcode_counts_2=$(head -n 4000000 "$TMP_R2" | awk -F'[: ]' 'NR % 4 == 1 { gsub(/\+/, "", $NF); count[$NF]++ } 
			END {
				for (b in count) print count[b], b
			}' | sort -nr)

		# Log the barcodes
		echo $barcode_counts_2 > "${OUTPUT_DIR}/${1}/barcode_counts_${RUN}_2.txt"
		
		# Get the most common barcode
		barcode2=$(echo "$barcode_counts_2" | head -n 1 | awk '{print $2}')

		# compare barcodes
		if [[ "$barcode1" != "$barcode2" ]]; then
			echo "WARNING: Barcodes do not match.
			$barcode1
			$TMP_R1
			$barcode2
			$TMP_R2" > ${OUTPUT_DIR}/${1}/barcode_error.txt
		fi
				
		## Assume everything is sequenced on the same flowcell unit
		UNIT_R1=$(head -n 1 "$TMP_R1" | cut -d':' -f1 | sed 's/^@//')
		UNIT_R2=$(head -n 1 "$TMP_R2" | cut -d':' -f1 | sed 's/^@//')

		## remove unsplit fastq files
		rm $TMP_R1 $TMP_R2

		for file in ${TMP_R1_L}[0-9].fastq; do
			LANE="${file#$TMP_R1_L}"
			LANE="${LANE%.fastq}"
			echo "Aligning Lane $LANE"
			# Defining read group information and Alignment
			printf "@RG\tID:${1}.${RUN}.${LANE}\tLB:${barcode1}\tPL:AVITI_TRINITY\tSM:${1}\tPU:${UNIT_R1}.${barcode1}.${LANE}\n"
			bwa-mem2 mem -v 1 -t $THREADS -R "@RG\tID:${1}.${RUN}.${LANE}\tLB:${barcode1}\tPL:AVITI_TRINITY\tSM:${1}\tPU:${UNIT_R1}.${barcode1}.${LANE}" "$REF_GENOME" "${file}" "${TMP_R2_L}${LANE}.fastq" > "${ALIGNED}_L${LANE}.sam"
			
			# remove fastq
			rm "${TMP_R2_L}${LANE}.fastq" $file
			# Convert to sorted BAM
		    echo "Converting to BAM and sorting."
			samtools view -@ $THREADS -Sb "${ALIGNED}_L${LANE}.sam" | samtools sort -@ $THREADS -o "${ALIGNED}_L${LANE}.bam"
			# remove sam
			rm ${ALIGNED}_L${LANE}.sam

		done

		# Combine BAM files
		echo "Merging"
		samtools merge -@ $THREADS -o "${ALIGNED}.bam" ${ALIGNED}_L[0-9].bam
		# samtools merge -@ $THREADS -o "${ALIGNED}.bam" ${ALIGNED}_L[0-9]_RG.bam

		## remove unmerged bams
		rm ${ALIGNED}_L[0-9].bam
		# rm ${ALIGNED}_L[0-9]_RG.bam

		# Index BAM
		echo "Indexing"
		samtools index -@ $THREADS "${ALIGNED}.bam"

    done
done
	
SNAME=$(basename $OUT_DIR)
OUTFOLDER="${OUTPUT_DIR}/$SNAME/"

# Combine BAM files
echo "Merging all"
	
samtools merge -@ $THREADS -o "${OUTFOLDER}${SNAME}.bam" $OUTFOLDER*_${SNAME}.bam

## remove unmerged
rm $OUTFOLDER*_${SNAME}.bam*

gatk MarkDuplicates -I "${OUTFOLDER}${SNAME}.bam" -O "${OUTFOLDER}${SNAME}_nodup.bam" -M "${OUTFOLDER}${SNAME}_mkdup_metrics.txt" --TMP_DIR $OUTFOLDER

## remove unmarked
rm ${OUTFOLDER}${SNAME}.bam

gatk BaseRecalibrator -I ${OUTFOLDER}${SNAME}_nodup.bam -R $REF_GENOME --known-sites $GERMLINE -O ${OUTFOLDER}${SNAME}_recal_data.table

gatk ApplyBQSR -I ${OUTFOLDER}${SNAME}_nodup.bam -R $REF_GENOME --bqsr-recal-file ${OUTFOLDER}${SNAME}_recal_data.table -O ${OUTFOLDER}${SNAME}_nodup_recal.bam	

rm ${OUTFOLDER}${SNAME}_nodup.bam

# Cleanup
# rm -r ${OUTFOLDER}RUN* ${OUTFOLDER}${SNAME}_nodup.bam ${OUTFOLDER}${SNAME}.bam
echo "Finished sample: $SNAME"