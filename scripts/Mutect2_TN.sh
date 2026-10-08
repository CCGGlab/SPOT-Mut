#!/bin/bash

## Usage: ./Mutect2_TN.sh <Tumor> <Normal>
##
## or (if you have GNU parallel installed)
##
## Usage: parallel --link ./Mutect2_TN.sh ::: <Tumor1> <Tumor2> ... ::: <Normal1> <Normal2> ...
##
## or (if you have a file with each line being a sample name)
##
## Usage: parallel --link -a Tumor_file.txt -a Normal_file.txt ./Mutect2_TN.sh


## <Tumor> points to the sample name of the tumor, assumed to be stored in the BAMFOLDERS variable; This folder should contain subfolders named after <Tumor>, within this folder, there should be a BAM file named <Tumor>_nodup_recal.bam (and corresponding index)
## <Normal>, same as <Tumor> but for the normal reference sample.
## When using GNU parallel, I reccomend using --dry-run to test the tumor/normal pairs

set -euo pipefail

REF_GENOME="hg38.fa" # indexed UCSC hg38 reference genome file
GERMLINE="somatic_hg38/af-only-gnomad.hg38.vcf.gz"
## Download from:
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
##
PON="somatic_hg38/1000g_pon.hg38.vcf.gz"
## Download from:
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
## https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
##
TARGETS="data/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed" # WES target regions file
BAMFOLDERS=""

# Run Mutect2
if [ $# -eq 0 ]; then
    echo "No inputs given"
    exit
fi
## copy BAM file to sample folder
BAMF=${BAMFOLDERS}${1}/

echo $1
if [ $# -eq 1 ]; then
    ## No matched normal
    gatk Mutect2 -I ${BAMF}${1}_nodup_recal.bam --germline-resource $GERMLINE --panel-of-normals $PON -O ${BAMF}${1}_unfiltered_no_normal.vcf -R $REF_GENOME --f1r2-tar-gz ${BAMF}f1r2.tar.gz --callable-depth 20
else
    BAM_N=${BAMFOLDERS}${2}/
    gatk Mutect2 -I ${BAMF}${1}_nodup_recal.bam -I ${BAM_N}${2}_nodup_recal.bam -normal $2 --germline-resource $GERMLINE --panel-of-normals $PON -O ${BAMF}${1}_unfiltered_normal.vcf -R $REF_GENOME --f1r2-tar-gz ${BAMF}f1r2.tar.gz --callable-depth 20
fi
# ## Prepare Filter Mutect Calls
gatk LearnReadOrientationModel -I ${BAMF}f1r2.tar.gz -O ${BAMF}read-orientation-model.tar.gz
gatk GetPileupSummaries -I ${BAMF}${1}_nodup_recal.bam -V $GERMLINE -L $TARGETS -O ${BAMF}pileups.table
rm ${BAMF}${1}_nodup_recal.ba*
rm ${BAMFOLDERS}${2}/${2}_nodup_recal.ba*
gatk CalculateContamination -I ${BAMF}pileups.table -O ${BAMF}contamination.table

## Filter
gatk FilterMutectCalls -V ${BAMF}${1}_unfiltered_normal.vcf -O ${BAMF}${1}_filtered.vcf --ob-priors ${BAMF}read-orientation-model.tar.gz --contamination-table ${BAMF}contamination.table -R $REF_GENOME --filtering-stats ${BAMF}filtering_stats.txt
grep -P '(\tPASS\t|##|#CHROM)' ${BAMF}${1}_filtered.vcf > ${BAMF}${1}_filtered_removed.vcf
		
## Annotate
## !! preconfigure Funcotator data sources before running 
## https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial
## Here, data source v1.7.20200521s was used (https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s)
gatk Funcotator --variant ${BAMF}${1}_filtered_removed.vcf --reference $REF_GENOME --ref-version hg38 --data-sources-path path/to/funcotator_dataSources.v1.7.20200521s/ --output ${BAMF}${1}_filtered_removed.maf --output-file-format MAF

## cleanup
rm ${BAMF}${1}_unfiltered_normal.vcf

gzip -f ${BAMF}${1}_filtered_removed.vcf
gzip -f ${BAMF}${1}_filtered_removed.maf
gzip -f ${BAMF}${1}_filtered.vcf
gzip -f ${BAMF}${1}_filtered.vcf.idx
gzip ${BAMF}${1}_unfiltered_normal.vcf.idx

echo "All samples processed."
