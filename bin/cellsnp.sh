#!/bin/bash

set -euo pipefail

##Mandatory inputs
sample_id=$1 #sample id (just for naming output folder)
barcodes_file=$2 #filtered barcodes file generated from mapping data to a reference genome
bam_file=$3 #bam file generated from mapping data to a reference genome
snp_vcf=$4 #SNP vcf, cellSNP has advice here: https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html

echo "Number of barcodes: $( cat ${barcodes_file} | wc -l )"

mkdir "${sample_id}-cellsnp"
echo "cellsnp-lite -s ${bam_file} -b ${barcodes_file} -O ${sample_id}-cellsnp -R ${snp_vcf} --gzip" > "${sample_id}-cellsnp/cmd.txt"

cellsnp-lite -s ${bam_file} -b ${barcodes_file} -O ${sample_id}-cellsnp -R ${snp_vcf} --gzip --nproc 16
