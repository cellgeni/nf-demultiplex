#!/bin/bash

set -euo pipefail

##Mandatory inputs
sample_id=$1 #sample id (just for naming output folder)
barcodes_file=$2 #filtered barcodes file generated from mapping data to a reference genome
bam_file=$3 #bam file generated from mapping data to a reference genome
known_genotypes=$4 #boolean variable saying whether to use known genotypes option or not
soc_vcf=$5 #the vcf file, either a common variants vcf downloaded using: https://github.com/wheaton5/souporcell/blob/master/README.md#souporcell or a known genotypes vcf 
soc_fasta=$6 #the reference fasta file for the organism the sample was aligned to
k_value=$7 #the number of donors to be demultiplexed
skip_remap=$8 #boolean value telling souporcell whether to skip remapping step or not
no_umi=$9 #boolean value telling souporcell whether bam file contains umis or not
ignore=False # souporcell checks first 10e5 reads and dies if less than 25% of them have barcodes not in barcode list. Usually it is ok, but one can change it to True if it fails (should it be default?).


echo "Number of barcodes: $( cat ${barcodes_file} | wc -l )"

common_or_known="--common_variants"
if "${known_genotypes}"; then
   common_or_known="--known_genotypes"
fi
common_or_known="${common_or_known} ${soc_vcf}"

mkdir "${sample_id}"

echo "souporcell_pipeline.py -i ${bam_file} -b ${barcodes_file} -f ${soc_fasta} -k ${k_value} ${common_or_known} -t 8 -o ${sample_id} --skip_remap ${skip_remap} --no_umi ${no_umi}" > "${sample_id}/cmd.txt"

souporcell_pipeline.py                  \
  -i "${bam_file}"                      \
  -b "${barcodes_file}"                 \
  -f "${soc_fasta}"                     \
  -k "${k_value}"                       \
  $common_or_known                      \
  --threads 16                          \
  -o "${sample_id}"                     \
  --skip_remap "${skip_remap}"          \
  --no_umi "${no_umi}"                  \
  --ignore ${ignore}
