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
if [[ "${barcodes_file}" == *.gz ]]; then
  barcode_count=$(gzip -cd -- "${barcodes_file}" | wc -l)
else
  barcode_count=$(wc -l < "${barcodes_file}")
fi
echo "Number of barcodes: ${barcode_count}"

common_or_known=(--common_variants "${soc_vcf}")
if "${known_genotypes}"; then
   common_or_known=(--known_genotypes "${soc_vcf}")
fi

mkdir "${sample_id}"

souporcell_args=(
  souporcell_pipeline.py
  -i "${bam_file}"
  -b "${barcodes_file}"
  -f "${soc_fasta}"
  -k "${k_value}"
  "${common_or_known[@]}"
  --threads 16
  -o "${sample_id}"
)

# Souporcell 2.5 parses --skip_remap with argparse type=bool, for which even
# the string "false" is truthy. Omit the option to request its false default.
if "${skip_remap}"; then
  souporcell_args+=(--skip_remap True)
fi

# --no_umi is parsed as a string and recognises only the exact value "True".
if "${no_umi}"; then
  souporcell_args+=(--no_umi True)
else
  souporcell_args+=(--no_umi False)
fi

printf '%q ' "${souporcell_args[@]}" > "${sample_id}/cmd.txt"
printf '\n' >> "${sample_id}/cmd.txt"

"${souporcell_args[@]}"
