#!/bin/bash

set -euo pipefail

##Mandatory inputs
sample_id=$1 #sample id (just for naming output folder)
snp_folder=$2 #the folder contains the output from cellSNP
k_value=$3 #the number of donors to be demultiplexed
known_genotypes=$4 #boolean variable saying whether to use known genotypes option or not
snp_vcf=$5 #SNP vcf, cellSNP has advice here: https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html
genoTag=$6

common_or_known=""
if "${known_genotypes}"; then
   common_or_known="--donorFile ${snp_vcf} --genoTag ${genoTag}"
fi



mkdir "${sample_id}-vireo"

cmd="vireo -c ${snp_folder} -N ${k_value} -o ${sample_id}-vireo $common_or_known --nproc 16"
echo $cmd  > "${sample_id}-vireo/cmd.txt"

$cmd
