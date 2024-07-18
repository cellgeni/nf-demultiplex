#!/bin/bash

set -euo pipefail

##Mandatory inputs
sample_id=$1 #sample id (just for naming output folder)
snp_folder=$2 #the folder contains the output from cellSNP
k_value=$3 #the number of donors to be demultiplexed

mkdir "${sample_id}-vireo"
echo "vireo -c ${snp_folder} -N ${k_value} -o ${sample_id}-vireo" > "${sample_id}-vireo/cmd.txt"

vireo -c ${snp_folder} -N ${k_value} -o ${sample_id}-vireo --nproc 16
