#!/bin/bash

set -euo pipefail

##Mandatory inputs
samplefile=$1 #the list of samples to have souporcell outputs identified for matching donors
k_value=$2 #the number of donors to be demultiplexed

cut -f 1 "${samplefile}" | while read s1; do
   cut -f 1 "${samplefile}" | while read s2; do
       shared_samples.py -1 $s1 -2 $s2 -n ${k_value} 1> "map${k_value}.${s1}-${s2}" 2> "err${k_value}.${s1}-${s2}"
   done
done
