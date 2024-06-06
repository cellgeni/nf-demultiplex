#!/bin/bash

set -euo pipefail

##Mandatory inputs
samplefile=$1 #the list of samples to have souporcell outputs identified for matching donors

# -n 100000 it seems to control only number of pairs to print. See no reasons why not to print all of them. or at least up to 100000
cut -f 1 "${samplefile}" | while read s1; do
   cut -f 1 "${samplefile}" | while read s2; do
       shared_samples.py -1 $s1 -2 $s2 -n 100000 1> "map.${s1}-${s2}" 2> "err.${s1}-${s2}"
   done
done
