#!/bin/bash

set -euo pipefail

##Mandatory inputs
samplefile=$1 #the list of samples to have souporcell outputs identified for matching donors

# -n 100000 it seems to control only number of pairs to print. See no reasons why not to print all of them. or at least up to 100000
# use tr-cut to make it work for both tab and space separated files
cat ${samplefile} | tr "\t" " " | cut -f1 -d ' ' | while read s1; do
   cat ${samplefile} | tr "\t" " " | cut -f1 -d ' ' | while read s2; do
       shared_samples.py -1 $s1 -2 $s2 -n 100000 1> "map.${s1}-${s2}" 2> "err.${s1}-${s2}"
   done
done
