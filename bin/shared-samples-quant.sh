#!/bin/bash

set -euo pipefail

IMAGE="/nfs/cellgeni/singularity/images/anndata_08.sif"
PYSCRIPT="/nfs/cellgeni/tickets/tic-2362/actions/shared-samples-quant.py"

SHARED_SAMPLES_DIR="/lustre/scratch127/cellgen/cellgeni/tickets/tic-2362/souporcell_results/simon_stuff"

cd $SHARED_SAMPLES_DIR
mkdir -p quantification/clusters
for shared in out*; do
 #Get sample id from shared sample file name
 sample=`echo $shared | cut -f 2 -d .`
 #Extract cluster comparison tsv from shared samples file
 tail -n +6 $shared > "quantification/clusters/${sample}.cluster"
done

cd quantification
mkdir -p loss-tables
for comparison in clusters/*.cluster; do
 sample=`echo $comparison | cut -f 2 -d / | cut -f 1 -d .`
 /software/singularity/v3.10.0/bin/singularity exec -B /lustre,/nfs $IMAGE python $PYSCRIPT --input $comparison --sample $sample
done

echo -e "experiment1_cluster\texperiment2_cluster\tloss\tsample_id" > total.tsv
cat loss-tables/* >> total.tsv
