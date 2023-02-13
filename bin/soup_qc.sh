#!/bin/bash

set -euo pipefail

sample=$1

##Get ambient RNA percentage
ambient=$(cat $sample/ambient_rna.txt | cut -d " " -f 5 | cut -c 1-4)
echo "${ambient}% ambient RNA" > "${sample}/qc.txt"
##Get number of doublets in sample
doublet=$(grep "removing .* as doublet" $sample/doublets.err | wc -l)
echo "${doublet} doublets removed" >> "${sample}/qc.txt"
##Get number of cells, doublets, singlets and unassigned for each cluster in sample
clusters=($(tail -n +2 $sample/clusters.tsv | cut -f 3 | sort | uniq))
printf "%-20s %-20s %-20s %-20s %-20s\n" "Cluster" "Total Cells" "Doublets" "Singlets" "Unassigned Cells" >> "${sample}/qc.txt"
for cluster in "${clusters[@]}" ; do
    count=$(grep $cluster $sample/clusters.tsv | wc -l);
    d_count=$(grep $cluster $sample/clusters.tsv | { grep doublet || true; } | wc -l);
    s_count=$(grep $cluster $sample/clusters.tsv | { grep singlet || true; } | wc -l);
    u_count=$(grep $cluster $sample/clusters.tsv | { grep unassigned || true; } | wc -l);
    printf "%-20s %-20s %-20s %-20s %-20s\n" "${cluster}" "${count}" "${d_count}" "${s_count}" "${u_count}" >> "${sample}/qc.txt"; 
done
