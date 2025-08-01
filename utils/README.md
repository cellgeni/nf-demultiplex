Here are scripts that can  be used to merge multiple bams together to ran demultiplication on on multiple sumples simultaneoulsy.
The procedure is following:

1. Create `actions/samples.tsv` with two columns: sample_id and irods paths
2. Download bam files into `data` folder
```
cd data
while read n p
do
 echo $n
 iget ${p}/possorted_genome_bam.bam ${n}.bam
 iget ${p}/possorted_genome_bam.bam.bai ${n}.bai
 iget ${p}/filtered_feature_bc_matrix/barcodes.tsv.gz ${n}_barcodes.tsv.gz
done < ../actions/samples.tsv
```
3. Add sample_id to cell barcodes in bam file:
```
cd ..
actions/nf-demultiplex/utils/bsub_flag_bam.sh
```
4. merge bsms
```
cd data
../actions/nf-demultiplex/utils/bsub_merge_bams.sh
```

5.  make barcode file
```
for s in `cut -d ' ' -f1 ../actions/samples.tsv`
do 
 zcat ${s}_barcodes.tsv.gz | sed s/^/${s}-/
done > merged_barcodes_all.tsv
 
gzip merged_barcodes_all.tsv
```
6. make sample file, set N to number of expected genomes
```
N=...
p=`pwd`
echo "all ${p}/data/merged_all.bam ${p}/data/merged_barcodes_all.tsv.gz ${N}" > actions/samples_merged.tsv
```
7. start pipeline
```
nextflow run actions/nf-demultiplex/main.nf \
 -entry all \
 --SAMPLEFILE actions/samples_merged.tsv \
 --outdir demultiplex-results_merged \
 --ngenomes 2 \
 -resume
```