MEM=10000
GROUP="cellgeni"
QUE='normal'

# loop through sample_ids
for i in `cut -d' ' -f1 actions/samples.tsv`
do
 echo $i
 bsub -n 2 -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE \
  -o logs/ooo.flag.%J.txt -e logs/eee.flag.%J.txt \
  "conda activate pysam-env; python actions/nf-demultiplex/utils/flag_bam.py data/${i}.bam data/${i}_flagged.bam ${i}"
done
