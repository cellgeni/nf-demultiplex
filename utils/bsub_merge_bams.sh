MEM=10000
GROUP="cellgeni"
QUE='normal'


ls -1 *flagged.bam > bam_all.list
bsub -n 16 -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE \
 -o logs/ooo.merge.%J.txt -e logs/eee.merge.%J.txt \
 "samtools merge -b bam_all.list  --write-index --threads 16 -o merged_all.bam; \
 samtools index merged_all.bam -@ 16"
