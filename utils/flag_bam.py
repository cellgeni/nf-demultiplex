import pysam
import sys

# this script was originally written by kp9

bampath = sys.argv[1]
bamout = sys.argv[2]
sample = sys.argv[3]

#open up the original bam, and the new one which will have CB flagged
bamfile = pysam.AlignmentFile(bampath, "rb")
tweaked = pysam.AlignmentFile(bamout, "wb", template=bamfile)
for read in bamfile.fetch():
    #do we have a CB tag? if so, stick the sample on as a prefix. write out to flagged
    if read.has_tag('CB'): 
        read.set_tag('CB',sample+'-'+read.get_tag('CB'))
    tweaked.write(read)

tweaked.close()
bamfile.close()
