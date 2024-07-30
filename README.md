# nf-demultiplex
Our [Souporcell repo](https://github.com/cellgeni/souporcell) combined with Vireo into a Nextflow pipeline.

There are two branches:

`main` - this branch contains the script for running demultiplexing on the FARM using Nextflow command line.

`nextflow-tower` - this branch contains the script for running demultiplexing on the FARM using Nextflow Tower.

## Quick start
First `git clone https://github.com/cellgeni/nf-demultiplex.git` into actions. Then create sample file `actions/samples.tsv`:
<pre>
name1	/path/to/1/bam	/path/to/1/barcodes.tsv.gz  K1
name2	/path/to/2/bam	/path/to/2/barcodes.tsv.gz  K2
</pre>

Then run from ticket folder:
<pre>
screen -S tic-N
fash 1 oversubscribed 1
nextflow run actions/nf-demultiplex/main.nf \
 -entry souporcell \
 --SAMPLEFILE actions/samples.tsv \
 -resume
</pre>


If you know number of expected genotypes you can following script to group souporcell sample-slusters:
<pre>
Rscript actions/nf-demultiplex/bin/group_genotypes.R actions/samples.tsv demultiplex-results/souporcell number_of_genotypes
</pre>
Script will produce three output files:
* `shared_samples_loss_heatmap.pdf` heatmap of loss function and cluster grouping by sample (columns) and genotypes (rows)
* `shared_samples_loss.csv` matrix of pairwise losses
* `shared_samples_clusters.csv` cluster information, including cell counts and genotype group

## Contents of Repo:
* `main.nf` - the Nextflow pipeline that executes demultiplexing.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `examples/samples.txt` - samplefile tsv containing 3 fields: sampleID, path to BAM file, path to barcodes tsv.gz file (The order of these files is important!). These paths can be IRODs paths or local paths.
* `examples/RESUME-demultiplex` - an example run script that executes the pipeline it has 1 hardcoded argument: `/path/to/sample/file` that needs to be changed based on your local set up.
* `bin/cellsnp.sh` - a bash script that runs cellsnp inside the pipeline.
* `bin/shared-samples-quant.py` - a python script that produces a tsv of top clusters shared between samples.
* `bin/shared_samples.sh` - a bash script that runs souporcell's shared samples funcitonality inside the pipeline.
* `bin/soup_qc.sh` - a quick qc script that enables quick sanity checks that souporcell worked correctly.
* `bin/souporcell.sh` - a bash script that runs souporcell inside the pipeline.
* `bin/vireo.sh` - a nash script that runs vireo inside the pipeline.
* `dockerfiles/Dockerfile-souporcell` - a dockerfile to reproduce the environment used to run souporcell in the pipeline.
* `dockerfiles/Dockerfile-vireo` - a dockerfile to reproduce the environment used to run vireo in the pipeline.
* `dockerfiles/Dockerfile-shared-samples-quantification` - a dockerfile to reproduce the environment used to run shared samples quantificaiton.

## Pipeline Arguments:
* `-entry` - The entrypoint to specify which determines whether souporcell or vireo or both.
* `--SAMPLEFILE` - The path to the sample file provided to the pipeline. This is a tab-separated file with one sample per line. Each line should contain a sample id, path to bam file, path to gzipped barcodes file, K - number of genotypes in this sample (in that order!).
* `--outdir` - The path to where the results will be saved.
* `--K` - The number of donors in the samples. All samples must contain the same number of donors.
* `--barcodes_on_irods` - Tells pipeline whether to look for the gzipped barcodes file on IRODS or the FARM (default false means look locally).
* `--bam_on_irods` - Tells pipeline whether to look for the bam file on IRODS or the FARM (default false means look locally).
#### Souporcell:
* `--soc_vcf` - Path to vcf file used for souporcell (default is 2p 1k genome with chr nomenclature). This default argument is hardcoded and needs to be changed to your local path to the file. 
* `--soc_fasta` - Path to  genome fasta for pipeline to use (by default GRCh38 2020A is used). This default argument is hardcoded and needs to be changed to your local path to the file. 
* `--known_genotypes` - Whether to use the `known_genotypes` option. If true is used then the `--soc_vcf` option needs to be provided with a path to the known genotypes VCF file. The number of genotypes in the known genotypes vcf must match the number of donors supplied in the K parameter (Defualt false means not to use known_genotypes).
* `--skip_remap` - Whether to skip remapping in souporcell pipeline (default true means skip remapping).
* `--no_umi` - Tells the pipeline whether the BAM files have a UMI tag (default false means BAM file has UMI tag).
#### Vireo
* `--snp_vcf` - The gzipped VCF file to provide to cellSNP which is ran to generate input for vireo (default genome1K.phase3.SNP_AF5e2.chr1toX.hg38). This default argument is hardcoded and needs to be changed to your local path to the file. 
