# nf-souporcell
Our Souporcell repo but for Nextflow Tower

There are two branches:

`main` - this branch contains the script for running STARsolo on the FARM using Nextflow command line

`nextflow-tower` - this branch conrains the script for running STARsolo on the FARM using Nextflow Tower

## Contents of Repo:
* `main.nf` - the Nextflow pipeline that executes souporcell.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `examples/samples.txt` - samplefile tsv containing 3 fields: sampleID, path to BAM file, path to barcodes tsv.gz file (The order of these files is important!). These paths can be IRODs paths or local paths.
* `examples/RESUME-souporcell` - an example run script that executes the pipeline it has 2 hardcoded arguments: `/path/to/sample/file` and `/path/to/config/file` that need to be changed based on your local set up.
* `bin/soup_qc.sh` - a quick qc script that enables quick sanity checks that souporcell worked correctly.
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Pipeline Arguments:
* `--SAMPLEFILE` - The path to the sample file provided to the pipeline. This is a tab-separated file with one sample per line. Each line should contain a sample id, path to bam file, path to barcodes file (in that order!).
* `--outdir` - The path to where the results will be saved.
* `--K` - The number of donors in the samples. All samples must contain the same number of donors.
* `--known_genotypes` - "Whether to use the known_genotypes option. The requires a VCF with the known genotypes in it. The number of genotypes must match the number of donors supplied in the K parameter (Defualt "no" means not to use known_genotypes).
* `--barcodes_on_irods` - Tells pipeline whether to look for the gzipped barcodes file on IRODS or the FARM (default yes means look on IRODS).
* `--bam_on_irods` - Tells pipeline whether to look for the bam file on IRODS or the FARM (default yes means look on IRODS).
* `--reference_fasta` - Path to  genome fasta for pipeline to use (by default GRCh38 2020A is used). This argument is hardcoded and needs to be changed to your local path to the fasta file. 
* `--vcf` - Path to vcf file used for souporcell (default is 2p 1k genome with chr nomenclature). This argument is hardcoded and needs to be changed to your local path to the vcf file.
* `--skip_remap` - Whether to skip remapping in souporcell pipeline (default true means skip remapping).
* `--no_umi` - Tells the pipeline whether the BAM files have a UMI tag (default false means BAM file has UMI tag).
