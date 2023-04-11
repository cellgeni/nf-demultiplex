#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    ===================
    souporcell pipeline
    ===================
    This pipeline runs Souporcell.
    The parameters you need to input are:
      --SAMPLEFILE /full/path/to/sample/file 
      -K k-value
    This file should be a tsv with 3 columns: SAMPLEID\t/PATH/TO/BAM/t/PATH/TO/BARCODES
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/nf-souporcell/blob/main/examples/example.txt
    The pipeline uses a common variants file by default, to use a known genotypes vcf please use:
    --known_genotypes true
    The default reference fasta used is: GRCh38_v32_modified (2020A)
    The default common variants vcf used is: filtered_2p_1kgenomes_GRCh38
    To change these defaults input:
      --reference_fasta /path/to/reference/fasta
      --vcf /path/to/vcf
    """.stripIndent()
}

if (params.HELP) {
  helpMessage()
  exit 0
}

def errorMessage() {
    log.info"""
    ================
    souporcell error
    ================
    You failed to provide the --SAMPLEFILE parameter
    Please provide this parameter as follows:
      --SAMPLEFILE /full/path/to/sample/file
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

process get_data {

  input:
  val(sample) 

  output:
  env(NAME), emit: sampleid 
  path('*barcodes.tsv'), emit: barcode_file 
  path('*bam'), emit: bam_file 
  path('*bam.bai'), emit: index_file 

  shell:
  '''
  NAME=`echo !{sample} | cut -f 1 -d " "`
  bam_path=`echo !{sample} | cut -f 2 -d " "`
  barcodes_path=`echo !{sample} | cut -f 3 -d " "`
  
  if [[ "!{params.bam_on_irods}" == "no" ]]; then
    cp "${bam_path}" "${NAME}.bam"
    cp "${bam_path}.bai" "${NAME}.bam.bai"
  elif [[ "!{params.bam_on_irods}" == "yes" ]]; then
    iget -f -v -K "${bam_path}" "${NAME}.bam"
    iget -f -v -K "${bam_path}.bai" "${NAME}.bam.bai"
  else
    echo "incorrect bam option"
    exit 1
  fi
  if [[ "!{params.barcodes_on_irods}" == "no" ]]; then
    cp "${barcodes_path}" "${NAME}.barcodes.tsv.gz"
    gunzip -f "${NAME}.barcodes.tsv.gz"
  elif [[ "!{params.barcodes_on_irods}" == "yes" ]]; then
    iget -f -v -K "${barcodes_path}" "${NAME}.barcodes.tsv.gz"
    gunzip -f "${NAME}.barcodes.tsv.gz"
  else
    echo "incorrect barcodes option"
    exit 1
  fi
  '''
}


process run_souporcell {

  publishDir "${params.outdir}", mode: 'copy'

  input:
  val(name)
  path(barcodes)
  path(bam)
  path(index)

  output:
  path(name), emit: output 

  shell:
  '''
  common_or_known="--common_variants"
  if [[ "!{params.known_genotypes}" == "yes" ]]; then
    common_or_known="--known_genotypes"
  fi
  common_or_known="${common_or_known} !{params.vcf}"
  mkdir "!{name}"
  echo "souporcell_pipeline.py -i !{bam} -b !{barcodes} -f !{params.reference_fasta} -k !{params.K} ${common_or_known} -t 8 -o !{name} --skip_remap !{params.skip_remap} --no_umi !{params.no_umi}" > !{name}/cmd.txt
  souporcell_pipeline.py                \
  -i "!{bam}"                           \
  -b "!{barcodes}"                      \
  -f "!{params.reference_fasta}"        \
  -k "!{params.K}"                      \
  $common_or_known                      \
  -t 8                                  \
  -o "!{name}"                          \
  --skip_remap "!{params.skip_remap}"   \
  --no_umi "!{params.no_umi}"
  !{projectDir}/bin/soup_qc.sh !{name}
  '''
}

process run_shared_samples {
  
  publishDir "${params.outdir}/shared_samples", mode: 'copy'

  input:
  path(name) 
  
  output:
  path('map*')

  shell:
  '''
  cut -f 1 "!{launchDir}/!{params.SAMPLEFILE}" | while read s1; do 
    cut -f 1 "!{launchDir}/!{params.SAMPLEFILE}" | while read s2; do 
      shared_samples.py -1 $s1 -2 $s2 -n !{params.K} > "map!{params.K}.${s1}-${s2}" 2> "err!{params.K}.${s1}-${s2}"
    done 
  done
  '''
}

workflow {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | get_data | run_souporcell | collect | run_shared_samples
  }
}
