#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =======================
    demultiplexing pipeline
    =======================
    This pipeline runs Souporcell/Vireo.
    The parameters you need to input are:
      --SAMPLEFILE /full/path/to/sample/file 
      -K k-value
    This file should be a tsv with 3 columns: SAMPLEID\t/PATH/TO/BAM/t/PATH/TO/BARCODES
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/nf-souporcell/blob/main/examples/example.txt
    ----------
    souporcell
    ----------
    Souporcell uses a common variants file by default, to use a known genotypes vcf please use:
    --known_genotypes true
    The default reference fasta used is: GRCh38_v32_modified (2020A)
    The default common variants vcf used is: filtered_2p_1kgenomes_GRCh38
    To change these defaults input:
      --soc_fasta /path/to/reference/fasta
      --soc_vcf /path/to/vcf
    -----
    vireo
    -----
    The default reference vcf is: genome1K.phase3.SNP_AF5e2.chr1toX.hg38
    To change this default input:
      --vireo_vcf /path/to/vcf
    """.stripIndent()
}


def errorMessage() {
    log.info"""
    ====================
    demultiplexing error
    ====================
    You failed to provide the --SAMPLEFILE or --K parameter
    Please provide this parameter as follows:
      --SAMPLEFILE /full/path/to/sample/file
      --K donor-number
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

process get_data {

  input:
  tuple val(id), val(bam_path), val(barcodes_path)

  output:
  val(id)
  path('*barcodes.tsv')
  path('*bam')
  path('*bam.bai')

  shell:
  '''
  if "!{params.bam_on_irods}"; then
    iget -f -v -K "!{bam_path}" "!{id}.bam"
    iget -f -v -K "!{bam_path}.bai" "!{id}.bam.bai"
  else
    cp "!{bam_path}" "!{id}.bam"
    cp "!{bam_path}.bai" "!{id}.bam.bai"
  fi
  if "!{params.barcodes_on_irods}"; then
    iget -f -v -K "!{barcodes_path}" "!{id}.barcodes.tsv.gz"
    gunzip -f "!{id}.barcodes.tsv.gz"
  else
    cp "!{barcodes_path}" "!{id}.barcodes.tsv.gz"
    gunzip -f "!{id}.barcodes.tsv.gz"
  fi
  '''
}

process run_cellsnp {

  label 'vireo'

  publishDir "${params.outdir}/vireo/cellsnp", mode: 'copy'

  input:
  val(name)
  path(barcodes)
  path(bam)
  path(index)

  output:
  val(name)
  path("*cellsnp")

  shell:
  '''
  cellsnp-lite -s !{bam} -b !{barcodes} -O "!{name}-cellsnp" -R !{params.vireo_vcf} --gzip
  echo "cellsnp-lite -s !{bam} -b !{barcodes} -O !{name}-cellsnp -R !{params.vireo_vcf} --gzip" > "!{name}-cellsnp/cmd.txt"
  '''
}

process run_vireo {

  label 'vireo'

  publishDir "${params.outdir}/vireo", mode: 'copy'

  input:
  val(name)
  path(cellsnp)

  output:
  path("*vireo")

  shell:
  '''
  vireo -c !{cellsnp} -N !{params.K} -o "!{name}-vireo"
  echo "vireo -c !{cellsnp} -N !{params.K} -o !{name}-vireo" > "!{name}-vireo/cmd.txt"
  '''
}

process run_souporcell {

  publishDir "${params.outdir}/souporcell", mode: 'copy'

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
  common_or_known="${common_or_known} !{params.soc_vcf}"
  mkdir "!{name}"
  echo "souporcell_pipeline.py -i !{bam} -b !{barcodes} -f !{params.soc_fasta} -k !{params.K} ${common_or_known} -t 8 -o !{name} --skip_remap !{params.skip_remap} --no_umi !{params.no_umi}" > !{name}/cmd.txt
  souporcell_pipeline.py                \
  -i "!{bam}"                           \
  -b "!{barcodes}"                      \
  -f "!{params.soc_fasta}"              \
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
  
  publishDir "${params.outdir}/souporcell/shared_samples", mode: 'copy'

  input:
  path(name) 
  path(samplefile)
  
  output:
  path('map*')

  shell:
  '''
  cut -f 1 "!{samplefile}" | while read s1; do 
    cut -f 1 "!{samplefile}" | while read s2; do 
      shared_samples.py -1 $s1 -2 $s2 -n !{params.K} 1> "map!{params.K}.${s1}-${s2}" 2> "err!{params.K}.${s1}-${s2}"
    done 
  done
  '''
}

workflow vireo {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    params.K != null ?: errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2] ] } | get_data | run_cellsnp | run_vireo
  }
}

workflow souporcell {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    params.K != null ?: errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2] ] } | get_data | run_souporcell
    ch_soc = run_souporcell.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list)
  }
}

workflow all {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    params.K != null ?: errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2] ] } | get_data | set { ch_data }
    run_cellsnp(ch_data) | run_vireo
    run_souporcell(ch_data)
    ch_soc = run_souporcell.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list)
  }
}
