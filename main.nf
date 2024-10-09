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
    This file should be a tsv with 3 columns: SAMPLEID\t/PATH/TO/BAM/\t/PATH/TO/BARCODES
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/nf-demultiplex/blob/main/examples/samples.txt
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
      --snp_vcf /path/to/vcf
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
  label 'get_data'
  
  input:
  tuple val(id), val(bam_path), val(barcodes_path), val(K)

  output:
  tuple val(id), path('*barcodes.tsv'), path('*bam'), path('*bam.bai'), val(K)

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
  tuple val(name), path(barcodes), path(bam), path(index), val(K)

  output:
  val(name)
  path("*cellsnp")
  val(K)

  shell:
  '''
  !{projectDir}/bin/cellsnp.sh !{name} !{barcodes} !{bam} !{params.snp_vcf}
  '''
}

process run_vireo {

  label 'vireo'

  publishDir "${params.outdir}/vireo", mode: 'copy'

  input:
  val(name)
  path(cellsnp)
  val(K)

  output:
  path("*vireo"), emit: output

  shell:
  '''
  !{projectDir}/bin/vireo.sh !{name} !{cellsnp} !{K}
  '''
}

process run_souporcell {

  publishDir "${params.outdir}/souporcell", mode: 'copy'

  input:
  tuple val(name), path(barcodes), path(bam), path(index), val(K)

  output:
  path(name), emit: output

  shell:
  '''
  !{projectDir}/bin/souporcell.sh !{name} !{barcodes} !{bam} !{params.known_genotypes} !{params.soc_vcf} !{params.soc_fasta} !{K} !{params.skip_remap} !{params.no_umi}
  !{projectDir}/bin/soup_qc.sh !{name}
  '''
}

process run_shared_samples {
  
  publishDir "${params.outdir}/souporcell/shared_samples", mode: 'copy'

  input:
  path(name) 
  path(samplefile)
  
  output:
  path('map*'), emit: mapping

  shell:
  '''
  !{projectDir}/bin/shared_samples.sh !{samplefile}
  '''
}

process quantify_shared_samples {

  publishDir "${params.outdir}/souporcell/shared_samples", mode: 'copy'

  input:
  path(mapping)

  output:
  path('quantification')

  shell:
  '''
  mkdir -p quantification/clusters
  for shared in map*; do
    #Get sample id from shared sample file name
    sample=`echo $shared | cut -f 2 -d .`
    #Extract cluster comparison tsv from shared samples file
    tail -n +6 $shared > "quantification/clusters/${sample}.cluster"
  done
  cd quantification
  mkdir -p loss-tables
  for comparison in clusters/*.cluster; do
    sample=`echo $comparison | cut -f 2 -d / | cut -f 1 -d .`
    python !{projectDir}/bin/shared-samples-quant.py --input $comparison --sample $sample
  done
  echo -e "experiment1_cluster\texperiment2_cluster\tloss\tsample_id" > total.tsv
  cat loss-tables/* >> total.tsv
  '''
}

process group_shared_samples {
  label 'rscript'
  publishDir "${params.outdir}/souporcell", mode: 'copy'

  input:
  path(shared_samples, stageAs: "shared_samples/*") 
  path(souporcell, stageAs: "souporcell/*") 
  path(samplefile)

  output:
  path('group_samples')

  shell:
  '''
  Rscript !{projectDir}/bin/group_genotypes.R shared_samples souporcell !{samplefile} !{params.ngenomes}
  '''
}

process compare_methods {
  label 'rscript'
  
  publishDir "${params.outdir}/compare_methods", mode: 'copy'
  
  input:
  path(souporcell, stageAs: "souporcell/*") 
  path(vireo, stageAs: "vireo/*")
  path(samplefile)
  
  output:
  path('soc2vir*')

  shell:
  '''
  Rscript !{projectDir}/bin/compare_methods.R !{samplefile}
  '''
}

process get_gex_data {
  label 'get_data'
  
  input:
  path(samplefile)
  
  output:
  path('out/*')

  shell:
  '''
  mkdir out
  # load expression data 
  if "!{params.barcodes_on_irods}"; then
    cmd='iget'
  else
    cmd='cp'
  fi
  
  while read name bam b n; 
  do 
   ${cmd} -r $(dirname $bam)/filtered_feature_bc_matrix out/${name}
  done < !{samplefile}
  '''
}

process check_sex {
  label 'rscript'
  
  publishDir "${params.outdir}/", mode: 'copy'
  
  input:
  path(souporcell, stageAs: "souporcell/*") 
  path(vireo, stageAs: "vireo/*")
  path(gex, stageAs: "gex/*")
  path(samplefile)
  
  output:
  path('sex')

  shell:
  '''
  Rscript !{projectDir}/bin/check_sex.R !{samplefile}
  '''
}

workflow vireo {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2], it.split()[3]] } | get_data | run_cellsnp | run_vireo
  }
}

workflow souporcell {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2], it.split()[3] ] } | get_data | run_souporcell
    ch_soc = run_souporcell.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list)
    group_shared_samples(run_shared_samples.out.mapping,ch_soc,ch_sample_list)
  }
}

workflow all {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2] , it.split()[3] ] } | get_data | set { ch_data }
	  ch_data_vireo = ch_data.filter { it[4].toInteger() > 1}
    run_cellsnp(ch_data_vireo) | run_vireo
    run_souporcell(ch_data)
    ch_soc = run_souporcell.out.output | collect
    ch_vir = run_vireo.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list)
    // quantify_shared_samples is older version of group_shared_samples and probably should be depricated
    quantify_shared_samples(run_shared_samples.out.mapping)
    group_shared_samples(run_shared_samples.out.mapping,ch_soc,ch_sample_list)
    compare_methods(ch_soc, ch_vir, ch_sample_list)
    get_gex_data(ch_sample_list)
    if(params.check_sex)
      check_sex(ch_soc, ch_vir, get_gex_data.out, ch_sample_list)
  }
}
