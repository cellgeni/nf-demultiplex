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
      --known_genotypes false - change to true to use known donor genotypes
    This file should be a tsv with 3 columns: SAMPLEID\t/PATH/TO/BAM/\t/PATH/TO/BARCODES
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/nf-demultiplex/blob/main/examples/samples.txt
    ----------
    souporcell
    ----------
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
    ln -s "!{bam_path}" "!{id}.bam"
    ln -s "!{bam_path}.bai" "!{id}.bam.bai"
  fi
  if "!{params.barcodes_on_irods}"; then
    iget -f -v -K "!{barcodes_path}" "!{id}.barcodes.tsv.gz"
    gunzip -f "!{id}.barcodes.tsv.gz"
  else
    ln -s "!{barcodes_path}" "!{id}.barcodes.tsv.gz"
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

process prepare_bam_for_merge {
  input:
  tuple val(name), path(barcodes), path(bam), path(index), val(K)

  output:
  tuple path("${name}.flagged.bam"), path("${name}.prefixed.barcodes.tsv"), val(K), emit: prepared

  shell:
  '''
  awk -v prefix="!{name}-" '{ print prefix $0 }' "!{barcodes}" > "!{name}.prefixed.barcodes.tsv"

  samtools view -h "!{bam}" \
    | awk -v prefix="!{name}-" 'BEGIN { OFS="\t" } /^@/ { print; next } { for(i=12;i<=NF;i++) { if($i ~ /^CB:Z:/) { sub(/^CB:Z:/, "", $i); $i = "CB:Z:" prefix $i; break } } print }' \
    | samtools view -b -o "!{name}.flagged.bam" -
  '''
}

process merge_bams {
  input:
  path(flagged_bams)
  path(prefixed_barcodes)
  val(k_values)

  output:
  tuple val('all'), path('merged_barcodes_all.tsv'), path('merged_all.bam'), path('merged_all.bam.bai'), path('merged_k.txt'), emit: merged_data
  path('samples_merged.tsv'), emit: samplefile

  script:
  def uniqueK = k_values.collect { it.toString() }.unique()
  if (uniqueK.size() != 1) {
    error "merge_bams=true requires all K values in SAMPLEFILE to be identical."
  }
  def merged_k = uniqueK[0]
  """
  ls -1 *.flagged.bam | sort > bam_all.list
  samtools merge -b bam_all.list --threads ${task.cpus} merged_all.bam
  samtools index merged_all.bam -@ ${task.cpus}

  for f in \$(ls -1 *.prefixed.barcodes.tsv | sort); do
    cat "\$f"
  done > merged_barcodes_all.tsv
  gzip -c merged_barcodes_all.tsv > merged_barcodes_all.tsv.gz

  printf "%s\\n" "${merged_k}" > merged_k.txt

  printf "%s\\t%s\\t%s\\t%s\\n" "all" "\$PWD/merged_all.bam" "\$PWD/merged_barcodes_all.tsv.gz" "${merged_k}" > samples_merged.tsv
  """
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
  !{projectDir}/bin/vireo.sh !{name} !{cellsnp} !{K} !{params.known_genotypes} !{params.snp_vcf} !{params.vireo_genoTag}
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
  
  while read -r name bam b n || [ -n "$name" ];
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

workflow resolve_inputs {
  take:
  ch_sample_list

  main:
  ch_data_raw = ch_sample_list | flatMap { it.readLines() } | map { it -> [it.split()[0], it.split()[1], it.split()[2], it.split()[3]] } | get_data

  if (params.merge_bams) {
    prepare_bam_for_merge(ch_data_raw)
    ch_merge_inputs = prepare_bam_for_merge.out.prepared
    merge_bams(
      ch_merge_inputs.map { it[0] }.collect(),
      ch_merge_inputs.map { it[1] }.collect(),
      ch_merge_inputs.map { it[2] }.collect()
    )
    ch_data = merge_bams.out.merged_data.map { id, barcodes, bam, bai, k_file -> [id, barcodes, bam, bai, k_file.text.trim()] }
    ch_sample_list_effective = merge_bams.out.samplefile
  } else {
    ch_data = ch_data_raw
    ch_sample_list_effective = ch_sample_list
  }

  emit:
  data = ch_data
  samplefile = ch_sample_list_effective
}

workflow vireo {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    resolve_inputs(ch_sample_list)
    ch_data = resolve_inputs.out.data

    ch_data_vireo = ch_data.filter { it[4].toInteger() > 1}
    run_cellsnp(ch_data_vireo) | run_vireo
  }
}

workflow souporcell {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    resolve_inputs(ch_sample_list)
    ch_data = resolve_inputs.out.data
    ch_sample_list_effective = resolve_inputs.out.samplefile

    run_souporcell(ch_data)
    ch_soc = run_souporcell.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list_effective)
    group_shared_samples(run_shared_samples.out.mapping,ch_soc,ch_sample_list_effective)
  }
}

workflow all {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    resolve_inputs(ch_sample_list)
    ch_data = resolve_inputs.out.data
    ch_sample_list_effective = resolve_inputs.out.samplefile

    ch_data_vireo = ch_data.filter { it[4].toInteger() > 1}
    run_cellsnp(ch_data_vireo) | run_vireo
    run_souporcell(ch_data)
    ch_soc = run_souporcell.out.output | collect
    ch_vir = run_vireo.out.output | collect
    run_shared_samples(ch_soc, ch_sample_list_effective)
    // quantify_shared_samples is older version of group_shared_samples and probably should be depricated
    quantify_shared_samples(run_shared_samples.out.mapping)
    group_shared_samples(run_shared_samples.out.mapping,ch_soc,ch_sample_list_effective)
    compare_methods(ch_soc, ch_vir, ch_sample_list_effective)
    if(params.check_sex){
      if (params.merge_bams) {
        log.warn "merge_bams=true with check_sex=true is not supported; skipping check_sex."
      } else {
        get_gex_data(ch_sample_list_effective)
        check_sex(ch_soc, ch_vir, get_gex_data.out, ch_sample_list_effective)
      }
    }
  }
}
