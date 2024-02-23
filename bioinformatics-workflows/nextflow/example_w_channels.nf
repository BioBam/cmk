#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define pipeline input parameters 
params.ref = 's3://gwfcore-runworkflows/bioproject/PRJNA419302/ref/GGAL01.1.fsa_nt.gz'
params.reads = [
  ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312174_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312174_2.fastq.gz'],
  // ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312175_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312175_2.fastq.gz'],
  // ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312181_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312181_2.fastq.gz'],
  // ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312182_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312182_2.fastq.gz'],
  // ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312187_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312187_2.fastq.gz'],
  ['s3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312190_1.fastq.gz', 's3://gwfcore-runworkflows/bioproject/PRJNA419302/samples/SRR6312190_2.fastq.gz'],
  // Add more entries if needed...
]
params.outdir = 'results'

Channel
  .from(params.reads)
  .map { reads -> tuple("group", tuple(reads[0], reads[1])) }
  .set { reads_ch }

// Create reference transcriptome index using Salmon
process SALMON_INDEX {
  containerOptions "--env C_TOOL=SALMON_INDEX --env C_DATASET=$ref --env C_WMS=nextflow"
  cpus 2

  input: 
    path ref
  output:
    path 'index'


  """
    salmon index -t '${ref}' -i index
  """
}

// Transcriptome alignment and quantification using Salmon
process SALMON_ALIGN_QUANT {
  containerOptions "--env C_TOOL=SALMON_ALIGN_QUANT --env C_DATASET=${reads[0]} --env C_WMS=nextflow"
  cpus 2
  
  // publishDir params.outdir
  publishDir "${params.outdir}/${task.process}", mode: 'copy'

  input:
    path index
    tuple val(group), path(reads)
  output:
    path 'quant'

  """
    salmon quant -i $index -l A -1 '${reads[0]}' -2 '${reads[1]}' --validateMappings -o quant
  """
}

// FastQC
process FASTQC {
  containerOptions "--env C_TOOL=FASTQC --env C_DATASET=${reads[0]} --env C_WMS=nextflow"
  cpus 2
  publishDir params.outdir

  input:
    tuple val(group), path(reads)
  output:
    path 'qc'

  """
    mkdir qc && fastqc --quiet '${reads[0]}' '${reads[1]}' --outdir qc
  """
}

workflow {
  index_ch = SALMON_INDEX(params.ref)

  SALMON_ALIGN_QUANT(index_ch, reads_ch)
  FASTQC(reads_ch)
}