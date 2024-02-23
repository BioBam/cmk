#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define pipeline input parameters - require the use of absolute paths
params.ref = '/path/to/ref.fasta'
params.outdir = 'results'


// Create reference transcriptome index using Salmom
process SALMON_INDEX {
  containerOptions "--env C_TOOL=SALMON_INDEX --env C_DATASET=${params.ref} --env C_WMS=nextflow"

  input: 
    path ref
  output:
    path index

  """
    salmon index -t '${ref}' -i index
  """
}

workflow {
  SALMON_INDEX(params.ref)
}
