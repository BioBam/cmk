#!/usr/bin/env bash

function run_workflow {
    nextflow run example.nf \
        -c nextflow-aws.config \
        -with-docker docker.io/robertbio/wf \
        --ref $1 \
        --left $2 \
        --right $3 \
        --outdir $4/nextflow_results/ \
        -bucket-dir $4/nextflow_bucket_dir/
}

function stage_input_files {

  # Create a list of tuples with the input files
  # The first element of the tuple is the path to the left reads
  # The second element of the tuple is the path to the right reads
  declare -a input_files = (
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/004/SRR6312174/SRR6312174_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/004/SRR6312174/SRR6312174_2.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/005/SRR6312175/SRR6312175_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/005/SRR6312175/SRR6312175_2.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/000/SRR6312190/SRR6312190_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/000/SRR6312190/SRR6312190_2.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/007/SRR6312187/SRR6312187_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/007/SRR6312187/SRR6312187_2.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/002/SRR6312182/SRR6312182_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/002/SRR6312182/SRR6312182_2.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/001/SRR6312181/SRR6312181_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/001/SRR6312181/SRR6312181_2.fastq.gz"
  )


}

#Define the main function that runs the workflow
function main {
    # Define the paths to the input files
    ref=https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GG/AL/GGAL01/GGAL01.1.fsa_nt.gz
    left=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/004/SRR6312174/SRR6312174_1.fastq.gz
    right=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/004/SRR6312174/SRR6312174_2.fastq.gz

    # Define the paths to the output files
    data_s3=$BUCKET/final

    # Run the workflow
    run_workflow $ref $left $right $data_s3
}

main
