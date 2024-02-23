workflow wdl_poc {
  input {
    Array[String] reads1_pairs
    Array[String] reads2_pairs
    String ref
    String output_s3_path
  }

  scatter (i in range(length(reads1_pairs))) {
    call fastqc {
      input:
        reads1 = reads1_pairs[i],
        reads2 = reads2_pairs[i],
        output_s3_path = output_s3_path,
    }
  }
  
  call salmonindex {
    input:
      ref = ref,
      output_s3_path = output_s3_path,
  }

  scatter (i in range(length(reads1_pairs))) {
    call salmonalignquant {
      input:
        reads1 = reads1_pairs[i],
        reads2 = reads2_pairs[i],
        index = salmonindex.index,
        output_s3_path = output_s3_path,
    }
  }
}



task fastqc {
  input {
     String reads1
     String reads2
     String output_s3_path
     env String C_TOOL = "FASTQC"
     env String C_WMS = "miniwdl"
     env String C_DATASET = basename(reads1)
  }

  command <<<
    aws s3 cp "~{reads1}" /tmp/input/ \
    && aws s3 cp "~{reads2}" /tmp/input/ \
    && mkdir /tmp/fastqc_res \
    && fastqc --quiet "/tmp/input/$(basename ~{reads1})" "/tmp/input/$(basename ~{reads2})" --outdir /tmp/fastqc_res \
    && aws s3 sync /tmp/fastqc_res "~{output_s3_path}fastqc_res"
  >>>

  runtime {
    cpu: 2
    docker: "docker.io/robertbio/wf"
  }

  output {
	 String fastqc_res = "~{output_s3_path}fastqc_res"
  }
}

task salmonindex {
  input {
     String ref
     String output_s3_path
     env String C_TOOL = "SALMON_INDEX"
     env String C_WMS = "miniwdl"
     env String C_DATASET = basename(ref)
  }

  command {
     aws s3 cp "~{ref}" /tmp/input/ && salmon index -t /tmp/input/$(basename ~{ref}) -i /tmp/index && aws s3 sync /tmp/index "~{output_s3_path}index"
  }

  runtime {
    cpu: 2
    docker: "docker.io/robertbio/wf"
  }

  output {
     String index = "~{output_s3_path}index"
  }
}

task salmonalignquant {
  input {
	 String reads1
	 String reads2
   String index
   String output_s3_path
   env String C_TOOL = "SALMON_ALIGN_QUANT"
   env String C_WMS = "miniwdl"
   env String C_DATASET = basename(reads1)
  }

  command <<<
    mkdir -p /tmp/input/index

    aws s3 sync "~{index}" /tmp/input/index && \
    aws s3 cp "~{reads1}" /tmp/input/ && \
    aws s3 cp "~{reads2}" /tmp/input/ && \
    salmon quant -i /tmp/input/index -l A -1 /tmp/input/$(basename ~{reads1}) -2 /tmp/input/$(basename ~{reads2}) --validateMappings -o /tmp/quant && \
    aws s3 sync /tmp/quant "~{output_s3_path}quant"
  >>>

  runtime {
    cpu: 2
    docker: "docker.io/robertbio/wf"
  }

  output {
	 String quant = "~{output_s3_path}quant"
  }
}



