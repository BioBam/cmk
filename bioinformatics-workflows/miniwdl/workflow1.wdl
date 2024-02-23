version 1.1

workflow wdl_poc {
   input {
     String ref
     String output_s3_path
   }
   call salmonindex {
     input:
        ref = ref,
        output_s3_path = output_s3_path
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
     aws s3 cp "~{ref}" /tmp/inputs/ && salmon index -t /tmp/inputs/$(basename ~{ref}) -i /tmp/index && aws s3 sync /tmp/index "~{output_s3_path}index"
  }

  runtime {
    docker: "docker.io/robertbio/wf"
  }

  output {
     String index = "~{output_s3_path}index"
  }
}


