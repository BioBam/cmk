version 1.1

workflow wdl_poc {
   input {
     File reads1
     File reads2
     File ref_txome
   }
   call fastqc {
      input:
         reads1 = reads1,
         reads2 = reads2,
   }
   call salmonindex {
     input:
        ref_txome = ref_txome,
  }
   call salmonalignquant {
     input:
        reads1 = reads1,
		reads2 = reads2,
		index = salmonindex.index
  }
}

task fastqc {
  input {
     File reads1
     File reads2
     env String C_USER = "rnica"
  }

  command {
     mkdir fastqc_res; fastqc --quiet "${reads1}" "${reads2}" --outdir fastqc_res
  }

  runtime {
    docker: "quay.io/nextflow/rnaseq-nf:v1.0"
  }

  output {
	 File fastqc_res = "fastqc_res"
  }
}

task salmonindex {
  input {
     File ref_txome
  }

  command {
     salmon index -t "${ref_txome}" -i index
  }

  runtime {
    docker: "combinelab/salmon:latest"
  }

  output {
	 File index = "index"
  }
}

task salmonalignquant {
  input {
	 File reads1
	 File reads2
     File index
  }

  command {
     salmon quant -i "${index}" -l A -1 "${reads1}" -2 "${reads2}" --validateMappings -o quant
  }

  runtime {
    docker: "combinelab/salmon:latest"
  }

  output {
	 File quant = "quant"
  }
}


