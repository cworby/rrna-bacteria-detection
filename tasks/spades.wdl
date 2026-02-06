task metaspades {
  input {
    File reads1
    File reads2
    Int threads = 16
    Int memory_gb
  }

  command <<<
    spades.py \
      --meta \
      -1 ${reads1} \
      -2 ${reads2} \
      -o metaspades_out \
      -t ${threads} \
      -m ${memory_gb}
  >>>

  output {
    File contigs = "metaspades_out/contigs.fasta"
  }

  runtime {
    docker: "staphb/spades:3.15.5"
    cpu: threads
    memory: "${memory_gb}G"
    preemptible: 0
  }
}

task rnaspades {
  input {
    File reads1
    File reads2
    Int threads = 16
    Int memory_gb
  }

  command <<<
    spades.py \
      --rna \
      -1 ${reads1} \
      -2 ${reads2} \
      -o rnaspades_out \
      -t ${threads} \
      -m ${memory_gb}
  >>>

  output {
    File contigs = "rnaspades_out/transcripts.fasta"
  }

  runtime {
    docker: "staphb/spades:3.15.5"
    cpu: threads
    memory: "${memory_gb}G"
    preemptible: 0
  }
}