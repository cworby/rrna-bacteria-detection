version 1.0

task metaspades {
  input {
    File reads1
    File reads2
    Int cpu
    Int mem
    Int disk
    Int preemptible
  }

  command <<<
    spades.py \
      --meta \
      -1 ~{reads1} \
      -2 ~{reads2} \
      -o metaspades_out \
      -t ~{cpu} \
      -m ~{mem}
  >>>

  output {
    File contigs = "metaspades_out/contigs.fasta"
  }

  runtime {
    docker: "staphb/spades:3.15.5"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: preemptible
  }
}

task rnaspades {
  input {
    File reads1
    File reads2
    Int cpu
    Int mem
    Int disk
    Int preemptible
  }

  command <<<
    spades.py \
      --rna \
      -1 ~{reads1} \
      -2 ~{reads2} \
      -o rnaspades_out \
      -t ~{cpu} \
      -m ~{mem}
  >>>

  output {
    File contigs = "rnaspades_out/transcripts.fasta"
  }

  runtime {
    docker: "staphb/spades:3.15.5"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: preemptible
  }
}