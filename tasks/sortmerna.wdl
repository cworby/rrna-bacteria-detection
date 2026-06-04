version 1.0

task extract_rrna {
  input {
    File fastq1
    File fastq2
    Array[File] rrna_db
    Int cpu
    Int mem
    Int disk
    Int preemptible
  }

  command <<<
    sortmerna \
      ~{sep=" " prefix("--ref ", rrna_db)} \
      --reads ~{fastq1} \
      --reads ~{fastq2} \
      --paired_in \
      --workdir . \
      --fastx \
      --out2 \
      --aligned rrna \
      --threads ~{cpu}
  >>>

  output {
    File rrna_1 = "rrna_fwd.fq.gz"
    File rrna_2 = "rrna_rev.fq.gz"
    File extract_rrna_summary = "rrna.log"
  }

  runtime {
    docker: "cworby/sortmerna:4.3.7"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: preemptible
  }
}