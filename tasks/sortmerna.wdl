task extract_rrna {
  input {
    File fastq1
    File fastq2
    File rrna_db_lsu
    File rrna_db_ssu
    Int threads = 8
  }

  command <<<
    sortmerna \
      --ref ${rrna_db_lsu} \
      --ref ${rrna_db_ssu} \
      --reads ${fastq1} \
      --reads ${fastq2} \
      --paired_in \
      --workdir . \
      --fastx \
      --out2 \
      --aligned rrna \
      --threads ${threads}
  >>>

  output {
    File rrna_1 = "rrna_fwd.fq.gz"
    File rrna_2 = "rrna_rev.fq.gz"
  }

  runtime {
    docker: "nanozoo/sortmerna:4.3.4--7b48a67"
    cpu: threads
    memory: "32G"
    preemptible: 3
  }
}