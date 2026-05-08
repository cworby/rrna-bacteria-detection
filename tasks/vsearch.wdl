version 1.0

task chimera_filter {
  input {
    File contigs
    File chimera_ref_db
    Int cpu
    Int mem
    Int preemptible
  }

  command <<<
    vsearch \
      --uchime_ref ~{contigs} \
      --db ~{chimera_ref_db} \
      --nonchimeras contigs_filtered.fasta \
      --chimeras contigs_chimeras.fasta \
      --uchimeout chimera_stats.txt \
      --threads ~{cpu}
  >>>

  output {
    File contigs_filtered = "contigs_filtered.fasta"
    File contigs_chimeras = "contigs_chimeras.fasta"
    File chimera_stats    = "chimera_stats.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/vsearch:2.28.1--h6a68c12_0"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk 50 HDD"
    preemptible: preemptible
  }
}
