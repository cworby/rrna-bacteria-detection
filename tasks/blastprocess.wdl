version 1.0

task contig_tax {

  input {
    File blast_results
    File taxmap_lsu
    File taxmap_ssu
    Int cpu
    Int mem
    Int preemptible
  }

  command <<<
    blast_taxonomy.r \
      --blast ~{blast_results} \
      --ssu ~{taxmap_ssu} \
      --lsu ~{taxmap_lsu} \
      --out contig_taxonomy.tsv
  >>>

  output {
    File summary = "contig_taxonomy.tsv"
  }

  runtime {
    docker: "cworby/parse_taxonomy:v1.0.1"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }
}