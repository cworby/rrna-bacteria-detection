
task contig_tax {

  input {
    File blast_results
    File taxmap_lsu
    File taxmap_ssu
  }

  command <<<
    blast_taxonomy.r \
      --blast ${blast_results} \
      --ssu ${taxmap_ssu} \
      --lsu ${taxmap_lsu} \
      --out contig_taxonomy.tsv
  >>>

  output {
    File summary = "contig_taxonomy.tsv"
  }

  runtime {
    docker: "gcr.io/gcid-bacterial/gcid-bacterial/parse_taxonomy:v1.0.1"
    cpu: 1
    memory: "4G"
  }
}