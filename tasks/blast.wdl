task contigs_silva {
  input {
    File contigs
    File blast_db
    String blast_db_name
    Int threads = 8
  }

  command <<<
    mkdir silva_blast_db
    tar -xzvf ${blast_db} -C silva_blast_db --strip-components 1

    blastn \
      -query ${contigs} \
      -db silva_blast_db/${blast_db_name} \
      -out blastresults.tsv \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
      -evalue 1e-05 \
      -num_threads ${threads}
  >>>

  output {
    File blast_results = "blastresults.tsv"
  }

  runtime {
    docker: "ncbi/blast:2.14.1"
    cpu: threads
    memory: "16G"
  }
}