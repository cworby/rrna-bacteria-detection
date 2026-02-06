task bam_to_fastq {
  input {
    File bamfile
  }

  command <<<
    samtools fastq ${bamfile} \
        -1 reads_R1.fastq.gz \
        -2 reads_R2.fastq.gz
  >>>

  output {
    File fq1 = "reads_R1.fastq.gz"
    File fq2 = "reads_R2.fastq.gz"
  }

  runtime {
    docker: "staphb/samtools:1.21"
    memory: "8G"
    preemptible: 3
  }
}