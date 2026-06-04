version 1.0

task bam_to_fastq {
  input {
    File bamfile
    Int cpu
    Int mem
    Int disk
    Int preemptible
  }

  command <<<
    READ_COUNT=$(samtools view -c ~{bamfile})
    if [ "$READ_COUNT" -eq 0 ]; then
      echo "ERROR: BAM file contains 0 reads: ~{bamfile}" >&2
      exit 1
    fi
    samtools fastq ~{bamfile} \
        -1 reads_R1.fastq.gz \
        -2 reads_R2.fastq.gz
  >>>

  output {
    File fq1 = "reads_R1.fastq.gz"
    File fq2 = "reads_R2.fastq.gz"
  }

  runtime {
    docker: "staphb/samtools:1.21"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
}
