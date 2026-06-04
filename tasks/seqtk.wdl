version 1.0

task downsample_fastq {
  input {
    File fastq1
    File fastq2
    Int target_reads
    Int cpu
    Int mem
    Int disk
    Int preemptible
  }

  command <<<
    seqtk sample -s 42 ~{fastq1} ~{target_reads} | gzip > downsampled_R1.fq.gz
    seqtk sample -s 42 ~{fastq2} ~{target_reads} | gzip > downsampled_R2.fq.gz
  >>>

  output {
    File fq1 = "downsampled_R1.fq.gz"
    File fq2 = "downsampled_R2.fq.gz"
    Int reads_used = target_reads
  }

  runtime {
    docker: "staphb/seqtk:1.4"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: preemptible
  }
}