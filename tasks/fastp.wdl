task qc {
  input {
    File fastq1
    File fastq2
  }

  command <<<
    fastp \
      --in1 ${fastq1} \
      --in2 ${fastq2} \
      --out1 trimmed_R1.fq.gz \
      --out2 trimmed_R2.fq.gz \
      --detect_adapter_for_pe \
      --json qc_summary.json \
      --html qc_summary.html
  >>>

  output {
    File fq1 = "trimmed_R1.fq.gz"
    File fq2 = "trimmed_R2.fq.gz"
    File qc_summary = "qc_summary.json"
  }

  runtime {
    docker: "staphb/fastp:0.23.4"
    memory: "8G"
    preemptible: 3
  }
}