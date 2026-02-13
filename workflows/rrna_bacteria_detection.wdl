version 1.0

import "../tasks/samtools.wdl" as samtools
import "../tasks/fastp.wdl" as fastp
import "../tasks/sortmerna.wdl" as sortmerna
import "../tasks/spades.wdl" as spades
import "../tasks/blast.wdl" as blast
import "../tasks/blastprocess.wdl" as blastproc

workflow rrna_bacteria_detection {
  meta {
    description: "Workflow for identifying bacterial taxa from total RNAseq data"
  }
  input {
    Int base_cpu = 4
    Int base_mem = 8
    Int base_preempt = 3
    Int sortmerna_mem = 32
    Int blast_mem = 16

    File bamfile
    Array[File] rrna_db = ["gs://gcid-bacterial-public/rrna_databases/smr_v4.3_default_db_bacteria.fasta"]
    File blast_db = "gs://gcid-bacterial-public/rrna_databases/silva_blast_db.tar.gz"
    String blast_db_name = "blast_db"
    File taxmap_lsu = "gs://gcid-bacterial-public/rrna_databases/taxmap_slv_lsu_ref_nr_138.1.txt"
    File taxmap_ssu = "gs://gcid-bacterial-public/rrna_databases/taxmap_slv_ssu_ref_nr_138.1.txt"
  }

  call samtools.bam_to_fastq {
    input:
      bamfile = bamfile,
      cpu = base_cpu,
      mem = base_mem,
      preemptible = base_preempt
  }

  call fastp.qc {
    input:
      fastq1 = bam_to_fastq.fq1,
      fastq2 = bam_to_fastq.fq2,
      cpu = base_cpu,
      mem = base_mem,
      preemptible = base_preempt
  }

  Int total_reads_mb = ceil(size(qc.fq1, "MB") + size(qc.fq2, "MB"))
  Int sortmerna_cpu = if total_reads_mb > 5000 then 16 else if total_reads_mb > 1000 then 8 else 4
  Int sortmerna_disk = if total_reads_mb > 5000 then 500 else if total_reads_mb > 1000 then 300 else 100

  call sortmerna.extract_rrna {
    input:
      fastq1 = qc.fq1,
      fastq2 = qc.fq2,
      rrna_db = rrna_db,
      cpu = sortmerna_cpu,
      mem = sortmerna_mem,
      disk = sortmerna_disk,
      preemptible = base_preempt
  }

  Int rrna_reads_mb = ceil(size(extract_rrna.rrna_1, "MB") + size(extract_rrna.rrna_2, "MB"))
  Int spades_cpu = if rrna_reads_mb > 500 then 16 else if rrna_reads_mb > 100 then 8 else 4
  Int spades_mem = if rrna_reads_mb > 500 then 128 else if rrna_reads_mb > 100 then 72 else 16
  Int spades_disk = if rrna_reads_mb > 500 then 500 else if rrna_reads_mb > 100 then 300 else 100

  call spades.metaspades {
    input:
      reads1 = extract_rrna.rrna_1,
      reads2 = extract_rrna.rrna_2,
      cpu = spades_cpu,
      mem = spades_mem,
      disk = spades_disk,
      preemptible = base_preempt
  }

  Boolean successful_assembly = size(metaspades.contigs) > 0

  if (successful_assembly) {
    call blast.contigs_silva {
      input:
        contigs = metaspades.contigs,
        blast_db = blast_db,
        blast_db_name = blast_db_name,
        cpu = base_cpu,
        mem = blast_mem,
        preemptible = base_preempt
    }

    call blastproc.contig_tax {
      input:
        taxmap_lsu = taxmap_lsu,
        taxmap_ssu = taxmap_ssu,
        blast_results = contigs_silva.blast_results,
        cpu = base_cpu,
        mem = base_mem,
        preemptible = base_preempt
    }
  }
  
  output {
    File rrna_filtering_summary = extract_rrna.extract_rrna_summary
    File rrna_contigs = metaspades.contigs
    File? rrna_blast_results = contigs_silva.blast_results
    File? rrna_tax_summary = contig_tax.summary
  }
}