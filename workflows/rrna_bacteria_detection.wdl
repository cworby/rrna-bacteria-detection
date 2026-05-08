version 1.0

import "../tasks/samtools.wdl" as samtools
import "../tasks/fastp.wdl" as fastp
import "../tasks/sortmerna.wdl" as sortmerna
import "../tasks/spades.wdl" as spades
import "../tasks/vsearch.wdl" as vsearch
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
    Int blast_mem = 16
    Int vsearch_mem = 16
    Boolean run_chimera_filter = true

    File bamfile
    Array[File] rrna_db = ["gs://gcid-bacterial-public/rrna_databases/smr_v4.3_default_db_bacteria.fasta"]
    File blast_db = "gs://gcid-bacterial-public/rrna_databases/silva_blast_db.tar.gz"
    String blast_db_name = "blast_db"
    File taxmap_lsu = "gs://gcid-bacterial-public/rrna_databases/taxmap_slv_lsu_ref_nr_138.1.txt"
    File taxmap_ssu = "gs://gcid-bacterial-public/rrna_databases/taxmap_slv_ssu_ref_nr_138.1.txt"
    File chimera_ref_db = "gs://gcid-bacterial-public/rrna_databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta"
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
  Int sortmerna_cpu = if total_reads_mb > 3000 then 16 else if total_reads_mb > 1000 then 8 else 4
  Int sortmerna_mem = if total_reads_mb > 3000 then 64 else if total_reads_mb > 1000 then 48 else 32
  Int sortmerna_disk = if total_reads_mb > 3000 then 500 else if total_reads_mb > 1000 then 300 else 100

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
  Int spades_cpu = if rrna_reads_mb > 1000 then 32 else if rrna_reads_mb > 500 then 16 else if rrna_reads_mb > 100 then 8 else 4
  Int spades_mem = if rrna_reads_mb > 1000 then 256 else if rrna_reads_mb > 500 then 192 else if rrna_reads_mb > 100 then 128 else 32
  Int spades_disk = if rrna_reads_mb > 1000 then 750 else if rrna_reads_mb > 500 then 500 else if rrna_reads_mb > 100 then 300 else 100

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
    if (run_chimera_filter) {
      call vsearch.chimera_filter {
        input:
          contigs        = metaspades.contigs,
          chimera_ref_db = chimera_ref_db,
          cpu            = base_cpu,
          mem            = vsearch_mem,
          preemptible    = base_preempt
      }

      call blast.contigs_silva as contigs_silva_filtered {
        input:
          contigs       = chimera_filter.contigs_filtered,
          blast_db      = blast_db,
          blast_db_name = blast_db_name,
          cpu           = base_cpu,
          mem           = blast_mem,
          preemptible   = base_preempt
      }

      call blastproc.contig_tax as contig_tax_filtered {
        input:
          taxmap_lsu    = taxmap_lsu,
          taxmap_ssu    = taxmap_ssu,
          blast_results = contigs_silva_filtered.blast_results,
          cpu           = base_cpu,
          mem           = base_mem,
          preemptible   = base_preempt
      }
    }

    if (!run_chimera_filter) {
      call blast.contigs_silva {
        input:
          contigs       = metaspades.contigs,
          blast_db      = blast_db,
          blast_db_name = blast_db_name,
          cpu           = base_cpu,
          mem           = blast_mem,
          preemptible   = base_preempt
      }

      call blastproc.contig_tax {
        input:
          taxmap_lsu    = taxmap_lsu,
          taxmap_ssu    = taxmap_ssu,
          blast_results = contigs_silva.blast_results,
          cpu           = base_cpu,
          mem           = base_mem,
          preemptible   = base_preempt
      }
    }
  }

  output {
    File  rrna_filtering_summary      = extract_rrna.extract_rrna_summary
    File  rrna_contigs                = metaspades.contigs
    File? rrna_contigs_filtered       = chimera_filter.contigs_filtered
    File? rrna_chimera_stats          = chimera_filter.chimera_stats
    File? rrna_blast_results_filtered = contigs_silva_filtered.blast_results
    File? rrna_tax_summary_filtered   = contig_tax_filtered.summary
    File? rrna_blast_results          = contigs_silva.blast_results
    File? rrna_tax_summary            = contig_tax.summary
  }
}