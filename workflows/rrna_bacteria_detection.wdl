
import "../tasks/samtools.wdl" as samtools
import "../tasks/fastp.wdl" as fastp
import "../tasks/sortmerna.wdl" as sortmerna
import "../tasks/spades.wdl" as spades
import "../tasks/blast.wdl" as blast
import "../tasks/blastprocess.wdl" as blastproc

workflow rrna_bacteria_detection {
  input {
    File bamfile
    File rrna_db_lsu = "resources/silva-bac-16s-id90.fasta"
    File rrna_db_ssu = "resources/silva-bac-23s-id98.fasta"
    File blast_db = "resources/silva_blast_db.tar.gz"
    File taxmap_lsu = "resources/taxmap_slv_lsu_ref_nr_138.1.txt"
    File taxmap_ssu = "resources/taxmap_slv_ssu_ref_nr_138.1.txt"
    String blast_db_name = "blast_db"
    Int spades_mem = 64
  }

  call samtools.bam_to_fastq {
    input:
      bamfile = bamfile
  }

  call fastp.qc {
    input:
      fastq1 = bam_to_fastq.fq1,
      fastq2 = bam_to_fastq.fq2
  }

  call sortmerna.extract_rrna {
    input:
      fastq1 = qc.fq1,
      fastq2 = qc.fq2,
      rrna_db_lsu = rrna_db_lsu,
      rrna_db_ssu = rrna_db_ssu
  }

  call spades.metaspades {
    input:
      reads1 = extract_rrna.rrna_1,
      reads2 = extract_rrna.rrna_2,
      memory_gb = spades_mem
  }

  call blast.contigs_silva {
    input:
      contigs = metaspades.contigs,
      blast_db = blast_db,
      blast_db_name = blast_db_name
  }

  call blastproc.contig_tax {
    input:
      taxmap_lsu = taxmap_lsu,
      taxmap_ssu = taxmap_ssu,
      blast_results = contigs_silva.blast_results
  }
  
  output {
    File contigs = metaspades.contigs
    File blast_results = contigs_silva.blast_results
    File blast_summary = contig_tax.summary
  }
}