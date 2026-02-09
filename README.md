# rrna-bacteria-detection

Pipeline to detect bacterial organisms from total RNAseq data based on rRNA assembly.

1. BAM -> FASTQ conversion (Samtools)
2. Read filtering (fastp)
3. rRNA read classification (SortMeRNA)
4. Assembly (MetaSPAdes)
5. Contig profiling (BLAST vs. SILVA database)
6. Contig taxonomic identification (LCA of best hits per contig)