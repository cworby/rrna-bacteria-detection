#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-b", "--blast"), type="character", default=NULL, 
              help="Path to BLAST results", metavar="character"),
  make_option(c("-s", "--ssu"), type="character", default=NULL, 
              help="Path to SSU mapping file", metavar="character"),
  make_option(c("-l", "--lsu"), type="character", default=NULL, 
              help="Path to LSU mapping file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="summary.tsv", 
              help="Output filename [default= %default]", metavar="character")
)

# 2. Parse the arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

blast_out <- opt$blast
ssu_map <- opt$ssu
lsu_map <- opt$lsu
summary_table <- opt$out

#### Define functions
# String partitioning
breakup <- function(string, sep, elements=NULL, 
                    fromend=NULL, join=NULL, seq=FALSE) {
  if (seq) {
    ll <- lapply(strsplit(string, sep[1]), function(x){x[elements[1]]})
  } else {
    ll <- strsplit(string, sep[1])    
  }
  
  if (length(sep)>1) {
    for (i in 2:length(sep)) {
      if (seq) {
        ll <- lapply(ll, function(x){unlist(strsplit(x, sep[i]))[elements[i]]})
      } else {
        ll <- lapply(ll, function(x){unlist(strsplit(x, sep[i]))})
      }
    }
  }
  if (seq) {
    return(sapply(ll, function(x){paste(x,collapse=join)}))
  } else {
    return(sapply(ll, function(x){paste(x[c(elements, length(x)-fromend+1)],collapse=join)}))
  }
}

# Find lowest common ancestor
# Given a list of semicolon-separated taxonomic ranks, return the longest common path
tax_path_reducer <- function(x, separator=";") {
  ln <- min(unique(sapply(strsplit(x,separator),length)))
  lvl <- 1
  lca <- NA
  # Cycle down ranks from highest
  while(TRUE) {
    # multiple taxa called at this level?
    K <- unique(breakup(x,separator,1:lvl, join=separator))
    if (length(K)>1) {
      break
    } else {
      # only a single taxon called - go down to next level if available
      lvl <- lvl+1
      lca <- K
      if (lvl>ln) {
        break
      }
    }
  }
  return(lca)
}

# Assign taxonomy to contigs based on BLAST output
taxonomy_from_blast <- function(blasttable, # blast outfmt 6
                                tax_map, # list of tables (e.g. 16S, 23S) mapping blast DB IDs to names ('organism_name')
                                qnamecol="acc_full", # name of column in tax_map tables with blast DB IDs
                                eval_col=12, # BLAST table evalue column
                                idcol=3, # BLAST table ID column
                                matchlencol=4, # BLAST table match length column
                                subs=NULL) { # substitutions to make in taxonomic assignments (list of pairs); e.g. renamed taxa
  
  uq_contigs <- unique(blasttable[,1])
  tax_assignment_org <- list()
  tax_assignment_path <- list()  
  tax_rrna <- list()
  tax_maxID <- numeric(length(uq_contigs))
  tax_maxlen <- numeric(length(uq_contigs))  
  # For each contig, find best BLAST hits by e-value
  for (j in 1:length(uq_contigs)) {
    cat("Contig ", j, "/", length(uq_contigs), "\r", sep="")
    # subset the BLAST results pertaining to this contig
    B <- blasttable[which(blasttable[,1]==uq_contigs[j]),]
    # find BLAST hits with maximum e-value
    B <- B[which(B[,12]==max(B[,eval_col])),]
    contig_hits <- unique(B[,2])
    org <- NULL
    taxpath <- NULL
    rrna <- NULL
    maxID <- NULL
    # For each hit, extract taxonomy
    for (k in 1:length(contig_hits)) {
      # In which mapping table does the hit belong?
      db_hits <- which(unlist(lapply(tax_map, function(x)contig_hits[k]%in%x$acc_full)))
      for (db in db_hits) {
        # Organism
        org <- c(org,tax_map[[db]]$organism_name[which(tax_map[[db]][,qnamecol]==contig_hits[k])[1]])
        # Taxonomic path
        taxpath <- c(taxpath,tax_map[[db]]$path[which(tax_map[[db]][,qnamecol]==contig_hits[k])[1]])
        # Mapping table name
        rrna <- c(rrna, paste0(names(tax_map)[db],collapse=";"))
      }
    }
    # Make taxonomic substitutions, if any
    if (!is.null(subs)) {
      for (k in 1:length(subs)) {
        org <- gsub(subs[[k]][1], subs[[k]][2], org)
      }
    }
    tax_assignment_org[[j]] <- unique(org)
    tax_assignment_path[[j]] <- unique(taxpath)
    tax_rrna[[j]] <- unique(rrna)
    tax_maxID[j] <- max(B[,idcol])
    tax_maxlen[j] <- max(B[,matchlencol])      
  }
  # Semicolon separated best hits
  summary_org <- sapply(lapply(tax_assignment_org,function(x)unique(breakup(x, " ", 1:2, join=" "))),paste0,collapse=";")
  # LCA
  summary_path <- sapply(tax_assignment_path, tax_path_reducer)
  # DB hits
  summary_rrna <- sapply(tax_rrna, paste0, collapse=";")
  
  # taxonomic ranks of LCAs
  sp_hits <- as.numeric(sapply(summary_org,function(x)length(unlist(strsplit(x,";")))))
  # Best hit taxon for anything above kingdom level
  best_hit <- summary_org
  best_hit[which(sp_hits>1)] <- breakup(summary_path[which(sp_hits>1)],";",fromend=1)
  best_hit[grep("uncultured",best_hit)] <- breakup(summary_path[grep("uncultured",best_hit)],";",fromend=1)
  best_hit[grep("metagenome",best_hit)] <- breakup(summary_path[grep("metagenome",best_hit)],";",fromend=1)
  best_hit <- gsub(" sp[.]","",best_hit)
  # Extract length of contig (assuming SPAdes contig names)
  contiglen <- as.numeric(breakup(uq_contigs,"_",4))
  blast_tax_table <- cbind(contig=uq_contigs, length=contiglen, best_hit=best_hit, maxID=tax_maxID,
                           max_matchlen=round(100*tax_maxlen/contiglen,2), rrna=summary_rrna, path=summary_path, sp_hits=summary_org)
  return(data.frame(blast_tax_table))
}
  
# Read in SILVA mapping tables
tax_map_ssu <- read.table(ssu_map,
                          sep="\t", header = TRUE, quote = "%", comment.char = "")
tax_map_ssu$acc_full <- paste(tax_map_ssu$primaryAccession, tax_map_ssu$start, tax_map_ssu$stop, sep=".")
tax_map_lsu <- read.table(lsu_map, 
                          sep="\t", header = TRUE, quote = "%", comment.char = "")
tax_map_lsu$acc_full <- paste(tax_map_lsu$primaryAccession, tax_map_lsu$start, tax_map_lsu$stop, sep=".")

# Name mapping tables
tax_map <- list(`SSU`=tax_map_ssu, `LSU`=tax_map_lsu)

if (file.info(blast_out)$size == 0) {
  blast_tax_table <- data.frame(contig=character(), length=character(), best_hit=character(),
                                maxID=numeric(), max_matchlen=numeric(), rrna=character(),
                                path=character(), sp_hits=character(), blast_bacteria=logical(),
                                blast_best_lvl=character())
  write.table(blast_tax_table, file=summary_table, sep="\t", quote=FALSE, row.names=FALSE)
} else {
  blasttable <- read.table(blast_out, sep="\t", header=FALSE)
  blast_tax_table <- taxonomy_from_blast(blasttable, tax_map, qnamecol="acc_full",
                                        eval_col=12, idcol=3, matchlencol=4, subs=list(c("Mycoplasmoides", "Mycoplasma")))

  blast_tax_table$path <- gsub(";uncultured","",blast_tax_table$path)
  blast_tax_table$blast_bacteria <- breakup(blast_tax_table$path, ";", 1)=="Bacteria"
  lvl <- sapply(strsplit(blast_tax_table$path, ";"),length)
  blast_tax_table$blast_best_lvl <- c("K", "P", "C", "O", "F", "G")[lvl]
  blast_tax_table$blast_best_lvl[which(lvl==6 & grepl(" ", blast_tax_table$best_hit))] <- "S"

  write.table(blast_tax_table, file=summary_table, sep="\t", quote=FALSE, row.names=FALSE)
}

