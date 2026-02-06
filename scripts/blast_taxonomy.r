#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-b", "--blast"), type="character", default=NULL, 
              help="Path to BLAST results", metavar="character"),
  make_option(c("-s", "--ssu"), type="character", default=NULL, 
              help="Path to SSU mapping file", metavar="character"),
  make_option(c("-l", "--lsu"), type="character", default=NULL, 
              help="Path to LSU mapping file", metavar="character"),
  make_option(c("-i", "--idcol"), type="numeric", default=3, 
              help="Column in BLAST file for identity", metavar="NUM"),
  make_option(c("-g", "--lencol"), type="numeric", default=4, 
              help="Path to LSU mapping file", metavar="NUM"),
  make_option(c("-e", "--evalcol"), type="numeric", default=12, 
              help="Path to LSU mapping file", metavar="NUM"),
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
idcol <- opt$idcol
lencol <- opt$lencol
evalcol <- opt$evalcol

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
tax_path_reducer <- function(x, separator=";") {
  ln <- min(unique(sapply(strsplit(x,separator),length)))
  lvl <- 1
  lca <- NA
  while(TRUE) {
    K <- unique(breakup(x,separator,1:lvl, join=separator))
    if (length(K)>1) {
      break
    } else {
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
taxonomy_from_blast <- function(blasttable, tax_map, verbose=TRUE, qnamecol="acc_full", 
                                eval_col=12, idcol=3, matchlencol=4, subs=NULL) {
  
  uq_contigs <- unique(blasttable[,1])
  tax_assignment_org <- list()
  tax_assignment_path <- list()  
  tax_rrna <- list()
  tax_maxID <- numeric(length(uq_contigs))
  tax_maxlen <- numeric(length(uq_contigs))    
  for (j in 1:length(uq_contigs)) {
    if (verbose) cat(j, "/", length(uq_contigs), "\r")
    # subset the BLAST results pertaining to this contig
    B <- blasttable[which(blasttable[,1]==uq_contigs[j]),]
    # find BLAST hits with maximum e-value
    B <- B[which(B[,12]==max(B[,eval_col])),]
    contig_hits <- unique(B[,2])
    org <- NULL
    taxpath <- NULL
    rrna <- NULL
    maxID <- NULL
    for (k in 1:length(contig_hits)) {
      db_hits <- which(unlist(lapply(tax_map, function(x)contig_hits[k]%in%x$acc_full)))
      for (db in db_hits) {
        org <- c(org,tax_map[[db]]$organism_name[which(tax_map[[db]][,qnamecol]==contig_hits[k])[1]])
        taxpath <- c(taxpath,tax_map[[db]]$path[which(tax_map[[db]][,qnamecol]==contig_hits[k])[1]])
        rrna <- c(rrna, paste0(names(tax_map)[db],collapse=";"))
      }
    }
    
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
  summary_org <- sapply(lapply(tax_assignment_org,function(x)unique(breakup(x, " ", 1:2, join=" "))),paste0,collapse=";")
  summary_path <- sapply(tax_assignment_path, tax_path_reducer)
  summary_rrna <- sapply(tax_rrna, paste0, collapse=";")
  
  
  sp_hits <- as.numeric(sapply(summary_org,function(x)length(unlist(strsplit(x,";")))))
  best_hit <- summary_org
  best_hit[which(sp_hits>1)] <- breakup(summary_path[which(sp_hits>1)],";",fromend=1)
  best_hit[grep("uncultured",best_hit)] <- breakup(summary_path[grep("uncultured",best_hit)],";",fromend=1)
  best_hit[grep("metagenome",best_hit)] <- breakup(summary_path[grep("metagenome",best_hit)],";",fromend=1)
  best_hit <- gsub(" sp.","",best_hit)
  contiglen <- as.numeric(breakup(uq_contigs,"_",4))
  blast_tax_table <- cbind(contig=uq_contigs, length=contiglen, best_hit=best_hit, maxID=tax_maxID,
                           max_matchlen=round(100*tax_maxlen/contiglen,2), rrna=summary_rrna, path=summary_path, sp_hits=summary_org)
  return(blast_tax_table)
}
  
tax_map_ssu <- read.table(ssu_map,
                          sep="\t", header = TRUE, quote = "%", comment.char = "")
tax_map_ssu$acc_full <- paste(tax_map_ssu$primaryAccession, tax_map_ssu$start, tax_map_ssu$stop, sep=".")
tax_map_lsu <- read.table(lsu_map, 
                          sep="\t", header = TRUE, quote = "%", comment.char = "")
tax_map_lsu$acc_full <- paste(tax_map_lsu$primaryAccession, tax_map_lsu$start, tax_map_lsu$stop, sep=".")

tax_map <- list(`16S`=tax_map_ssu, `23S`=tax_map_lsu)

blasttable <- read.table(blast_out, sep="\t", header=FALSE)
blast_tax_table <- taxonomy_from_blast(blasttable, tax_map, verbose=TRUE, qnamecol="acc_full", 
                                       eval_col=evalcol, idcol=idcol, matchlencol=lencol, subs=list(c("Mycoplasmoides", "Mycoplasma")))
write.table(blast_tax_table, file=summary_table, sep="\t", quote=FALSE, row.names=FALSE)


