#' @name readMapper
#' @title Mapping reads
#'
#' @description Mapping reads to an RMS database using VSEARCH
#'
#' @param sample.files A vector of filenames, see details below.
#' @param centroids.file A filename, see details below.
#' @param identity Identity (numeric, 0.0-1.0) threshold for mapping.
#' @param threads Number of threads to be used by VSEARCH.
#' @param min.length Minimum fragment length (bases).
#' @param unmapped.dir If supplied, a directory to put non-mapped reads.
#' @param verbose Logical, if TRUE text is written to the Console during computations.
#'
#' @details The \code{sample.files} must be a vector of (full path) names of FASTA-files containing
#' the reads from the samples. The \code{centroids.file} is the name of a FASTA-file with the
#' database centroids to map against. Thus, the reads from each sample will be mapped to these
#' sequences, using the \code{identity} threshold. This results in a table of read-counts. The \code{sample.files}
#' are used as column names in this table. The first token in the Headers of the \code{centroids.file} are
#' used as row names.
#'
#' @return A
#'
#' @author Lars Snipen.
#'
#' @seealso This is version latest
#'
#' @importFrom microseq gregexpr
#' @importFrom data.table fread
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom stringr word str_c
#' @importFrom dplyr mutate filter select
#' @importFrom tibble tibble
#'
#' @examples more here.
#'
#' @export readMapper
#'
readMapper <- function(sample.files, centroids.file, identity, tmp.dir = ".",
                       threads = 1, min.length = 30, unmapped.dir = NULL,
                       verbose = TRUE){
  sample.id <- gsub("^.+/", "", gsub("\\.[a-z]+$", "", sample.files))
  ctr <- readFasta(centroids.file)
  tags <- word(ctr$Header, 1, 1, sep = ";")
  RMS.counts <- matrix(0, nrow = length(tags), ncol = length(sample.id))
  rownames(RMS.counts) <- tags
  colnames(RMS.counts) <- sample.id
  unmap <- ""
  quiet <- "--quiet"
  if(verbose) quiet <- ""
  for(i in 1:length(sample.id)){
    if(!is.null(unmapped.dir)){
      unmap <- paste("--notmatched", file.path(unmapped.dir, paste0(sample.id[i], ".fasta")))
    }
    if(verbose) cat("Mapping reads from sample", sample.id[i], "...\n")
    cmd <- paste("vsearch",
                 quiet,
                 "--threads", threads,
                 "--usearch_global", sample.files[i],
                 "--db", centroids.file,
                 "--id", identity,
                 "--iddef", "2",
                 "--minseqlength", min.length,
                 "--sizein --sizeout",
                 unmap,
                 "--otutabout", file.path(tmp.dir, "rmstab.txt"))
    system(cmd)
    rms.tbl <- suppressMessages(read_delim(file.path(tmp.dir, "rmstab.txt"), delim = "\t"))
    idx <- match(rms.tbl[[1]], tags)
    RMS.counts[idx,i] <- rms.tbl[[2]]
  }
  return(RMS.counts)
}



readcountSummary <- function(readcounts, amplicon.dbase){
  ntax <- ncol(amplicon.dbase$Cpn)
  nsamp <- ncol(readcounts)
  tab <- data.frame(Taxon=colnames(amplicon.dbase$Cpn),
                    N.amplicons=rep(0,ntax),
                    stringsAsFactors=F)
  stats <- matrix(0, nrow=ntax, ncol=nsamp*2)
  colnames(stats) <- c(paste0("RpC.", colnames(readcounts)),
                       paste0("ZeroFrac.", colnames(readcounts)))
  for(i in 1:ntax){
    idx <- which(amplicon.dbase$Cpn[,i]>0)
    tab$N.amplicons[i] <- length(idx)
    rcpm <- colMeans(readcounts[idx,]/amplicon.dbase$Cpn[idx,i])
    zerof <- apply(readcounts[idx,], 2, function(x){sum(x==0)})/length(idx)
    stats[i,] <- c(rcpm, zerof)
  }
  tab <- cbind(tab, stats)
  return(tab)
}
