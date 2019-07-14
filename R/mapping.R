#' @name readMapper
#' @title Mapping reads
#'
#' @description Mapping reads to an RMS database using VSEARCH
#'
#' @param sample.files A vector of FASTA-filenames (text).
#' @param centroids.file A FASTA-filename (text).
#' @param identity Identity threshold for mapping (0.0-1.0).
#' @param threads Number of threads to be used by vsearch (integer).
#' @param min.length Minimum fragment length used by vsearch (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param tmp.dir If supplied, the output directory used by vsearch.
#' @param unmapped.dir If supplied, a directory to output unmapped reads.
#' 
#'
#' @details The \code{sample.files} must be a vector of names of FASTA-files containing
#' the reads from the samples. The \code{centroids.file} is the name of a FASTA-file with the
#' database centroids to map against. Thus, the reads from each sample will be mapped to these
#' sequences, using the \code{identity} threshold.
#' 
#' This results in a matrix of read-counts. The column names of this matrix are the filenames in
#' \code{sample.files}, after the removal of the path and the file extension. The first token in
#' the Headers of the \code{centroids.file} are used as row names.
#'
#' @return A matrix of read-counts, with one row for each fragment cluster and one column for each
#' sample file.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSdbase}}, \code{\link{rmscols}}.
#'
#' @importFrom microseq readFasta
#' @importFrom readr read_delim
#' @importFrom stringr word
#'
#' @examples
#' \dontrun{
#' xpth <- file.path(path.package("microrms"),"extdata")
#' smp.files <- file.path(xpth, list.files(xpth, pattern = "Sample"))
#' ctr.file <- file.path(xpth, "centroids.fasta")
#' 
#' read.counts <- readMapper(smp.files, ctr.file)
#' }
#'
#' @export readMapper
#'
readMapper <- function(sample.files, centroids.file, identity = 0.99, 
                       threads = 1, min.length = 30, verbose = TRUE,
                       tmp.dir = NULL, unmapped.dir = NULL){
  sample.id <- gsub("^.+/", "", gsub("\\.[a-z]+$", "", sample.files))
  ctr <- readFasta(centroids.file)
  tags <- word(ctr$Header, 1, 1, sep = ";")
  RMS.counts <- matrix(0, nrow = length(tags), ncol = length(sample.id))
  rownames(RMS.counts) <- tags
  colnames(RMS.counts) <- sample.id
  unmap <- ""
  quiet <- "--quiet"
  if(verbose) quiet <- ""
  if(is.null(tmp.dir)) tab.file <- tempfile(pattern = "rmstab", fileext = ".txt")
  else tab.file <- file.path(tmp.dir, "rmstab.txt")
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
                 "--otutabout", tab.file)
    system(cmd)
    rms.tbl <- suppressMessages(read_delim(tab.file, delim = "\t"))
    idx <- match(rms.tbl[[1]], tags)
    RMS.counts[idx,i] <- rms.tbl[[2]]
  }
  if(is.null(tmp.dir)) ok <- file.remove(tab.file)
  return(RMS.counts)
}



# readcountSummary <- function(readcounts, amplicon.dbase){
#   ntax <- ncol(amplicon.dbase$Cpn)
#   nsamp <- ncol(readcounts)
#   tab <- data.frame(Taxon=colnames(amplicon.dbase$Cpn),
#                     N.amplicons=rep(0,ntax),
#                     stringsAsFactors=F)
#   stats <- matrix(0, nrow=ntax, ncol=nsamp*2)
#   colnames(stats) <- c(paste0("RpC.", colnames(readcounts)),
#                        paste0("ZeroFrac.", colnames(readcounts)))
#   for(i in 1:ntax){
#     idx <- which(amplicon.dbase$Cpn[,i]>0)
#     tab$N.amplicons[i] <- length(idx)
#     rcpm <- colMeans(readcounts[idx,]/amplicon.dbase$Cpn[idx,i])
#     zerof <- apply(readcounts[idx,], 2, function(x){sum(x==0)})/length(idx)
#     stats[i,] <- c(rcpm, zerof)
#   }
#   tab <- cbind(tab, stats)
#   return(tab)
# }
