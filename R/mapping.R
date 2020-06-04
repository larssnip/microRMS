#' @name readMapper
#' @title Mapping reads
#'
#' @description Mapping reads to clusters in an RMS object using VSEARCH
#'
#' @param rms.obj A \code{list} with the tables \code{Sample.tbl} and \code{Cluster.tbl}, see details below.
#' @param fa.dir A path to where the fasta files with reads are located.
#' @param identity Identity threshold for mapping (0.0-1.0).
#' @param threads Number of threads to be used by vsearch (integer).
#' @param min.length Minimum fragment length used by vsearch (integer).
#' @param verbose Turn on/off output text during processing (logical).
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
#' @seealso \code{\link{RMSobject}}, \code{\link{rmscols}}.
#'
#' @importFrom microseq readFasta
#' @importFrom readr read_delim
#' @importFrom stringr word
#'
#' @examples
#'
#' @export readMapper
#'
readMapper <- function(rms.obj, fa.dir, identity = 0.99,
                       threads = 1, min.length = 30, verbose = TRUE){
  if(!exists("Sample.tbl", rms.obj)) stop("The rms.obj must contain a Sample.tbl")
  if(length(grep("sample_id", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column sample_id")
  if(length(grep("reads_file", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column reads_file")
  if(!exists("Cluster.tbl", rms.obj)) stop("The rms.obj must contain a Cluster.tbl")
  if(length(grep("Header", colnames(rms.obj$Cluster.tbl))) == 0) stop("The rms.obj$Cluster.tbl must contain a column Header")
  if(length(grep("Sequence", colnames(rms.obj$Cluster.tbl))) == 0) stop("The rms.obj$Cluster.tbl must contain a column Sequence")
  ok <- available.external("vsearch")
  tags <- word(rms.obj$Cluster.tbl$Header, 1, 1, sep = ";")
  RMS.counts <- matrix(0, nrow = length(tags), ncol = nrow(rms.obj$Sample.tbl))
  rownames(RMS.counts) <- tags
  colnames(RMS.counts) <- rms.obj$Sample.tbl$sample_id
  centroids.file <- tempfile(pattern = "centroid", fileext = ".fasta")
  writeFasta(rms.obj$Cluster.tbl, out.file = centroids.file)
  tab.file <- tempfile(pattern = "rmstab", fileext = ".txt")
  for(i in 1:nrow(rms.obj$Sample.tbl)){
    if(verbose) cat("Mapping reads from sample", rms.obj$Sample.tbl$sample_id[i], "...\n")
    cmd <- paste("vsearch",
                 "--threads", threads,
                 "--usearch_global", file.path(fa.dir, rms.obj$Sample.tbl$reads_file[i]),
                 "--db", centroids.file,
                 "--id", identity,
                 "--iddef", "2",
                 "--minseqlength", min.length,
                 "--sizein --sizeout",
                 "--otutabout", tab.file)
    system(cmd)
    rms.tbl <- suppressMessages(read_delim(tab.file, delim = "\t"))
    idx <- match(rms.tbl[[1]], tags)
    RMS.counts[idx,i] <- rms.tbl[[2]]
  }
  ok <- file.remove(tab.file, centroids.file)
  return(c(rms.obj, list(Readcount.mat = RMS.counts)))
}


#' @name addSampleTable
#' @title Adds sample table
#'
#' @description Adds a sample table to an rms object
#'
#' @param rms.obj An \code{list} with RMS data structures, see \code{\link{RMSobject}}.
#' @param sample.tbl A table with sample metadata (data.frame or tibble).
#'
#' @details This small function just adds a table with metadata about samples to an 
#' already existing \code{rms.obj}. The latter is a \code{list}, and this function only
#' ensures the added element is correctly named \code{Sample.tbl}.
#' 
#' The \code{sample.tbl} must contain at least the two columns \code{sample_id} and 
#' \code{reads_file}. The first is a unique text to identify each sample, the latter is
#' the name of the fasta file with processed reads, see \code{\link{readMapper}}.
#' 
#' @return A \code{list} similar to the input \code{rms.obj}, but with one new element added.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}.
#'
#' @examples
#'
#' @export addSampleTable
#'
addSampleTable <- function(rms.obj, sample.tbl){
  return(c(rms.obj, list(Sample.tbl = sample.tbl)))
}
