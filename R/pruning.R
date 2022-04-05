#' @name pruneSamples
#' @title Discarding samples
#'
#' @description Discarding samples having too few reads mapped.
#'
#' @param rms.obj A \code{list} with the tables \code{Sample.tbl}, see details below.
#' @param mapping.threshold Minimum fraction of reads being mapped.
#' @param verbose Turn on/off output text during processing (logical).
#' 
#' @details The \code{rms.obj} must be a list with the required data structures. It must have a 
#' \code{Sample.tbl}, which is added to such an object with \code{\link{addSampleTable}}. It must
#' also have a \code{Readcount.mat}, which is added by the \code{\link{readMapper}}.
#' 
#' The pruning means sample where the fraction of mapped reads to the total number of reads is below
#' the threshold set by \code{mapping.threshold}. If very small fraction of reads are mapped, it probably indicates
#' none of the genomes are present in the sample, and even if abundance estimates can be made they are
#' also very unreliable in such cases.
#' 
#' @return An RMS object where the \code{Sample.tbl} and \code{Readcount.mat} have been pruned, i.e.
#' some samples have been discarded.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}, \code{\link{readMapper}}, \code{\link{rmscols}}.
#'
#' @importFrom dplyr slice
#'
#' @examples
#'
#' @export pruneSamples
#'
pruneSamples <- function(rms.obj, mapping.threshold = 0, verbose = TRUE){
  if(!exists("Sample.tbl", where = rms.obj)) stop("The rms.obj must contain a Sample.tbl")
  if(!exists("Readcount.mat", where = rms.obj)) stop("The rms.obj must contain a Readcount.mat")
  if(length(grep("reads_total", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column reads_total")
  if(length(grep("reads_mapped", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column reads_mapped")
  idx <- which((rms.obj$Sample.tbl$reads_mapped / rms.obj$Sample.tbl$reads_total) > mapping.threshold)
  if(length(idx) > 0){
    rms.obj$Sample.tbl <- slice(rms.obj$Sample.tbl, idx)
    rms.obj$Readcount.mat <- rms.obj$Readcount.mat[,idx,drop = FALSE]
  } else {
    warning("The mapping.threshold is too large, leaving no samples after pruning, skips pruning")
  }
  return(rms.obj)
}
