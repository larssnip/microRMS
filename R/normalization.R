#' @name normLength
#' @title Length-normalization
#'
#' @description Reducing the RMS readcount bias that is due to fragment length.
#'
#' @param rms.obj A \code{list} with RMS data structures, see \code{\link{RMSobject}}.
#' @param max.train Integer indicating maximum number of data points used to train the loess model.
#' @param readcount.low Lowest readcount considered as valid when training loess model.
#' @param verbose Logical, if TRUE text is written to the Console during computations.
#'
#' @details RMS amplicon signals (readcounts) usually have some bias due to varying lengths. This
#' function attempts to remove this bias. It has been described in detail in Snipen et al, 2019.
#'
#' @return An rms object similar to \code{rms.obj}, but where the \code{Readcount.mat}
#' matrix has been normalized, columnwise.
#'
#' @author Lars Snipen.
#'
#' @importFrom microseq readFasta gff2fasta
#' @importFrom stringr str_c
#' @importFrom dplyr mutate filter
#'
#' @examples
#'
#' @export normLength
#'
normLength <- function(rms.obj, max.train = 1000, readcount.low = 1, verbose = TRUE){
  frg.length <- rms.obj$Cluster.tbl$Length
  for(k in 1:ncol(rms.obj$Readcount.mat)){
    if(verbose) cat("Normalizing", colnames(rms.obj$Readcount.mat)[k], "\r")
    idx <- which(rms.obj$Readcount.mat[,k] >= readcount.low)
    idx <- idx[sample(idx, size = min(max.train, length(idx)))]
    mod <- loess(log2(rms.obj$Readcount.mat[idx,k]) ~ frg.length[idx],
                 control = loess.control(surface = "direct"))
    v.hat <- predict(mod, newdata = frg.length)
    ff <- 2^(max(v.hat) - v.hat)
    rms.obj$Readcount.mat[,k] <- rms.obj$Readcount.mat[,k] * ff
  }
  return(rms.obj)
}
