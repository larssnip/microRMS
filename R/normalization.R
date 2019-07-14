#' @name normLength
#' @title Length-normalization
#'
#' @description Reducing the RMS readcount bias that is due to fragment length.
#'
#' @param Y Matrix of readcounts, one row for each fragment cluster, one column for each sample.
#' @param L Fragment lengths.
#' @param max.train Integer indicating maximum number of data used to train the loess model
#' @param trunc Lowest readcount considered as valid when training loess model
#' @param verbose Logical, if TRUE text is written to the Console during computations.
#'
#' @details
#'
#' @return A matrix same size as \code{Y} with corrected readcounts.
#'
#' @author Lars Snipen.
#'
#' @seealso more here.
#'
#' @importFrom microseq readFasta gff2fasta
#' @importFrom stringr str_c
#' @importFrom dplyr mutate filter
#'
#' @examples more here.
#'
#' @export normLength
#'
normLength <- function(Y, L, max.train = 1000, trunc = 1, verbose = TRUE){
  for(k in 1:ncol(Y)){
    if(verbose) cat("Normalizing sample", k, "\r")
    idx <- which(Y[,k] >= trunc)
    idx <- idx[sample(idx, size = min(max.train, length(idx)))]
    mod <- loess(log2(Y[idx,k]) ~ L[idx], control = loess.control(surface = "direct"))
    v.hat <- predict(mod, newdata = L)
    ff <- 2^(max(v.hat) - v.hat)
    Y[,k] <- Y[,k] * ff
  }
  return(Y)
}





# normLibSize <- function(readcount.matrix, N.scale=FALSE){
#   library.size <- colSums(readcount.matrix, na.rm=T)
#   if(N.scale){
#     N.mark <- colSums(readcount.matrix>0, na.rm=T)
#     library.size <- library.size/N.mark
#   }
#   rrc <- t(t(readcount.matrix)/library.size)
#   return(rrc)
# }
