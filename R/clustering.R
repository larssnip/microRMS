#' @name genomeClustering
#' @title Clustering genomes based on shared RMS fragments
#'
#' @description Alters an RMS-object by clustering genomes who are too similar to distinguish.
#'
#' @param rms.obj A \code{list} with RMS data structures, see \code{\link{RMSobject}}.
#' @param max.cond Maximum condition value tolerated, see Details.
#'
#' @details This function will cluster genomes based on how similar they are in RMS fragment content.
#' The input is an RMS object (see \code{\link{RMSobject}}), with a \code{Genome.tbl} listing the genomes
#' and a \code{Cpn.mat} holding the information about RMS fragments in each genome. The latter has a 
#' column for genome.
#' 
#' The distance metric is 1 minus the correlation between the columns of \code{Cpn.mat}, i.e. if all RMS
#' fragments found in genome A are exactly those found in genome B, the distance between them is 0. Conversely,
#' if none of the fragments from genome A are found in genome B, the distance becomes 2 (maximum).
#' 
#' Clustering is done by hierchical clustering with complete linkage. Cutting the dendrogram tree at the various
#' heights produce a unique clustering of the genomes. Each of the produces a new \code{Cpn.mat}, with one column for each
#' genome cluster. The closer the genomes are to each other, the larger the condition number of this matrix. Thus, the
#' argument \code{max.cond} specifies the maximum tolerated condition value. A smaller value means harder clustering (fewer
#' genome clusters) but also more stable estimates of the abundance of each genome cluster. Conversely, a larger tolerated
#' condition value results in more genome clusters, but also more uncertain abundance estimates of each. 
#' 
#' More technical details: The problem of too similar genomes is reflected in a close to singular covariance matrix when
#' de-convolving the genome abundances. Clustering by correlation distance is, in theory, no guarantee against this. 
#' Two (or more) fairly uncorrelated genomes may still combine into something very correlated with a third genome. However,
#' in reality this rarely happens. Use the \code{\link{conditionValue}} function on the resulting \code{Cpn.mat} to
#' see if the clustering resulted in a lowering of the condition value to a tolerable size (e.g. around 1000 or less). 
#'
#' @return A updated RMS object, where \code{Genome.tbl} has an additional column named \code{members_genome_id}. This \code{Genome.tbl} 
#' typically has fewer rows than the original, one for each cluster, and this new column indicates which of the original
#' genomes are grouped into each cluster. Each cluster is represented by one of the original genomes (the cluster medoide).
#' The other objects inside the RMS object have also been updated accordingly.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{conditionValue}}.
#'
#' @importFrom Matrix crossprod Matrix
#'
#' @examples
#'
#' @export genomeClustering
#'
genomeClustering <- function(rms.obj, max.cond = 1000){
  D <- corrDist(rms.obj$Cpn.mat)
  tree <- hclust(as.dist(D), method = "complete")
  heights <- unique(c(0, sort(tree$height)))
  hh.lo <- heights[1]
  cc.lo <- log10(conditionValue(rms.obj$Cpn.mat))
  if(cc.lo <= log10(max.cond)) return(rms.obj)
  
  hh.up <- heights[length(heights) - 1]
  idx <- medoids(D, cutree(tree, h = hh.up))
  X <- rms.obj$Cpn.mat[,idx, drop = F]
  cc.up <- log10(conditionValue(X))
  if(cc.up > log10(max.cond)){
    cat("Impossible to reach condition value", max.cond, "with these genomes\n")
    return(rms.obj)
  }
  # we have: cc.up < max.cond < cc.lo
  
  # Search for optimal clustering
  idx2 <- 0
  idx1 <- 1
  while(idx1 != idx2){
    idx2 <- idx1
    dd <- abs(heights - (hh.lo + hh.up) / 2)
    idx1 <- which(dd == min(dd))[1]
    hh <- heights[idx1]
    med.idx <- medoids(D, cutree(tree, h = hh))
    cc <- log10(conditionValue(rms.obj$Cpn.mat[,med.idx, drop = F]))
    if(cc < log10(max.cond)){
      hh.up <- hh
      cc.up <- cc
    } else {
      hh.lo <- hh
      cc.lo <- cc
    }
    if(idx1 == idx2){
      idd <- max(which(c(cc.up, cc, cc.lo) <= log10(max.cond)))
      hh <- c(hh.up, hh, hh.lo)[idd[1]]
    }
  }
  
  # Finding cluster medoids and members
  clst <- cutree(tree, h = hh)
  med.idx <- medoids(D, clst)
  clst.tbl <- tibble(genome_id = rownames(D)[med.idx],
                     members_genome_id = tapply(rownames(D), clst, str_c, collapse = ","))

  # Pruning the data structure
  clst.tbl %>%
    left_join(rms.obj$Genome.tbl, by = "genome_id") -> rms.obj$Genome.tbl
  idx <- match(rms.obj$Genome.tbl$genome_id, colnames(rms.obj$Cpn.mat))
  rms.obj$Cpn.mat <- rms.obj$Cpn.mat[,idx]
  rms.obj$Cpn.mat <- rms.obj$Cpn.mat[which(rowSums(rms.obj$Cpn.mat) > 0),]
  rms.obj$Cluster.tbl %>%
    filter(Cluster %in% rownames(rms.obj$Cpn.mat)) -> rms.obj$Cluster.tbl

  # Updating N.unique
  pa <- (rms.obj$Cpn.mat > 0)
  pau <- pa[Matrix::rowSums(pa) == 1,]
  rms.obj$Genome.tbl %>%
    mutate(N_clusters = Matrix::colSums(pa)) %>%
    mutate(N_unique = Matrix::colSums(pau)) -> rms.obj$Genome.tbl

  return(rms.obj)
}




#' @name corrDist
#' @title Correlation distances in a sparse matrix
#'
#' @description Computes the correlation distances between pairs of columns of a sparse matrix.
#'
#' @param X A sparse matrix (Matrix).
#'
#' @details The correlation distance between two columns of a matrix is simply 1 minus the
#' correlation between them.
#'
#' @return An nxn distance matrix, when the input sparse matrix has n columns
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{genomeClustering}}.
#'
#' @importFrom Matrix crossprod colMeans
#'
#' @examples
#'
#' @export corrDist
#'
corrDist <- function(X){
  csd <- apply(X, 2, sd)
  if(sum(csd == 0) > 0){
    stop("Columns", which(csd == 0), "have no variance! Cannot compute correlation distances")
  }
  cm <- colMeans(X)
  n <- nrow(X)
  XX <- as.matrix(crossprod(X))
  CM <- cm %*% t(cm)
  CSD <- csd %*% t(csd)
  D <- 1 - (XX - n * CM) / ((n - 1) * CSD)
  return(D)
}


#' @name conditionValue
#' @title Condition value for a sparse matrix
#'
#' @description Computes the condition number for the cross-product of a sparse matrix to itself
#'
#' @param X A sparse matrix (Matrix).
#'
#' @details The condition number for X'X reflects how close this matrix is to being singular,
#' i.e. impossible to invert. The larger condition value, the closer it is to being singular.
#'
#' @return A positive numeric value
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{genomeClustering}}.
#'
#' @importFrom Matrix crossprod
#'
#' @examples
#'
#' @export conditionValue
#'
conditionValue <- function(X){
  lst <- svd(crossprod(X))
  cond <- max(lst$d)/min(lst$d)
  return(cond)
}

medoids <- function(D, clst){
  uclst <- unique(clst)
  med.idx <- rep(0, length(uclst))
  for(i in 1:length(uclst)){
    idx <- which(clst == uclst[i])
    d <- D[idx, idx, drop = F]
    rs <- rowSums(d)
    idd <- which(rs == min(rs))[1]
    med.idx[i] <- idx[idd]
  }
  return(med.idx)
}

