#' @name genomeClustering
#' @title Clustering genomes based on shared RMS fragments
#'
#' @description Alters an RMS-object by clustering genomes who are too similar to distinguish.
#'
#' @param rms.obj A \code{list} with RMS data structures, see \code{\link{RMSobject}}.
#' @param max.corr Maximum correlation between genomes, see Details.
#' @param verbose Logical, turning on/off screen report on progress during clustering.
#'
#' @details This function will cluster genomes based on how similar they are in 
#' RMS fragment content. If two genomes are very similar with respect to RMS 
#' fragment content, their corresponding columns in \code{rms.obj$Cpn.mat} are 
#' highly correlated. If genomes are too correlated it is impossible to estimate
#' their abundance separately, thus such genomes must be seen as a cluster, and 
#' we only estimate the abundance of this cluster.
#' 
#' Genomes are first represented as a graph, where two genomes are connected 
#' with an edge if they share at least 1 RMS fragment. For many highly unrelated
#' genomes this step will result in a disconnected graph, i.e. groups of genomes 
#' not sharing any fragments between them. These graph components 
#' are the first grouping of the genomes. Genomes sharing no fragments will 
#' always end up in different clusters anyway. This step saves a lot of memory
#' when the RMS object contain many and unrelated genomes, since all these
#' computations can be done on a sparse \code{Matrix}.
#'
#' Next, clustering within each graph component is done by hierarchical 
#' clustering with complete linkage. The distance metric is 1 minus correlation, 
#' i.e. genomes with a large correlation close to 1.0) has a small distance 
#' (close to 0.0). The \code{max.corr} argument indicates where to cut the 
#' dendrogram tree to group the genomes. With \code{max.corr = 0.8} we cut the 
#' dendrogram at distance \code{0.2}. Note that the distance matrix cannot be
#' a sparse \code{Matrix}, and if too many genomes in the \code{rms.obj} are
#' in the same graph component, you may run into memory problems.
#' 
#' More technical details: The problem of too similar genomes will be reflected in 
#' a close to singular covariance matrix when de-convolving the genome 
#' abundances. Clustering by correlation distance is, in theory, no guarantee 
#' against this. Two (or more) fairly uncorrelated genomes may still combine 
#' into something very correlated with a third genome. However, in reality this 
#' rarely happens with RMS data. Use the \code{\link{conditionValue}} function on the 
#' resulting \code{rms.obj$Cpn.mat} to see if the clustering resulted in a fairly
#' low condition value to a tolerable size (e.g. around 1e+3 to 1e+4 or less). 
#'
#' @return An updated RMS object, where \code{Genome.tbl} has an additional 
#' column named \code{members_genome_id}. This \code{Genome.tbl} typically has 
#' fewer rows than the original, one for each cluster, and this new column
#' indicates which of the original genomes are grouped into each cluster. Each
#' cluster is represented by one of the original genomes (the cluster medoide).
#' The other objects inside the RMS object have also been updated accordingly.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{corrDist}}, \code{\link{conditionValue}}.
#'
#' @importFrom Matrix crossprod Matrix
#' @importFrom igraph graph_from_adjacency_matrix clusters
#'
#' @examples
#'
#' @export genomeClustering
#'
genomeClustering <- function(rms.obj, max.corr = 0.80, verbose = TRUE){
  genome_id <- colnames(rms.obj$Cpn.mat)
  n.genomes <- length(genome_id)
  if(verbose) cat("genomeClustering:\n...starts with", n.genomes, "genomes...\n")
  
  if(verbose) cat("The genome graph...")
  adj.mat <- crossprod(rms.obj$Cpn.mat)                  # works on a sparse Matrix
  genome.graph <- graph_from_adjacency_matrix(adj.mat)   # works on a sparse Matrix
  graph.lst <- clusters(genome.graph)
  if(verbose) cat("has", graph.lst$no, "disconnected components\n")
  
  clustering <- graph.lst$membership * 10^ceiling(log10(max(graph.lst$csize)))
  is.medoide <- rep(FALSE, n.genomes)
  for(i in 1:graph.lst$no){
    if(verbose) cat("   clustering component", i, "by correlation distance...")
    idx <- which(graph.lst$membership == i)
    D <- corrDist(rms.obj$Cpn.mat[,idx])
    tree <- hclust(as.dist(D), method = "complete")
    clst <- cutree(tree, h = (1 - max.corr))
    if(verbose) cat("resulted in", length(unique(clst)), "clusters\n")
    med.idx <- medoids(D, clst)
    clustering[idx] <- clustering[idx] + clst
    is.medoide[idx[med.idx]] <- TRUE
    if(verbose) cat("resulted in", length(unique(clst)), "clusters\n")
  }
  clst.tbl <- tibble(genome_id = genome_id[is.medoide],
                     members_genome_id = tapply(genome_id, clustering, str_c, collapse = ","))
  
  if(verbose) cat("...pruning the data structure...\n")
  rms.obj$Genome.tbl <- clst.tbl %>%
    left_join(rms.obj$Genome.tbl, by = "genome_id")
  idx <- match(rms.obj$Genome.tbl$genome_id, colnames(rms.obj$Cpn.mat))
  rms.obj$Cpn.mat <- rms.obj$Cpn.mat[,idx]
  rms.obj$Cpn.mat <- rms.obj$Cpn.mat[which(rowSums(rms.obj$Cpn.mat) > 0),]
  rms.obj$Cluster.tbl <- rms.obj$Cluster.tbl %>%
    filter(Cluster %in% rownames(rms.obj$Cpn.mat))
  
  if(verbose) cat("...updating N_unique...\n")
  pa <- (rms.obj$Cpn.mat > 0)
  pau <- pa[Matrix::rowSums(pa) == 1,]
  rms.obj$Genome.tbl %>%
    mutate(N_clusters = Matrix::colSums(pa)) %>%
    mutate(N_unique = Matrix::colSums(pau)) -> rms.obj$Genome.tbl
  
  if(verbose) cat("...done!\n")
  return(rms.obj)
  
  
  
  
  
  
  
  
  
  
  # if(verbose) cat("...computing correlation distances...\n")
  # D <- corrDist(rms.obj$Cpn.mat)
  # tree <- hclust(as.dist(D), method = "complete")
  # 
  # if(verbose) cat("...finding maximum condition value = ")
  # heights <- unique(c(0, sort(tree$height)))
  # hh.lo <- heights[1]
  # cc.lo <- log10(conditionValue(rms.obj$Cpn.mat))
  # if(verbose) cat(10^cc.lo, "...\n")
  # if(cc.lo <= log10(max.cond)) return(rms.obj)
  # 
  # if(verbose) cat("...finding minimum condition value = ")
  # hh.up <- heights[length(heights) - 1]
  # idx <- medoids(D, cutree(tree, h = hh.up))
  # X <- rms.obj$Cpn.mat[,idx, drop = F]
  # cc.up <- log10(conditionValue(X))
  # if(verbose) cat(10^cc.up, "...\n")
  # if(cc.up > log10(max.cond)){
  #   cat("Impossible to reach condition value", max.cond, "with these genomes\n")
  #   return(rms.obj)
  # }
  # # we have: cc.up < max.cond < cc.lo
  
  # if(verbose) cat("...searching for optimal clustering...\n")
  # idx2 <- 0
  # idx1 <- 1
  # while(idx1 != idx2){
  #   idx2 <- idx1
  #   dd <- abs(heights - (hh.lo + hh.up) / 2)
  #   idx1 <- which(dd == min(dd))[1]
  #   hh <- heights[idx1]
  #   med.idx <- medoids(D, cutree(tree, h = hh))
  #   if(verbose) cat("...", length(med.idx), "clusters...\n")
  #   cc <- log10(conditionValue(rms.obj$Cpn.mat[,med.idx, drop = F]))
  #   if(cc < log10(max.cond)){
  #     hh.up <- hh
  #     cc.up <- cc
  #   } else {
  #     hh.lo <- hh
  #     cc.lo <- cc
  #   }
  #   if(idx1 == idx2){
  #     idd <- max(which(c(cc.up, cc, cc.lo) <= log10(max.cond)))
  #     hh <- c(hh.up, hh, hh.lo)[idd[1]]
  #   }
  # }
  
  # if(verbose) cat("...finding cluster members and medoides...\n")
  # clst <- cutree(tree, h = hh)
  # med.idx <- medoids(D, clst)
  # clst.tbl <- tibble(genome_id = rownames(D)[med.idx],
  #                    members_genome_id = tapply(rownames(D), clst, str_c, collapse = ","))
  # 
  # if(verbose) cat("...pruning the data structure...\n")
  # clst.tbl %>%
  #   left_join(rms.obj$Genome.tbl, by = "genome_id") -> rms.obj$Genome.tbl
  # idx <- match(rms.obj$Genome.tbl$genome_id, colnames(rms.obj$Cpn.mat))
  # rms.obj$Cpn.mat <- rms.obj$Cpn.mat[,idx]
  # rms.obj$Cpn.mat <- rms.obj$Cpn.mat[which(rowSums(rms.obj$Cpn.mat) > 0),]
  # rms.obj$Cluster.tbl %>%
  #   filter(Cluster %in% rownames(rms.obj$Cpn.mat)) -> rms.obj$Cluster.tbl
  # 
  # if(verbose) cat("...updating N_unique...\n")
  # pa <- (rms.obj$Cpn.mat > 0)
  # pau <- pa[Matrix::rowSums(pa) == 1,]
  # rms.obj$Genome.tbl %>%
  #   mutate(N_clusters = Matrix::colSums(pa)) %>%
  #   mutate(N_unique = Matrix::colSums(pau)) -> rms.obj$Genome.tbl
  # 
  # if(verbose) cat("...done!\n")
  # return(rms.obj)
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
  CSD <- csd %*% t(csd)
  cm <- colMeans(X)
  CM <- cm %*% t(cm)
  n <- nrow(X)
  XX <- crossprod(X)
  D <- 1 - (as.matrix(XX) - n * CM) / ((n - 1) * CSD)
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

