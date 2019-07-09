






pruneUnique <- function(rmsdb.obj){
  idx <- which(rmsdb.obj$Cluster.tbl$N.genomes == 1)
  rmsdb.obj$Cluster.tbl <- rmsdb.obj$Cluster.tbl[idx,]
  rmsdb.obj$Cpn <- rmsdb.obj$Cpn[idx,]
  pa <- (rmsdb.obj$Cpn > 0)
  pau <- pa[rowSums(pa) == 1,]
  rmsdb.obj$Genome.tbl$N.clusters <- colSums(pa)
  rmsdb.obj$Genome.tbl$N.unique <- colSums(pau)
  return(rmsdb.obj)
}

jaccard <- function(rmsdb.obj){
  pa <- rmsdb.obj$Cpn > 0
  ng <- ncol(pa)
  D <- matrix(0, nrow = ng, ncol = ng)
  rownames(D) <- colnames(pa)
  colnames(D) <- colnames(pa)
  for(i in 1:(ng-1)){
    for(j in (i+1):ng){
      D[i,j] <- D[j,i] <- 1 - (sum(pa[,i] & pa[,j]) / sum(pa[,i] | pa[,j]))
    }
  }
  return(as.dist(D))
}




# pruneAmplicons <- function(rmsdb.obj, min.amplicons=1){
#   pa <- rmsdb.obj$Cpn>0
#
#   # re-ordering rows to have the most shared amplicons first
#   n.share <- rowSums(pa)
#   idx <- order(n.share, decreasing=T)
#   pas <- pa[idx,]
#
#   idx <- which(rowSums(pas)==1)
#   pas.unik <- pas[idx,]
#   pas.share <- pas[-idx,]
#   N.unique <- colSums(pas.unik)
#   N.share <- colSums(pas.share)
#
#   # the pruning
#   is.pruned <- rep(F,nrow(pas.share))
#   names(is.pruned) <- rownames(pas.share)
#   keep.pruning <- T
#   while(keep.pruning){
#     N.share <- colSums(pas.share[which(!is.pruned),,drop=F])
#     idd.in <- which(N.unique+N.share>min.amplicons & N.share>0)
#     idd.out <- which(N.unique+N.share<=min.amplicons)
#     if(length(idd.in)>0){
#       idx <- which(rowSums(pas.share[,idd.in,drop=F])>0 &
#                      rowSums(pas.share[,idd.out,drop=F])==0 &
#                      !is.pruned)
#       if(length(idx)>0){
#         is.pruned[idx[1]] <- T
#       } else {
#         keep.pruning <- F
#       }
#     } else {
#       # no more genomes to prune
#       keep.pruning <- F
#     }
#     cat("Is pruned:", sum(is.pruned), "\n")
#   }
#   pruned <- names(is.pruned)[which(is.pruned)]
#   idx <- which(!(rownames(pa) %in% pruned))
#   rmsdb.obj$Meta <- rmsdb.obj$Meta[idx,]
#   rmsdb.obj$Cpn  <- rmsdb.obj$Cpn[idx,]
#   rmsdb.obj$Centroids <- rmsdb.obj$Centroids[idx,]
#   return(rmsdb.obj)
# }
