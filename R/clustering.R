#' @name genomeClustering
#' @title Clustering genomes based on shared amplicon patterns
#'
#' @description Clusters genomes to obtain a copy number matrix with a maximum
#' condition number in order to estimate composition.
#'
#' @param Cpn A matrix of amplicon cluster copy numbers, a column for each genome (cluster).
#' @param steps The correlation distances to visit during clustering.
#' @param cond.max The maximum condition number tolerated.
#'
#' @details more here.
#'
#' @return A new copy number matrix, with fewer columns than the original if clustering was needed
#' to get the condition number below \code{cond.max}.
#'
#' @author Lars Snipen.
#'
#' @seealso more here.
#'
#' @importFrom Matrix crossprod Matrix
#'
#' @examples more here.
#'
#' @export genomeClustering
#'
genomeClustering <- function(Cpn, steps=cumsum(0.01*1.1^(0:15)), cond.max=1000){
  require(Matrix)
  N <- nrow(Cpn)
  P <- ncol(Cpn)
  cat("genomeClustering: Starts with", P, "genomes\n")
  cat("   computing condition number...")
  lst <- svd(crossprod(Cpn))
  cond <- max(lst$d)/min(lst$d)
  if(cond <= cond.max){
    cat("No need for clustering, well conditioned matrix\n")
    return(Cpn)
  } else {
    cat("condition number is", cond, "\n")
  }
  cn <- colnames(Cpn)
  D <- corrDist(Cpn)
  tree <- hclust(D, method="single")
  cc <- 1
  while(cond>cond.max & cc<=length(steps)){
    cat("   iteration", cc, " ")
    clst <- cutree(tree, h=steps[cc])
    tc <- table(clst)
    cat("results in", length(tc), "clusters")
    idx <- which(tc==1)
    if(length(idx)>0){
      cls1 <- as.integer(names(tc[idx]))
      idd <- which(clst %in% cls1)
      Z1 <- Cpn[,idd]
      atr1 <- cn[idd]
    } else {
      Z1 <- Matrix(0, nrow=N, ncol=0)
      atr1 <- character(0)
    }
    idx <- which(tc>1)
    uclst <- as.integer(names(tc[idx]))
    nu <- length(uclst)
    Z2 <- Matrix(0, nrow=N, ncol=nu)
    atr2 <- character(nu)
    for(i in 1:nu){
      idd <- which(clst==uclst[i])
      Z2[,i] <- rowMeans(Cpn[,idd])
      atr2[i] <- paste(cn[idd], collapse="|")
      if((i/10)==round(i/10)) cat(".")
    }
    Z <- cbind(Z2,Z1)
    atr <- c(atr2, atr1)

    lst <- svd(crossprod(Z))
    cond <- max(lst$d)/min(lst$d)
    cat("condition number =", cond, "\n")
    cc <- cc +1
  }

  cat("New copy-number matrix has", ncol(Z), "columns and condition number:", cond, "\n")
  rownames(Z) <- rownames(Cpn)
  colnames(Z) <- atr
  return(Z)
}

### Local functions
corrDist <- function(X){
  require(Matrix)
  csd <- apply(X,2,sd)
  if(sum(csd==0)>0){
    stop("Columns", which(csd==0), "have no variance! Cannot compute correlation distances")
  }
  cm <- colMeans(X)
  n <- nrow(X)
  XX <- as.matrix(crossprod(X))
  CM <- cm %*% t(cm)
  CSD <- csd %*% t(csd)
  D <- 1-abs((XX-n*CM)/((n-1)*CSD))
  return(as.dist(D))
}
