#' @name rmscols
#' @title Estimating microbiota composition
#'
#' @description Estimates the fraction of each genome in a sample, based on read counts and copy
#' numbers for each amplicon cluster.
#'
#' @param Y A matrix of read counts, a column for each sample.
#' @param X A matrix of amplicon cluster copy numbers, a column for each genome (cluster).
#' @param trim Fraction of extreme readcounts to discard when fitting linear model.
#' @param reltol Relative stopping tolerance for the iterative constrained leats square search
#' @param verbose Logical, if TRUE text is written to the Console during computations.
#'
#' @details more here.
#'
#' @return A matrix of the fraction for each genome.
#'
#' @author Lars Snipen.
#'
#' @seealso more here.
#'
#' @importFrom Matrix crossprod Diagonal
#' @importFrom stats optim
#'
#' @examples more here.
#'
#' @export rmscols
#'
rmscols <- function(Y, X, trim = 0, reltol = 1e-6, verbose = TRUE){
  J <- ncol(Y)
  C <- nrow(X)
  G <- ncol(X)
  beta.matrix <- matrix(0, nrow = G, ncol = J)
  rownames(beta.matrix) <- colnames(X)
  colnames(beta.matrix) <- colnames(Y)
  if(trim > 0) residuals <- matrix(0, nrow = C, ncol = J)
  for(j in 1:J){
    if(verbose) cat("Deconvolving sample", j, "...\n")
    w <- rep(1,C)
    s.hat <- constrLS(Y[,j], X, w, reltol, verbose)
    
    if(trim > 0){
      if(verbose) cat("Trimming residuals...\n")
      residuals[,j] <- as.numeric(Y[,j] - X %*% s.hat)
      rsd <- sd(residuals[,j])
      for(g in 1:G){
        idx <- which(X[,g] > 0)
        qq <- quantile(residuals[idx,j], c(trim/2, 1-trim/2))
        idd <- which(residuals[idx,j] < qq[1] | residuals[idx,j] > qq[2])
        w[idx[idd]] <- 0
      }
      s.hat <- constrLS(Y[,j], X, w, reltol, verbose)
    }
    beta.hat <- s.hat/sum(s.hat)
    beta.matrix[,j] <- as.numeric(beta.hat)
  }
  return(beta.matrix)
}

### Local functions
constrLS <- function(y, X, w, reltol = 1e-6, verbose = FALSE){
  require(Matrix)
  p <- ncol(X)

  # Initial estimate
  if(verbose) cat("   initial estimate...\n")
  W <- Diagonal(x = w)
  theta0 <- Matrix::solve(crossprod(X, crossprod(W, X)), crossprod(X, y))
  theta0 <- pmax(1e-10, theta0)
  names(theta0) <- colnames(X)
  
  # Constrained estimate
  ctl <- list(factr = reltol/.Machine$double.eps)
  if(verbose){
    cat("   constrained optimization...\n")
    ctl <- list(trace = 1, REPORT = 1, factr = reltol/.Machine$double.eps)
  }
  lst <- optim(theta0, fn = objectFun, gr = NULL, y, X, w,
               method = "L-BFGS-B", lower = rep(0, p),
               control = ctl)
  return(lst$par)
}
objectFun <- function(b, y, X, w){
  require(Matrix)
  r <- w * (y - (X %*% b))
  return(sum(r^2))
}

