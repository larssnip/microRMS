#' @name rmscols
#' @title Estimating microbiota composition
#'
#' @description Estimates the fraction of each genome in a sample, based on read counts and copy
#' numbers for each amplicon cluster.
#'
#' @param rms.obj A \code{list} with the matrices \code{Readcount.mat} and \code{Cpn.mat}, see details below.
#' @param trim Fraction of extreme readcounts to discard when fitting linear model.
#' @param reltol Relative stopping tolerance for the iterative constrained least square search.
#' @param verbose Logical, if TRUE text is written to the Console during computations.
#'
#' @details The \code{rms.obj} must be a list with the required data structures for performing a 
#' Constrained Ordinary Least Square estimation of abundances.
#' 
#' The \code{rms.obj} is typically constructed by the use of \code{\link{RMSobject}}. In this step
#' the copy number matrix \code{Cpn.mat} is constructed, based on RMS fragments in a selection
#' (database) of genomes.
#' 
#' In addition, a matrix of readcounts from one or more samples is required to be found in
#' \code{rms.obj$Readcount.mat}, see \code{\link{readMapper}} and tutorial for more details on this.
#'
#' @return A matrix with one row for each genome in \code{rms.obj$Cpn.mat} and one column for each 
#' sample in \code{rms.obj$Readcount.mat}. Each column contains the estimated relative abundance for 
#' all genomes in the corresponding sample. 
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}, \code{\link{readMapper}}.
#'
#' @importFrom Matrix crossprod Diagonal
#' @importFrom stats optim
#'
#' @examples
#'
#' @export rmscols
#'
rmscols <- function(rms.obj, trim = 0, reltol = 1e-6, verbose = TRUE){
  if(!exists("Readcount.mat", rms.obj)) stop("The rms.obj must contain a Readcount.mat")
  if(!exists("Cpn.mat", rms.obj)) stop("The rms.obj must contain a Cpn.mat")
  J <- ncol(rms.obj$Readcount.mat)
  C <- nrow(rms.obj$Cpn.mat)
  G <- ncol(rms.obj$Cpn.mat)
  beta.matrix <- matrix(0, nrow = G, ncol = J)
  rownames(beta.matrix) <- colnames(rms.obj$Cpn.mat)
  colnames(beta.matrix) <- colnames(rms.obj$Readcount.mat)
  if(trim > 0) residuals <- matrix(0, nrow = C, ncol = J)
  for(j in 1:J){
    tot.reads <- sum(rms.obj$Readcount.mat[,j])
    nc <- sum(rms.obj$Readcount.mat[,j] > 0)
    if(verbose) cat("De-convolving sample", rms.obj$Sample.tbl$sample_id[j],
                    "having", tot.reads, "reads mapped to", nc, "clusters...\n")
    if(tot.reads > 0){
      w <- rep(1,C)
      s.hat <- constrLS(rms.obj$Readcount.mat[,j], rms.obj$Cpn.mat, w, reltol, verbose)
      
      if(trim > 0){
        if(verbose) cat("Trimming residuals...\n")
        residuals[,j] <- as.numeric(rms.obj$Readcount.mat[,j] - rms.obj$Cpn.mat %*% s.hat)
        rsd <- sd(residuals[,j])
        for(g in 1:G){
          idx <- which(rms.obj$Cpn.mat[,g] > 0)
          qq <- quantile(residuals[idx,j], c(trim/2, 1-trim/2))
          idd <- which(residuals[idx,j] < qq[1] | residuals[idx,j] > qq[2])
          w[idx[idd]] <- 0
        }
        s.hat <- constrLS(rms.obj$Readcount.mat[,j], rms.obj$Cpn.mat, w, reltol, verbose)
      }
      beta.hat <- s.hat/sum(s.hat)
      beta.matrix[,j] <- as.numeric(beta.hat)
    }
  }
  return(beta.matrix)
}

### Local functions
constrLS <- function(y, X, w, reltol = 1e-6, verbose = FALSE){
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
  r <- w * (y - (X %*% b))
  return(sum(r^2))
}

