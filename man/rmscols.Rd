% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composition.R
\name{rmscols}
\alias{rmscols}
\title{Estimating microbiota composition}
\usage{
rmscols(rms.obj, trim = 0, reltol = 1e-06, verbose = TRUE)
}
\arguments{
\item{rms.obj}{A \code{list} with the matrices \code{Readcount.mat} and \code{Cpn.mat}, see details below.}

\item{trim}{Fraction of extreme readcounts to discard when fitting linear model.}

\item{reltol}{Relative stopping tolerance for the iterative constrained least square search.}

\item{verbose}{Logical, if TRUE text is written to the Console during computations.}
}
\value{
A matrix with one row for each genome in \code{rms.obj$Cpn.mat} and one column for each 
sample in \code{rms.obj$Readcount.mat}. Each column contains the estimated relative abundance for 
all genomes in the corresponding sample.
}
\description{
Estimates the fraction of each genome in a sample, based on read counts and copy
numbers for each amplicon cluster.
}
\details{
The \code{rms.obj} must be a list with the required data structures for performing a 
Constrained Ordinary Least Square estimation of abundances.

The \code{rms.obj} is typically constructed by the use of \code{\link{RMSobject}}. In this step
the copy number matrix \code{Cpn.mat} is constructed, based on RMS fragments in a selection
(database) of genomes.

In addition, a matrix of readcounts from one or more samples is required to be found in
\code{rms.obj$Readcount.mat}, see \code{\link{readMapper}} and tutorial for more details on this.
}
\seealso{
\code{\link{RMSobject}}, \code{\link{readMapper}}.
}
\author{
Lars Snipen.
}
