% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{genomeClustering}
\alias{genomeClustering}
\title{Clustering genomes based on shared amplicon patterns}
\usage{
genomeClustering(Cpn, steps = cumsum(0.01 * 1.1^(0:15)),
  cond.max = 1000)
}
\arguments{
\item{Cpn}{A matrix of amplicon cluster copy numbers, a column for each genome (cluster).}

\item{steps}{The correlation distances to visit during clustering.}

\item{cond.max}{The maximum condition number tolerated.}
}
\value{
A new copy number matrix, with fewer columns than the original if clustering was needed
to get the condition number below \code{cond.max}.
}
\description{
Clusters genomes to obtain a copy number matrix with a maximum
condition number in order to estimate composition.
}
\details{
more here.
}
\examples{
more here.

}
\seealso{
more here.
}
\author{
Lars Snipen.
}