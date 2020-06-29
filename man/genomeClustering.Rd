% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{genomeClustering}
\alias{genomeClustering}
\title{Clustering genomes based on shared RMS fragments}
\usage{
genomeClustering(rms.obj, max.cond = 1000)
}
\arguments{
\item{rms.obj}{A \code{list} with RMS data structures, see \code{\link{RMSobject}}.}

\item{max.cond}{Maximum condition value tolerated, see Details.}
}
\value{
A updated RMS object, where \code{Genome.tbl} has an additional column named \code{members_genome_id}. This \code{Genome.tbl} 
typically has fewer rows than the original, one for each cluster, and this new column indicates which of the original
genomes are grouped into each cluster. Each cluster is represented by one of the original genomes (the cluster medoide).
The other objects inside the RMS object have also been updated accordingly.
}
\description{
Alters an RMS-object by clustering genomes who are too similar to distinguish.
}
\details{
This function will cluster genomes based on how similar they are in RMS fragment content.
The input is an RMS object (see \code{\link{RMSobject}}), with a \code{Genome.tbl} listing the genomes
and a \code{Cpn.mat} holding the information about RMS fragments in each genome. The latter has a 
column for genome.

The distance metric is 1 minus the correlation between the columns of \code{Cpn.mat}, i.e. if all RMS
fragments found in genome A are exactly those found in genome B, the distance between them is 0. Conversely,
if none of the fragments from genome A are found in genome B, the distance becomes 2 (maximum).

Clustering is done by hierchical clustering with complete linkage. Cutting the dendrogram tree at the various
heights produce a unique clustering of the genomes. Each of the produces a new \code{Cpn.mat}, with one column for each
genome cluster. The closer the genomes are to each other, the larger the condition number of this matrix. Thus, the
argument \code{max.cond} specifies the maximum tolerated condition value. A smaller value means harder clustering (fewer
genome clusters) but also more stable estimates of the abundance of each genome cluster. Conversely, a larger tolerated
condition value results in more genome clusters, but also more uncertain abundance estimates of each. 

More technical details: The problem of too similar genomes is reflected in a close to singular covariance matrix when
de-convolving the genome abundances. Clustering by correlation distance is, in theory, no guarantee against this. 
Two (or more) fairly uncorrelated genomes may still combine into something very correlated with a third genome. However,
in reality this rarely happens. Use the \code{\link{conditionValue}} function on the resulting \code{Cpn.mat} to
see if the clustering resulted in a lowering of the condition value to a tolerable size (e.g. around 1000 or less).
}
\examples{

}
\seealso{
\code{\link{conditionValue}}.
}
\author{
Lars Snipen.
}