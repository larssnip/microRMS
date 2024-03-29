% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbase.R
\name{RMSobject}
\alias{RMSobject}
\title{Constructing an RMS object}
\usage{
RMSobject(
  genome.tbl,
  frg.dir,
  vsearch.exe = "vsearch",
  identity = 0.99,
  min.length = 30,
  max.length = 500,
  verbose = TRUE,
  threads = 1,
  tmp.dir = "tmp"
)
}
\arguments{
\item{genome.tbl}{A table (data.frame or tibble) with genome information, see below.}

\item{frg.dir}{Path to folder with fragment fasta files.}

\item{vsearch.exe}{Text with the VSEARCH executable command.}

\item{identity}{The sequence identity for clustering fragments (0.0-1.0).}

\item{min.length}{Minimum fragment length (integer).}

\item{max.length}{Maximum fragment length (integer).}

\item{verbose}{Turn on/off output text during processing (logical).}

\item{threads}{Number of threads to be used by \code{vsearch} (integer).}

\item{tmp.dir}{Name of folder for temporary output, will be created if not already existing.}
}
\value{
A list with the following objects: \code{Cluster.tbl}, \code{Cpn.mat}
and \code{Genome.tbl}.

The \code{Cluster.tbl} is a \code{\link{tibble}} with data about all fragment clusters.
It contains columns with data about each cluster, including the centroid \code{Sequence}
and its \code{Header}, making it possible to write the table to a fasta file using
\code{\link{writeFasta}}.

The \code{Cpn.mat} is the copy number matrix, implemented as a sparse dgeMatrix from the
\code{\link{Matrix}} package. It has one row for each fragment cluster and one column
for each genome. This is the central data structure for de-convolving the genome
content from read-count data, see \code{\link{rmscols}}.

The \code{Genome.tbl} is a copy of the argument \code{genome.tbl}, but with columns 
\code{N_cluster} and \code{N_unique} added, containing the
number of clusters and the number of unique fragment clusters to each genome.
}
\description{
Constructs an RMS object with information about a set of genomes.
}
\details{
The \code{genome.tbl} has a row for each genome to include in the RMS database.
There must be a column named \code{genome_file}, containing fasta filenames. These must be the
names of the fasta files containing the RMS fragments from each genome. Use \code{\link{getRMSfragments}}
to create these fasta files, ensuring the fasta headers follow the pattern
<genome.ID>_RMSx, where <genome.ID> is some text unique to each genome and x is some integer.
The \code{genome.tbl} may contain other columns as well, but \code{genome_file} is required.

The \code{vsearch.exe} is the exact command to invoke the VSEARCH software. This is normally just "vsearch", 
but if you run this as a singularity container (or any other container) it may be something like
"srun singularity exec <container_name> vsearch".
}
\seealso{
\code{\link{getRMSfragments}}.
}
\author{
Lars Snipen.
}
