#' @name readMapper
#' @title Mapping reads
#'
#' @description Mapping reads to clusters in an RMS object using VSEARCH.
#'
#' @param rms.obj A \code{list} with the tables \code{Sample.tbl} and \code{Cluster.tbl}, see details below.
#' @param fa.dir A path to where the fasta files with reads are located.
#' @param vsearch.exe Text with the VSEARCH executable command.
#' @param identity Identity threshold for mapping (0.0-1.0).
#' @param threads Number of threads to be used by vsearch (integer).
#' @param min.length Minimum fragment length used by vsearch (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param tmp.dir Folder for temporary files.
#' 
#' @details The \code{rms.obj} must be a list with the required data structures for mapping
#' reads, i.e. it must contain a \code{Sample.tbl}, see \code{\link{addSampleTable}}. Note
#' that this table must have a column named \code{fasta_file} naming the fasta files with
#' the reads to be mapped, one file for each sample. The argument \code{fa.dir} is used to
#' specify the path to these files.
#' 
#' The \code{rms.obj} must also contain a \code{Cluster.tbl}, see \code{\link{RMSobject}}.
#' The reads from each sample will be mapped to the cluster
#' sequences, using the \code{identity} threshold.
#' 
#' This results in a matrix of read counts, with one row for each fragment cluster and one column
#' for each sample. The column names of this matrix are the \code{sample_id}
#' texts in the \code{Sample.tbl}. The first token in the \code{Header} of the \code{Cluster.tbl} 
#' are used as row names.
#'
#' The \code{vsearch.exe} is the exact command to invoke the VSEARCH software. This is normally just "vsearch", 
#' but if you run this as a singularity container (or any other container) it may be something like
#' "singularity exec <container_name> vsearch".
#'
#' @return An RMS object with the matrix \code{Readcount.mat} added. Also, two new columns, 
#' \code{reads_total} and \code{reads_mapped}, have been added to the \code{Sample.tbl}.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}, \code{\link{rmscols}}.
#'
#' @importFrom microseq readFasta
#' @importFrom stringr word
#'
#' @examples
#'
#' @export readMapper
#'
readMapper <- function(rms.obj, fa.dir, vsearch.exe = "vsearch", identity = 0.99,
                       threads = 1, min.length = 30, verbose = TRUE,
                       tmp.dir = "tmp"){
  if(!exists("Sample.tbl", where = rms.obj)) stop("The rms.obj must contain a Sample.tbl")
  if(length(grep("sample_id", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column sample_id")
  if(length(grep("fasta_file", colnames(rms.obj$Sample.tbl))) == 0) stop("The rms.obj$Sample.tbl must contain a column fasta_file")
  fa.files <- list.files(fa.dir)
  idx <- which(!(rms.obj$Sample.tbl$fasta_file %in% fa.files))
  if(length(idx) > 0) stop("Missing fasta files:", str_c(rms.obj$Sample.tbl$fasta_file[idx], sep = ","))
  if(!exists("Cluster.tbl", where = rms.obj)) stop("The rms.obj must contain a Cluster.tbl")
  if(length(grep("Header", colnames(rms.obj$Cluster.tbl))) == 0) stop("The rms.obj$Cluster.tbl must contain a column Header")
  if(length(grep("Sequence", colnames(rms.obj$Cluster.tbl))) == 0) stop("The rms.obj$Cluster.tbl must contain a column Sequence")
  ok <- available.external(vsearch.exe)
  tags <- word(rms.obj$Cluster.tbl$Header, 1, 1, sep = ";")
  RMS.counts <- matrix(0, nrow = length(tags), ncol = nrow(rms.obj$Sample.tbl))
  rownames(RMS.counts) <- tags
  colnames(RMS.counts) <- rms.obj$Sample.tbl$sample_id
  centroids.file <- file.path(tmp.dir, "centroids.fasta")
  writeFasta(rms.obj$Cluster.tbl, out.file = centroids.file)
  tab.file <- file.path(tmp.dir, "rmstab.txt")
  tot <- numeric(nrow(rms.obj$Sample.tbl))
  for(i in 1:nrow(rms.obj$Sample.tbl)){
    if(verbose) cat("Mapping reads from sample", rms.obj$Sample.tbl$sample_id[i], "...\n")
    cmd <- paste(vsearch.exe,
                 "--threads", threads,
                 "--usearch_global", file.path(fa.dir, rms.obj$Sample.tbl$fasta_file[i]),
                 "--db", centroids.file,
                 "--id", identity,
                 "--iddef", "2",
                 "--minseqlength", min.length,
                 "--sizein --sizeout",
                 "--otutabout", tab.file)
    system(cmd)
    lines <- readLines(tab.file)
    if(length(lines) > 1){
      M <- str_split(lines[-1], pattern = "\t", simplify = T)
      idx <- match(M[,1], tags)
      RMS.counts[idx,i] <- as.numeric(M[,2])
    }
    readFasta(file.path(fa.dir, rms.obj$Sample.tbl$fasta_file[i])) %>% 
      mutate(size = as.numeric(str_remove(str_extract(Header, pattern = "size=[0-9]+"), "size="))) %>% 
      summarize(size_sum = sum(size)) -> stb
    tot[i] <- stb$size_sum
  }
  ok <- file.remove(tab.file, centroids.file)
  rms.obj$Sample.tbl %>% 
    mutate(reads_total = tot, reads_mapped = colSums(RMS.counts)) -> rms.obj$Sample.tbl
  file.remove(centroids.file, tab.file)
  
  return(c(rms.obj, list(Readcount.mat = RMS.counts)))
}


#' @name addSampleTable
#' @title Adds sample table
#'
#' @description Adds a sample table to an rms object
#'
#' @param rms.obj An \code{list} with RMS data structures, see \code{\link{RMSobject}}.
#' @param sample.tbl A table with sample metadata (data.frame or tibble).
#'
#' @details This small function just adds a table with metadata about samples to an 
#' already existing \code{rms.obj}. The latter is a \code{list}, and this function only
#' ensures the added element is correctly named \code{Sample.tbl}.
#' 
#' If the supplied \code{rms.obj} already contains a \code{Sample.tbl}, it is simply replaced by the new \code{sample.tbl}.
#' 
#' The \code{sample.tbl} must contain at least the two columns \code{sample_id} and 
#' \code{fasta_file}. The first is a unique text to identify each sample, the latter is
#' the name of the fasta file with processed reads, see \code{\link{readMapper}}.
#' 
#' @return A \code{list} similar to the input \code{rms.obj}, but with one new element added.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}.
#'
#' @examples
#'
#' @export addSampleTable
#'
addSampleTable <- function(rms.obj, sample.tbl){
  if(exists("Sample.tbl", where = rms.obj)){
    rms.obj$Sample.tbl <- sample.tbl
  } else {
    rms.obj <- c(rms.obj, list(Sample.tbl = sample.tbl))
  }
  return(rms.obj)
}
