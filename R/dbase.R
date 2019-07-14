#' @name RMSdbase
#' @title Constructing RMS database
#'
#' @description Constructs a database of information about a set of RMS fragments.
#'
#' @param fragment.file A FASTA-file with all the RMS fragments (text).
#' @param identity The sequence identity for clustering fragments (0.0-1.0).
#' @param min.length Minimum fragment length (integer).
#' @param max.length Maximum fragment length (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param threads Number of threads to be used by \code{vsearch} (integer).
#'
#' @details The \code{fragment.file} must be a FASTA-file containing all fragments to be clustered.
#' The first token of each Header-line must indicate the genome of its origin. It must follow the
#' pattern <genome.ID>_RMSx, where <genome.ID> is some text unique to each genome, and x
#' is some integer. See \code{\link{getRMSfragments}} for how to make such FASTA-files.
#'
#' @return A list with the following objects: \code{Cluster.tbl}, \code{Cpn.mat}
#' and \code{Genome.tbl}.
#'
#' The \code{Cluster.tbl} is a \code{\link{tibble}} with data about all fragment clusters.
#' It contains columns with data about each cluster, including the centroid \code{Sequence}
#' and its \code{Header}, making it possible to write the table to a FASTA-file using
#' \code{\link{writeFasta}}.
#'
#' The \code{Cpn.mat} is the copy number matrix, implemented as a sparse dgeMatrix from the
#' \code{\link{Matrix}} package. It has one row for each fragment cluster and one column
#' for each genome. This is the central data structure for de-convolving the genome
#' content from read-count data, see \code{\link{rmscols}}.
#'
#' The \code{Genome.tbl} is a \code{\link{tibble}} with a row for each genome, their
#' number of clusters and the number of unique fragment clusters to each genome.
#'
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{getRMSfragments}}.
#'
#' @importFrom microseq readFasta
#' @importFrom data.table fread
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom stringr word str_length str_remove str_c
#' @importFrom dplyr mutate filter select group_by slice summarise
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' # Some fragment files in this package
#' xpth <- file.path(path.package("microrms"),"extdata")
#' frg.files <- file.path(xpth, list.files(xpth, pattern = ".frg"))
#' 
#' # Merging all fragments into one file
#' tmp.frg.file <- tempfile(pattern = "all", fileext = ".frg")
#' ok <- file.append(tmp.frg.file, frg.files)
#' 
#' # RMSdbase needs the external software vsearch
#' rms.db <- RMSdbase(tmp.frg.file)
#' print(rms.db$Genome.tbl)
#' 
#' # clean up
#' ok <- file.remove(tmp.frg.file)
#' }
#' 
#' @export RMSdbase
#'
RMSdbase <- function(fragment.file, identity = 0.99, min.length = 30, max.length = 500,
                     verbose = TRUE, threads = 1){
  quiet <- "--quiet"
  if(verbose) quiet <- ""

  ### The VSEARCH clustering
  if(verbose) cat("VSEARCH clustering...\n")
  ctr.file <- tempfile(pattern = "centroide", fileext = ".fasta")
  uc.file <- tempfile(pattern = "uc", fileext = ".txt")
  cmd <- paste("vsearch", quiet,
               "--threads", threads,
               "--cluster_fast", fragment.file,
               "--id", identity,
               "--minseqlength", min.length,
               "--maxseqlength", max.length,
               "--iddef", "2",
               "--qmask", "none",
               "--strand plus --sizeout",
               "--relabel CLST",
               "--relabel_keep",
               "--uc", uc.file,
               "--centroids", ctr.file)
  system(cmd)
  readFasta(ctr.file) %>%
    mutate(Cluster = word(Header, 1, 1, sep = ";")) -> centroids
  if(verbose) cat("...produced", nrow(centroids), "clusters\n...the cluster table...")
  fread(uc.file, sep = "\t", header = F, drop = c(3,4,5,6,7,8,10), data.table = F) %>%
    filter(V1 != "C") %>%
    select(Cluster = V2, Tag = V9) %>%
    mutate(Cluster = Cluster + 1) %>%
    mutate(Cluster = str_remove(str_c("CLST", format(Cluster, scientific = F)), " +")) %>%
    mutate(Genome.ID = str_remove(Tag, "_RMS[0-9]+$")) -> uc.tbl

  uc.tbl %>%
    group_by(Cluster) %>%
    summarise(Members = str_c(Tag, collapse = ";"), N.genomes = length(unique(Genome.ID))) %>%
    slice(match(centroids$Cluster, Cluster)) %>%
    mutate(Header = centroids$Header, Sequence = centroids$Sequence) %>%
    mutate(Length = str_length(Sequence)) %>%
    mutate(GC = baseCount(Sequence, c("C", "G")) / Length) %>%
    select(Cluster, Length, GC, N.genomes, Members, Header, Sequence) -> cluster.tbl
  if(verbose) cat("done\n")

  if(verbose) cat("...the copy number matrix\n")
  ug <- unique(uc.tbl$Genome.ID)
  cpn <- Matrix(0, nrow = nrow(cluster.tbl), ncol= length(ug))
  colnames(cpn) <- ug
  rownames(cpn) <- cluster.tbl$Cluster
  for(j in 1:length(ug)){
    if(verbose) cat("genome", j, "/", length(ug), "\r")
    uc.tbl %>%
      filter(Genome.ID == ug[j]) %>%
      group_by(Cluster) %>%
      summarise(Count = n()) -> tbl
    cpn[match(tbl$Cluster, cluster.tbl$Cluster),j] <- tbl$Count
  }

  if(verbose) cat("\n...the genome table...")
  pa <- (cpn > 0)
  pau <- pa[Matrix::rowSums(pa) == 1,]
  genome.tbl <- tibble(Genome.ID = colnames(cpn),
                       N.clusters = Matrix::colSums(pa),
                       N.unique = Matrix::colSums(pau))

  ### Cleaning up
  file.remove(ctr.file, uc.file)
  rmsdb.obj <- list(Cluster.tbl = cluster.tbl,
                    Cpn.mat = cpn,
                    Genome.tbl = genome.tbl)
  if(verbose) cat("done\n")
  return(rmsdb.obj)
}

# local function
baseCount <- function(seq, bases){
  return(sapply(strsplit(seq, split = ""), function(x){sum(x %in% bases)}))
}




