#' @name RMSobject
#' @title Constructing an RMS object
#'
#' @description Constructs an RMS object with information about a set of genomes.
#'
#' @param genome.tbl A table (data.frame or tibble) with genome information, see below.
#' @param frg.dir Path to folder with fragment fasta files.
#' @param identity The sequence identity for clustering fragments (0.0-1.0).
#' @param min.length Minimum fragment length (integer).
#' @param max.length Maximum fragment length (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param threads Number of threads to be used by \code{vsearch} (integer).
#'
#' @details The \code{genome.tbl} has a row for each genome to include in the RMS database.
#' There must be a column named \code{genome_file}, containing fasta filenames. These must be the
#' names of the fasta files containing the RMS fragments from each genome. Use \code{\link{getRMSfragments}}
#' to create these fasta files, ensuring the fasta headers follow the pattern
#' <genome.ID>_RMSx, where <genome.ID> is some text unique to each genome and x is some integer.
#' The \code{genome.tbl} may contain other columns as well, but \code{genome_file} is required.
#'
#' @return A list with the following objects: \code{Cluster.tbl}, \code{Cpn.mat}
#' and \code{Genome.tbl}.
#'
#' The \code{Cluster.tbl} is a \code{\link{tibble}} with data about all fragment clusters.
#' It contains columns with data about each cluster, including the centroid \code{Sequence}
#' and its \code{Header}, making it possible to write the table to a fasta file using
#' \code{\link{writeFasta}}.
#'
#' The \code{Cpn.mat} is the copy number matrix, implemented as a sparse dgeMatrix from the
#' \code{\link{Matrix}} package. It has one row for each fragment cluster and one column
#' for each genome. This is the central data structure for de-convolving the genome
#' content from read-count data, see \code{\link{rmscols}}.
#'
#' The \code{Genome.tbl} is a copy of the argument \code{genome.tbl}, but with columns 
#' \code{N_cluster} and \code{N_unique} added, containing the
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
#' 
#' @export RMSobject
#'
RMSobject <- function(genome.tbl, frg.dir, identity = 0.99, min.length = 30, max.length = 500, verbose = TRUE, threads = 1){
  if(length(grep("genome_id", colnames(genome.tbl))) == 0) stop("The genome.tbl must contain a column 'genome_id' with unique texts")
  if(length(genome.tbl$genome_id) != length(unique(genome.tbl$genome_id))) stop("The genome_id's must be unique for each genome (row)")
  if(length(grep("genome_file", colnames(genome.tbl))) == 0) stop("The genome.tbl must contain a column 'genome_file' with filenames")
  frg_files <- normalizePath(file.path(frg.dir, genome.tbl$genome_file))
  ok <- file.exists(frg_files)
  idx <- which(!ok)
  if(length(idx) > 0) stop("The file(s)", frg_files[idx], "does not exist")
  ok <- available.external("vsearch")
  all.frg <- tempfile(pattern = "all_frg", fileext = ".fasta")
  ok <- file.append(all.frg, frg_files)
  if(min(ok) == 0) stop("Could not copy all fragment fasta files from", frg.dir)

  ### The VSEARCH clustering
  if(verbose) cat("VSEARCH clustering of RMS fragments...\n")
  ctr.file <- tempfile(pattern = "centroide", fileext = ".fasta")
  uc.file <- tempfile(pattern = "uc", fileext = ".txt")
  cmd <- paste("vsearch",
               "--threads", threads,
               "--cluster_fast", all.frg,
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
    mutate(Genome.id = str_remove(Tag, "_RMS[0-9]+$")) -> uc.tbl
  
  uc.tbl %>%
    group_by(Cluster) %>%
    summarise(Members = str_c(Tag, collapse = ";"), N.genomes = length(unique(Genome.id))) %>% 
    right_join(centroids, by = "Cluster") %>% 
    mutate(Length = str_length(Sequence)) %>%
    mutate(GC = baseCount(Sequence, c("C", "G")) / Length) %>%
    select(Cluster, Length, GC, N.genomes, Members, Header, Sequence) -> cluster.tbl
  if(verbose) cat("done\n")
  
  if(verbose) cat("...the copy number matrix\n")
  cpn <- Matrix(0, nrow = nrow(cluster.tbl), ncol= nrow(genome.tbl))
  colnames(cpn) <- genome.tbl$genome_id
  rownames(cpn) <- cluster.tbl$Cluster
  for(j in 1:nrow(genome.tbl)){
    if(verbose) cat("genome", j, "/", nrow(genome.tbl), "\r")
    uc.tbl %>%
      filter(Genome.id == genome.tbl$genome_id[j]) %>%
      group_by(Cluster) %>%
      summarise(Count = n()) -> tbl
    cpn[match(tbl$Cluster, cluster.tbl$Cluster),j] <- tbl$Count
  }

  if(verbose) cat("\n...the genome table...")
  pa <- (cpn > 0)
  pau <- pa[Matrix::rowSums(pa) == 1,]
  genome.tbl %>% 
    mutate(N_clusters = Matrix::colSums(pa)) %>% 
    mutate(N_unique = Matrix::colSums(pau)) -> genome.tbl

  ### Cleaning up
  file.remove(all.frg, ctr.file, uc.file)
  rms.obj <- list(Cluster.tbl = cluster.tbl,
                  Cpn.mat = cpn,
                  Genome.tbl = genome.tbl)
  if(verbose) cat("done\n")
  return(rms.obj)
}

# local function
baseCount <- function(seq, bases){
  return(sapply(strsplit(seq, split = ""), function(x){sum(x %in% bases)}))
}

