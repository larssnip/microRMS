#' @name getRMSfragments
#' @title Retrieving RMS fragments from genomes
#'
#' @description Retrieves a set of fragments from a genome, given restriction enzyme cutting motifs.
#'
#' @param genome A table (fasta object) with genome data.
#' @param genome.id Unique identifier for each genome, will be added to FASTA-headers (text).
#' @param min.length Minimum fragment length (integer).
#' @param max.length Maximum fragment length (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param left Text with first, long, restriction enzyme cut motif (text).
#' @param right Text with second, short, restriction enzyme cut motif (text).
#'
#' @details This function is used to find and retrieve all RMS fragments from a \code{genome}.
#' This is a \code{\link{tibble}} with sequence data in FASTA format, see \code{\link{readFasta}}.
#' In addition, a \code{genome.id} is required,
#' which is a text unique to each genome to be analyzed. This \code{genome.id} will be added to the
#' fasta headers of the output, and all headers start with the token <genome.id>_RMSx, where x is
#' an integer (1,2,...,). This first token is followed by a blank. This ensures that all first tokens
#' are unique and that the genome of its origin is indicated for all fragments
#'
#' The default restriction enzymes are EcoRI and MseI, with cutting motifs \code{"G|AATTC"} and
#' \code{"T|TAA"}, respectively. The vertical bar indicates the cut site in the motif. 
#'
#' @return A \code{\link{tibble}} with with all fragment sequences (5'-3') in FASTA format.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSobject}}.
#'
#' @importFrom microseq readFasta gff2fasta
#' @importFrom stringr str_c str_sub str_length str_remove str_remove_all
#' @importFrom dplyr mutate filter
#'
#' @examples
#' # A small genome in this package
#' xpth <- file.path(path.package("microrms"),"extdata")
#' genome.file <- file.path(xpth,"GCF_000009605.1_ASM960v1_genomic.fna")
#' 
#' # Read genome, find fragments
#' gnm <- readFasta(genome.file)
#' frg <- getRMSfragments(gnm, "genome_1")
#'
#' # Write to file with writeFasta(frg, out.file = <filename>)
#' 
#' @export getRMSfragments
#'
getRMSfragments <- function(genome, genome.id, min.length = 30, max.length = 500,
                            left = "G|AATTC", right = "T|TAA", verbose = TRUE){
  if(verbose) cat("Genome", genome.id, ": ")
  lft <- str_remove_all(left, "\\|")
  rght <- str_remove_all(right, "\\|")
  gff <- getRMS(genome, lft, rght)
  if(nrow(gff)>0){
    ct.lft <- str_length(str_remove(left, "\\|.+"))
    ct.rght <- str_length(str_remove(right, "\\|.+"))
    fsa <- gff2fasta(gff, genome) %>% 
      mutate(Sequence = str_sub(Sequence, ct.lft+1, -(ct.rght+1))) %>% 
      mutate(Length = str_length(Sequence)) %>%
      filter((Length >= min.length) & (Length <= max.length)) %>% 
      mutate(Header = str_c(str_c(genome.id, str_c("RMS", 1:n()), sep = "_"), Header, sep = " "))
    if(nrow(fsa) == 0){
      if(verbose) cat("found no RMS-fragments within min and max length!\n")
      fsa <- tibble(Header = NULL, Sequence = NULL)
    } else {
      if(verbose) cat("found", nrow(fsa), "RMS-fragments\n")
    }
  } else {
    if(verbose) cat("found no RMS-fragments!\n")
    fsa <- tibble(Header = NULL, Sequence = NULL)
  }
  return(fsa)
}



### Local function
getRMS <- function(genome, left, right){
  gff <- tibble(Seqid      = NULL,
                Source     = NULL,
                Type       = NULL,
                Start      = NULL,
                End        = NULL,
                Score      = NULL,
                Strand     = NULL,
                Phase      = NULL,
                Attributes = NULL)
  ### looping over genome-sequences
  for(i in 1:nrow(genome)){
    left.m <- str_locate_all(genome$Sequence[i], pattern = left)[[1]]
    right.m <- str_locate_all(genome$Sequence[i], pattern = right)[[1]]
    if(length(left.m) > 0 & length(right.m) > 0){
      ### positive strand
      lft <- left.m[,1]
      rght <- right.m[,2]
      lft.rght <- t(sapply(lft, function(ll){
        ll2 <- suppressWarnings(min(lft[lft > ll]))
        ss <- suppressWarnings(min(rght[rght > ll]))
        return(c(ll, ifelse(ss < ll2, ss, Inf)))
      }))
      nf <- nrow(lft.rght)
      if(nf > 0){
        gff <- rbind(gff, tibble(Seqid      = rep(word(genome$Header[i], 1, 1), nf),
                                 Source     = NA,
                                 Type       = rep("RMS_fragment", nf),
                                 Start      = lft.rght[,1],
                                 End        = lft.rght[,2],
                                 Score      = NA,
                                 Strand     = rep("+", nf),
                                 Phase      = NA,
                                 Attributes = NA))
      }
      ### negative strand
      lft <- left.m[,2]
      rght <- right.m[,1]
      lft.rght <- t(sapply(lft, function(ll){
        ll2 <- suppressWarnings(max(lft[lft < ll]))
        ss <- suppressWarnings(max(rght[rght < ll]))
        return(c(ifelse(ss > ll2, ss, -Inf), ll))
      }))
      nn <- nrow(lft.rght)
      if(nn > 0){
        gff <- rbind(gff, tibble(Seqid      = rep(word(genome$Header[i], 1, 1), nn),
                                 Source     = NA,
                                 Type       = rep("RMS_fragment", nf),
                                 Start      = lft.rght[,1],
                                 End        = lft.rght[,2],
                                 Score      = NA,
                                 Strand     = rep("-", nf),
                                 Phase      = NA,
                                 Attributes = NA))
      } # end if
    } # end if
  }  # end for
  idx <- which(is.finite(gff$Start) & is.finite(gff$End))
  return(gff[idx,])
}

